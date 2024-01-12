classdef NmpcControl < handle
    
    properties
        solver
        nx, nu, N
        nlp_x0
        nlp_lbx, nlp_ubx
        nlp_lbg, nlp_ubg
        nlp_p
        
        T_opt
        sol
        idx
        
        % Delay compensation
        rocket
        expected_delay
        mem_u
        
        % Warmstart
        nlp_lam_x0
        nlp_lam_g0
    end
    
    methods
        function obj = NmpcControl(rocket, tf, expected_delay)
            if nargin < 3, expected_delay = 0; end
            
            import casadi.*
            
            N_segs = ceil(tf/rocket.Ts); % MPC horizon
            nx = 12; % Number of states
            nu = 4;  % Number of inputs
            
            % Decision variables (symbolic)
            N = N_segs + 1; % Index of last point
            X_sym = SX.sym('X_sym', nx, N); % state trajectory
            U_sym = SX.sym('U_sym', nu, N-1); % control trajectory)
            
            % Parameters (symbolic)
            x0_sym  = SX.sym('x0_sym', nx, 1);  % initial state
            ref_sym = SX.sym('ref_sym', 4, 1);  % target position
            
            % Default state and input constraints
            ubx = inf(nx, 1);
            lbx = -inf(nx, 1);
            ubu = inf(nu, 1);
            lbu = -inf(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Cost
            cost = 0;

            % Equality constraints (Casadi SX), each entry == 0
            eq_constr = [ ; ];

            % Inequality constraints (Casadi SX), each entry <= 0
            ineq_constr = [ ; ];

            % --- Model Parameters ---
            
            % Define the state weighting matrix Q, emphasizing orientation and position
            Q = diag([ ...
                        1.5 1.5 0.75  ...          % Angular velocities wx, wy, wz
                        0 0 3  ...                 % Orientation angles alpha, beta, gamma
                        20 20 40  ...              % Linear velocities vx, vy, vz
                        70 70 70 ...               % Position coordinates x, y, z
                    ]);
            
            % Define the control input weighting matrix R
            R = diag([ ...
                        0.5 0.5 ...               % Control inputs d1, d2
                        0.001 0.001 ...           % Power inputs Pavg, Pdiff
                    ]);
            
            % --- System Dynamics ---
            
            % Define the system dynamics function
            f = @(x, u) rocket.f(x, u);
            % Define the function to compute the next state using the RK4 integrator
            next_X = @(x, u) RK4(x, u, rocket.Ts, f);
            
            % --- Reference Tracking ---
            
            % Initialize the tracking matrix T to transform states into the reference frame
            T = zeros(nx, 4);
            T(10, 1) = 1; % Mapping for x coordinate
            T(11, 2) = 1; % Mapping for y coordinate
            T(12, 3) = 1; % Mapping for z coordinate
            T( 6, 4) = 1; % Mapping for gamma angle
            
            % --- Terminal Cost ---
            
            % Linearize the system around the trim condition for terminal cost computation
            [xs, us] = rocket.trim();
            linearSys = rocket.linearize(xs, us);
            % Discretize the linear system with the specified sampling time
            linearSys = c2d(linearSys,rocket.Ts);
            % Compute the terminal cost matrix Qf using the DLQR method
            [~, Qf, ~] = dlqr(linearSys.A, linearSys.B, Q, R);
            
            % Set the initial equality constraint
            eq_constr = [eq_constr; X_sym(:, 1) - x0_sym];
            % Build the cost and dynamics equality constraints over the prediction horizon
            for i=1:N-1
                % Add control input cost
                cost = cost + (U_sym(:, i) - us)' * R * (U_sym(:, i) - us);
                % Add state deviation cost
                cost = cost + (X_sym(:, i) - T * ref_sym)' * Q * (X_sym(:, i) - T * ref_sym);
                % Add dynamics
                eq_constr = [eq_constr; X_sym(:, i + 1) - next_X(X_sym(:, i), U_sym(:, i))];
            end
            % Add terminal state cost, transformed by T
            cost = cost + (X_sym(:,end) - T * ref_sym)' * Qf * (X_sym(:,end) - T * ref_sym);
            
            % --- Constraints ---
            
            % Specify state and input bounds, overwriting any previous definitions
            
            % State bounds for beta angle
            lbx(5, 1) = -deg2rad(75);
            ubx(5, 1) =  deg2rad(75);            
            
            % Input bounds for control and power inputs
            lbu(:, 1) = [-deg2rad(15), -deg2rad(15), 50, -20];
            ubu(:, 1) = [ deg2rad(15),  deg2rad(15), 80,  20];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ---- Assemble NLP ------
            nlp_x = [X_sym(:); U_sym(:)];
            nlp_p = [x0_sym; ref_sym];
            nlp_f = cost;
            nlp_g = [eq_constr; ineq_constr];
            
            nlp = struct('x', nlp_x, 'p', nlp_p, 'f', nlp_f, 'g', nlp_g);
            
            % ---- Setup solver ------
            opts = struct('ipopt', struct('print_level', 0), 'print_time', false);
            obj.solver = nlpsol('solver', 'ipopt', nlp, opts);
            
            % ---- Assemble NLP bounds ----
            obj.nlp_x0  = zeros(size(nlp_x));
            
            obj.nlp_ubx = [repmat(ubx, N, 1); repmat(ubu, (N-1), 1)];
            obj.nlp_lbx = [repmat(lbx, N, 1); repmat(lbu, (N-1), 1)];
            
            obj.nlp_ubg = [zeros(size(eq_constr)); zeros(size(ineq_constr))];
            obj.nlp_lbg = [zeros(size(eq_constr)); -inf(size(ineq_constr))];
            
            obj.nlp_p = [zeros(size(x0_sym)); zeros(size(ref_sym))];
            
            obj.nlp_lam_x0 = [];
            obj.nlp_lam_g0 = [];
            
            obj.nx = nx;
            obj.nu = nu;
            obj.N = N;
            obj.T_opt = linspace(0, N * rocket.Ts, N);
            
            obj.idx.X = [1, obj.N * obj.nx];
            obj.idx.U = obj.idx.X(2) + [1, (obj.N-1) * obj.nu];
            obj.idx.u0 = obj.idx.U(1) + [0, obj.nu-1];
            
            % Members for delay compensation
            obj.rocket = rocket;
            obj.expected_delay = expected_delay;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            [~, u_init] = rocket.trim();
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.mem_u = repmat(u_init, 1, expected_delay);
        end
        
        function [u, T_opt, X_opt, U_opt] = get_u(obj, x0, ref)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            delay = obj.expected_delay;
            mem_u = obj.mem_u;
            
            % Delay compensation: Predict x0 delay timesteps later.
            % Simulate x_ for 'delay' timesteps
            x_ = x0;
            % ...
       
            x0 = x_;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute solution from x0
            obj.solve(x0, ref);
            
            % Evaluate u0
            nlp_x = obj.sol.x;
            id = obj.idx.u0;
            u = full( nlp_x(id(1):id(2)) );      
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Delay compensation: Save current u
            if obj.expected_delay > 0
               % obj.mem_u = ...
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargout > 1, T_opt = obj.get_T_opt(); end
            if nargout > 2, X_opt = obj.get_X_opt(); end
            if nargout > 3, U_opt = obj.get_U_opt(); end
            return
            
            % Additional evaluation
            % Complete trajectory
            % % X_opt = full(reshape(nlp_x(idx_X(1):idx_X(2)), obj.nx, obj.N));
            % % U_opt = full(reshape(nlp_x(idx_U(1):idx_U(2)), obj.nu, obj.N - 1));
            % %
            % % cost_opt = full(sol.f);
            % % constr_opt = full(sol.g);
            % %
            % % stats = obj.solver.stats;
        end
        
        function solve(obj, x0, ref)
            
            % ---- Set the initial state and reference ----
            obj.nlp_p = [x0; ref];     % Initial condition
            obj.nlp_x0(1:obj.nx) = x0; % Initial guess consistent
            
            % ---- Solve the optimization problem ----
            args = {'x0', obj.nlp_x0, ...
                'lbg', obj.nlp_lbg, ...
                'ubg', obj.nlp_ubg, ...
                'lbx', obj.nlp_lbx, ...
                'ubx', obj.nlp_ubx, ...
                'p', obj.nlp_p, ...
                %                 'lam_x0', obj.nlp_lam_x0, ...
                %                 'lam_g0', obj.nlp_lam_g0
                };
            
            obj.sol = obj.solver(args{:});
            if obj.solver.stats.success ~= true
                solve_status_str = obj.solver.stats.return_status;
                fprintf([' [' class(obj) ': ' solve_status_str '] ']);
                obj.sol.x(obj.idx.u0) = nan;
            end
            
            % Use the current solution to speed up the next optimization
            obj.nlp_x0 = obj.sol.x;
            obj.nlp_lam_x0 = obj.sol.lam_x;
            obj.nlp_lam_g0 = obj.sol.lam_g;
        end
        function T_opt = get_T_opt(obj)
            T_opt = obj.T_opt;
        end
        function X_opt = get_X_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.X;
            X_opt = full(reshape(nlp_x(id(1):id(2)), obj.nx, obj.N));
        end
        function U_opt = get_U_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.U;
            U_opt = full(reshape(nlp_x(id(1):id(2)), obj.nu, obj.N - 1));
        end
    end
end