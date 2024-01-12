classdef MpcControl_roll < MpcControlBase
    properties
        % Define the cost parameters
        Q = diag([350,1500]);
        R = 0.01;
        
        % Define the constraints
        F = nan;
        f = nan;
        M = [1;-1];
        m = [20;20];
        Xmax = [Inf,Inf];
        Xmin = [-Inf,-Inf];
        Umax = 20;
        Umin = -20;
    end
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing
            
            [nx, nu] = size(mpc.B);
            
            % Steady-state targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Parameters definition
            useSlacks = 0;
            stateConstraints = ~all(all(isnan(mpc.F))) && ~all(isnan(mpc.f));

            % Define the slack variables
            useSlacks = stateConstraints && useSlacks;
            if useSlacks; eps = sdpvar(size(mpc.f,1), N); end

            % Define constraints and objectives
            con = [];
            obj = 0;


            % Set the constraints for t=1...N-1
            for i=1:N-1
                con = [con, (X(:,i+1)-x_ref) == mpc.A * (X(:,i)-x_ref) + mpc.B * (U(:,i)-u_ref)];
                if stateConstraints
                    if useSlacks
                        con = [con,mpc.F*(X(:,i)-x_ref) <= mpc.f-mpc.F*x_ref-eps(:,i)]; 
                    else
                        con = [con,mpc.F*(X(:,i)-x_ref) <= mpc.f-mpc.F*x_ref];
                    end
                end
                con = [con,mpc.M*(U(:,i)-u_ref) <= mpc.m-mpc.M*u_ref];
                if useSlacks
                    obj = obj + (X(:,i)-x_ref)'*mpc.Q*(X(:,i)-x_ref) + (U(:,i)-u_ref)'*mpc.R*(U(:,i)-u_ref) + eps(:,i)'*mpc.S*eps(:,i);
                else
                    obj = obj + (X(:,i)-x_ref)'*mpc.Q*(X(:,i)-x_ref) + (U(:,i)-u_ref)'*mpc.R*(U(:,i)-u_ref);
                end
            end

            % Set the constraints for t=N
            if stateConstraints
                if useSlacks
                    con = [con,mpc.F*(X(:,N)-x_ref) <= mpc.f-mpc.F*x_ref-eps(:,N)]; 
                else
                    con = [con,mpc.F*(X(:,N)-x_ref) <= mpc.f-mpc.F*x_ref];
                end
            end
            if useSlacks
                obj = obj + (X(:,N)-x_ref)'*mpc.Q*(X(:,N)-x_ref) + eps(:,N)'*mpc.S*eps(:,N);
            else
                obj = obj + (X(:,N)-x_ref)'*mpc.Q*(X(:,N)-x_ref);
            end
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, {U(:,1), X, U});
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state targets
            nx = size(mpc.A, 1);
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            % Initialize parameters
            Rs = 1;
            [ny, nx] = size(mpc.C);

            % Define constraints and cost function
            obj = us'*Rs*us;
            con = [eye(nx,nx) - mpc.A, -mpc.B ; mpc.C, 0]*[xs; us] == [zeros(nx,ny);ref];
            if ~all(all(isnan(mpc.F))) && ~all(all(isnan(mpc.f))); con = [con, mpc.F*xs <= mpc.f]; end
            if ~all(isnan(mpc.M)) && ~all(isnan(mpc.m)); con = [con, mpc.M*us <= mpc.m]; end
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
