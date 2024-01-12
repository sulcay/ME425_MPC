addpath(fullfile('..', 'src'));
addpath(fullfile('..', 'Deliverable_4_1')); % For the linear controllers

close all
clear
clc

% Create rocket
Ts          = 1/20;
rocket      = Rocket(Ts);
rocket.anim_rate = 4;
x0          = zeros(12,1);

% Create MPC controller5
H           = 1.0; % Horizon length in seconds
nmpc        = NmpcControl(rocket, H);

%% Open loop visualization
% Reference
ref         = [0.5, 0, 1, deg2rad(65)]';

% Evaluate once and plot optimal open−loop trajectory,
% pad last input to get consistent size with time and state
disp('Open Loop')
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref);
U_opt(:,end+1) = nan;
ph          = rocket.plotvis(T_opt, X_opt, U_opt, ref);

%% Closed loop simulation
Tf          = 30;

% Reference with normal roll
disp('Normal Roll')
roll_max    = deg2rad(15);
ref         = @(t, x) ref_TVC(t, roll_max);
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);


% Reference with large roll
disp('Large Roll')
roll_max    = deg2rad(50);
ref         = @(t, x) ref_TVC(t, roll_max);
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);


%% Linear Controller to Compare

disp('Linear Controller')

% Initialize the system
[xs, us]    = rocket.trim(); % Find one equilibrium point
sys         = rocket.linearize(xs, us); % Linearize the system in the equilibrium point


[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); % Decompose the system in four indipendent systems

% Designed the controllers for each sub-system
H           =   5.0; % Horizon length     [s]

% Apply the controllers designed for each sub-system
mpc_x       =   MpcControl_x(sys_x, Ts, H);
mpc_y       =   MpcControl_y(sys_y, Ts, H);
mpc_z       =   MpcControl_z(sys_z, Ts, H);
mpc_roll    =   MpcControl_roll(sys_roll, Ts, H);

% Merge the four sub−system controllers into one full−system controller
lmpc        = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Simulate
% Reference with large roll
Tf          = 30;
roll_max    = deg2rad(50);
ref         = @(t, x) ref_TVC(t, roll_max);
[T, X, U, Ref] = rocket.simulate(x0, Tf, @lmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);

