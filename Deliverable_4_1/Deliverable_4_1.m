%% Problem definition
addpath(fullfile('..', 'src')); % Add the path with all the informations about the rocket 
import gurobi.*;
clc; clear; close all;

% General settings
H           = 8; % Horizon length   [s]
Ts          = 1/20; % Sample time   [s]
Tf          = 30; % Simulation time [s]
simAnimationTime = 20; %            [s]
closedLoop  = false;

% Initialize the system
x0          = zeros(12,1); %vector intial condition

rocket      = Rocket(Ts); % Create the rocket
[xs,us]     = rocket.trim(); % Find one equilibrium point
sys         = rocket.linearize(xs, us); % Linearize the system in the equilibrium point

[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); % Decompose the system in four indipendent systems

% Setup the controllers
mpc_x       = MpcControl_x(sys_x, Ts, H); % Controller for x state
mpc_y       = MpcControl_y(sys_y, Ts, H); % Contcroller for y state
mpc_z       = MpcControl_z(sys_z, Ts, H); % Controller for z state
mpc_roll    = MpcControl_roll(sys_roll, Ts, H); % Controller for roll state

lmpc        = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Open Loop
ref         = [2 2 2 deg2rad(40)]'; %when the ref is a point
[u, T_opt, X_opt, U_opt] = lmpc.get_u(x0, ref);
U_opt(:,end+1) = nan;
RockPlot    = rocket.plotvis(T_opt, X_opt, U_opt, ref); % Plot as usual

% Closed Loop
% ref = @(t_, x_) ref_TVC(t_); %ref is a trajectory
% [T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
% rocket.anim_rate = simAnimationTime; % Increase this to make the animation faster
% ph = rocket.plotvis(T, X, U, Ref);
% ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title
% 
% if closedLoop == true
%     ref = @(t_, x_) ref_TVC(t_); %ref is a trajectory
%     [T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
%     rocket.anim_rate = simAnimationTime; % Increase this to make the animation faster
%     ph = rocket.plotvis(T, X, U, Ref);
%     ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title
% else
%     ref = [2 2 2 deg2rad(40)]'; %when the ref is a point
%     [u, T_opt, X_opt, U_opt] = lmpc.get_u(x0, ref);
%     U_opt(:,end+1) = nan;
%     ph = rocket.plotvis(T_opt, X_opt, U_opt, ref); % Plot as usual
% end