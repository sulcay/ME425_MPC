addpath(fullfile('..', 'src'));

close all
clear
clc

Tf          = 2.5;

x0          = zeros(12, 1);

ref         = [0.5, 0, 1, deg2rad(65)]'; % Constant reference

H           = 1.0; % Horizon length in seconds

% Create rocket
Ts          = 1/40;
rocket      = Rocket(Ts);
rocket.anim_rate = 1;
rocket.mass = 1.75;

%% Without compensation

d_comp      = 0;
d_actual    = 2;


% Create MPC controller
nmpc        = NmpcControl(rocket, H, d_comp);

% Add delay
rocket.delay = d_actual; % at 2 drop, 3 fail

% Simulate and Plot
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);


%% With partial compensation

d_comp      = 3;
d_actual    = 4;


% Create MPC controller
nmpc        = NmpcControl(rocket, H, d_comp);

% Add delay
rocket.delay = d_actual;

% Simulate and Plot
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);


%% With full compensation

d_comp      = 4;
d_actual    = 4;


% Create MPC controller
nmpc        = NmpcControl(rocket, H, d_comp);

% Add delay
rocket.delay = d_actual;

% Simulate and Plot
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph          = rocket.plotvis(T, X, U, Ref);