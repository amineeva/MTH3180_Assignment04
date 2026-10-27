%%% Day 2 testing on planetary system


% Orbit params
orbit_params.G = 4*pi^2;
orbit_params.m_sun = 1; % holy sun mass
orbit_params.m_planet = 3.003e-6; % holy earth mass

% initial conditions -> circular orbit with r = 1AU
X0 = [1; 0; 0; 2*pi];
tspan = [0, 1]; %one year
h_ref = 1/12; % initial step size in years
p = 5; %error-estimate scaling component
error_desired = 1e-8; % error tolerance


% Chosen rate function + solution function pair
DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
1/5, 0, 0, 0,0,0,0;...
3/40, 9/40, 0, 0, 0, 0,0;...
44/45, -56/15, 32/9, 0, 0, 0,0;...
19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];



% wrapping gravity_rate_func to add orbit_params
rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);


[t_list,X_list,h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,DormandPrince,p,error_desired);

