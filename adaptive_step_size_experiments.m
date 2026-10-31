% Day 3 Plotting

% using the orbital ODE %

%% setting up
% Orbit params
orbit_params.G = 4*pi^2;
orbit_params.m_sun = 1; % holy sun mass
orbit_params.m_planet = 3.003e-6; % holy earth mass

% initial conditions
X0 = [1; 0; 0; 2*pi]; % circular orbit
X0 = [1; 0; 0; 4.5]; % eccentric orbit
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


%% running integration
% need a reference solution
X_ref_final = compute_planetary_motion(tspan(2), X0, orbit_params);


% vary error tolerance for adaptive integrator
error_desired_list = logspace(-4, -10, 7);
num_tol = length(error_desired_list);
gte_adaptive = zeros(1,num_tol);
h_avg_adaptive = zeros(1,num_tol);
evals_adaptive = zeros(1,num_tol);
failure_rate = zeros(1,num_tol);

for i = 1:num_tol
    error_desired = error_desired_list(i);
    [t_list,X_list,h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration ...
    (rate_func_in,tspan,X0,h_ref,DormandPrince,p,error_desired);

    % find GTE at final timestep
    gte_adaptive(i) = norm(X_list(end,:)' - X_ref_final);

    % store vals
    h_avg_adaptive(i) = h_avg;
    evals_adaptive(i) = num_evals;
    failure_rate(i) = step_failure_rate;
end


% F=vary step size for fixed-step integrator
DormanPrince_fixed = DormandPrince;
DormanPrince_fixed.B = DormandPrince.B(1,:)'; %getting the 5th order RK weights used by adaptive method
h_fixed_list = logspace(log10(1/10), log10(1/1000), 8);
num_h = length(h_fixed_list);
gte_fixed = zeros(1,num_h);
h_avg_fixed = zeros(1,num_h);
evals_fixed = zeros(1,num_h);

for i = 1:num_h
    h_ref = h_fixed_list(i);
    [t_fixed,X_fixed,h_avg, num_evals] = explicit_RK_fixed_step_integration_orbital ...
    (@gravity_rate_func,tspan,X0,h_ref,DormanPrince_fixed, orbit_params);

    % find GTE at final timestep
    gte_fixed(i) = norm(X_fixed(end,:)' - X_ref_final);

    % store vals
    h_avg_fixed(i) = h_avg;
    evals_fixed(i) = num_evals;

end
%% plotting

% GTE vs. average step size
% getting fit lines
% Fit lines in log-log space
fit_fixed = polyfit(log10(h_avg_fixed), log10(gte_fixed), 1);
fit_adaptive = polyfit(log10(h_avg_adaptive), log10(gte_adaptive), 1);
h_fit_fixed = logspace(log10(min(h_avg_fixed)), log10(max(h_avg_fixed)), 200);
h_fit_adaptive = logspace(log10(min(h_avg_adaptive)), log10(max(h_avg_adaptive)), 200);
gte_fit_fixed = 10.^(polyval(fit_fixed, log10(h_fit_fixed)));
gte_fit_adaptive = 10.^(polyval(fit_adaptive, log10(h_fit_adaptive)));

figure();
hold on;
loglog(h_avg_fixed,gte_fixed,'bo', 'LineWidth',1.5,'MarkerFaceColor','b', 'DisplayName','Fixed step');
loglog(h_avg_adaptive,gte_adaptive,'ro', 'LineWidth',1.5,'MarkerFaceColor','r', 'DisplayName','Adaptive step');
loglog(h_fit_fixed, gte_fit_fixed, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fixed step fit');
loglog(h_fit_adaptive, gte_fit_adaptive, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Adaptive step fit');

xlabel('Average Step Size ($h$) (years)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global truncation error ($\epsilon$)','Interpreter', 'latex', 'FontSize', 18);
title('Global Truncation Error vs. Average Step Size', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 15, 'TickLabelInterpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);


grid on;
box on;
hold off;
% exportgraphics(gcf, 'gte_avg_step_size_RK_comparison.png', 'Resolution', 600);



% GTE vs. num of rate function calls during integration
% getting fit lines
% Fit lines in log-log space
fit_fixed = polyfit(log10(evals_fixed), log10(gte_fixed), 1);
fit_adaptive = polyfit(log10(evals_adaptive), log10(gte_adaptive), 1);
evals_fit_fixed = logspace(log10(min(evals_fixed)), log10(max(evals_fixed)), 200);
evals_fit_adaptive = logspace(log10(min(evals_adaptive)), log10(max(evals_adaptive)), 200);
gte_fit_fixed = 10.^(polyval(fit_fixed, log10(evals_fit_fixed)));
gte_fit_adaptive = 10.^(polyval(fit_adaptive, log10(evals_fit_adaptive)));

figure();
hold on;
loglog(evals_fixed,gte_fixed,'bo', 'LineWidth',1.5,'MarkerFaceColor','b', 'DisplayName','Fixed step');
loglog(evals_adaptive,gte_adaptive,'ro', 'LineWidth',1.5,'MarkerFaceColor','r', 'DisplayName','Adaptive step');
loglog(evals_fit_fixed, gte_fit_fixed, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fixed step fit');
loglog(evals_fit_adaptive, gte_fit_adaptive, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Adaptive step fit');

xlabel('Number of Function Evaluations (-)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global truncation error ($\epsilon$)','Interpreter', 'latex', 'FontSize', 18);
title('Global Truncation Error vs. Function Evaluations', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 15, 'TickLabelInterpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);


grid on;
box on;
hold off;
% exportgraphics(gcf, 'gte_num_function_evals_size_RK_comparison.png', 'Resolution', 600);

% failure rate vs. average step size
figure();
hold on;
loglog(h_avg_adaptive,failure_rate,'ro-', 'LineWidth',1.5,'MarkerFaceColor','r', 'DisplayName','Adaptive step');

xlabel('Average Step Size ($h$) (years)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Step Failure Rate (-)','Interpreter', 'latex', 'FontSize', 18);
title('Failure Rate vs. Average Step Size', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 15, 'TickLabelInterpreter', 'latex');
xlim([1e-3 1e-1]);
ylim([1e-3 1e0]);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);


grid on;
box on;
hold off;
% exportgraphics(gcf, 'failure_rate_avg_step_size_RK_comparison.png', 'Resolution', 600);
