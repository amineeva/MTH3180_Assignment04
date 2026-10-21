%%% Day 1 RK Global truncation plotting


%%% Global truncation error: Runge-Kutta %%%
clear;
clc;
close all;

%% Heun's Method
BT_struct_H = struct(...
    'A', [0, 0; 1, 0],...
    'B', [1/2; 1/2],...
    'C', [0; 1]);

%% Ralston's Method
BT_struct_R = struct(...
    'A', [0, 0; 2/3, 0],...
    'B', [1/4; 3/4],...
    'C', [0; 2/3]);

%% Kutta's Third Order Method
BT_struct_K = struct(...
    'A', [0, 0, 0; 1/2, 0, 0; -1, 2, 0],...
    'B', [1/6; 2/3; 1/6],...
    'C', [0; 1/2; 1]);


%% setup
% time interval
t0 = 0;
tf = 1;

% getting the exact solution for comparison
X_exact = @solution01;
X0 = solution01(t0);

% step sizes 
% varying N instead of h_ref case bc # of steps impacts step sizes
N_refs = unique(round(logspace(1, 4, 50))); % number of time steps
h_refs = (tf - t0) ./ N_refs; % step sizes 

num_h = length(h_refs);

%% setup arrays to store vals
gte_heun   = zeros(size(h_refs));
gte_ralston = zeros(size(h_refs));
gte_kutta  = zeros(size(h_refs));


%% Calculate global truncation error

for i = 1:num_h

    h = h_refs(i);
    N = N_refs(i);

    % initial conditions for each method
    XH = X0;
    XR = X0;
    XK = X0;

    for n = 1:N

        t = t0 + (n-1)*h;

        %%% solution 01
        XH = explicit_RK_step(@rate_func01, t, XH, h, BT_struct_H);
        XR = explicit_RK_step(@rate_func01, t, XR, h, BT_struct_R);
        XK = explicit_RK_step(@rate_func01, t, XK, h, BT_struct_K);

    end

    % Global error at the final time
    Xf_exact = X_exact(tf);

    gte_heun(i) = norm(XH - Xf_exact);
    gte_ralston(i) = norm(XR - Xf_exact);
    gte_kutta(i) = norm(XK - Xf_exact);

end

%% fit lines

error_floor = 1e-14;
h_max_fit = 1e-2;
valid_h = h_refs <= h_max_fit;

valid_heun = valid_h & (gte_heun > error_floor);
valid_ralston = valid_h & (gte_ralston > error_floor);
valid_kutta = valid_h & (gte_kutta > error_floor);

fit_heun = polyfit(log10(h_refs(valid_heun)), log10(gte_heun(valid_heun)), 1);
fit_ralston = polyfit(log10(h_refs(valid_ralston)), log10(gte_ralston(valid_ralston)), 1);
fit_kutta = polyfit(log10(h_refs(valid_kutta)), log10(gte_kutta(valid_kutta)), 1);

p_heun = fit_heun(1);
p_ralston = fit_ralston(1);
p_kutta = fit_kutta(1);

h_fit_heun = logspace(log10(min(h_refs(valid_heun))), log10(max(h_refs(valid_heun))), 200);
y_fit_heun = 10.^polyval(fit_heun, log10(h_fit_heun));

h_fit_ralston = logspace(log10(min(h_refs(valid_ralston))), log10(max(h_refs(valid_ralston))), 200);
y_fit_ralston = 10.^polyval(fit_ralston, log10(h_fit_ralston));

h_fit_kutta = logspace(log10(min(h_refs(valid_kutta))), log10(max(h_refs(valid_kutta))), 200);
y_fit_kutta = 10.^polyval(fit_kutta, log10(h_fit_kutta));


%% Plot global truncation error vs. step size

figure;
hold on;

% measured data points
loglog(h_refs, gte_heun, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'LineStyle', 'none', 'DisplayName', 'Heun Measured')

loglog(h_refs, gte_ralston, 'ro','MarkerFaceColor', 'r', 'MarkerSize', 4, 'LineStyle', 'none', 'DisplayName', 'Ralston Measured');
loglog(h_refs, gte_kutta, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineStyle', 'none', 'DisplayName', 'Kutta Measured');

% fitted lines
loglog(h_fit_heun, y_fit_heun, 'b-', 'LineWidth', 1.8, 'DisplayName', sprintf('Heun fit: p = %.2f', p_heun));
loglog(h_fit_ralston, y_fit_ralston, 'r-','LineWidth', 1.8, 'DisplayName', sprintf('Ralston fit: p = %.2f', p_ralston));
loglog(h_fit_kutta, y_fit_kutta, 'k-', 'LineWidth', 1.8,'DisplayName', sprintf('Kutta fit: p = %.2f', p_kutta));

xlabel('Step size ($h$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Global truncation error ($\epsilon$)','Interpreter', 'latex', 'FontSize', 18);
title('Global Truncation Error vs. Step Size', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 15, 'TickLabelInterpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
box on;
hold off;