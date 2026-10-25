%%% Day 2 - LTE log-log plot

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

% Setup
t_ref = 0.492;
h_refs = logspace(-5, 1, 1000); % 

% Test function 1: X(t) = cos(t)
X_ref = solution01(t_ref); % change this based on what to calculate

% Preallocate error arrays
lte_XB1 = zeros(size(h_refs));
lte_XB2 = zeros(size(h_refs));
estimate_diff = zeros(size(h_refs));
rate_diff = zeros(size(h_refs));

% for tracking evals
num_evals = zeros(size(h_refs));



% Calculate local truncation error for each step size

%%% solution 01
% Calculate local truncation error for each step size
for i = 1:length(h_refs)

    h = h_refs(i);

    % Exact solution after one step
    X_exact_next = solution01(t_ref + h);

    % embedded RK estimates
    [XB1, XB2, num_evals(i)] = RK_step_embedded(@rate_func01,t_ref,X_ref,h,DormandPrince);

    % calculate ltes
    lte_XB1(i) = norm(XB1 - X_exact_next);
    lte_XB2(i) = norm(XB2 - X_exact_next);

    % difference between estimates
    estimate_diff(i) = norm(XB1 - XB2);

    % difference between rate fun and exact solution
    f_start = rate_func01(t_ref, X_ref);
    f_end = rate_func01(t_ref + h, X_exact_next);
    rate_diff(i) = norm(f_end - f_start);

end


% log-log plot of LTE for XB1 and XB2 as a function of h
figure;
hold on;

loglog(h_refs, lte_XB1, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'LineStyle', 'none', 'DisplayName', 'LTE of XB1 (order 5)');
loglog(h_refs, lte_XB2, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'LineStyle', 'none', 'DisplayName', 'LTE of XB2 (order 4)');
loglog(h_refs, estimate_diff, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'g', 'LineStyle', 'none', 'DisplayName', '$|XB1 - XB2|$');
loglog(h_refs, rate_diff, 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.5 0 0.7], 'LineStyle', 'none', 'DisplayName', '$|f(t+h) - f(t)|$');
xlabel('Step size ($h$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Local Truncation Error ($\epsilon$)', 'Interpreter', 'latex', 'FontSize', 18);
title('Embedded RK: Errors vs. Step Size', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14, 'TickLabelInterpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 11);
grid on;
box on;
hold off;
