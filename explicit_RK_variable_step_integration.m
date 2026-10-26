%%% Day 2 - variable step RK integration
%%% Day 3 addition -> calculating step failure rate

%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
%step_failure_rate: # of failed steps / # of attempted steps
function [t_list,X_list,h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    %your code here

    % get the start and end times
    t_start = tspan(1);
    t_end = tspan(2);

    % setup
    t_n = t_start;
    X_n = X0(:);
    h = h_ref;
    t_list = t_start;
    X_list = X_n(:)';
    num_evals = 0;
    total_step_size = 0;
    accepted_steps = 0;
    attempted_steps = 0;
    total_failed_steps = 0;

    % run loop
    while t_n < t_end
        % orion's safeguard -> make sure the h size never goes beyond t_end
        h = min(h, t_end - t_n);

        % variable RK step
        [XB, step_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in,t,XA,h,BT_struct,p,error_desired);
        
        % add to counter
        num_evals = num_evals + step_evals;

        if redo % need to retry with a different step size
            h = h_next;
            total_failed_steps = total_failed_steps + 1;
        else % step is fine, updating happens
            t_n = t_n + h;
            X_n = XB(:);
            t_list(end+1, 1) = t_n;
            X_list(end+1, :) = X_n(:)';
            total_step_size = total_step_size + h;
            accepted_steps = accepted_steps + 1;
            
            h = h_next;

        end
        attempted_steps = attempted_steps + 1;
    end

    % calculate avg accepted step size
    h_avg = total_step_size / accepted_steps;

    % calculate step failure rate
    step_failure_rate = total_failed_steps / attempted_steps;
end
