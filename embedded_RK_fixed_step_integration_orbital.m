%%% Day 1 RK implementation, file 2

%Runs numerical integration arbitrary RK method
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
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, h_list, local_error_list] = embedded_RK_fixed_step_integration_orbital ...
(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired, orbit_params)
    % get the start and end times
    t_start = tspan(1);
    t_end = tspan(2);

    % find N
    N = ceil((t_end - t_start)/h_ref);

    % set vectors
    t_vector = [t_start]'; % needs to be a column vector
    X_vector = zeros(N+1, length(X0));

    % initial conditions
    t_vector(1) = t_start;
    X_vector(1, :) = X0(:)'; % rows are time steps, columns are components of X

    % num evals
    num_evals = 0;
    h_avg = (t_end - t_start)/N; % same as h_ref?
   
    alpha=4;
    
h_list=[h_avg];
local_error_list=[0];

t_n=t_start;
    % loop for N values
    while t_n<t_end
        %List of step size history
        
        h_list(end+1)=h_avg;
        t_n = t_vector(end);
        X_n = X_vector(end,:)'; % X_n is a column vector
        % run one forward euler step     
        [XB1, XB2, step_evals,redo] = RK_step_embedded(rate_func_in,t_n,X_n,h_avg,BT_struct,error_desired, orbit_params);
        % update XA and time
        if redo==0
        t_vector(end+1) = t_n + h_avg;
        X_vector(end+1,:) = XB2(:)'; % XB is a column vector 
        h_avg = min(.9*(error_desired/abs(XB1-XB2))^(1/p),alpha);
        %Condition just in case algorithm creates step size exceeding t_end
        if t_n+h_avg>t_end
            t_vector(end+1) = t_end;
            h_avg=t_end-t_n;
        end
        local_error_list(end+1)=abs(XB1-XB2);
        % update num_evals
        num_evals = num_evals + step_evals;
       
        else
           h_avg = min(.9*(error_desired/abs(XB1-XB2))^(1/p),alpha);
        end

    end
    % outputs
    t_list = t_vector;
    X_list = X_vector;
end