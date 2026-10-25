%%% Day 1 RK Step, file 1

%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals, redo] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct, e_d,orbit_params)
    %your code here

    % initialize m x s  matrix of zeros for K
    m = length(XA); % height of X (height)
    s = length(BT_struct.B); % number of stages for this method (width)
    K = zeros(m, s);

    % populate K matrix
    for i = 1:s
        % sum of dot product of vector of k's
        sum_val1 = K*(BT_struct.A(i,:)');

        % current stage
        K(:,i) = rate_func_in(t + BT_struct.C(i)*h, XA + h*sum_val1, orbit_params);

    end

    % doc product of K and B
    sum_val2_1 = K*BT_struct.B(1);
    sum_val2_2 = K*BT_struct.B(2);
    % finding XB:
    XB1 = XA + h*sum_val2_1;
    XB2 = XA + h*sum_val2_2;

    if abs(XB1-XB2)>e_d
        redo=1;
    else
        redo=0;
    end
    

    % how can you figure out what num_evals should be?
    num_evals = s;
end