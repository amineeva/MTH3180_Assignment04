%The test

%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant

% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
clc
clf
clear
%% HUEN'S METHOD!!!
BT_struct_H = struct();
BT_struct_H.A = [0 0; 1 0];
BT_struct_H.B = [.5; .5];
BT_struct_H.C = [0; 1];

%% VAN DER HOUWEN'S/WRAY'S METHOD!!!
BT_struct_V = struct();
BT_struct_V.A = [0 0 0; 8/15 0 0; 1/4 5/12 0];
BT_struct_V.B = [1/4; 0; 3/4];
BT_struct_V.C = [0; 8/15; 2/3];

%% KUTTA'S METHOD!!!
BT_struct_K = struct();
BT_struct_K.A = [0 0 0; .5 0 0; -1 2 0];
BT_struct_K.B = [1/6; 2/3; 1/6];
BT_struct_K.C = [0; 1/2; 1];

%% Fehlberg's Embedded Method
BT_struct_F = struct();
BT_struct_F.A = [0,0,0,0,0,0;...
1/4, 0,0,0,0,0;...
3/32, 9/32, 0,0,0,0;...
1932/2197, -7200/2197, 7296/2197, 0,0,0;...
439/216, -8, 3680/513, -845/4104, 0,0;...
-8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
BT_struct_F.B =  [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
25/216, 0, 1408/2565, 2197/4104, -1/5, 0];

BT_struct_F.C = [0, 1/4, 3/8, 12/13, 1, 1/2];

%P, Orders of the method
p=length(BT_struct_F.C);
%εd, desired error
error_desired=10^-10;

orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_stepcount = 800;
    t_max = 60;
    t_range = linspace(0,t_max,t_stepcount);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');
    hold off

    figure;
subplot(2, 1, 1);
hold on
    plot(t_range, V_list(:,1), 'k-',  'DisplayName', "Exact Solution", LineWidth=2 )
    xlabel('Time (-)', Interpreter='latex', FontSize=16)
    ylabel('X Position (-)', Interpreter='latex', FontSize=16)
    title('X Position of Planet vs Time', Interpreter='latex', FontSize=20)
    ylim([-3,9])

 subplot(2, 1, 2)
hold on
    plot(t_range,V_list(:,2), 'k-', LineWidth=2)
    title('Y Position of Planet vs Time', Interpreter='latex', FontSize=20)
    xlabel('Time (-)', Interpreter='latex',FontSize=16)
    ylabel('Y Position (-)', Interpreter='latex', FontSize=16)

tspan = [0, t_max];
% %HUEN'S METHOD!
% 
% subplot(2, 1, 1);
% hold on
% xlim(tspan);
% h_ref=(t_max+1)/t_stepcount;
% 
% 
%     x0 = 8;
%     y0 = 0;
%     dxdt0 = 0;
%     dydt0 = 1.5;
% 
%     V0 = [x0;y0;dxdt0;dydt0];
% 
%    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_orbital(@gravity_rate_func,tspan,V0,h_ref,BT_struct_H, orbit_params);
% 
%     plot(t_list, X_list(:,1), 'r--', 'LineWidth', 1.5, 'MarkerSize', 1.5, 'DisplayName', "Huen's Method");
% 
% %VAN DER HOUWEN'S/WRAY'S METHOD!
% 
% subplot(2, 1, 1);
% hold on
% xlim(tspan);
% h_ref=(t_max+1)/t_stepcount;
% 
% 
%     x0 = 8;
%     y0 = 0;
%     dxdt0 = 0;
%     dydt0 = 1.5;
% 
%     V0 = [x0;y0;dxdt0;dydt0];
% 
%    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_orbital(@gravity_rate_func,tspan,V0,h_ref,BT_struct_V, orbit_params);
% 
%     plot(t_list, X_list(:,1), 'g--', 'LineWidth', 1.5, 'MarkerSize', 1.5, 'DisplayName', "Wray's Method");

%% EMBEDDED FEHLBERG'S METHOD!

subplot(2, 1, 1);
hold on
xlim(tspan);
h_ref=(t_max+1)/t_stepcount;


    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];

   [t_list,X_list,h_avg, num_evals, h_list, local_error_list] = embedded_RK_fixed_step_integration_orbital(@gravity_rate_func,tspan,V0,h_ref,BT_struct_V, p,error_desired, orbit_params);

    plot(t_list, X_list(:,1), 'g--', 'LineWidth', 1.5, 'MarkerSize', 1.5, 'DisplayName', "Wray's Method");

%     %KUTTA'S METHOD!
% 
% subplot(2, 1, 1);
% hold on
% xlim(tspan);
% h_ref=(t_max+1)/t_stepcount;
% 
% 
%     x0 = 8;
%     y0 = 0;
%     dxdt0 = 0;
%     dydt0 = 1.5;
% 
%     V0 = [x0;y0;dxdt0;dydt0];
% 
%    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration_orbital(@gravity_rate_func,tspan,V0,h_ref,BT_struct_K, orbit_params);
% 
%     plot(t_list, X_list(:,1), 'b--', 'LineWidth', 1.5, 'MarkerSize', 1.5, 'DisplayName',"Kutta's Method");
%     legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
