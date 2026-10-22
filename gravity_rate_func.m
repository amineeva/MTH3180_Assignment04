%%% Day 1 Testing Runge-Kutta


%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
x_p=V(1);
y_p=V(2);
dxdt_p=V(3);
dydt_p=V(4);
d2xdt2_p = x_p*(-orbit_params.m_sun*orbit_params.G/((x_p^2+y_p^2)^1.5));
d2ydt2_p = y_p*(-orbit_params.m_sun*orbit_params.G/((x_p^2+y_p^2)^1.5));
dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p];
end