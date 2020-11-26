function [min_val,max_val] = phi_d(x)
syms u;
tm = x(6);
u_v_y = [3 * u^2,4 * u^3,5 * u^4] / tm;
u_v_x = [1,2 * u,3 * u^2,4 * u^3,5 * u^4] / tm;
yu_a = [10 -15 6]';
xu_a = [1 0 0 0 0 0
        0 0 0 0 0 0
        -6 -4 0 0 10 0
        8 7 0 0 -15 0
        -3 -3 0 0 6 0];
YB = 3.75;
yud = u_v_y * yu_a * YB;
xud = u_v_x * xu_a * x;
% xudd = diff(xud,u) / tm;
% yudd = diff(yud,u) / tm;
phi = atan(yud / xud);
phi_d = diff(phi,u) / tm;
% phi_d = (xud * yudd - xudd * yud) / (xud^2 + yud^2);
phi_fun = matlabFunction(phi_d);
y = [0:0.01:1];
phi_val = phi_fun(y);
min_val = min(phi_val);
max_val = max(phi_val);
end