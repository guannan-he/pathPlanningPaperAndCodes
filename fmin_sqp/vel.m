function [min_val,max_val] = vel(x)
syms u;
tm = x(6);
% u_v_y = [3 * u^2,4 * u^3,5 * u^4] / tm;
u_v_x = [1,2 * u,3 * u^2,4 * u^3,5 * u^4] / tm;
% yu_a = [10 -15 6]';
xu_a = [1 0 0 0 0 0
        0 0 0 0 0 0
        -6 -4 0 0 10 0
        8 7 0 0 -15 0
        -3 -3 0 0 6 0];
% YB = 3.75;
% yud = u_v_y * yu_a * YB;
xud = u_v_x * xu_a * x;
% spd = sqrt(xud^2 + yud^2);
spd = xud ;
y = [0:0.01:1];
spd_fun = matlabFunction(spd);
spd_val = spd_fun(y);
max_val = max(spd_val);
min_val = min(spd_val);
end