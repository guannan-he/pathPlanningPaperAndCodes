function [min_val,max_val] = yd(x)
syms u;
tm = x(6);
u_v = [3 * u^2,4 * u^3,5 * u^4] / tm;
yu_a = [10 -15 6]';
YB = 3.75;
%% y
yud = u_v * yu_a * YB;
yd_fun = matlabFunction(yud);
y = [0:0.01:1];
yud_val = yd_fun(y);
min_val = min(yud_val);
max_val = max(yud_val);
end