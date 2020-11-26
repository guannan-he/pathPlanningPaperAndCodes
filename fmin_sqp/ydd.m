function [min_val,max_val] = ydd(x)
syms u;
tm = x(6);
u_v = [6 * u,12 * u^2,20 * u^3] / tm^2;
yu_a = [10 -15 6]';
YB = 3.75;
%% y
yudd = u_v * yu_a * YB;
ydd_fun = matlabFunction(yudd);
y = [0:0.01:1];
yudd_val = ydd_fun(y);
min_val = min(yudd_val);
max_val = max(yudd_val);
end