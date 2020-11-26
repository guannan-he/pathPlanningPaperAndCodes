function [min_k,max_k] = get_k(x)
syms u;
c1 = x(1);%d1
c2 = x(2);%d4
c3 = x(3);%x2
ki = x(7);
xt = x(8);
yt = x(9);
pt = x(10);
p0x = 0;                    p0y = 0;
p1x = c1;                   p1y = 0;
p2x = c3;                   p2y =4 * ki * c1^2 / 3;
p3x = xt - c2 * cos(pt);    p3y = yt - c2 * sin(pt);
p4x = xt;                   p4y = yt;
xu = p0x * (1 - u)^4 ...
    + 4 * p1x * (1 - u)^3 * u ...
    + 6 * p2x * (1 - u)^2 * u^2 ...
    + 4 * p3x * (1 - u) * u^3 ...
    + p4x * u^4;
yu = p0y * (1 - u)^4 ...
    + 4 * p1y * (1 - u)^3 * u ...
    + 6 * p2y * (1 - u)^2 * u^2 ...
    + 4 * p3y * (1 - u) * u^3 ...
    + p4y * u^4;
xud = diff(xu,u);
xudd = diff(xud,u);
yud = diff(yu,u);
yudd = diff(yud,u);
ku = (xud * yudd - yud * xudd) / (xud^2 + yud^2)^(1.5);
% kud = diff(ku,u);
func_ku = matlabFunction(ku);
u_ran = 0:0.005:1;
u_res = func_ku(u_ran);
min_k = min(u_res);
max_k = max(u_res);
% func_kud = matlabFunction(kud);
end