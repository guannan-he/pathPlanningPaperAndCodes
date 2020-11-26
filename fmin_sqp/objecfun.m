function f = objecfun(x)
syms u;
u_v = [1 u u^2 u^3 u^4 u^5];
xu_a = [0 0 0 0 0 0
        1 0 0 0 0 0
        0 0 0 0 0 0
        -6 -4 0 0 10 0
        8 7 0 0 -15 0
        -3 -3 0 0 6 0];
yu_a = [0 0 0 10 -15 6]';
YB = 3.75;
%% x
xu = u_v * xu_a * x;
xud = diff(xu,u);
xudd = diff(xud,u);
%% y
yu = u_v * yu_a * YB;
yud = diff(yu,u);
yudd = diff(yud,u);
%% 计算代价函数
ds = sqrt(xud^2 + yud^2);
dk = abs((xud * yudd - xudd * yud) / (xud^2 + yud^2)^(3/2));% * ds;
dk_f = matlabFunction(dk);
ds_f = matlabFunction(ds);
s = integral(ds_f,0,1);
k = integral(dk_f,0,1);% / s;
%% 输出
w_1 = 1000;
w_2 = 1;
f = w_1 * k +w_2 * s;
% toc;
end