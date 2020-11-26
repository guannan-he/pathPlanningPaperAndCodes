clear;
clc;
%% 始末状态设置
KI = 0.1;%k start
xt = 15;
yt = -2;
phit = 0.05;
%% 生成 等式约束
XI = [0,0,0,KI];%start
XT = [xt,yt,phit];%final
c = [zeros(7,3),eye(7)];%eq matrix
b = [XI,XT];%eq val
%% 生成线性不等式约束
c_ieq_l = [-1 * eye(2),zeros(2,8)];
b_ieq_l = [0;0];
%% 求解
x_start = [0.5,0.5,0.5 * xt,XI,XT];
options=optimoptions('fmincon','Display','iter','Algorithm','sqp');
out = fmincon(@cost_func,x_start,c_ieq_l,b_ieq_l,c,b,[],[],@ieq_constr,options);
%% 绘制
c1 = out(1);%d1
c2 = out(2);%d4
c3 = out(3);%x2
syms u;
p0x = 0;                      p0y = 0;
p1x = c1;                     p1y = 0;
p2x = c3;                     p2y = 4 * KI * c1^2 / 3;
p3x = xt - c2 * cos(phit);    p3y = yt - c2 * sin(phit);
p4x = xt;                     p4y = yt;
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
xu_func = matlabFunction(xu);
yu_func = matlabFunction(yu);
u_ran = 0:0.005:1;
x_res = xu_func(u_ran);
y_res = yu_func(u_ran);
figure(1);
clf;
plot(x_res,y_res);
title('path_gen');
xlabel('x/m');
ylabel('y/m');