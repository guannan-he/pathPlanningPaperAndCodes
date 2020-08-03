clear;
clc;
load 'track_352.mat';%原始赛道信息
load 'body.mat';%车辆信息
%% 平滑插值赛道
tolerance = 5;%先线性插值，去掉大间隙，步长不要大于下一步的二倍（经验）
[track_x,track_y] = track_smooth_linear(track_x,track_y,tolerance);
tolerance = 3;%然后三次样条线插值，达到最终步长
[track_x,track_y] = track_smooth_spline(track_x,track_y,tolerance);
%% 求得与赛道长度对应的样条线参数
[c,~] = size(track_x);
[spline_A] = spline_matrix_gen(c);%样条线求解矩阵，要留在内存中，加速计算
para_x = spline_A * track_x;
para_y = spline_A * track_y;

%% 单位切向量
[vector] = vector_gen(para_x,para_y,c);

%% 绘制中心线、左右边界及障碍物
draw_track(track_x,track_y,vector,track_w,obs_x,obs_y,c);
clear x_l y_l x_r y_r;

%% 第一次QP优化
alpha = zeros(c,1);
[H,F] = HF_gen(spline_A,para_x,para_y,track_x,track_y,vector,c);
[lb,ub] = update_lub(track_w,body.Wb,alpha,c);
[Aeq,beq] = eq_gen(c);
[Aieq,bieq] = ieq_gen(track_x,track_y,obs_x,obs_y,track_w,body.Wb,obs_w,alpha,c);
alpha_last = alpha;
alpha = quadprog(H,F,Aieq,bieq,Aeq,beq,lb,ub);
draw_optimised_patn(track_x,track_y,vector,alpha,'m');
%% 第二次QP优化
alpha = alpha_last * 1 / 3 + alpha * 2 / 3;%更新alpha
[track_x,track_y] = path_gen(track_x,track_y,vector,alpha);%更新中心线
para_x = spline_A * track_x;
para_y = spline_A * track_y;
[vector] = vector_gen(para_x,para_y,c);%更新向量
[lb,ub] = update_lub(track_w,body.Wb,alpha,c);
[Aieq,bieq] = ieq_gen(track_x,track_y,obs_x,obs_y,track_w,body.Wb,obs_w,alpha,c);
alpha_last = alpha;
alpha = quadprog(H,F,Aieq,bieq,Aeq,beq,lb,ub);
draw_optimised_patn(track_x,track_y,vector,alpha,'k');
%% 绘制单位向量
figure(2)
clf;
plot(vector(:,1));
hold on;
plot(vector(:,2));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%子函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 生成样条线求解矩阵
function [spline_A] = spline_matrix_gen(c) %[a;b;c;d] = A * y
%% 输入D
D = 3 * [eye(c - 2),zeros(c - 2,2)] -6 * [zeros(c - 2,1),eye(c - 2),zeros(c - 2,1)] + 3 * [zeros(c - 2,2),eye(c - 2)];
tmp = [-3,3,zeros(1,c - 4),3,-3];
D = [tmp;D;tmp];
%% c
A = eye(c) + 0.5 * [zeros(c - 1,1),eye(c - 1);zeros(1,c)];
A = A + A';
A(1,c - 1) = 0.5;
A(c,2) = 0.5;
matrix_c = 2 * A \ D;
%% a
matrix_a = eye(c);
%% b
Ab = eye(c) / 6 + [zeros(c - 1,1),eye(c - 1) / 3;zeros(1,c)];
Ab(c,2) = 1 / 3;
Aby = [-eye(c - 1),zeros(c - 1,1)] + [zeros(c - 1,1),eye(c - 1)];
Aby = [Aby;Aby(1,:)];
matrix_b = Ab / A * D + Aby;
%% d
Ad = eye(c) / 6 - [zeros(1,c);[eye(c - 1) / 6,zeros(c - 1,1)]];
Ad(1,c - 1) = -1 / 6;
matrix_d = Ad / A * D;
%% 整合
zero = zeros(c);
A_sig = [matrix_a,zero,zero,zero
    zero,matrix_b,zero,zero
    zero,zero,matrix_c,zero
    zero,zero,zero,matrix_d];
Dy = [eye(c);eye(c);eye(c);eye(c)];
spline_A = A_sig * Dy;
end

%% 生成切向量
function [vector] = vector_gen(para_x,para_y,c)
select_b = [zeros(c),eye(c),zeros(c),zeros(c)];
Dx = select_b * para_x;
Dy = select_b * para_y;
vector = zeros(c,2);%该点处切向量
for i = 1:c%切向量
    len = sqrt(Dx(i)^2 + Dy(i)^2);
    vector(i,1) = Dx(i) / len;
    vector(i,2) = Dy(i) / len;
    if isnan(vector(i,1))%滤掉NaN
        vector(i,1) = vector(i - 1,1);
    end
    if isnan(vector(i,2))
        vector(i,2) = vector(i - 1,2);
    end
    if (i>1)&&(vector(i,1) - vector(i - 1,1) > 1.8)%滤掉过大跳变
        vector(i,1) = vector(i - 1,1);
    end
    if (i>1)&&(vector(i,2) - vector(i - 1,2) > 1.8)
        vector(i,2) = vector(i - 1,2);
    end
end
end

%% 绘制原始赛道
function [] = draw_track(track_x,track_y,vector,track_w,obs_x,obs_y,c)
figure(1);
clf;
alpha_r = track_w / 2 *ones(c,1);
alpha_l = -alpha_r;
[x_l,y_l] = path_gen(track_x,track_y,vector,alpha_l);
[x_r,y_r] = path_gen(track_x,track_y,vector,alpha_r);
plot(track_x,track_y,'k--');
hold on;
plot(x_l,y_l,'r');
hold on;
plot(x_r,y_r,'b');
hold on;
plot(track_x(1),track_y(1),'r*');
hold on;
plot(obs_x,obs_y,'r*');
end

%% 生成HF矩阵
function [H,F] = HF_gen(spline_A,para_x,para_y,track_x,track_y,vector,c)
select_b = [zeros(c),eye(c),zeros(c),zeros(c)];
select_c = [zeros(c),zeros(c),eye(c),zeros(c)];
A = select_c * spline_A;
Dx = select_b * para_x;
Dy = select_b * para_y;
base = (Dx.^2 + Dy.^2).^3;
Pxx = diag(Dy.^2 ./ base);
Pxy = diag(-2 * Dx .* Dy ./ base);
Pyy = diag(Dx.^2 ./ base);
Vx = diag(vector(:,2));%alpah左负右正
Vy = diag(-vector(:,1));
Tc = 2 * A;
Tnx = 2 * A * Vx;
Tny = 2 * A * Vy;
H = Tnx' * Pxx * Tnx ...
    + Tny' * Pxy * Tnx + ...
    Tny' * Pyy * Tny;
H = H + H';
F = 2 * Tnx' * Pxx' * Tc * track_x ...
    + Tny' * Pxy' * Tc * track_x + Tnx' * Pxy' * Tc * track_y ...
    + 2 * Tny' * Pyy' * Tc * track_y;
F = F';
end

%% 更新alpha上下限
function [lb,ub] = update_lub(track_w,wb,alpha,c)
hard_bound = track_w / 2 - wb / 2000;
lb = -hard_bound * ones(c,1) - alpha;
ub = hard_bound * ones(c,1) - alpha;
end

%% 生成不等式约束
function [Aieq,bieq] = ieq_gen(track_x,track_y,obs_x,obs_y,track_w,wb,obs_w,alpha,c)
hard_bound = track_w / 2 - wb / 2000;
[num_obs,~] = size(obs_x);
distance = zeros(c,num_obs);
for i = 1:c
    for k = 1:num_obs
        distance(i,k) = (track_x(i) - obs_x(k))^2 + (track_y(i) - obs_y(k))^2;
    end
end
pos = zeros(num_obs,1);
for k = 1:num_obs
    pos(k) = find(distance(:,k)==min(min(distance(:,k))));
end
Aieq = zeros(num_obs,c);
bieq = zeros(num_obs,1);
for k = 1:num_obs
    if (track_x(pos(k)) - obs_x(k) + track_y(pos(k)) - obs_y(k)) >=0
        Aieq(k,pos(k)) = -1;
        bieq(k) = hard_bound + alpha(pos(k)) - sqrt(distance(pos(k),k));
    else
        Aieq(k,pos(k)) = 1;
        bieq(k) = hard_bound - alpha(pos(k)) - sqrt(distance(pos(k),k));
    end
end
end

%% 生成等式约束
function [Aeq,beq] = eq_gen(c) %保证始末点连续
Aeq = zeros(1,c);
Aeq(1) = 1;
Aeq(c) = -1;
beq = [0];
end

%% 生成路径
function [p_x,p_y] = path_gen(track_x,track_y,vector,alpha)
p_x = track_x + alpha .* vector(:,2);
p_y = track_y - alpha .* vector(:,1);
end

%% 绘制优化后路径
function [] = draw_optimised_patn(track_x,track_y,vector,alpha,color_inpt)
[p_x,p_y] = path_gen(track_x,track_y,vector,alpha);
plot(p_x,p_y,color_inpt);
end

%% 线性插值赛道 前处理
function [x_out,y_out] = track_smooth_linear(track_x,track_y,tolerance)
[c,~] = size(track_x);
x_out = [];
y_out = [];
for i = 1:c - 1
    len = sqrt((track_x(i + 1) - track_x(i))^2 + (track_y(i + 1) - track_y(i))^2);
    if len >= tolerance
        n = ceil(len / tolerance);
        y_ran = track_y(i + 1) - track_y(i);
        x_ran = track_x(i + 1) - track_x(i);
        temp_x = zeros(n,1);
        temp_y = zeros(n,1);
        for k = 1:n
            val = (k - 1) / n;
            temp_x(k) = track_x(i) + val * x_ran;
            temp_y(k) = track_y(i) + val * y_ran;
        end
        x_out = [x_out;temp_x];
        y_out = [y_out;temp_y];
    else
        x_out = [x_out;track_x(i)];
        y_out = [y_out;track_y(i)];
    end
end
x_out = [x_out;track_x(i + 1)];
y_out = [y_out;track_y(i + 1)];
end

%% 三次样条插值赛道 第二次处理
function [x_out,y_out] = track_smooth_spline(track_x,track_y,tolerance)
[c,~] = size(track_x);
[spline_A] = spline_matrix_gen(c);
para_x = spline_A * track_x;
para_y = spline_A * track_y;
select_c = 2 * [zeros(c),zeros(c),eye(c),zeros(c)];
Mx = select_c * para_x;
My = select_c * para_y;
x_out = [];
y_out = [];
for i = 1:c - 1
    len = sqrt((track_x(i + 1) - track_x(i))^2 + (track_y(i + 1) - track_y(i))^2);
    if len >= tolerance
        n = ceil(len / tolerance);
        temp_x = zeros(n,1);
        temp_y = zeros(n,1);
        for k = 1:n
            val = (k - 1) / n;
            temp_x(k) = Mx(i) * (1 - val)^3 / 6 ...
                      + Mx(i + 1) * (val - 0)^3 / 6 ...
                      + (track_x(i) - Mx(i) / 6) * (1 - val) ...
                      + (track_x(i + 1) - Mx(i + 1) / 6) * (val - 0);
            temp_y(k) = My(i) * (1 - val)^3 / 6 ...
                      + My(i + 1) * (val - 0)^3 / 6 ...
                      + (track_y(i) - My(i) / 6) * (1 - val) ...
                      + (track_y(i + 1) - My(i + 1) / 6) * (val - 0);
        end
        x_out = [x_out;temp_x];
        y_out = [y_out;temp_y];
    else
        x_out = [x_out;track_x(i)];
        y_out = [y_out;track_y(i)];
    end
end
x_out = [x_out;track_x(i + 1)];
y_out = [y_out;track_y(i + 1)];
end
%%
%%
%%
%%
%%
%%
%%