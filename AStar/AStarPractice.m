%%
clear;
clc;
figure(1);
clf;
cmap = [
    0 0 0; % 0 - black - 障碍
    1 1 1; % 1 - white - 空地
    1 0 0; % 2 - red - 已搜索过的地方
    0 0 1; % 3 - blue - 下次搜索备选中心
    0 1 0; % 4 - green - 起始点/终点
    1 1 0];% 5 - yellow -  到目标点的路径
colormap(cmap);
fid = "map.bmp";
map = int8(imbinarize(imread(fid)));
image(map);
%% 参数设置
startPt = [5, 5];
endPt = [48, 48];
startIndex = rolCol2Index(startPt(1), startPt(2), 50);
endIndex  = rolCol2Index(endPt(1), endPt(2), 50);
G = inf(size(map));% 代价函数
H = inf(size(map));% 启发函数
curMap = zeros(size(map));% 指针路径
G(startIndex) = 0;% 初始点置零
for i = 1: 2500
    cur = index2RolCol(i, 50);
    H(i) = 2 * (abs(cur(1) - 48) + abs(cur(2) - 48));
    G(i) = abs(cur(1) - 5) + abs(cur(2) - 5);
end

map(startIndex) = 4;%标记起点
map(endIndex) = 4;%标记终点
image(map);
print(1,'-dbmp',sprintf('image/%d',1));
fastClose = 1;% 快速结束（终点是全局最小就停止更新）
%% 主循环 BFS
direction = [
    -1, 0;
    0, 1;
    1, 0;
    0, -1];
queue = [startIndex, H(startIndex) + G(startIndex)];%出发点ID入队列
tic;
while size(queue, 1) ~= 0
    currentIndex = queue(1, 1);
    if fastClose
        if currentIndex == endIndex
            tTotal = toc;
            break;
        end
    end
    queue = queue(2:end, :);%当前点出队列
    if (map(currentIndex) == 2)% 被其他点完成搜索
        continue;
    end
    map(currentIndex) = 2;
    currentCur = index2RolCol(currentIndex, 50);
    for i = 1:4
        tmpCur = currentCur + direction(i, :);
        if (tmpCur(1) < 1) || (tmpCur(1) > 50) || (tmpCur(2) < 1) || (tmpCur(2) > 50)%边界限制
            continue;
        end
        tmpIndex = rolCol2Index(tmpCur(1), tmpCur(2), 50);
        if (map(tmpIndex) == 0) || (map(tmpIndex) == 2) || (map(tmpIndex) == 3)
            continue;
        end
        map(tmpIndex) = 3;
        tmpVal = H(tmpIndex) + G(tmpIndex);
        queue = [queue; tmpIndex, tmpVal];%符号要求节点入队列
        curMap(tmpIndex) = currentIndex;
    end
    queue = sortrows(queue, 2); % 对queue排序
    image(map);
    pause(0.0001);
end
tTotal = toc;
print(1,'-dbmp',sprintf('image/%d',2));
%% 绘制轨迹
currentIndex = endIndex;
while currentIndex ~= 0
    map(currentIndex) = 5;
    currentIndex = curMap(currentIndex);
    image(map);
    pause(0.0001);
end
%% 后处理杂项
map(startIndex) = 4;%标记起点
map(endIndex) = 4;%标记终点
image(map);
print(1,'-dbmp',sprintf('image/%d',3));
fprintf("Done!\n");
fprintf("路径长度%d\n", H(endIndex));
fprintf("耗时%f秒\n", tTotal);
%% 子函数
% 返回线性序列
function index = rolCol2Index(rolCur, colCur, rol)
index = (colCur - 1) * rol + rolCur;
end
% 线性序列转回坐标序列
function cur = index2RolCol(index, rol)
cur = [0, 0];
cur(1) = mod(index, rol);
cur(2) = floor(index / rol) + 1;
end