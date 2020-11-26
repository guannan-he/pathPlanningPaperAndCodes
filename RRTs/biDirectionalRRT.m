% https://www.cnblogs.com/21207-iHome/p/7210543.html
% http://rkala.in/codes/RRT.zip
% 该方法是概率完备且不最优的
% 两颗生长树，分别从起点和终点生长，最后生长的点互为目标点
%% preset
% close all;
figure(1);
clf;
clear;
clc;
map = imbinarize(imread('RRT_reference\map2.bmp'));
source = [10, 10];
goal = [490, 490];
stepLen = 20;
threshold = 20;
maxFailCnt = 10000;
display = true;
RRTreeSource = double([source, -1]);
RRTreeGoal = double([goal, -1]);
pathFound = [];%中间连接点
pathPoint = [];%最终路径
tree1ExpansionFail = false; % sets to true if expansion after set number of attempts fails
tree2ExpansionFail = false; % sets to true if expansion after set number of attempts fails
if (display)
    imshow(map);
    rectangle('position',[1 1 size(map)-1],'edgecolor','k');
    hold on;
    plot(source(1), source(2), 'b*');
    plot(goal(1), goal(2), 'r*');
end

%% mainProg
tic;
% 两侧都找不到点就结束
while (~tree1ExpansionFail || ~tree2ExpansionFail)
    % 树1生长
    if (~tree1ExpansionFail)
        [RRTreeSource, pathFound, tree1ExpansionFail] = rrtExtend(RRTreeSource, RRTreeGoal, goal, stepLen, maxFailCnt, threshold, map);
        if (~tree1ExpansionFail && isempty(pathFound) && display)
            line([RRTreeSource(end,2);RRTreeSource(RRTreeSource(end,3),2)],[RRTreeSource(end,1);RRTreeSource(RRTreeSource(end,3),1)],'color','b');
        end
    end
    if (~isempty(pathFound))
        break;
    end
    % 树2生长
    if (~tree2ExpansionFail)
        [RRTreeGoal, pathFound, tree2ExpansionFail] = rrtExtend(RRTreeGoal, RRTreeSource, source, stepLen, maxFailCnt, threshold, map);
        if (~isempty(pathFound))
            pathFound(3 : 4) = pathFound(4 : -1 : 3);%确保pathFound（3）指向1，pathFound（4）指向4
            break;
        end
        if (~tree1ExpansionFail  && display)
            line([RRTreeGoal(end,2);RRTreeGoal(RRTreeGoal(end,3),2)],[RRTreeGoal(end,1);RRTreeGoal(RRTreeGoal(end,3),1)],'color','r');
        end
    end
end
if (~isempty(pathFound))
    % 画中间连接线
    line([pathFound(2), RRTreeSource(pathFound(3), 2)], [pathFound(1), RRTreeSource(pathFound(3), 1)], 'color', 'g');
    line([pathFound(2), RRTreeGoal(pathFound(4), 2)], [pathFound(1), RRTreeGoal(pathFound(4), 1)], 'color', 'g');
    pathLen = 0;
    nextCur = pathFound(3);
    lastPoint = RRTreeSource(nextCur, 1:2);
    % 从尾到头连接路径，起点是连接点
    while (nextCur ~= -1)
        currentPoint = RRTreeSource(nextCur, 1:2);
        pathLen = pathLen + distanceCost(lastPoint, currentPoint);
        pathPoint = [currentPoint; pathPoint];
        nextCur = RRTreeSource(nextCur, 3);
    end
    % 添加连接点
    pathPoint = [pathPoint; pathFound(1), pathFound(2)];
    nextCur = pathFound(4);
    lastPoint = RRTreeGoal(nextCur, 1:2);
    % 从头到尾连接路径，起点是连接点
    while (nextCur ~= -1)
        currentPoint = RRTreeGoal(nextCur, 1:2);
        pathLen = pathLen + distanceCost(lastPoint, currentPoint);
        pathPoint = [pathPoint; currentPoint];
        nextCur = RRTreeGoal(nextCur, 3);
    end
    fprintf("找到路径\n长度%f像素点\n共%d个点\n", pathLen, size(pathPoint, 1));
    if (display)
        for i = 1: size(pathPoint, 1) - 1
            line([pathPoint(i, 2), pathPoint(i + 1, 2)], [pathPoint(i, 1), pathPoint(i + 1, 1)], 'Color', 'c', 'LineWidth', 2);
        end
    end
else
    fprintf("未找到路径");
end

t = toc;
fprintf("耗时%f秒\nDone\n", t);
%% subFunctions
% 检查点是不是在图像范围内且不在障碍物上
function feasible = checkPoint(point, map)
feasible = false;
if (point(1) > 0 && point(1) < size(map, 1) && point(2) > 0 && ...
        point(2) < size(map, 2) && map(point(1), point(2)) == 1)
    feasible = true;
end
end

% 计算点欧式距离
function h = distanceCost(a,b)
h = sqrt(sum((a - b).^2, 2));
end

% 检测新点形成的路径是否合法
function feasible = checkPath(n,newPos,map)
feasible = true;
dir=atan2(newPos(1) - n(1), newPos(2) - n(2));
if (~checkPoint(newPos, map))
    feasible = false;
    return;
end
for r=0: 0.5: sqrt(sum((n - newPos).^2))
    posCheck = n + r .* [sin(dir), cos(dir)];
    if (~(checkPoint(ceil(posCheck), map) && checkPoint(floor(posCheck), map) && ...
            checkPoint([ceil(posCheck(1)), floor(posCheck(2))],map) && ...
            checkPoint([floor(posCheck(1)), ceil(posCheck(2))],map)))
        feasible = false;
        return;
    end
end
end

% 将树扩展函数封装起来
function [RRTree1,pathFound,extendFail] = rrtExtend(RRTree1,RRTree2,goal,stepsize,maxFailCnt,threshold,map)
pathFound=[];
failedAttempts=0;
extendFail = true;
while (failedAttempts <= maxFailCnt)
    if (rand < 0.5)
        sample = rand(1,2) .* size(map);
    else
        sample = goal;
    end
    % 找最近点
    [~, I] = min( distanceCost(RRTree1(:,1:2),sample) ,[],1);
    closestNode = RRTree1(I(1),:);
    % 生成新点
    theta = atan2((sample(1)-closestNode(1)),(sample(2)-closestNode(2)));
    newPoint = double(int32(closestNode(1:2) + stepsize * [sin(theta), cos(theta)]));
    % 新点是否符合要求
    if (~checkPath(closestNode(1:2), newPoint, map))
        failedAttempts = failedAttempts + 1;
        continue;
    end
    % 另一棵树中最近点
    [A, I2] = min(distanceCost(RRTree2(:, 1:2),newPoint), [], 1); % find closest in the second tree
    % 新点与另一棵树最近点也在阈值内
    if (A < threshold)
        % 一个下标指向tree1， 另一个指向tree2
        pathFound = [newPoint I(1) I2(1)];
        extendFail=false;
        break;
    end
    A = min(distanceCost(RRTree1(:, 1:2),newPoint), [], 1); % check if new node is not already pre-existing in the tree
    if (A < threshold)
        failedAttempts = failedAttempts + 1;
        continue;
    end
    RRTree1 = [RRTree1;newPoint I(1)];
    extendFail = false;
    % 找到一个就可以停止
    break;
end
end