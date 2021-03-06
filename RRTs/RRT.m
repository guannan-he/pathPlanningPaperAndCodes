% https://www.cnblogs.com/21207-iHome/p/7210543.html
% http://rkala.in/codes/RRT.zip
% 该方法是概率完备且不最优的
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
failCnt = 0;
pathFound = false;
display = true;
RRTree = double([source, -1]);
pathPoint = [goal];
if (display)
    imshow(map);
    rectangle('position',[1 1 size(map)-1],'edgecolor','k');
    hold on;
    plot(source(1), source(2), 'b*');
    plot(goal(1), goal(2), 'r*');
end

%% mainProg
tic;
while (failCnt < maxFailCnt)
    % 有二分之一概率被拉向结束点
    if (rand < 0.5)
        sample = rand(1, 2) .* size(map);
    else
        sample = goal;
    end
    % 取得树上最近点
    [A, I] = min( distanceCost(RRTree(:,1:2),sample) ,[],1); % find the minimum value of each column
    closestNode = RRTree(I(1),1:2);
    % 生成新点
    theta = atan2(sample(1)-closestNode(1),sample(2)-closestNode(2));
    newPoint = double(int32(closestNode(1:2) + stepLen * [sin(theta)  cos(theta)]));
    if (~checkPath(closestNode, newPoint, map))
        failCnt = failCnt + 1;
        continue;
    end
    % 找到路径
    if (distanceCost(newPoint, goal) < threshold)
        pathFound = true;
        RRTree = [RRTree; newPoint, I(1)];
        break;
    end
    % 确保不重复出现在附近区域
    if (min(distanceCost(RRTree(:, 1:2), newPoint)) < threshold)
        failCnt = failCnt + 1;
        continue;
    end
    RRTree = [RRTree; newPoint, I(1)];
    failCnt = 0;
    if (display)
        line([closestNode(2);newPoint(2)],[closestNode(1);newPoint(1)]);
    end
end
if (pathFound)
    pathLen = 0;
    nextCur = size(RRTree, 1);
    lastPoint = RRTree(nextCur, 1:2);
    while (nextCur ~= -1)
        currentPoint = RRTree(nextCur, 1:2);
        pathLen = pathLen + distanceCost(lastPoint, currentPoint);
        pathPoint = [currentPoint; pathPoint];
        nextCur = RRTree(nextCur, 3);
    end
    fprintf("找到路径\n长度%f像素点\n共%d个点\n", pathLen, size(pathPoint, 1));
    if (display)
        for i = 1: size(pathPoint, 1) - 1
            line([pathPoint(i, 2), pathPoint(i + 1, 2)], [pathPoint(i, 1), pathPoint(i + 1, 1)], 'Color', 'r', 'LineWidth', 2);
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