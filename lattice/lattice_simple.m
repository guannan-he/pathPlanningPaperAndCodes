% 基于状态空间的采样
%% preProcess
clear;
clc;
figure(1);
clf;
figName = './roadmap.bmp';
road = imread(figName);
[yDirLim, xDirLim] = size(road);
xDirStep = 20;
yDirStep = 25;
xPos = 0;
xPosSet = 0: 1: xDirStep;
yPos = yDirLim / 2;
yPosSet = 0:yDirStep:yDirLim;
yPosSet(1) = 1;
[~, yPosCnt] = size(yPosSet);
imshow(road);
%% mainProgress
while xPos < xDirLim
    % 添加有效终止状态点到集合中
    effectiveY =[];
    for i = 1:yPosCnt
        if noColi(xPos + xDirStep, yPosSet(i), road)
            effectiveY = [effectiveY,yPosSet(i)];
        end
    end
    [~, effectCnt] = size(effectiveY);
    yEval = [];
    xLoca = xPosSet + xPos;
    minCur = 1;
    minVal = inf;
    plotCnt = 0;
    yLocaSet = [];
    for i = 1:effectCnt
        yPoly = getYPoly(xDirStep, yPos, effectiveY(i));% 计算y参数
        yLoca = curveGen(xPosSet, yPoly);% 求y点坐标
        if polynoCoil(xLoca, yLoca, road)
            plotCnt = plotCnt + 1;
            yLocaSet = [yLocaSet; yLoca];
            tmp = evaluate(xLoca, yLoca, road, xDirLim);
            if tmp < minVal
                minVal = tmp;
                minCur = plotCnt;
            end
        end
    end
    for i = 1:plotCnt
        if i == minCur
            hold on;
            plot(xLoca, yLocaSet(i, :), 'r');
            continue;
        end
        hold on;
        plot(xLoca, yLocaSet(i, :), 'k');
    end
    yPos = floor(yLocaSet(minCur, end));
    xPos = xPos + xDirStep;
end
print(1,'-dpng','.\result.png');
%% subFunctions
% 终点是否在有效范围内
function ret = noColi(yPos, xPos, road)
if xPos > 250 || yPos > 1600
    ret = false;
    return;
end
if yPos == 0
    yPos = 1;
end
if xPos == 0
    xPos = 1;
end
ret = road(xPos, yPos) == 1;
end
% 获取多项式参数
function yPoly = getYPoly(x, yPosStart, yPosEnd)
A = [x^5, x^4, x^3, x^2, x, 1
    5 * x^4, 4 * x^3, 3 * x^2, 2 * x, 1, 0
    20 * x^3, 12 * x^2, 6 * x, 2, 0, 0
    0, 0, 0, 0, 0, 1
    0, 0, 0, 0, 1, 0
    0, 0, 0, 2, 0, 0];
b = [yPosEnd; 0; 0; yPosStart; 0; 0];
yPoly = A \ b;
end
% 获取路径坐标
function yLoca = curveGen(xPosSet, yPoly)
yLoca = (((((yPoly(1) * xPosSet) + yPoly(2)) .* xPosSet + yPoly(3)) .* xPosSet + yPoly(4)) .* xPosSet + yPoly(5)) .* xPosSet + yPoly(6);
end
% 检测路径有无碰撞
function ret = polynoCoil(xLoca, yLoca, road)
[~, r] = size(xLoca);
for i = 1: r
    if noColi(xLoca(i), floor(yLoca(i)), road)
        continue
    end
    ret = false;
    return;
end
ret = true;
end
% 计算代价函数
% 评估方法
function ret = evaluate(xLoca, yLoca, road, xDirLim)
xEnd = floor(xLoca(end));
yEnd = floor(yLoca(end));
while xEnd <= xDirLim && noColi(xEnd, yEnd, road)
    xEnd = xEnd + 1;
end
ret = - xEnd + abs(yEnd - yLoca(1));
end