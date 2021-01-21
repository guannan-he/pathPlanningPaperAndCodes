%%
clear;
clc;
figure(1);
clf;
cmap = [
    0 0 0; % 0 - black - �ϰ�
    1 1 1; % 1 - white - �յ�
    1 0 0; % 2 - red - ���������ĵط�
    0 0 1; % 3 - blue - �´�������ѡ����
    0 1 0; % 4 - green - ��ʼ��/�յ�
    1 1 0];% 5 - yellow -  ��Ŀ����·��
colormap(cmap);
fid = "map.bmp";
map = int8(imbinarize(imread(fid)));
image(map);
%% ��������
startPt = [5, 5];
endPt = [48, 48];
startIndex = rolCol2Index(startPt(1), startPt(2), 50);
endIndex  = rolCol2Index(endPt(1), endPt(2), 50);
G = inf(size(map));% ���ۺ���
H = inf(size(map));% ��������
curMap = zeros(size(map));% ָ��·��
G(startIndex) = 0;% ��ʼ������
for i = 1: 2500
    cur = index2RolCol(i, 50);
    H(i) = 2 * (abs(cur(1) - 48) + abs(cur(2) - 48));
    G(i) = abs(cur(1) - 5) + abs(cur(2) - 5);
end

map(startIndex) = 4;%������
map(endIndex) = 4;%����յ�
image(map);
print(1,'-dbmp',sprintf('image/%d',1));
fastClose = 1;% ���ٽ������յ���ȫ����С��ֹͣ���£�
%% ��ѭ�� BFS
direction = [
    -1, 0;
    0, 1;
    1, 0;
    0, -1];
queue = [startIndex, H(startIndex) + G(startIndex)];%������ID�����
tic;
while size(queue, 1) ~= 0
    currentIndex = queue(1, 1);
    if fastClose
        if currentIndex == endIndex
            tTotal = toc;
            break;
        end
    end
    queue = queue(2:end, :);%��ǰ�������
    if (map(currentIndex) == 2)% ���������������
        continue;
    end
    map(currentIndex) = 2;
    currentCur = index2RolCol(currentIndex, 50);
    for i = 1:4
        tmpCur = currentCur + direction(i, :);
        if (tmpCur(1) < 1) || (tmpCur(1) > 50) || (tmpCur(2) < 1) || (tmpCur(2) > 50)%�߽�����
            continue;
        end
        tmpIndex = rolCol2Index(tmpCur(1), tmpCur(2), 50);
        if (map(tmpIndex) == 0) || (map(tmpIndex) == 2) || (map(tmpIndex) == 3)
            continue;
        end
        map(tmpIndex) = 3;
        tmpVal = H(tmpIndex) + G(tmpIndex);
        queue = [queue; tmpIndex, tmpVal];%����Ҫ��ڵ������
        curMap(tmpIndex) = currentIndex;
    end
    queue = sortrows(queue, 2); % ��queue����
    image(map);
    pause(0.0001);
end
tTotal = toc;
print(1,'-dbmp',sprintf('image/%d',2));
%% ���ƹ켣
currentIndex = endIndex;
while currentIndex ~= 0
    map(currentIndex) = 5;
    currentIndex = curMap(currentIndex);
    image(map);
    pause(0.0001);
end
%% ��������
map(startIndex) = 4;%������
map(endIndex) = 4;%����յ�
image(map);
print(1,'-dbmp',sprintf('image/%d',3));
fprintf("Done!\n");
fprintf("·������%d\n", H(endIndex));
fprintf("��ʱ%f��\n", tTotal);
%% �Ӻ���
% ������������
function index = rolCol2Index(rolCur, colCur, rol)
index = (colCur - 1) * rol + rolCur;
end
% ��������ת����������
function cur = index2RolCol(index, rol)
cur = [0, 0];
cur(1) = mod(index, rol);
cur(2) = floor(index / rol) + 1;
end