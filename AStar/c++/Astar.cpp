#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <fstream>

#include <io.h>
#include <opencv2/opencv.hpp>

using namespace std;

class commonMethod {
public:
	pair<int, int> getPointCoord(string& str) {//字符串转坐标
		int x = 0, y = 0;
		bool sw = true;
		for (char& ch : str) {
			if (ch < '0' || ch > '9') {
				sw = false;
			}
			else if (sw) {
				x *= 10;
				x += ch - '0';
			}
			else {
				y *= 10;
				y += ch - '0';
			}
		}
		return make_pair(x, y);
	}
};
int main(int argc, char* argv[]) {
	//异常处理
	if (argc < 4) {
		cout << "more arguments" << endl;
		return -1;
	}
	if (argv[1] == "-h" || argv[1] == "-help") {
		cout << "[source] + [result] + [point txt]" << endl;
		return 0;
	}
	//前处理及定义变量
	cv::Mat img = cv::imread(argv[1], 0);
	size_t rolSize = img.rows, colSize = img.cols;
	commonMethod dijkstraMethod;
	vector<vector<int>> dir = { {0, 1}, {1, 0}, {-1, 0}, {0, -1} };
	vector<vector<int>> H(rolSize, vector<int>(colSize, INT_MAX)), G(rolSize, vector<int>(colSize, INT_MAX));
	vector<vector<pair<int, int>>> curMap(rolSize, vector<pair<int, int>>(colSize, { 0, 0 }));
	//读取bmp文件到矩阵
	vector<vector<uchar>> graph(rolSize, vector<uchar>(colSize));
	for (int i = 0; i < rolSize; i++) {
		for (int j = 0; j < colSize; j++) {
			graph[i][j] = img.at<uchar>(j, i);
		}
	}
	/*
	* 读取到的数值中：
	* 0是黑块
	* 50是起始点和终点
	* 100是搜索过的点
	* 150是队列中的点
	* 200是找到的路径
	* 255是空地
	*/
	//读取起始点和终点
	ifstream fin(argv[3]);
	string startStr, endStr;
	getline(fin, startStr);
	getline(fin, endStr);
	pair<int, int> startPoint = dijkstraMethod.getPointCoord(startStr);
	pair<int, int> endPoint = dijkstraMethod.getPointCoord(endStr);
	//初始化矩阵
	for (int i = 0; i < rolSize; i++) {
		for (int j = 0; j < colSize; j++) {
			H[i][j] = 2 * (abs(i - endPoint.first)) + 2 * abs(j - endPoint.second);
			G[i][j] = abs(i - startPoint.first) + abs(j - startPoint.second);
		}
	}
	//主算法
	if (graph[startPoint.first][startPoint.second] == 0 || graph[endPoint.first][endPoint.second] == 0) {
		cout << "invalid points" << endl;
		return -1;
	}
	int tmpVal = H[startPoint.first][startPoint.second] + G[startPoint.first][startPoint.second];
	priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> que;
	que.push(make_pair(tmpVal, startPoint.first * colSize + startPoint.second));
	while (!que.empty()) {
		int pixelCur = que.top().second;
		pair<int, int> cur = make_pair(pixelCur / colSize, pixelCur % colSize);
		que.pop();
		if (graph[cur.first][cur.second] == 100) {
			continue;
		}
		graph[cur.first][cur.second] = 100;
		for (int i = 0; i < 4; i++) {
			int tmpRol = cur.first + dir[i][0];
			int tmpCol = cur.second + dir[i][1];
			if (tmpRol < 0 || tmpRol >= rolSize || tmpCol < 0 || tmpCol >= colSize) {
				continue;
			}
			if (graph[tmpRol][tmpCol] == 0 || graph[tmpRol][tmpCol] == 100 || graph[tmpRol][tmpCol] == 150) {
				continue;
			}
			graph[tmpRol][tmpCol] == 150;
			tmpVal = H[tmpRol][tmpCol] + G[tmpRol][tmpCol];
			que.push(make_pair(tmpVal, tmpRol * colSize + tmpCol));
			curMap[tmpRol][tmpCol] = cur;
		}
		if (cur == endPoint) {
			break;
		}
	}
	//追回最优路径
	pair<int, int> cur = endPoint;
	while (cur != startPoint) {
		graph[cur.first][cur.second] = 200;
		cur = curMap[cur.first][cur.second];
	}
	graph[startPoint.first][startPoint.second] = 50;
	graph[endPoint.first][endPoint.second] = 50;
	//显示或保存图像
	cv::Mat res(rolSize, colSize, CV_8UC3, { 0, 0, 0 });
	for (int i = 0; i < rolSize; i++) {
		for (int j = 0; j < colSize; j++) {
			int graphVal = graph[i][j];
			if (graphVal == 0) {
				continue;
			}
			if (graphVal == 50) {
				res.at<cv::Vec3b>(j, i) = { 0, 255, 0 };
			}
			else if (graphVal == 100) {
				res.at<cv::Vec3b>(j, i) = { 127, 127, 127 };
			}
			else if (graphVal == 150) {
				res.at<cv::Vec3b>(j, i) = { 0, 0, 128 };
			}
			else if (graphVal == 200) {
				res.at<cv::Vec3b>(j, i) = { 255, 255, 0 };
			}
			else {
				res.at<cv::Vec3b>(j, i) = { 255, 255, 255 };
			}
		}
	}
	/*string s = "res";
	cv::imshow(s, res);
	cv::waitKey(0);*/
	cv::imwrite(argv[2], res);
	return 0;
}