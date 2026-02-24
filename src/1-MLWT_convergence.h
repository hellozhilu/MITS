#ifndef TEST_MLWT_V2_0_0_H_
#define TEST_MLWT_V2_0_0_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <libgen.h>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <filesystem>

#include "tools.h"

using namespace std;

#define PARENT(idx) ((idx-1)/2)
#define MAX_VALUE 99999999.0

#define DEBUG_PRT1 0				// 输出禁忌搜索的优化过程
#define DEBUG_PRT2 0				// 输出ILS的优化过程
#define DEBUG_White 0				// 测试逻辑分支覆盖
#define DEBUG_VRF 0					// 验证测试
#define DEBUG_VRF1 0				// 验证测试
#define DEBUG_ils 0					// 验证测试ils
#define DEBUG_cost_tsso 0
#define DEBUG_cost_ts 0
#define USE_NBH0 1
#define USE_NBH1 1
#define USE_NBH2 1
#define USE_NBH3 1
#define DEBUGx 0


// 固定参数
const double epsilon = 0.000001;

int nexit = 0;
vector<double> rcd_cost;
vector<double> rcd_time;


//const double epsilon = 0.001;

class Vertex
{
public:
	int level;							// 当前层级
	int id;								// 在该层级的id 
	vector<double> attrs;				// 在各属性上的权重
	vector<int> pre_vts;				// 序列（原始点中的原始点序列是自己）

	Vertex() {}

	/*
	 * 生成原始点的构造函数
	 */
	Vertex(int id, vector<double> &attrs) : level(0), id(id), attrs(attrs)
	{
		pre_vts.push_back(id);
	}

	/*
	 * 两个上一级点合并成一个点的构造函数
	 */

	Vertex(int level, int id, Vertex &vtx1, Vertex &vtx2) : level(level), id(id), pre_vts({ vtx1.id, vtx2.id })
	{
		// 更新属性权重
		int tau = vtx1.attrs.size();
		this->attrs.resize(tau);
		for (int i = 0; i < tau; i++)
			this->attrs[i] = vtx1.attrs[i] + vtx2.attrs[i];
	}

	/*
	 * 一个点直接变成后一级的压缩点
	 * paramaters:
	 * 		level:新点所在层级
	 * 		id	 :新点的编号（从0开始）
	 * 		vtx1 :要加入新点的上一级（压缩）点
	 */
	Vertex(int level, int id, Vertex &vtx1) : level(level), id(id), attrs(vtx1.attrs), pre_vts({ vtx1.id }) {}

};

class Gain_node {
public:
	int elem1;						// 要移动的点
	int elem2;						// 移动到的分区
	int elem3;						// new pid
	double ifsb_mg;				// infeasible_move_gain or cost_move_gain
	double cost_mg;
	int type;						// type=0表示是移动，type=1表示是交换，type=2表示是push

	Gain_node() : elem1(-1), elem2(-1), elem3(-1), ifsb_mg(MAX_VALUE),  cost_mg(MAX_VALUE), type(-1){}
	Gain_node(int id, int cluster, int elem, double ifsb_mg, double cost_mg, int type) : elem1(id), elem2(cluster), elem3(elem), ifsb_mg(ifsb_mg), cost_mg(cost_mg), type(type){}

	void clear()
	{
		elem1 = -1;
		elem2 = -1;
		ifsb_mg = MAX_VALUE;
		cost_mg = MAX_VALUE;
		type = -1;
	}


	~Gain_node() {}
};

class EdgeNode
{
public:
	int vtx1;
	int vtx2;
	double dist;

	EdgeNode(int vtx1, int vtx2, double dist) : vtx1(vtx1), vtx2(vtx2), dist(dist) {}

	void cpy(EdgeNode &en2){
		vtx1 = en2.vtx1;
		vtx2 = en2.vtx2;
		dist = myRound(en2.dist);
	}
};

class Graph
{
public:
	int n;							// （压缩）点的个数
	int k;							// 分区个数
	int tau;						// 属性个数τ:num_attrs
	int cur_level;					// 当前压缩图的级别
	double cg_dist;					// 被折叠的点之间的距离之和
//	double sum_wavg;				// 计算【每个属性的平均权重】之和
	vector<Vertex> vertices;		// 每个顶点的属性
	vector<double> w_upper;			// W_upper[i]表示每个分区在属性i上的上界
	vector<double> w_lower;			// W_lower[i]表示每个分区在属性i上的下界
	vector<vector<double>> edge;	// 完全图使用邻接矩阵来存储边权重
	vector<vector<vector<double>>> vtx_attr_diff;		// -i+j vtx_attr_diff[i][j][t]表示点i和点j在属性t上的差距，vtx_attr_diff[i][j][t] = attr[j][t] - attr[i][t]
	vector<EdgeNode> ord_edge;		// 顺序存储所有边结点
//	vector<int> heap_pos;			// heap_pos[i] = p表示边界点i在heap中的位置是p

	Graph() {}

	int part(int low, int high) {
	    // 选择枢轴（这里选择最后一个元素作为枢轴）
	    double pivot = ord_edge[high].dist;
	    int i = low - 1;  // 小于枢轴的元素的最后一个索引

	    for (int j = low; j < high; j++) {
	        // 如果当前元素小于或等于枢轴
	        if (ord_edge[j].dist <= pivot) {
	            i++;  // 增加小于枢轴的元素的索引
	            EdgeNode temp_en = ord_edge[i];
	            ord_edge[i].cpy(ord_edge[j]);
	            ord_edge[j].cpy(temp_en);
	        }
	    }

	    // 将枢轴元素放到正确的位置
        EdgeNode temp_en = ord_edge[i+1];
        ord_edge[i+1].cpy(ord_edge[high]);
        ord_edge[high].cpy(temp_en);

        return i + 1;  // 返回枢轴元素的索引
	}

	void quickSort(int low, int high) {
	    if (low < high) {
	        // 获取分区索引
	        int pi = part(low, high);

	        // 递归排序分区
	        quickSort(low, pi - 1);
	        quickSort(pi + 1, high);
	    }
	}

	/*
	 * 构造函数：根据上一级图和新的压缩点集合构造下一级的图
	 */
	Graph(Graph &pre_graph, vector<Vertex> &vertices, double cg_dist)
		: n(vertices.size()),		// n指的是压缩点的个数
		k(pre_graph.k),
		tau(pre_graph.tau),
		cur_level(pre_graph.cur_level + 1),
		cg_dist(cg_dist),
		vertices(vertices),
		w_upper(pre_graph.w_upper),
		w_lower(pre_graph.w_lower)
	{
		try {
			edge.resize(n, vector<double>(n));						// 分配空间，重设大小
		} catch (const bad_alloc& e) {
		    // 内存分配失败，处理异常
		    cerr << "edge Memory allocation failed: " << e.what() << endl;
		}
//		edge.resize(n, vector<double>(n));						// 分配空间，重设大小

		try {
			ord_edge.reserve(n*(n-1)/2);						// 为边堆分配空间
		} catch (const bad_alloc& e) {
		    // 内存分配失败，处理异常
		    cerr << "Memory allocation failed: " << e.what() << endl;
		}

		// 更新edge
		for (int i = 0; i < n; i++)
		{
			int size_i = vertices[i].pre_vts.size();
			for (int j = i + 1; j < n; j++)
			{
				int size_j = vertices[j].pre_vts.size();

				// 遍历所有上一层级点，计算更新后的边（最多2*2）
				double sum_d = 0.0;
				for (int ii = 0; ii < size_i; ii++)				// 点i中的上级点
				{
					int idx_i = vertices[i].pre_vts[ii];
					for (int jj = 0; jj < size_j; jj++)			// 点j中的上级点
					{
						int idx_j = vertices[j].pre_vts[jj];
						sum_d = myRound(sum_d + pre_graph.edge[idx_i][idx_j]);
					}
				}

				// 更新邻接边（无向图）
				edge[i][j] = edge[j][i] = sum_d;

				EdgeNode en(i, j, edge[i][j]);
				ord_edge.push_back(en);
			}
		}

		// 更新vtx_attr_diff
		vtx_attr_diff.resize(n, vector<vector<double>>(n, vector<double>(tau)));
		for(int i = 0; i < n - 1; i++)
		{
			for(int j = i + 1; j < n; j++)
			{
				for(int t = 0; t < tau; t++)
				{
					vtx_attr_diff[i][j][t] = vertices[j].attrs[t] - vertices[i].attrs[t];
					vtx_attr_diff[j][i][t] = -vtx_attr_diff[i][j][t];
				}
			}
		}

		// 更新ord_edge
		int esize = ord_edge.size();
		quickSort(0, esize-1);
	}

	/*
	 * function:读入原始图
	 */
	void read_file(string filename) {
//		cout << "\nfilename = " << filename << endl;

		ifstream fin(filename);		// 打开文件
		// 判断是否成功打开文件
		if (!fin.is_open())
		{
			cerr << "Can not open the file!" << endl;
			exit(-11);
		}

		// 检查文件流是否处于错误状态
		if (fin.fail())
		{
			cerr << "Error occurred during file operation." << filename << endl;
			exit(-12);
		}

		// 判断文件是否为空
		if (fin.eof())
		{
			cerr << "Empty file" << filename << endl;
			exit(-13);
		}

		// 开始读图
		int cur_id = 0;										// 记录当前点id
		cur_level = 0;										// 初始化当前级别为0
		cg_dist = 0.0;

		fin >> n;
		fin >> k;
		tau = 2;

		// 分配内存，并调整大小
		w_upper.resize(tau);
		w_lower.resize(tau);

		//		edge.resize(n);
		vertices.resize(n);									// 分配内存，改变大小
		edge.resize(n, vector<double>(n));
		try {
			ord_edge.reserve(n*(n-1)/2);					// 分配空间，重设大小
		} catch (const bad_alloc& e) {
			// 内存分配失败，处理异常
			cerr << "ord_edge Memory allocation failed-1: " << e.what() << endl;
		}

		// 开始读每一行的信息
		string line;
		while (getline(fin, line))
		{
			if (line.empty()) continue;  // 跳过空行
			istringstream iss(line);		//  转为文件流格式(转换后用的是iss)

			// 每次while循环会读取一整行的信息
			vector<double> attrs(tau);

			for (int i = 0; i < tau; i++) {					// 读每个点在所有属性上的权重
				iss >> attrs[i];							// 输入到数组中
				w_upper[i] = myRound(w_upper[i] + attrs[i]);						// 将W_upper临时用来计算所有点在属性i上的累计权重
			}

			Vertex vtx = Vertex(cur_id, attrs);				// 生成点
			vertices[cur_id] = vtx;							// 将点加入到点集中

			int j = 0;
			while (!iss.eof()) {							// 读邻接边
				iss >> edge[cur_id][j++];
				// 以下两句保证在Visual Studio中也能运行
				if (j == n)
					break;
			}

			cur_id++;
		}

		fin.close();

		// 计算平均值和标准差
		vector<double> avg(tau), dev(tau);
		for (int i = 0; i < tau; i++)
		{
			avg[i] = w_upper[i] / n;
		}

		for (int i = 0; i < tau; i++)
		{
			for (int j = 0; j < n; j++)
			{
				dev[i] += pow(vertices[j].attrs[i] - avg[i], 2) / n;
			}
			dev[i] = sqrt(dev[i]);
		}

		// 计算W_l 和W_u
		for (int i = 0; i < tau; i++)
		{
			w_upper[i] = avg[i] * n / k + 2 * dev[i];
			w_lower[i] = avg[i] * n / k - 2 * dev[i];
		}

		// 更新vtx_attr_diff
		vtx_attr_diff.resize(n, vector<vector<double>>(n, vector<double>(tau)));
		for(int i = 0; i < n - 1; i++)
		{
			for(int j = i + 1; j < n; j++)
			{
				EdgeNode en(i, j, edge[i][j]);
				ord_edge.push_back(en);
				for(int t = 0; t < tau; t++)
				{
					vtx_attr_diff[i][j][t] = vertices[j].attrs[t] - vertices[i].attrs[t];
					vtx_attr_diff[j][i][t] = -vtx_attr_diff[i][j][t];
				}
			}
		}

		// 更新ord_edge
		int esize = ord_edge.size();
		quickSort(0, esize-1);
	}

	/*
	 * 根据2024的图格式来读入图
	 */
	void read_file2(string filename) {
//		cout << "\nfilename = " << filename << endl;

		ifstream fin(filename);		// 打开文件
		// 判断是否成功打开文件
		if (!fin.is_open())
		{
			cerr << "Can not open the file! filename=" << filename <<endl;
			exit(-11);
		}

		// 检查文件流是否处于错误状态
		if (fin.fail())
		{
			cerr << "Error occurred during file operation." << filename << endl;
			exit(-12);
		}

		// 判断文件是否为空
		if (fin.eof())
		{
			cerr << "Empty file" << filename << endl;
			exit(-13);
		}

		// 开始读图
		string instancetype, ignore;
		fin >> instancetype;
		int cur_id = 0;										// 记录当前点id
		cur_level = 0;										// 初始化当前级别为0
		cg_dist = 0.0;

		fin >> n >> k >> tau;;								//

		// 分配内存，并调整大小
		w_upper.resize(tau);
		w_lower.resize(tau);
		vertices.resize(n);									// 分配内存，改变大小
		edge.resize(n, vector<double>(n));
		try {
			ord_edge.reserve(n*(n-1)/2);					// 分配空间，重设大小
		} catch (const bad_alloc& e) {
			// 内存分配失败，处理异常
			cerr << "ord_edge Memory allocation failed-1: " << e.what() << endl;
		}
//		ord_edge.reserve(n*(n-1)/2);							// 为边堆分配空间

		if(instancetype == "p")
		{
			// 开始读每一行的信息
			string line;
			while (getline(fin, line))
			{
				if (line.empty()) continue;  // 跳过空行
				istringstream iss(line);		//  转为文件流格式(转换后用的是iss)

				// 每次while循环会读取一整行的信息
				vector<double> attrs(tau);
				iss >> ignore;
				for (int i = 0; i < tau; i++) {					// 读每个点在所有属性上的权重
					iss >> attrs[i];							// 输入到数组中
					w_upper[i] = myRound(w_upper[i] + attrs[i]);						// 将W_upper临时用来计算所有点在属性i上的累计权重
				}

				Vertex vtx = Vertex(cur_id, attrs);				// 生成点
				vertices[cur_id] = vtx;							// 将点加入到点集中

				int j = 0;
				while (!iss.eof()) {							// 读邻接边
					iss >> edge[cur_id][j++];
					// 以下两句保证在Visual Studio中也能运行
					if (j == n)
						break;
				}

				cur_id++;
			}
		}
		else if(instancetype == "t")
		{
			// 开始读每一行的信息
			string line;
			int **coordinates = new int*[n];
			for(int i = 0; i < n; i++)
			{
				coordinates[i] = new int[2];
			}
			while (getline(fin, line))
			{
				if (line.empty()) continue;  // 跳过空行
				istringstream iss(line);		//  转为文件流格式(转换后用的是iss)

				// 每次while循环会读取一整行的信息
				// 读点的各属性权重
				vector<double> attrs(tau);
				iss >> ignore;									// 标识去掉
				for (int i = 0; i < tau; i++) {					// 读每个点在所有属性上的权重
					iss >> attrs[i];							// 输入到数组中
					w_upper[i] = myRound(w_upper[i] + attrs[i]);						// 将W_upper临时用来计算所有点在属性i上的累计权重

				}
				Vertex vtx = Vertex(cur_id, attrs);				// 生成点
				vertices[cur_id] = vtx;							// 将点加入到点集中

				// 读点的坐标
				iss >> coordinates[cur_id][0] >> coordinates[cur_id][1];
				cur_id++;
			}
			// 计算边权
			unsigned decimals = 6;
			double factor = pow(10, decimals);
			for(int i = 0; i < n; i++) {
				for(int j = i + 1; j < n; j++) {
					edge[i][j] = round(sqrt((pow(coordinates[i][0]-coordinates[j][0],2) + pow(coordinates[i][1]-coordinates[j][1],2)))*factor)/factor;
					edge[j][i] = edge[i][j];
				}
			}

			// 释放辅助数组内存
			for(int i = 0; i < n; i++)
			{
				delete[] coordinates[i];
				coordinates[i] = NULL;
			}
			delete[] coordinates;
			coordinates = NULL;
		}
		else
		{
			cout << "\n Wrong instance file.";
			cin.get();
			exit(1);
		}

		fin.close();

		// 计算平均值和标准差
		vector<double> avg(tau), dev(tau);
		for (int i = 0; i < tau; i++)
		{
			avg[i] = w_upper[i] / n;
		}

		for (int i = 0; i < tau; i++)
		{
			for (int j = 0; j < n; j++)
			{
				dev[i] += pow(vertices[j].attrs[i] - avg[i], 2) / n;
			}
			dev[i] = sqrt(dev[i]);
		}

		// 计算W_l 和W_u
		unsigned decimals = 2;
		double factor = pow(10, decimals);
		for (int i = 0; i < tau; i++)
		{
			w_upper[i] = round((avg[i]*((double)n/(double)k) + 2.0*dev[i])*factor)/factor;
			w_lower[i] = round((avg[i]*((double)n/(double)k) - 2.0*dev[i])*factor)/factor;
		}

		// 更新vtx_attr_diff
		vtx_attr_diff.resize(n, vector<vector<double>>(n, vector<double>(tau)));
		for(int i = 0; i < n - 1; i++)
		{
			for(int j = i + 1; j < n; j++)
			{
				EdgeNode en(i, j, edge[i][j]);
				ord_edge.push_back(en);
				for(int t = 0; t < tau; t++)
				{
					vtx_attr_diff[i][j][t] = vertices[j].attrs[t] - vertices[i].attrs[t];
					vtx_attr_diff[j][i][t] = -vtx_attr_diff[i][j][t];
				}
			}
		}

		// 更新ord_edge
		int esize = ord_edge.size();
		quickSort(0, esize-1);
	}


	void print_graph() {
		// 输出图基本信息
		cout << "n=" << n << ", k=" << k << ", tau=" << tau << endl;

		// 输出属性限制信息
		for (int i = 0; i < tau; i++)
			cout << "w_u_" << i << "=" << w_upper[i] << "w_l_" << i << "=" << w_lower[i] << endl;

		// 输出每行的信息（与原始图一致）
		for (int i = 0; i < n; i++)
		{
			cout << "---origin" << i << endl;

			// 输出点权重
			for (int j = 0; j < tau; j++)
				cout << vertices[i].attrs[j] << " ";

			cout << "||";
			// 输出边权重
			for (int j = 0; j < n; j++)
				cout << edge[i][j] << " ";
			cout << endl;
		}
	}

	void print_coarsen_graph() {
		// 输出图基本信息
		cout << "n=" << n << endl;

		// 输出点
		for (int i = 0; i < n; i++)
		{
			cout << "coarsen_vertex_" << i << ": ";
			int num_vtx = vertices[i].pre_vts.size();
			for (int j = 0; j < num_vtx; j++)
				cout << vertices[i].pre_vts[j] << " ";

			cout << "|| attrs = ";
			for (int t = 0; t < tau; t++)
				cout << vertices[i].attrs[t] << " ";
			cout << endl;
		}

		// 输出图（邻接矩阵）
		for (int i = 0; i < n; i++)
		{
			cout << i << "--";
			for (int j = 0; j < n; j++)
				cout << edge[i][j] << "(" << j << ")" << " ";
			cout << endl;
		}
	}

	~Graph()
	{

	}
};

class Solution
{
public:
	int level;
	double ifsb;
	double dist;

	short int *ptn;									// partition，sol[i] = n表示折叠点i被分到分区n 
	short int *cluster_size;												// cluster_size[k]=n指的是分区k中有n个点
	double **dist_vtx_to_cluster;											// dist_vtx_to_cluster[i][k]表示点i到分区k的距离
	double **cluster_attrs;													// cluster_weight[k][t]表示分区k在属性t上的总权重
	double *cluster_infeasible;												// cluster_infeasible[k]表示分区k的不可行程度
	double **ifsb_change;													// attrs_change[i][k]表示点i移入/移出分区k引起的不可行程度变化
	double **dist_change;													// dist_change[i][k]表示点i移入/移出分区k引起的dist变化

	Solution() {
		ifsb = MAX_VALUE;
		dist = MAX_VALUE;
	}

	/*
	 * 初始化：调整大小，分配内存
	 */
	Solution(Graph &graph) : level(graph.cur_level), ifsb(MAX_VALUE), dist(MAX_VALUE)
	{
		// 分配内存
		ptn = new short int[graph.n];
		memset(ptn, 0, sizeof(short int)*graph.n);

		cluster_size = new short int[graph.k];
		memset(cluster_size, 0, sizeof(short int)*graph.k);

		cluster_infeasible = new double[graph.k];
		memset(cluster_infeasible, 0, sizeof(double)*graph.k);

		dist_vtx_to_cluster = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			dist_vtx_to_cluster[i] = new double[graph.k];
			memset(dist_vtx_to_cluster[i], 0, sizeof(double)*graph.k);
		}

		cluster_attrs = new double *[graph.k];
		for (int i = 0; i < graph.k; i++)
		{
			cluster_attrs[i] = new double[graph.tau];
			memset(cluster_attrs[i], 0, sizeof(double)*graph.tau);
		}

		ifsb_change = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			ifsb_change[i] = new double[graph.k];
			memset(ifsb_change[i], 0, sizeof(double)*graph.k);
		}

		dist_change = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			dist_change[i] = new double[graph.k];
			memset(dist_change[i], 0, sizeof(double)*graph.k);
		}
	}


	Solution(Solution &sol1, Graph &graph) : level(graph.cur_level), ifsb(sol1.ifsb), dist(sol1.dist)
	{
		// 分配内存
		ptn = new short int[graph.n];
		memcpy(ptn, sol1.ptn, sizeof(short int)*graph.n);

		cluster_size = new short int[graph.k];
		memcpy(cluster_size, sol1.cluster_size, sizeof(short int)*graph.k);

		cluster_infeasible = new double[graph.k];
		memcpy(cluster_infeasible, sol1.cluster_infeasible, sizeof(double)*graph.k);

		dist_vtx_to_cluster = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			dist_vtx_to_cluster[i] = new double[graph.k];
			memcpy(dist_vtx_to_cluster[i], sol1.dist_vtx_to_cluster[i], sizeof(double)*graph.k);
		}

		cluster_attrs = new double *[graph.k];
		for (int i = 0; i < graph.k; i++)
		{
			cluster_attrs[i] = new double[graph.tau];
			memcpy(cluster_attrs[i], sol1.cluster_attrs[i], sizeof(double)*graph.tau);
		}

		ifsb_change = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			ifsb_change[i] = new double[graph.k];
			memcpy(ifsb_change[i], sol1.ifsb_change[i], sizeof(double)*graph.k);
		}

		dist_change = new double *[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			dist_change[i] = new double[graph.k];
			memcpy(dist_change[i], sol1.dist_change[i], sizeof(double)*graph.k);
		}
	}


	/*
	 * 验证不可解程度
	 */
	void verify_infeasible(Graph &graph, string error_info)
	{
		bool is_i_equal = false;

		// 初始化数组
		double **clst_attrs = new double *[graph.k];
		for (int i = 0; i < graph.k; i++)
		{
			clst_attrs[i] = new double[graph.tau];
			memset(clst_attrs[i], 0, sizeof(double)*graph.tau);
		}

		// 更新每个分区的各属性权重之和
		for (int i = 0; i < graph.n; i++)
		{
			int clst = ptn[i];									// 点i所在的分区clst

			// 各属性来说
			for (int t = 0; t < graph.tau; t++)
			{
				clst_attrs[clst][t] = myRound(clst_attrs[clst][t] + graph.vertices[i].attrs[t]);
			}
		}

		// 更新总的不可解程度
		double v_ifsb = 0.0;
		for (int i = 0; i < graph.k; i++)
		{
			double vcifsb = 0.0;
			for (int t = 0; t < graph.tau; t++)
			{
				if (clst_attrs[i][t] > graph.w_upper[t])
					vcifsb = myRound(vcifsb + clst_attrs[i][t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst_attrs[i][t] < graph.w_lower[t])
					vcifsb = myRound(vcifsb + graph.w_lower[t] - clst_attrs[i][t]);		// 下溢的部分
			}
			if(myRound(vcifsb) - myRound(cluster_infeasible[i]) > epsilon)
			{
				printf("分区 %d 的不可解程度错误\n", i);
				fflush(stdout);
			}
			v_ifsb = myRound(v_ifsb + vcifsb);
		}

		// 释放内存
		for (int i = 0; i < graph.k; i++)
		{
			delete[] clst_attrs[i];
		}
		delete[] clst_attrs;
		clst_attrs = NULL;

		if (abs(myRound(v_ifsb) - myRound(ifsb)) < epsilon)
			is_i_equal = true;

		// 输出错误信息
		if (is_i_equal == false)
		{
			cout << "【！！！不可解程度验证错误---】" << error_info << "出错:" << "verify_infeasible = " << v_ifsb << ", cur_infeasible = " << ifsb << endl;
			exit(-4);
		}
		else
			cout << " infeasible √ ";
	}

	/*
	 * 验证ifsb_change数组
	 */
	void verify_ifsb_change(Graph &graph, string error_info)
	{
		// 重新计算ifsb_change
		double **v_ifsb_change = new double*[graph.n];
		for (int i = 0; i < graph.n; i++)
		{
			v_ifsb_change[i] = new double[graph.k];
			memset(v_ifsb_change[i], 0, sizeof(double) * graph.k);
		}

		for (int i = 0; i < graph.n; i++)
		{
			int pre_clst = ptn[i];

			// k1=sol[i]中移除点i
			double quit_i = 0;

			for (int t = 0; t < graph.tau; t++)
			{
				if (cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
					quit_i = myRound(quit_i + cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
					quit_i = myRound(quit_i + graph.w_lower[t] - cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
			}
			v_ifsb_change[i][pre_clst] = quit_i - cluster_infeasible[pre_clst];						// 变化

			for (int next_clst = 0; next_clst < graph.k; next_clst++)
			{
				if (pre_clst != next_clst)
				{
					// k2中加入点i
					double join_i = 0;
					for (int t = 0; t < graph.tau; t++)
					{
						if (cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
							join_i = myRound(join_i + cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
						else if (cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
							join_i = myRound(join_i + graph.w_lower[t] - cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
					}
					v_ifsb_change[i][next_clst] = join_i - cluster_infeasible[next_clst];
				}
			}
		}

		// 对比
		bool flag = true;
		for (int i = 0; i < graph.n; i++)
		{
			for (int j = 0; j < graph.k; j++)
			{
				if (myRound(v_ifsb_change[i][j]) != myRound(ifsb_change[i][j]))
				{
					cout << "\n" << error_info << "!error--------ifsb_change_" << i << "_" << j << "，vic=" << v_ifsb_change[i][j] << "  ic=" << ifsb_change[i][j] << endl;
					flag = false;
				}

			}
		}

		if(flag)
			cout << " ifsb_change √ ";
		else
			exit(-8);

		// 回收内存
		for (int i = 0; i < graph.n; i++)
			delete[] v_ifsb_change[i];
		delete[] v_ifsb_change;
	}

	/*
	 * 验证dist
	 */
	void verify_dist(Graph &graph, string error_info)
	{
		bool is_i_equal = false;

		double vdist = myRound(graph.cg_dist);
		for (int i = 0; i < graph.n-1; i++)
		{
			for (int j = i + 1; j < graph.n; j++)
				if (ptn[i] == ptn[j])
					vdist = myRound(vdist + graph.edge[i][j]);
		}

		if (myRound(vdist) == myRound(dist))
			is_i_equal = true;

		// 输出错误信息
		if (is_i_equal == false)
		{
			cout << "【！！！目标函数值验证错误---】" << error_info << "出错:" << "verify_dist = " << vdist << ", dist = " << dist << endl;
			exit(-5);
		}
		else
			cout << error_info << " cost √ \n";
	}

	/*
	 * 验证每个cluster的属性权重
	 */
	void verify_cluster_attrs(Graph &graph, string error_info)
	{
		double **v_cluster_attrs = new double*[graph.k];
		for (int i = 0; i < graph.k; i++)
		{
			v_cluster_attrs[i] = new double[graph.tau];
			memset(v_cluster_attrs[i], 0, sizeof(double)*graph.tau);
		}
		for (int i = 0; i < graph.n; i++)						// 对于每个分区
		{
			int clst = ptn[i];
			// 更新分区的属性权重
			for (int j = 0; j < graph.tau; j++)
				v_cluster_attrs[clst][j] = myRound(v_cluster_attrs[clst][j] + graph.vertices[i].attrs[j]);
		}

		// 验证
		bool flag = true;
		for (int i = 0; i < graph.k; i++)
		{
			for (int t = 0; t < graph.tau; t++)
			{
				if (myRound(v_cluster_attrs[i][t]) != myRound(cluster_attrs[i][t]))
				{
					cout << "\n" << error_info << "!error--------cluster_attrs_" << i << "_" << t << ", vca=" << v_cluster_attrs[i][t] << "  ca=" << cluster_attrs[i][t] << endl;
					flag = false;
				}
			}
		}
		if(flag)
			cout << " cluster_attrs √ ";
		else
			exit(-7);

		// 清除空间
		for (int i = 0; i < graph.k; i++)
			delete[] v_cluster_attrs[i];
		delete[] v_cluster_attrs;
	}

	/*
	 * 输出解决方案
	 */
	void print_sol(Graph &graph)
	{
		// 输出ptn数组
		for (int i = 0; i < graph.n; i++)
			cout << "v_" << i << "(" << ptn[i] << ")" << " ";
		cout << endl;

		// 输出cluster_size数组
		for (int i = 0; i < graph.k; i++)
			cout << cluster_size[i] << " ";
		cout << endl;
	}

	/*
	 * 随机初始化解决方案
	 */
	void init_sol(Graph &graph)
	{
		int size = graph.n;

		// 随机初始化
		int *aux_vct = new int[size];
		for (int i = 0; i < graph.n; i++)
			aux_vct[i] = i;

		// 先给每个分区分一个
		for (int i = 0; i < graph.k; i++)			// 对于每个分区
		{
			int rnd = rand() % size;				// 随机选择一个点
			int idx = aux_vct[rnd];					// 确认点的下标
			ptn[idx] = i;							// 更新对应的解决方案

			// 更新分区的属性权重
			for (int j = 0; j < graph.tau; j++)
			{
				cluster_attrs[i][j] = myRound(cluster_attrs[i][j] + graph.vertices[idx].attrs[j]);
			}

			// 更新每个分区的大小
			cluster_size[i] = 1;
			size--;

			aux_vct[rnd] = aux_vct[size];
		}

		// 剩下的随机选分区
		for (int i = 0; i < size; i++)						// 对于剩余的每个点
		{
			int idx = aux_vct[i];							// 当前点的编号
			int rnd = rand() % graph.k;						// 随机选一个分区

			ptn[idx] = rnd;									// 更新解决方案

			// 更新分区的属性权重
			for (int j = 0; j < graph.tau; j++)
			{
				cluster_attrs[rnd][j] = myRound(cluster_attrs[rnd][j] + graph.vertices[idx].attrs[j]);
			}

			cluster_size[rnd]++;
		}

//		print_sol(graph);

		// 释放空间
		delete[] aux_vct;
		aux_vct = NULL;
	}

	/*
	 * edge是当前level的压缩图的邻接矩阵
	 */
	void cal_dist(Graph &graph)
	{
		// 根据解决方案的内容计算总距离并初始化move_gain需要的东西

		// 遍历所有原始点
		for (int i = 0; i < graph.n; i++)
		{
			// 更新点i与所有点所在分区的相关边距离
			for (int j = i + 1; j < graph.n; j++)
			{
				// 距离相关
				dist_vtx_to_cluster[i][ptn[j]] = myRound(dist_vtx_to_cluster[i][ptn[j]] + graph.edge[i][j]);
				dist_vtx_to_cluster[j][ptn[i]] = myRound(dist_vtx_to_cluster[j][ptn[i]] + graph.edge[i][j]);

				if (ptn[i] == ptn[j])
					dist = myRound(dist + graph.edge[i][j]);
			}
		}
	}

	/*
	 * 计算不可行程度
	 */
	void cal_infeasible(Graph &graph)
	{
		// 遍历所有分区
		for (int i = 0; i < graph.k; i++)
		{
			for (int t = 0; t < graph.tau; t++)
			{
				// 上限
				if (cluster_attrs[i][t] > graph.w_upper[t])
					cluster_infeasible[i] = myRound(cluster_infeasible[i] + cluster_attrs[i][t] - graph.w_upper[t]);
				// 下限
				if (cluster_attrs[i][t] < graph.w_lower[t])
					cluster_infeasible[i] = myRound(cluster_infeasible[i] + graph.w_lower[t] - cluster_attrs[i][t]);

				// 加和
				ifsb = myRound(ifsb + cluster_infeasible[i]);
			}
		}
	}

	/*
	 * 将sol2中的关键信息复制到当前sol中
	 */
	void cpy(Solution &sol2, Graph &graph)
	{
		dist = myRound(sol2.dist);
		ifsb = myRound(sol2.ifsb);
		memcpy(ptn, sol2.ptn, sizeof(short int)*graph.n);
		memcpy(cluster_size, sol2.cluster_size, sizeof(short int)*graph.k);
		memcpy(cluster_infeasible, sol2.cluster_infeasible, sizeof(double)*graph.k);
		for (int i = 0; i < graph.n; i++)
		{
			memcpy(dist_vtx_to_cluster[i], sol2.dist_vtx_to_cluster[i], sizeof(double)*graph.k);
		}
		for (int i = 0; i < graph.k; i++)
		{
			memcpy(cluster_attrs[i], sol2.cluster_attrs[i], sizeof(double)*graph.tau);
		}
		for (int i = 0; i < graph.n; i++)
		{
			memcpy(ifsb_change[i], sol2.ifsb_change[i], sizeof(double)*graph.k);
		}
		for (int i = 0; i < graph.n; i++)
		{
			memcpy(dist_change[i], sol2.dist_change[i], sizeof(double)*graph.k);
		}
	}

	/*
	 * 释放内存
	 */
	void free_memory(Graph &graph)
	{
		delete[] ptn;
		delete[] cluster_size;
		delete[] cluster_infeasible;

		for (int i = 0; i < graph.n; i++)
			delete[] dist_vtx_to_cluster[i];
		delete[] dist_vtx_to_cluster;

		for (int i = 0; i < graph.k; i++)
			delete[] cluster_attrs[i];
		delete[] cluster_attrs;

		for (int i = 0; i < graph.n; i++)
			delete[] ifsb_change[i];
		delete[] ifsb_change;

		for (int i = 0; i < graph.n; i++)
			delete[] dist_change[i];
		delete[] dist_change;
	}
};

/* ============================================= 全局变量 ============================================= */
//vector<vector<vector<int>>> cg_vertices;				// cg_vertices[level][i]表示level级压缩图中no=i的压缩点中包含的原始点id数组
Gain_node min_i_gain;
vector<Gain_node> candidate_gain;										// 用于选相同move_gain的邻域移动
vector<Gain_node> tabu_candidate_gain;
Gain_node min_c_gain;

int cur_lvl;											// 当前层级
//int max_lvl = 15;										// 压缩到的最大层级最多15
int max_cv_count = 1;									// 每一层最多折叠的点数
//double pct_mcc = 0.1;
vector<int> cg_size;									// cg_size[level]=n说明level级压缩图有n个压缩点
//vector<vector<vector<int>>> cg_edge;					// cg_edge[level]表示level级压缩图的邻接矩阵
vector<bool> matched;									// matched[i] = true表示（压缩）点i已经被匹配
//vector<int> aux_match;								// 剩余未匹配的
vector<Graph> coarsen_graph;							// coarsen_graph[level]是层级level的压缩图
//double best_dist;										// 最好的dist，即最短的分区内距离之和

// 最优记录
Solution glb_best_sol;									// global best sol
Solution ml_best_sol;
Solution ils_best_sol;
Solution ts_best_sol;									// tabu search best sol
bool is_hit;

double glb_best_cost;
double ml_best_cost;
double ils_best_cost;
double ts_best_cost;
double ts_best_ifsb;
double target_cost;
double rbest_cost;			// 用于更新best time
double rbest_time;

// 过程记录相关
int runs = 1;
const int max_iter = 1;
double *each_run_rlt;
double *each_run_time;
//double *each_hit_time;
vector<Solution> each_best_sol;
clock_t start, end;
double hit_time = 0;
double hit_cost = 0;
double avg_cost, avg_time, avg_htime, std_dev;
double sum_avg_cost = 0, sum_avg_time = 0, sum_avg_htime = 0;
int hit;

// SOS
int iter_mult;
int MAXLEN;

// ils相关
const int max_frozen = 20;
//const double np_pct = 0.2;												// 扰动强度0.35
//int num_perturb;

// tabu相关
short int **tabu_list;													// 禁忌表
//const double tt_pct = 0.2;												// 0.16
//const double td_pct = 0.3;												// 0.4

//const double np_pct = 0.2;												// 扰动强度0.35
//const double mti_pct = 0.5;											// 禁忌搜索的最大迭代次数
int tabu_tenure;														// 禁忌期
int tabu_depth;															// 允许的最大连续未改进次数
int ts_iter;															// 禁忌搜索的当前迭代次数
//int num_perturb;

//int max_ts_iter;


// 战略震荡SO相关
//double so_rate = 0.2;
double all_ifsb_upper;

int seed = 14;
char filename[1000] = "./Instances/pollster/muestra1_30_4.txt";
string graph_name;
double tmax = 10;

/* =============================================调参的参数============================================= */
int max_lvl = 7;										// 压缩到的最大层级最多
//const double so_rate = 0.1;												// 战略震荡的范围
double pct_mcc = 0.04;											// 每一层合并的点数
double tt_pct = 0.9;													// tabu tenure
double td_pct = 0.80;													// tabu depth
//const int max_iter = 1;													// ils的重启次数
double np_pct = 0.1;												// 扰动强度
double param_beta = 1;													// 系数
int nrcd = 5;
double param_brate = 2.0;

/*
 * 读取参数
 */
void Read_Parameters(int argc, char **argv)
{
//	for(int i = 1; i < argc; i++)
//	{
//		cout << argv[i] << endl;
//	}

	for (int i = 1; i < argc; i += 2)							// [0]表示'-'，[1]表示类型，[2]表示数据
	{
		if (argv[i][0] != '-')
		{
//			Show_Usage();
			exit(0);
		}
		else if (argv[i][1] == 'i')  // The file name
		{
			strncpy(filename, argv[i + 1], 1000);
//		cerr << filename << endl;
		}
		else if (argv[i][1] == 's') // The maximum time
		{
			seed = atoi(argv[i + 1]);
		}
		else if (argv[i][1] == 'm') // The maximum time
			max_lvl = atoi(argv[i + 1]);			// ils 的迭代次数
		else if (argv[i][1] == 'c') //not used
			pct_mcc = atof(argv[i + 1]);			// ts的最大迭代次数比率
		else if (argv[i][1] == 't')
			tt_pct = atof(argv[i + 1]);				// 禁忌期比率
		else if (argv[i][1] == 'd')
			td_pct = atof(argv[i + 1]);				// 连续未改进次数比率
		else if (argv[i][1] == 'n') //not used
			np_pct = atof(argv[i + 1]);			// ts的最大迭代次数比率
		else if (argv[i][1] == 'b')
			param_beta = atof(argv[i + 1]);				// 连续未改进次数比率
		else if (argv[i][1] == 'r')
			param_brate = atof(argv[i + 1]);				// 禁忌期比率


	}
	// check parameters
//	if (strlen(filename) == 0)
//	{
//		cerr << "No input data" << endl;
//		exit(1);
//	}
}


/*
 * 分配内存
 */
void allocate_memery(Graph &graph)
{
	//	best_dist = -1.0;											
//	max_cv_count = ceil(graph.n * pct_mcc);

	// set cutoff time
	tmax = 10;
	if(graph.n >=100 && graph.n < 200)
		tmax = 60;
	else if(graph.n >= 200 && graph.n < 400)
		tmax = 300;
	else if(graph.n >= 400 && graph.n < 500)
		tmax = 1500;

	max_cv_count = int(graph.n * pct_mcc);

	coarsen_graph.resize(max_lvl + 1);							// 给 压缩图数组 预先分配空间

	each_run_rlt = new double[runs];
	memset(each_run_rlt, 0, sizeof(double)*runs);

	each_run_time = new double[runs];
	memset(each_run_time, 0, sizeof(double)*runs);

//	each_hit_time = new double[runs];
//	memset(each_hit_time, 0, sizeof(double)*runs);

	each_best_sol.resize(runs);

	// 禁忌搜索
	tabu_tenure = ceil(tt_pct * graph.n);
	tabu_depth = ceil(td_pct * graph.n);
//	tabu_depth = 20;
//	max_ts_iter = ceil(mti_pct * graph.n * graph.k);
	tabu_list = new short int *[graph.n];
	for(int i = 0; i < graph.n; i++)
	{
		tabu_list[i] = new short int[graph.k];
		memset(tabu_list[i], 0, sizeof(short int)*graph.k);
	}

	// 最优记录
	glb_best_sol = Solution(graph);
	ml_best_sol = Solution(graph);

	// 战略震荡
//	all_ifsb_upper = graph.sum_wavg * so_rate;

	rcd_cost.reserve(10000);
	rcd_time.reserve(10000);
}


/*
 * 释放内存
 */
void free_memory(Graph &graph)
{
	delete[] each_run_rlt;
	delete[] each_run_time;
//	delete[] each_hit_time;

	each_run_rlt = NULL;
	each_run_time = NULL;
//	each_hit_time = NULL;

	for(int i = 0; i <runs; i++)
		each_best_sol[i].free_memory(graph);

	// 禁忌搜索的参数
	for(int i = 0; i < graph.n; i++)
	{
		delete[] tabu_list[i];
		tabu_list[i] = NULL;
	}
	delete[] tabu_list;
	tabu_list = NULL;

	// 类中的空间释放
	glb_best_sol.free_memory(graph);
	ml_best_sol.free_memory(graph);
	coarsen_graph.clear();
}


/*
 * 从小跟堆中取一个符合条件的边,该边在ord_edge的位置为last_idx
 * ①该边的两个端点都没有被标记过 ②两个点合并之后属性权重不会超过上限 ③在以上限制中找距离最小的边
 */
void find_edge(Graph &graph, int &last_idx)
{
	/* 1.从小跟堆中找到符合条件的dist最小的边 */
	// 从last_idx开始找起
	int size = graph.ord_edge.size();
	while(last_idx < size)
	{
		if(matched[graph.ord_edge[last_idx].vtx1] == false && matched[graph.ord_edge[last_idx].vtx2] == false)
		{
			// 检查属性权重是否超过上限
			bool flag = true;
			for (int t = 0; t < graph.tau; t++)
				if (graph.vertices[graph.ord_edge[last_idx].vtx1].attrs[t] + graph.vertices[graph.ord_edge[last_idx].vtx2].attrs[t] > graph.w_upper[t])
					flag = false;

			// 如果满足条件②：两个点合并之后属性权重不会超过上限
			if(flag)
				break;
		}
		last_idx++;
	}

	// 如果没有符合条件的边
	if(last_idx >= size)
		last_idx = -1;

	/* 2.更新match相关 */
	if(last_idx != -1)
	{
		matched[graph.ord_edge[last_idx].vtx1] = true;
		matched[graph.ord_edge[last_idx].vtx2] = true;
	}
}


/*
 * 将原始图压缩到max_lvl
 */
void coarsen_process(Graph &graph)
{
	cur_lvl = 0;
	coarsen_graph[0] = graph;													// 将原始图添加到压缩图数组中

	double cg_dist = 0.0;														// 被压缩的点的权重（相当于小的图中的边权之和）
	while (cur_lvl < max_lvl)
	{
		/* 1.重置参数，初始化*/
		// 初始化
		int cur_cvid = 0;														// 用于分配id
		vector<Vertex> vts;														// 未知大小，最好使用reserve和push_back
		matched.assign(coarsen_graph[cur_lvl].n, 0);							// 全置为0
		vts.resize(coarsen_graph[cur_lvl].n);									// 最差情况下都不能折叠

		/* 2.确定并生成新的压缩点 */
		int cur_cvc = 0;															// 当前层级已经折叠的点数
		int last_idx = 0;

		/* 判断是否需要折叠*/
		if(max_cv_count == 0)
			break;

		// 开始折叠
		while(cur_cvc < max_cv_count)												// 需要折叠max_cv_count个点
		{
			find_edge(coarsen_graph[cur_lvl], last_idx);

			if (last_idx != -1)														//  存在匹配点
			{
				Vertex vtx1 = coarsen_graph[cur_lvl].vertices[coarsen_graph[cur_lvl].ord_edge[last_idx].vtx1];
				Vertex vtx2 = coarsen_graph[cur_lvl].vertices[coarsen_graph[cur_lvl].ord_edge[last_idx].vtx2];
				Vertex c_vtx = Vertex(cur_lvl + 1, cur_cvid, vtx1, vtx2);
				cg_dist = myRound(cg_dist + coarsen_graph[cur_lvl].ord_edge[last_idx].dist);
				vts[cur_cvid++] = c_vtx;
				cur_cvc++;														// 折叠成功的点
			}
			else
				break;																// 如果没有符合条件的就退出
			//			cout << "Debug" << cur_cvid << "------idx1=" << idx1 << ", idx2=" << idx2 << endl;	// To_debug
		}

		/* 3.剩余点直接进下一层图 */
		for(int i = 0; i < coarsen_graph[cur_lvl].n; i++)
		{
			if(matched[i] == false)													// 这里没必要改变matched的值
			{
				Vertex c_vtx = Vertex(cur_lvl + 1, cur_cvid, coarsen_graph[cur_lvl].vertices[i]);
				vts[cur_cvid++] = c_vtx;
//				cur_cvid;
			}
		}

		// 调整为真实大小
		vts.resize(cur_cvid);														// 这里不需要加一，也不需要减一，刚刚好

		// 判断是否需要继续折叠（保证每个分区都至少有一个折叠点）
		if (cur_cvid < 2*graph.k)														
			break;

		// 填充cg的内容
		Graph new_cg = Graph(coarsen_graph[cur_lvl], vts, cg_dist);

		cur_lvl++;
		coarsen_graph[cur_lvl] = new_cg;
		//		new_cg.print_coarsen_graph();												// 输出折叠点信息
	}
}


/*
 * 初始化
 */
void init_i_move_gain(Solution &sol, Graph &graph, Gain_node &min_i_gain)
{
	/* 1.重置参数 */
	sol.ifsb = 0.0;
	memset(sol.cluster_infeasible, 0, sizeof(double)*graph.k);
	for (int i = 0; i < graph.n; i++)
		memset(sol.ifsb_change[i], 0, sizeof(double)*graph.k);

	// 更新不可解程度
	for (int i = 0; i < graph.k; i++)
	{
		for (int j = 0; j < graph.tau; j++)
		{
			if (sol.cluster_attrs[i][j] > graph.w_upper[j])
				sol.cluster_infeasible[i] = myRound(sol.cluster_infeasible[i] + sol.cluster_attrs[i][j] - graph.w_upper[j]);		// 上溢的部分
			else if (sol.cluster_attrs[i][j] < graph.w_lower[j])
				sol.cluster_infeasible[i] = myRound(sol.cluster_infeasible[i] + graph.w_lower[j] - sol.cluster_attrs[i][j]);		// 下溢的部分
		}
		sol.ifsb = myRound(sol.ifsb + sol.cluster_infeasible[i]);
	}

	candidate_gain.reserve(50);													 
	candidate_gain.resize(0);													// 这里可以删掉，在这里重复添加是为了防止其他地方没有及时清空
	double min_gain = MAX_VALUE;
	min_i_gain.clear();

	// 计算move gain
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		// k1=sol[i]中移除点i
		double quit_i = 0;

		for (int t = 0; t < graph.tau; t++)
		{
			if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
				quit_i = myRound(quit_i + sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
				quit_i = myRound(quit_i + graph.w_lower[t] - sol.cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
		}
		sol.ifsb_change[i][pre_clst] = myRound(quit_i - sol.cluster_infeasible[pre_clst]);						// 变化

		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				// k2中加入点i
				double join_i = 0;
				for (int t = 0; t < graph.tau; t++)
				{
					if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
						join_i = myRound(join_i + sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
						join_i = myRound(join_i + graph.w_lower[t] - sol.cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
				}
				sol.ifsb_change[i][next_clst] = myRound(join_i - sol.cluster_infeasible[next_clst]);

				// infeasible相关计算和更新
				double i_delta = myRound(sol.ifsb_change[i][pre_clst] + sol.ifsb_change[i][next_clst]);
#if(USE_NBH0)
				if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) < myRound(min_gain))								// 能移动的才移动
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
					min_gain = i_delta;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
					candidate_gain.push_back(new_node);
				}
#endif
			}
		}
	}

	// 寻找最小邻域移动（exchange）
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if (myRound(ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
					min_gain = ifsb_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
					candidate_gain.push_back(new_node);
				}
			}
		}
	}
#endif
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int vtx1 = idx_list[i];
		int vtx2 = idx_list[(i+1) % graph.n];
		int vtx3 = idx_list[(i+2) % graph.n];
		int move_clst1 = sol.ptn[vtx1];
		int move_clst2 = sol.ptn[vtx2];
		int move_clst3 = sol.ptn[vtx3];

		// 探索 push
#if(USE_NBH2)
		if(move_clst1 != move_clst2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == move_clst1 ||  npid == move_clst2)
					continue;
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] - graph.vertices[vtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[vtx2][vtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[vtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				if (myRound(ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(vtx1, vtx2, npid, ifsb_change, MAX_VALUE, 2);
					min_gain = ifsb_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(vtx1, vtx2, npid, ifsb_change, MAX_VALUE, 2);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(move_clst1 != move_clst2 && move_clst1 != move_clst3 && move_clst2 != move_clst3)
		{
			double ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[vtx1][vtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[vtx2][vtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[move_clst3][t] + graph.vtx_attr_diff[vtx3][vtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2] - sol.cluster_infeasible[move_clst3]); 						// 如果是infeasible下降阶段就不是-0了

			if (myRound(ifsb_change) < myRound(min_gain))
			{
				Gain_node new_node = Gain_node(vtx1, vtx2, vtx3, ifsb_change, MAX_VALUE, 3);
				min_gain = ifsb_change;
				candidate_gain.resize(1);
				candidate_gain[0] = new_node;
			}
			else if (myRound(ifsb_change) == myRound(min_gain))
			{
				Gain_node new_node = Gain_node(vtx1, vtx2, vtx3, ifsb_change, MAX_VALUE, 3);
				candidate_gain.push_back(new_node);
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	idx_list = NULL;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;

	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_i_gain = candidate_gain[rnd];
	}
	else
		min_i_gain.clear();

	candidate_gain.resize(0);

	// 验证
//	sol.verify_ifsb_change(graph, "initiate infeasible之后");
}


/*
 * 初始化
 */
void init_i_move_gain2(Solution &sol, Graph &graph, Gain_node &min_i_gain)
{
	/* 1.重置参数 */
	sol.ifsb = 0.0;
	memset(sol.cluster_infeasible, 0, sizeof(double)*graph.k);
	for (int i = 0; i < graph.n; i++)
		memset(sol.ifsb_change[i], 0, sizeof(double)*graph.k);

	// 更新不可解程度
	candidate_gain.reserve(50);	
	candidate_gain.resize(0);
	double min_gain = MAX_VALUE;
	min_i_gain.clear();

	// 计算move gain
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		// k1=sol[i]中移除点i
		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				// k2中加入点i
				// infeasible相关计算和更新
				double i_delta = myRound(sol.ifsb_change[i][pre_clst] + sol.ifsb_change[i][next_clst]);

#if(USE_NBH0)
				if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) < myRound(min_gain))								// 能移动的才移动
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
					min_gain = i_delta;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
					candidate_gain.push_back(new_node);
				}
#endif
			}
		}
	}

	// 寻找最小邻域移动（exchange）
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了


				if (myRound(ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
					min_gain = ifsb_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
					candidate_gain.push_back(new_node);
				}
			}
		}
	}
#endif

	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];

	for(int i = 0; i < graph.n; i++)
	{
		int vtx1 = idx_list[i];
		int vtx2 = idx_list[(i+1) % graph.n];
		int vtx3 = idx_list[(i+2) % graph.n];
		int move_clst1 = sol.ptn[vtx1];
		int move_clst2 = sol.ptn[vtx2];
		int move_clst3 = sol.ptn[vtx3];

#if(USE_NBH2)
		// 探索 push
		if(move_clst1 != move_clst2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == move_clst1 ||  npid == move_clst2)
					continue;
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] - graph.vertices[vtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[vtx2][vtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[vtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - 0 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了

				if (myRound(ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(vtx1, vtx2, npid, ifsb_change, MAX_VALUE, 2);
					min_gain = ifsb_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(vtx1, vtx2, npid, ifsb_change, MAX_VALUE, 2);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(move_clst1 != move_clst2 && move_clst1 != move_clst3 && move_clst2 != move_clst3)
		{
			double ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[vtx1][vtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[vtx2][vtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[move_clst3][t] + graph.vtx_attr_diff[vtx3][vtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了

			if (myRound(ifsb_change) < myRound(min_gain))
			{
				Gain_node new_node = Gain_node(vtx1, vtx2, vtx3, ifsb_change, MAX_VALUE, 3);
				min_gain = ifsb_change;
				candidate_gain.resize(1);
				candidate_gain[0] = new_node;
			}
			else if (myRound(ifsb_change) == myRound(min_gain))
			{
				Gain_node new_node = Gain_node(vtx1, vtx2, vtx3, ifsb_change, MAX_VALUE, 3);
				candidate_gain.push_back(new_node);
			}
		}
#endif
	}


	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	idx_list = NULL;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;

	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_i_gain = candidate_gain[rnd];
	}
	else
		min_i_gain.clear();

	candidate_gain.resize(0);

	// 验证
//	sol.verify_ifsb_change(graph, "initiate infeasible之后");
}


/*
 * 初始化cost move gain
 */
void init_c_move_gain(Solution &sol, Graph &graph, Gain_node &min_c_gain)
{
	/* 1.重置参数 */
	candidate_gain.reserve(50);													 
	double min_gain = MAX_VALUE;

	sol.dist = myRound(graph.cg_dist);
	candidate_gain.resize(0);																		 
	min_c_gain.clear();

	for (int i = 0; i < graph.n; i++)
	{
		memset(sol.dist_vtx_to_cluster[i], 0, sizeof(double)*graph.k);
		memset(sol.dist_change[i], 0, sizeof(double)*graph.k);
	}

	/* 2.计算并更新 */
	// 更新dist_vtx_to_cluster，dist,寻找最小的邻域移动
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		// 更新点i与所有点所在分区的相关边距离
		for (int j = i + 1; j < graph.n; j++)
		{
			sol.dist_vtx_to_cluster[i][sol.ptn[j]] = myRound(sol.dist_vtx_to_cluster[i][sol.ptn[j]] + graph.edge[i][j]);
			sol.dist_vtx_to_cluster[j][sol.ptn[i]] = myRound(sol.dist_vtx_to_cluster[j][sol.ptn[i]] + graph.edge[i][j]);

			if (sol.ptn[i] == sol.ptn[j])
				sol.dist = myRound(sol.dist + graph.edge[i][j]);
		}

		// 计算并更新move gain，寻找最小邻域移动(insert)
		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				// infeasible相关
				double ifsb_change = myRound(sol.ifsb_change[i][next_clst] + sol.ifsb_change[i][pre_clst]);
				if(abs(ifsb_change) < epsilon)
				{
					// dist相关
					sol.dist_change[i][next_clst] = myRound(sol.dist_vtx_to_cluster[i][next_clst] - sol.dist_vtx_to_cluster[i][pre_clst]);
#if(USE_NBH0)
					if (myRound(sol.dist_change[i][next_clst]) < myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						min_gain = sol.dist_change[i][next_clst];
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if ((myRound(sol.dist_change[i][next_clst]) == myRound(min_gain)) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						candidate_gain.push_back(new_node);
					}
#endif
				}
			}
		}
	}

	// 寻找最小邻域移动（exchange）
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double dist_change, ifsb_change;
				dist_change = myRound(sol.dist_vtx_to_cluster[i][move_clst2] + sol.dist_vtx_to_cluster[j][move_clst1] - sol.dist_vtx_to_cluster[i][move_clst1]  - sol.dist_vtx_to_cluster[j][move_clst2] - 2.0*graph.edge[i][j]);

				// infeasible相关

				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					if (myRound(dist_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						min_gain = dist_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
	}
#endif


	/* 4.部分探索其他邻域 */
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];

		// 探索 push
#if(USE_NBH2)
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;

				double dist_change, ifsb_change;
				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					// 计算cost变化
					dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] + sol.dist_vtx_to_cluster[nvtx2][npid]
								- sol.dist_vtx_to_cluster[nvtx1][npid1] - sol.dist_vtx_to_cluster[nvtx2][npid2]
								- graph.edge[nvtx1][nvtx2]);

					// 更新move gain结点
					if (myRound(dist_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						min_gain = dist_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{
			double dist_change, ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid3]); 						// 如果是infeasible下降阶段就不是-0了


			if(abs(ifsb_change) < epsilon)
			{
				// 计算cost变化
				dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] - sol.dist_vtx_to_cluster[nvtx1][npid1]
							+ sol.dist_vtx_to_cluster[nvtx2][npid3] - sol.dist_vtx_to_cluster[nvtx2][npid2]
							+ sol.dist_vtx_to_cluster[nvtx3][npid1] - sol.dist_vtx_to_cluster[nvtx3][npid3]
							- graph.edge[nvtx1][nvtx2] - graph.edge[nvtx1][nvtx3] - graph.edge[nvtx2][nvtx3]);

				// 更新move gain结点
				if (myRound(dist_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					min_gain = dist_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;
	idx_list = NULL;

	/* 4.选择邻域移动 */
	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_c_gain = candidate_gain[rnd];
	}
	else
		min_c_gain.clear();

	candidate_gain.resize(0);

	//	verify_infeasible(cur_sol, "初始化dist");
}

/*
 * 初始化cost move gain，使用扩展的目标函数
 */
void init_c_move_gain_so(Solution &sol, Graph &graph, Gain_node &min_c_gain)
{
	/* 1.重置参数 */
	candidate_gain.reserve(50);													 
	double min_gain = MAX_VALUE;

	sol.dist = myRound(graph.cg_dist);
	candidate_gain.resize(0);																		 
	min_c_gain.clear();

	for (int i = 0; i < graph.n; i++)
	{
		memset(sol.dist_vtx_to_cluster[i], 0, sizeof(double)*graph.k);
		memset(sol.dist_change[i], 0, sizeof(double)*graph.k);
	}

	/* 2.计算并更新 */
	// 更新dist_vtx_to_cluster，dist,寻找最小的邻域移动
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		// 更新点i与所有点所在分区的相关边距离
		for (int j = i + 1; j < graph.n; j++)
		{
			sol.dist_vtx_to_cluster[i][sol.ptn[j]] = myRound(sol.dist_vtx_to_cluster[i][sol.ptn[j]] + graph.edge[i][j]);
			sol.dist_vtx_to_cluster[j][sol.ptn[i]] = myRound(sol.dist_vtx_to_cluster[j][sol.ptn[i]] + graph.edge[i][j]);

			if (sol.ptn[i] == sol.ptn[j])
				sol.dist = myRound(sol.dist + graph.edge[i][j]);
		}

		// 计算并更新move gain，寻找最小邻域移动(insert)
		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				// infeasible相关
				double ifsb_change = myRound(sol.ifsb_change[i][next_clst] + sol.ifsb_change[i][pre_clst]);
				// dist相关
				sol.dist_change[i][next_clst] = myRound(sol.dist_vtx_to_cluster[i][next_clst] - sol.dist_vtx_to_cluster[i][pre_clst]);
#if(USE_NBH0)
				if (myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) < myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
					min_gain = myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change);
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if ((myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) == myRound(min_gain)) && sol.cluster_size[pre_clst] > 1)
				{
					Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
					candidate_gain.push_back(new_node);
				}
#endif
			}
		}
	}

	// 寻找最小邻域移动（exchange）
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double dist_change, ifsb_change;
				dist_change = myRound(sol.dist_vtx_to_cluster[i][move_clst2] + sol.dist_vtx_to_cluster[j][move_clst1] - sol.dist_vtx_to_cluster[i][move_clst1]  - sol.dist_vtx_to_cluster[j][move_clst2] - 2.0*graph.edge[i][j]);

				// infeasible相关

				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
					min_gain = myRound(dist_change + param_beta*ifsb_change);
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
					candidate_gain.push_back(new_node);
				}
			}
		}
	}
#endif


	/* 4.部分探索其他邻域 */
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];

		// 探索 push
#if(USE_NBH2)
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;

				double dist_change, ifsb_change;
				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				// 计算cost变化
				dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] + sol.dist_vtx_to_cluster[nvtx2][npid]
							- sol.dist_vtx_to_cluster[nvtx1][npid1] - sol.dist_vtx_to_cluster[nvtx2][npid2]
							- graph.edge[nvtx1][nvtx2]);

				// 更新move gain结点
				if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
					min_gain = myRound(dist_change + param_beta*ifsb_change);
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{
			double dist_change, ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid3]); 						// 如果是infeasible下降阶段就不是-0了

			// 计算cost变化
			dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] - sol.dist_vtx_to_cluster[nvtx1][npid1]
						+ sol.dist_vtx_to_cluster[nvtx2][npid3] - sol.dist_vtx_to_cluster[nvtx2][npid2]
						+ sol.dist_vtx_to_cluster[nvtx3][npid1] - sol.dist_vtx_to_cluster[nvtx3][npid3]
						- graph.edge[nvtx1][nvtx2] - graph.edge[nvtx1][nvtx3] - graph.edge[nvtx2][nvtx3]);

			// 更新move gain结点
			if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
			{
				Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
				min_gain = myRound(dist_change + param_beta*ifsb_change);
				candidate_gain.resize(1);
				candidate_gain[0] = new_node;
			}
			else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
			{
				Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
				candidate_gain.push_back(new_node);
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;
	idx_list = NULL;

	/* 4.选择邻域移动 */
	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_c_gain = candidate_gain[rnd];
	}
	else
		min_c_gain.clear();

	candidate_gain.resize(0);

	//	verify_infeasible(cur_sol, "初始化dist");
}


/*
 * 只更新可行程度的gain相关
 */
void update_i_move_gain_ts(Solution &sol, Graph &graph, Gain_node cur_move)
{
//	sol.verify_ifsb_change(graph, "move gain");
#if(DEBUG_VRF)
	sol.verify_cluster_attrs(graph, "tabu更新之前");
	sol.verify_ifsb_change(graph, "tabu更新之前");
#endif
	int move_vtx1;
	int move_vtx2;
	int move_vtx3;
	int pid1;								// 点i所在的cluster
	int pid2;
	int pid3;

	/* 1.更新部分sol内容：*/
	if(cur_move.type == 0)
	{
		move_vtx1 = cur_move.elem1;			// 点1
		move_vtx2 = -1;
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];			// 点1所在的原分区
		pid2 = cur_move.elem2;				// 点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[cur_move.elem1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vertices[cur_move.elem1].attrs[t]);
		}
		sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.ifsb_change[cur_move.elem1][pid1]);
		sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.ifsb_change[cur_move.elem1][pid2]);

		// 更新move_gain值
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.cluster_size[pid2]++;
	}
	else if(cur_move.type == 1)
	{
		move_vtx1 = cur_move.elem1;			// 点1
		move_vtx2 = cur_move.elem2;			// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];			// 点1所在的原分区，点2要移向的新分区
		pid2 = sol.ptn[move_vtx2];			// 点2所在的原分区，点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx2][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分
		}
 

		// 更新move_gain值
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
//		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid1;
//		sol.cluster_size[pid2]++;
	}
	else if(cur_move.type == 2)
	{
		move_vtx1 = cur_move.elem1;
		move_vtx2 = cur_move.elem2;
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];					// 点1所在的原分区
		pid2 = sol.ptn[move_vtx2];					// 点2所在的原分区，点1要移向的新分区
		pid3 = cur_move.elem3;						// 点2要移向的新分区

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vertices[move_vtx2].attrs[t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}
 

		// 更新move_gain值
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.cluster_size[pid3]++;
	}
	else if(cur_move.type == 3)
	{
		move_vtx1 = cur_move.elem1;
		move_vtx2 = cur_move.elem2;
		move_vtx3 = cur_move.elem3;
		pid1 = sol.ptn[cur_move.elem1];				// 点1所在的原分区，点3要移向的新分区
		pid2 = sol.ptn[cur_move.elem2];				// 点2所在的原分区，点1要移向的新分区
		pid3 = sol.ptn[cur_move.elem3];				// 点3所在的原分区，点2要移向的新分区

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx3][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vtx_attr_diff[move_vtx3][move_vtx2][t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}
 

		// 更新move_gain值
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;
		tabu_list[move_vtx3][pid3] = ts_iter + tabu_tenure;


		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.ptn[move_vtx3] = pid1;
	}

	//	verify_infeasible(sol, "move gain");															// DEBUG

	// 寻找最小的移动（insert）
//	vector<Gain_node> candidate_gain;
	candidate_gain.resize(0);																		 
	tabu_candidate_gain.resize(0);
	double min_gain = MAX_VALUE;
	double tabu_min_gain = MAX_VALUE;
	min_i_gain.clear();
#if(DEBUG_White)
	printf("邻域遍历1之前 —— min_gain=%f\n", min_gain);
	fflush(stdout);
#endif
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		double quit_i = 0;																			// k1中移出点i引起的infeasible degree变化
		if ((pre_clst == pid1 || pre_clst == pid2 || pre_clst == pid3))		// 更新部分v_ifsb_change
		{
			for (int t = 0; t < graph.tau; t++)
			{
				if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
					quit_i = myRound(quit_i + sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
					quit_i = myRound(quit_i + graph.w_lower[t] - sol.cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
			}
			sol.ifsb_change[i][pre_clst] = myRound(quit_i - sol.cluster_infeasible[pre_clst]);						// 变化
		}

		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				double join_i = 0;																			// k2中加入点i引起的infeasible degree变化
				if ((next_clst == pid1 || next_clst == pid2 || next_clst == pid3))
				{
					for (int t = 0; t < graph.tau; t++)
					{
						if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
							join_i = myRound(join_i + sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
						else if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
							join_i = myRound(join_i + graph.w_lower[t] - sol.cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
					}
					sol.ifsb_change[i][next_clst] = myRound(join_i - sol.cluster_infeasible[next_clst]);
				}

				// infeasible相关计算和更新
				double i_delta = myRound(sol.ifsb_change[i][pre_clst] + sol.ifsb_change[i][next_clst]);
#if(USE_NBH0)
				if(tabu_list[i][next_clst] > ts_iter)
				{
					if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) < myRound(tabu_min_gain))								// 能移动的才移动
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
						tabu_min_gain = i_delta;
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else
				{
					if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) < myRound(min_gain))								// 能移动的才移动
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
						min_gain = i_delta;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (sol.cluster_size[sol.ptn[i]] > 1 && myRound(i_delta) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, i_delta, MAX_VALUE, 0);
						candidate_gain.push_back(new_node);
					}
				}
#endif
			}
		}
	}

#if(DEBUG_White)
	printf("邻域遍历2之前 —— min_gain=%f\n", min_gain);
	fflush(stdout);
#endif
	// 寻找最小邻域移动（exchange）
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];					// 点1所在的分区
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if(tabu_list[i][move_clst2] > ts_iter && tabu_list[j][move_clst1] > ts_iter)
				{
					if (myRound(ifsb_change) < myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
						tabu_min_gain = ifsb_change;
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (myRound(ifsb_change) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else
				{
					if (myRound(ifsb_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
						min_gain = ifsb_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(ifsb_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, MAX_VALUE, 1);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
	}
#endif
#if(DEBUG_White)
	printf("邻域遍历3/4之前 —— min_gain=%f\n", min_gain);
	fflush(stdout);
#endif
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];

#if(USE_NBH2)
		// 探索 push
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;
				double ifsb_change;

				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				if(tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid] > ts_iter)
				{
					if (myRound(ifsb_change) < myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, MAX_VALUE, 2);
						tabu_min_gain = ifsb_change;
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (myRound(ifsb_change) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, MAX_VALUE, 2);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else
				{
					if (myRound(ifsb_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, MAX_VALUE, 2);
						min_gain = ifsb_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(ifsb_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, MAX_VALUE, 2);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{

			double ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid3]); 						// 如果是infeasible下降阶段就不是-0了

			if(tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid3] > ts_iter && tabu_list[nvtx3][npid1] > ts_iter)
			{
				if (myRound(ifsb_change) < myRound(tabu_min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, MAX_VALUE, 3);
					tabu_min_gain = ifsb_change;
					tabu_candidate_gain.resize(1);
					tabu_candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(tabu_min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, MAX_VALUE, 3);
					tabu_candidate_gain.push_back(new_node);
				}
			}
			else
			{
				if (myRound(ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, MAX_VALUE, 3);
					min_gain = ifsb_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(ifsb_change) == myRound(min_gain))	
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, MAX_VALUE, 3);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
	}

#if(DEBUG_White)
	printf("邻域遍历3/4之后 —— min_gain=%f\n", min_gain);
	fflush(stdout);
#endif
	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	idx_list = NULL;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;

	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	int ts_size_cl = tabu_candidate_gain.size();
#if(DEBUG_White)
	printf("size_cl=%d\n", size_cl);
	fflush(stdout);
#endif
	if (ts_size_cl > 0 && (min_gain > tabu_min_gain && myRound(tabu_min_gain + sol.ifsb) == 0))
	{
		int rnd = rand() % ts_size_cl;
		min_i_gain = tabu_candidate_gain[rnd];
	}
	else if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_i_gain = candidate_gain[rnd];
	}
	else
		min_i_gain.clear();

	//return min_i_gain;
	//	verify_infeasible(sol, "邻域移动");															// DEBUG
}


/*
 * 更新纯下降搜索移动的gain相关:修改邻域结构，部分探索多邻域（不带SO）
 */
void update_c_move_gain2(Solution &sol, Graph &graph, Gain_node cur_move)
{
	int move_vtx1;
	int move_vtx2;
	int move_vtx3;
	int pid1;
	int pid2;
	int pid3;

	//	verify_ifsb_change(sol, "update_c_move gain");

	/* 1.更新部分sol内容：*/
	if(cur_move.type == 0)
	{
		move_vtx1 = cur_move.elem1;
		move_vtx2 = -1;
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];			// 点1所在的原分区
		pid2 = cur_move.elem2;				// 点1要移向的新分区
		pid3 = -1;
		// 更新权重,两个分区的不可解程度
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vertices[move_vtx1].attrs[t]);
		}

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i != cur_move.elem1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
		}

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[cur_move.elem1] = pid2;
		sol.cluster_size[pid2]++;
	}
	else if(cur_move.type == 1)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点2要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = -1;
		// 更新权重,两个分区的不可解程度
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx2][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);
		}

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
			}
		}

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid1;
	}
	else if(cur_move.type == 2)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];							// 点1所在的原分区
		pid2 = sol.ptn[move_vtx2];							// 点2所在的原分区， 点1要移向的新分区
		pid3 = cur_move.elem3;								// 点2要移向的新分区
//		move_to_clst1 = pid2;
//		move_to_clst2 = pid3;

		// 更新各属性上的权重
		for (int t = 0; t < graph.tau; t++)
		{

			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vertices[move_vtx2].attrs[t]);
		}
 

		// 更新cost值
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.cluster_size[pid3]++;
	}
	else if(cur_move.type == 3)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = cur_move.elem3;								// 点3
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点3要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = sol.ptn[move_vtx3];								// 点3所在的原分区，点2要移向的新分区

		// 更新属性权重
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx3][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vtx_attr_diff[move_vtx3][move_vtx2][t]);
		}
 

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i]);
			}
			else if(i == move_vtx3)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
		}

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.ptn[move_vtx3] = pid1;
	}

	// 重置参数
	candidate_gain.resize(0);																		 
	double min_gain = 0;
	min_c_gain.clear();																				// 重置最小的move gain

	/* 2.寻找最小的移动（insert），并更新相关参数*/

	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		double quit_i = 0;																			// k1中移出点i引起的infeasible degree变化
		if ((pre_clst == pid1 || pre_clst == pid2 || pre_clst == pid3))		// 更新部分v_ifsb_change
		{
			for (int t = 0; t < graph.tau; t++)
			{
				if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
					quit_i = myRound(quit_i + sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
					quit_i = myRound(quit_i + graph.w_lower[t] - sol.cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
			}
			sol.ifsb_change[i][pre_clst] = myRound(quit_i - sol.cluster_infeasible[pre_clst]);						// 变化
		}

		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				double join_i = 0;																			// k2中加入点i引起的infeasible degree变化
				if ((next_clst == pid1 || next_clst == pid2 || next_clst == pid3))
				{
					for (int t = 0; t < graph.tau; t++)
					{
						if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
							join_i = myRound(join_i + sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
						else if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
							join_i = myRound(join_i + graph.w_lower[t] - sol.cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
					}
					sol.ifsb_change[i][next_clst] = myRound(join_i - sol.cluster_infeasible[next_clst]);
				}

				double ifsb_change = myRound(sol.ifsb_change[i][next_clst] + sol.ifsb_change[i][pre_clst]);

				// 不会改变不可解程度才执行以下内容
				if(abs(ifsb_change) < epsilon)
				{
					// 计算并更新move gain
					sol.dist_change[i][next_clst] = myRound(sol.dist_vtx_to_cluster[i][next_clst] - sol.dist_vtx_to_cluster[i][pre_clst]);
#if(USE_NBH0)
					if (myRound(sol.dist_change[i][next_clst]) < myRound(min_gain)  && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						min_gain = sol.dist_change[i][next_clst];
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(sol.dist_change[i][next_clst]) == myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						candidate_gain.push_back(new_node);
					}
#endif
				}
			}
		}
	}

	/* 3.寻找最小邻域移动（exchange）*/
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double dist_change, ifsb_change;
				dist_change = myRound(sol.dist_vtx_to_cluster[i][move_clst2] + sol.dist_vtx_to_cluster[j][move_clst1] - sol.dist_vtx_to_cluster[i][move_clst1]  - sol.dist_vtx_to_cluster[j][move_clst2] -2*graph.edge[i][j]);

				// infeasible相关

				// 计算交换后的attr
				double ic1 = 0;						// exchange后clst1的infeasible
				double ic2 = 0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					if (myRound(dist_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						min_gain = dist_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
	}
#endif

	/* 4.部分探索其他邻域 */
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];

#if(USE_NBH2)
		// 探索 push
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;

				double dist_change, ifsb_change;
				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - 0 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					// 计算cost变化
					dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] + sol.dist_vtx_to_cluster[nvtx2][npid]
								- sol.dist_vtx_to_cluster[nvtx1][npid1] - sol.dist_vtx_to_cluster[nvtx2][npid2]
								- graph.edge[nvtx1][nvtx2]);

					// 更新move gain结点
					if (myRound(dist_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						min_gain = dist_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{
			double dist_change, ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - 0 - 0); 						// 如果是infeasible下降阶段就不是-0了

			if(abs(ifsb_change) < epsilon)
			{
				// 计算cost变化
				dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] - sol.dist_vtx_to_cluster[nvtx1][npid1]
							+ sol.dist_vtx_to_cluster[nvtx2][npid3] - sol.dist_vtx_to_cluster[nvtx2][npid2]
							+ sol.dist_vtx_to_cluster[nvtx3][npid1] - sol.dist_vtx_to_cluster[nvtx3][npid3]
							- graph.edge[nvtx1][nvtx2] - graph.edge[nvtx1][nvtx3] - graph.edge[nvtx2][nvtx3]);

				// 更新move gain结点
				if (myRound(dist_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					min_gain = dist_change;
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;
	idx_list = NULL;

	/* 4.选择邻域移动 */
	int size_cl = candidate_gain.size();

	if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_c_gain = candidate_gain[rnd];
	}
	else
		min_c_gain.clear();

//	sol.verify_clst_attrs(graph, "禁忌邻域移动");
//	sol.verify_ifsb_change(graph, "禁忌邻域移动");
//	sol.verify_infeasible(graph, "禁忌邻域移动");											// DEBUG
//	sol.verify_dist(graph, "禁忌邻域移动");
}

/*
 * 更新禁忌搜索移动的gain相关:修改邻域结构，部分探索多邻域（不带SO）,会更新move_flag数组
 */
void update_tabu_move_gain2(Solution &sol, Graph &graph, Gain_node cur_move, int *move_flag)
{
	int move_vtx1;
	int move_vtx2;
	int move_vtx3;
	int pid1;
	int pid2;
	int pid3;

	//	verify_ifsb_change(sol, "update_c_move gain");

	/* 1.更新部分sol内容：*/
	if(cur_move.type == 0)
	{
		move_vtx1 = cur_move.elem1;
		move_vtx2 = -1;
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];			// 点1所在的原分区
		pid2 = cur_move.elem2;				// 点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[cur_move.elem1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vertices[cur_move.elem1].attrs[t]);
		}
		sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.ifsb_change[cur_move.elem1][pid1]);
		sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.ifsb_change[cur_move.elem1][pid2]);


		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i != cur_move.elem1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]));
				sol.dist_vtx_to_cluster[i][pid2] = myRound(myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]));
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
	//	tabu_list[cur_move.vertex] += tabu_tenure;				// 不知道为什么，这个的效果更好

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.cluster_size[pid2]++;
		move_flag[move_vtx1]++;
	}
	else if(cur_move.type == 1)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点2要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx2][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分
		}

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid1;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
	}
	else if(cur_move.type == 2)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];							// 点1所在的原分区
		pid2 = sol.ptn[move_vtx2];							// 点2所在的原分区， 点1要移向的新分区
		pid3 = cur_move.elem3;								// 点2要移向的新分区
//		move_to_clst1 = pid2;
//		move_to_clst2 = pid3;

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vertices[move_vtx2].attrs[t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}
 

		// 更新cost值
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.cluster_size[pid3]++;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
	}
	else if(cur_move.type == 3)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = cur_move.elem3;								// 点3
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点3要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = sol.ptn[move_vtx3];								// 点3所在的原分区，点2要移向的新分区

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx3][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vtx_attr_diff[move_vtx3][move_vtx2][t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}
 

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i]);
			}
			else if(i == move_vtx3)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;
		tabu_list[move_vtx3][pid3] = ts_iter + tabu_tenure;

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.ptn[move_vtx3] = pid1;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
		move_flag[move_vtx3]++;
	}

#if(DEBUG_VRF)
	printf("type=%d\n", min_c_gain.type);
	fflush(stdout);
	sol.verify_cluster_attrs(graph, "tabu更新之后");
	sol.verify_ifsb_change(graph, "tabu更新之后");
	sol.verify_infeasible(graph, "更新之后5");
#endif
	// 重置参数
	candidate_gain.resize(0);																		 
	tabu_candidate_gain.resize(0);
	double min_gain = MAX_VALUE;
	double tabu_min_gain = MAX_VALUE;
	min_c_gain.clear();																				// 重置最小的move gain

	/* 2.寻找最小的移动（insert），并更新相关参数*/
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		double quit_i = 0;																			// k1中移出点i引起的infeasible degree变化
		if ((pre_clst == pid1 || pre_clst == pid2 || pre_clst == pid3))		// 更新部分v_ifsb_change
		{
			for (int t = 0; t < graph.tau; t++)
			{
				if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
					quit_i = myRound(quit_i + sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
					quit_i = myRound(quit_i + graph.w_lower[t] - sol.cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
			}
			sol.ifsb_change[i][pre_clst] = myRound(quit_i - sol.cluster_infeasible[pre_clst]);						// 变化
		}

		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				double join_i = 0;																			// k2中加入点i引起的infeasible degree变化
				if ((next_clst == pid1 || next_clst == pid2 || next_clst == pid3))
				{
					for (int t = 0; t < graph.tau; t++)
					{
						if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
							join_i = myRound(join_i + sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
						else if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
							join_i = myRound(join_i + graph.w_lower[t] - sol.cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
					}
					sol.ifsb_change[i][next_clst] = myRound(join_i - sol.cluster_infeasible[next_clst]);
				}

				// 不会改变不可解程度才执行以下内容
				double ifsb_change = myRound(sol.ifsb_change[i][next_clst] + sol.ifsb_change[i][pre_clst]);
				if(abs(ifsb_change) < epsilon)
				{
					// 计算并更新move gain
					sol.dist_change[i][next_clst] = myRound(sol.dist_vtx_to_cluster[i][next_clst] - sol.dist_vtx_to_cluster[i][pre_clst]);
#if(USE_NBH0)
					if (tabu_list[i][next_clst] > ts_iter)  // 禁忌点
					{
						if (myRound(sol.dist_change[i][next_clst]) < myRound(tabu_min_gain) && sol.cluster_size[pre_clst] > 1)
						{
							Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
							tabu_min_gain = sol.dist_change[i][next_clst];
							tabu_candidate_gain.resize(1);
							tabu_candidate_gain[0] = new_node;
						}
						else if (myRound(sol.dist_change[i][next_clst]) == myRound(tabu_min_gain) && sol.cluster_size[pre_clst] > 1)
						{
							Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
							tabu_candidate_gain.push_back(new_node);
						}
					}
					else  // 非禁忌点
					{
						if (myRound(sol.dist_change[i][next_clst]) < myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
						{
							Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
							min_gain = sol.dist_change[i][next_clst];
							candidate_gain.resize(1);
							candidate_gain[0] = new_node;
						}
						else if (myRound(sol.dist_change[i][next_clst]) == myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
						{
							Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
							candidate_gain.push_back(new_node);
						}
					}
#endif
				}
			}
		}
	}

	/* 3.寻找最小邻域移动（exchange）*/
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double dist_change, ifsb_change;
				dist_change = myRound(sol.dist_vtx_to_cluster[i][move_clst2] + sol.dist_vtx_to_cluster[j][move_clst1] - sol.dist_vtx_to_cluster[i][move_clst1]  - sol.dist_vtx_to_cluster[j][move_clst2] -2*graph.edge[i][j]);

				// infeasible相关

				// 计算交换后的attr
				double ic1 = 0;						// exchange后clst1的infeasible
				double ic2 = 0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t]> graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t]> graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					if (tabu_list[i][move_clst2] > ts_iter && tabu_list[j][move_clst1] > ts_iter)  // 禁忌点
					{
						if (myRound(dist_change) < myRound(tabu_min_gain))
						{
							Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
							tabu_min_gain = dist_change;
							tabu_candidate_gain.resize(1);
							tabu_candidate_gain[0] = new_node;
						}
						else if (myRound(dist_change) == myRound(tabu_min_gain))
						{
							Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
							tabu_candidate_gain.push_back(new_node);
						}
					}
					else  // 非禁忌点
					{
						if (myRound(dist_change) < myRound(min_gain))
						{
							Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
							min_gain = dist_change;
							candidate_gain.resize(1);
							candidate_gain[0] = new_node;
						}
						else if (myRound(dist_change) == myRound(min_gain))
						{
							Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
							candidate_gain.push_back(new_node);
						}
					}
				}
			}
		}
	}
#endif

	/* 4.部分探索其他邻域 */
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];

#if(USE_NBH2)
		// 探索 push
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;

				double dist_change, ifsb_change;
				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				if(abs(ifsb_change) < epsilon)
				{
					// 计算cost变化
					dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] + sol.dist_vtx_to_cluster[nvtx2][npid]
								- sol.dist_vtx_to_cluster[nvtx1][npid1] - sol.dist_vtx_to_cluster[nvtx2][npid2]
								- graph.edge[nvtx1][nvtx2]);

					// 更新move gain结点
					if (tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid] > ts_iter)  // 禁忌点
					{
						if (myRound(dist_change) < myRound(tabu_min_gain))
						{
							Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
							tabu_min_gain = dist_change;
							tabu_candidate_gain.resize(1);
							tabu_candidate_gain[0] = new_node;
						}
						else if (myRound(dist_change) == myRound(tabu_min_gain))
						{
							Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
							tabu_candidate_gain.push_back(new_node);
						}
					}
					else  // 非禁忌点
					{
						if (myRound(dist_change) < myRound(min_gain))
						{
							Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
							min_gain = dist_change;
							candidate_gain.resize(1);
							candidate_gain[0] = new_node;
						}
						else if (myRound(dist_change) == myRound(min_gain))
						{
							Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
							candidate_gain.push_back(new_node);
						}
					}
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{
			double dist_change, ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid3]); 						// 如果是infeasible下降阶段就不是-0了


			if(abs(ifsb_change) < epsilon)
			{
				// 计算cost变化
				dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] - sol.dist_vtx_to_cluster[nvtx1][npid1]
							+ sol.dist_vtx_to_cluster[nvtx2][npid3] - sol.dist_vtx_to_cluster[nvtx2][npid2]
							+ sol.dist_vtx_to_cluster[nvtx3][npid1] - sol.dist_vtx_to_cluster[nvtx3][npid3]
							- graph.edge[nvtx1][nvtx2] - graph.edge[nvtx1][nvtx3] - graph.edge[nvtx2][nvtx3]);

				// 更新move gain结点
				if (tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid3] > ts_iter && tabu_list[nvtx3][npid1] > ts_iter)  // 禁忌点
				{
					if (myRound(dist_change) < myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
						tabu_min_gain = dist_change;
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else  // 非禁忌点
				{
					if (myRound(dist_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
						min_gain = dist_change;
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;
	idx_list = NULL;

	/* 4.选择邻域移动 */
	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	int ts_size_cl = tabu_candidate_gain.size();

	// 特赦准则：①非禁忌点中没有能改进的，禁忌点中有②禁忌点中的最大改进比best_cost要好（换成ts_best_cost也行）
	if (ts_size_cl > 0 && min_gain > 0 && tabu_min_gain < 0 && myRound(sol.dist + tabu_min_gain) != myRound(ts_best_cost))
	{
		int rnd = rand() % ts_size_cl;
		min_c_gain = tabu_candidate_gain[rnd];
	}
	else if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_c_gain = candidate_gain[rnd];
	}
	else
		min_c_gain.clear();

//	sol.verify_clst_attrs(graph, "禁忌邻域移动");
//	sol.verify_ifsb_change(graph, "禁忌邻域移动");
//	sol.verify_infeasible(graph, "禁忌邻域移动");											// DEBUG
//	sol.verify_dist(graph, "禁忌邻域移动");
}


/*
 * move gain是param_beta*ifsb_mg + cost_mg，SO中使用
 * 记录move_flag
 * 有SO的禁忌搜索move gain更新（SO）
 */
void update_tabu_move_gain3_1(Solution &sol, Graph &graph, Gain_node cur_move, int *move_flag)
{
	int move_vtx1;
	int move_vtx2;
	int move_vtx3;
	int pid1;
	int pid2;
	int pid3;

	//	verify_ifsb_change(sol, "update_c_move gain");

	/* 1.更新部分sol内容：*/
	if(cur_move.type == 0)
	{
		move_vtx1 = cur_move.elem1;
		move_vtx2 = -1;
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];			// 点1所在的原分区
		pid2 = cur_move.elem2;				// 点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vertices[move_vtx1].attrs[t]);

		}

		// 更新分区的不可解程度
		sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.ifsb_change[cur_move.elem1][pid1]);
		sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.ifsb_change[cur_move.elem1][pid2]);

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i != cur_move.elem1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[cur_move.elem1] = pid2;
		sol.cluster_size[pid2]++;
		move_flag[move_vtx1]++;
	}
	else if(cur_move.type == 1)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点2要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = -1;

		// 更新权重,两个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx2][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分
		}

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid1;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
	}
	else if(cur_move.type == 2)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = -1;
		pid1 = sol.ptn[move_vtx1];							// 点1所在的原分区
		pid2 = sol.ptn[move_vtx2];							// 点2所在的原分区， 点1要移向的新分区
		pid3 = cur_move.elem3;								// 点2要移向的新分区
//		move_to_clst1 = pid2;
//		move_to_clst2 = pid3;

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] - graph.vertices[move_vtx1].attrs[t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);;
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vertices[move_vtx2].attrs[t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}

		// 更新cost值
		sol.dist = myRound(sol.dist + cur_move.cost_mg);
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;

		// 移动
		sol.cluster_size[pid1]--;
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.cluster_size[pid3]++;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
	}
	else if(cur_move.type == 3)
	{
		move_vtx1 = cur_move.elem1;								// 点1
		move_vtx2 = cur_move.elem2;								// 点2
		move_vtx3 = cur_move.elem3;								// 点3
		pid1 = sol.ptn[move_vtx1];								// 点1所在的原分区，点3要移向的新分区
		pid2 = sol.ptn[move_vtx2];								// 点2所在的原分区，点1要移向的新分区
		pid3 = sol.ptn[move_vtx3];								// 点3所在的原分区，点2要移向的新分区

		// 更新权重，三个分区的不可解程度
		sol.cluster_infeasible[pid1] = 0.0;
		sol.cluster_infeasible[pid2] = 0.0;
		sol.cluster_infeasible[pid3] = 0.0;

		for (int t = 0; t < graph.tau; t++)
		{
			sol.cluster_attrs[pid1][t] = myRound(sol.cluster_attrs[pid1][t] + graph.vtx_attr_diff[move_vtx1][move_vtx3][t]);
			sol.cluster_attrs[pid2][t] = myRound(sol.cluster_attrs[pid2][t] + graph.vtx_attr_diff[move_vtx2][move_vtx1][t]);
			sol.cluster_attrs[pid3][t] = myRound(sol.cluster_attrs[pid3][t] + graph.vtx_attr_diff[move_vtx3][move_vtx2][t]);

			if (sol.cluster_attrs[pid1][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + sol.cluster_attrs[pid1][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid1][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid1] = myRound(sol.cluster_infeasible[pid1] + graph.w_lower[t] - sol.cluster_attrs[pid1][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid2][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + sol.cluster_attrs[pid2][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid2][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid2] = myRound(sol.cluster_infeasible[pid2] + graph.w_lower[t] - sol.cluster_attrs[pid2][t]);		// 下溢的部分

			if (sol.cluster_attrs[pid3][t] > graph.w_upper[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + sol.cluster_attrs[pid3][t] - graph.w_upper[t]);		// 上溢的部分
			else if (sol.cluster_attrs[pid3][t] < graph.w_lower[t])
				sol.cluster_infeasible[pid3] = myRound(sol.cluster_infeasible[pid3] + graph.w_lower[t] - sol.cluster_attrs[pid3][t]);		// 下溢的部分
		}
 

		// 更新cost值（infeasible必然是0，所以没有更新）
		sol.dist = myRound(sol.dist + cur_move.cost_mg);
		sol.ifsb = myRound(sol.ifsb + cur_move.ifsb_mg);

		// 更新dist_vtx_to_cluster
		for (int i = 0; i < graph.n; i++)
		{
			if (i == move_vtx1)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
			else if(i == move_vtx2)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i]);
			}
			else if(i == move_vtx3)
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] + graph.edge[move_vtx2][i]);
			}
			else
			{
				sol.dist_vtx_to_cluster[i][pid1] = myRound(sol.dist_vtx_to_cluster[i][pid1] - graph.edge[move_vtx1][i] + graph.edge[move_vtx3][i]);
				sol.dist_vtx_to_cluster[i][pid2] = myRound(sol.dist_vtx_to_cluster[i][pid2] - graph.edge[move_vtx2][i] + graph.edge[move_vtx1][i]);
				sol.dist_vtx_to_cluster[i][pid3] = myRound(sol.dist_vtx_to_cluster[i][pid3] - graph.edge[move_vtx3][i] + graph.edge[move_vtx2][i]);
			}
		}

		// 更新禁忌表
		tabu_list[move_vtx1][pid1] = ts_iter + tabu_tenure;
		tabu_list[move_vtx2][pid2] = ts_iter + tabu_tenure;
		tabu_list[move_vtx3][pid3] = ts_iter + tabu_tenure;

		// 移动
		sol.ptn[move_vtx1] = pid2;
		sol.ptn[move_vtx2] = pid3;
		sol.ptn[move_vtx3] = pid1;
		move_flag[move_vtx1]++;
		move_flag[move_vtx2]++;
		move_flag[move_vtx3]++;
	}

#if(DEBUG_cost_tsso)
	printf("type=%d ", cur_move.type);
	fflush(stdout);
	sol.verify_cluster_attrs(graph, "tabu更新之后");
	sol.verify_infeasible(graph, "更新之后4");
#endif

	// 重置参数
	candidate_gain.resize(0);																		 
	tabu_candidate_gain.resize(0);
	double min_gain = MAX_VALUE;
	double tabu_min_gain = MAX_VALUE;
	min_c_gain.clear();																				// 重置最小的move gain

	/* 2.寻找最小的移动（insert），并更新相关参数*/
	for (int i = 0; i < graph.n; i++)
	{
		int pre_clst = sol.ptn[i];

		double quit_i = 0;																			// k1中移出点i引起的infeasible degree变化
		if ((pre_clst == pid1 || pre_clst == pid2 || pre_clst == pid3))		// 更新部分v_ifsb_change
		{
			for (int t = 0; t < graph.tau; t++)
			{
				if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] > graph.w_upper[t])
					quit_i = myRound(quit_i + sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (sol.cluster_attrs[pre_clst][t] - graph.vertices[i].attrs[t] < graph.w_lower[t])
					quit_i = myRound(quit_i + graph.w_lower[t] - sol.cluster_attrs[pre_clst][t] + graph.vertices[i].attrs[t]);		// 下溢的部分
			}
			sol.ifsb_change[i][pre_clst] = myRound(quit_i - sol.cluster_infeasible[pre_clst]);						// 变化
		}

		for (int next_clst = 0; next_clst < graph.k; next_clst++)
		{
			if (pre_clst != next_clst)
			{
				double join_i = 0;																			// k2中加入点i引起的infeasible degree变化
				if ((next_clst == pid1 || next_clst == pid2 || next_clst == pid3))
				{
					for (int t = 0; t < graph.tau; t++)
					{
						if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] > graph.w_upper[t])
							join_i = myRound(join_i + sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] - graph.w_upper[t]);		// 上溢的部分
						else if (sol.cluster_attrs[next_clst][t] + graph.vertices[i].attrs[t] < graph.w_lower[t])
							join_i = myRound(join_i + graph.w_lower[t] - sol.cluster_attrs[next_clst][t] - graph.vertices[i].attrs[t]);		// 下溢的部分
					}
					sol.ifsb_change[i][next_clst] = myRound(join_i - sol.cluster_infeasible[next_clst]);
				}

				double ifsb_change = myRound(sol.ifsb_change[i][next_clst] + sol.ifsb_change[i][pre_clst]);

				// 计算并更新move gain
				sol.dist_change[i][next_clst] = myRound(sol.dist_vtx_to_cluster[i][next_clst] - sol.dist_vtx_to_cluster[i][pre_clst]);

#if(USE_NBH0)
				if (tabu_list[i][next_clst] > ts_iter)  // 禁忌点
				{
					if (myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) < myRound(tabu_min_gain) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						tabu_min_gain = myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change);
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if ((myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) == myRound(tabu_min_gain)) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else  // 非禁忌点
				{
					if (myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) < myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						min_gain = myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change);
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(sol.dist_change[i][next_clst] + param_beta*ifsb_change) == myRound(min_gain) && sol.cluster_size[pre_clst] > 1)
					{
						Gain_node new_node = Gain_node(i, next_clst, -1, ifsb_change, sol.dist_change[i][next_clst], 0);
						candidate_gain.push_back(new_node);
					}
				}
#endif
			}
		}
	}
#if(DEBUG_VRF)
	sol.verify_ifsb_change(graph, "So第一轮move gain更新完之后");
#endif

	/* 3.寻找最小邻域移动（exchange）*/
	double *clst1_attrs = new double[graph.tau];
	double *clst2_attrs = new double[graph.tau];
#if(USE_NBH1)
	for(int i = 0; i < graph.n; i++)
	{
		int move_clst1 = sol.ptn[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.ptn[i] != sol.ptn[j])
			{
				int move_clst2 = sol.ptn[j];
				double dist_change, ifsb_change;
				dist_change = myRound(sol.dist_vtx_to_cluster[i][move_clst2] + sol.dist_vtx_to_cluster[j][move_clst1] - sol.dist_vtx_to_cluster[i][move_clst1]  - sol.dist_vtx_to_cluster[j][move_clst2] -2*graph.edge[i][j]);

				// infeasible相关

				// 计算交换后的attr
				double ic1 = 0;						// exchange后clst1的infeasible
				double ic2 = 0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[move_clst1][t] + graph.vtx_attr_diff[i][j][t];
					clst2_attrs[t] = sol.cluster_attrs[move_clst2][t] + graph.vtx_attr_diff[j][i][t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 - sol.cluster_infeasible[move_clst1] - sol.cluster_infeasible[move_clst2]); 						// 如果是infeasible下降阶段就不是-0了

				if (tabu_list[i][move_clst2] > ts_iter && tabu_list[j][move_clst1] > ts_iter)  // 禁忌点
				{
					if (myRound(dist_change + param_beta*ifsb_change) < myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						tabu_min_gain = myRound(dist_change + param_beta*ifsb_change);
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change + param_beta*ifsb_change) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else  // 非禁忌点
				{
					if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						min_gain = myRound(dist_change + param_beta*ifsb_change);
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(i, j, -1, ifsb_change, dist_change, 1);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
	}
#endif

	/* 4.部分探索其他邻域 */
	// 寻找最小邻域移动（其他部分探索）
	int *idx_list = new int[graph.n];
	Generate_Rand_List(idx_list, graph.n);
	double *clst3_attrs = new double[graph.tau];
	for(int i = 0; i < graph.n; i++)
	{
		int nvtx1 = idx_list[i];
		int nvtx2 = idx_list[(i+1) % graph.n];
		int nvtx3 = idx_list[(i+2) % graph.n];
		int npid1 = sol.ptn[nvtx1];
		int npid2 = sol.ptn[nvtx2];
		int npid3 = sol.ptn[nvtx3];
#if(USE_NBH2)
		// 探索 push
		if(npid1 != npid2)
		{
			for (int npid = 0; npid < graph.k; npid++)
			{
				if(npid == npid1 ||  npid == npid2)
					continue;

				double dist_change, ifsb_change;
				/* 计算infeasible相关 */
				// 计算交换后的attr
				double ic1 = 0.0;						// exchange后clst1的infeasible
				double ic2 = 0.0;						// exchange后clst2的infeasible
				double ic3 = 0.0;						// exchange后clst2的infeasible
				for(int t = 0; t < graph.tau; t++)
				{
					clst1_attrs[t] = sol.cluster_attrs[npid1][t] - graph.vertices[nvtx1].attrs[t];
					clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
					clst3_attrs[t] = sol.cluster_attrs[npid][t] + graph.vertices[nvtx2].attrs[t];

					if (clst1_attrs[t] > graph.w_upper[t])
						ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst1_attrs[t] < graph.w_lower[t])
						ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

					if (clst2_attrs[t] > graph.w_upper[t])
						ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst2_attrs[t] < graph.w_lower[t])
						ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

					if (clst3_attrs[t] > graph.w_upper[t])
						ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
					else if (clst3_attrs[t] < graph.w_lower[t])
						ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
				}
				ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid]); 						// 如果是infeasible下降阶段就不是-0了

				// 计算cost变化
				dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] + sol.dist_vtx_to_cluster[nvtx2][npid]
							- sol.dist_vtx_to_cluster[nvtx1][npid1] - sol.dist_vtx_to_cluster[nvtx2][npid2]
							- graph.edge[nvtx1][nvtx2]);

				// 更新move gain结点
				if (tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid] > ts_iter)  // 禁忌点
				{
					if (myRound(dist_change + param_beta*ifsb_change) < myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						tabu_min_gain = myRound(dist_change + param_beta*ifsb_change);
						tabu_candidate_gain.resize(1);
						tabu_candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change + param_beta*ifsb_change) == myRound(tabu_min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						tabu_candidate_gain.push_back(new_node);
					}
				}
				else  // 非禁忌点
				{
					if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						min_gain = myRound(dist_change + param_beta*ifsb_change);
						candidate_gain.resize(1);
						candidate_gain[0] = new_node;
					}
					else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
					{
						Gain_node new_node = Gain_node(nvtx1, nvtx2, npid, ifsb_change, dist_change, 2);
						candidate_gain.push_back(new_node);
					}
				}
			}
		}
#endif
#if(USE_NBH3)
		// 探索 relocation
		if(npid1 != npid2 && npid1 != npid3 && npid2 != npid3)
		{
			double dist_change, ifsb_change;

			/* 计算infeasible相关 */
			// 计算交换后的attr
			double ic1 = 0.0;						// exchange后clst1的infeasible
			double ic2 = 0.0;						// exchange后clst2的infeasible
			double ic3 = 0.0;						// exchange后clst2的infeasible
			for(int t = 0; t < graph.tau; t++)
			{
				clst1_attrs[t] = sol.cluster_attrs[npid1][t] + graph.vtx_attr_diff[nvtx1][nvtx3][t];
				clst2_attrs[t] = sol.cluster_attrs[npid2][t] + graph.vtx_attr_diff[nvtx2][nvtx1][t];
				clst3_attrs[t] = sol.cluster_attrs[npid3][t] + graph.vtx_attr_diff[nvtx3][nvtx2][t];

				if (clst1_attrs[t] > graph.w_upper[t])
					ic1 = myRound(ic1 + clst1_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst1_attrs[t] < graph.w_lower[t])
					ic1 = myRound(ic1 + graph.w_lower[t] - clst1_attrs[t]);		// 下溢的部分

				if (clst2_attrs[t] > graph.w_upper[t])
					ic2 = myRound(ic2 + clst2_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst2_attrs[t] < graph.w_lower[t])
					ic2 = myRound(ic2 + graph.w_lower[t] - clst2_attrs[t]);		// 下溢的部分

				if (clst3_attrs[t] > graph.w_upper[t])
					ic3 = myRound(ic3 + clst3_attrs[t] - graph.w_upper[t]);		// 上溢的部分
				else if (clst3_attrs[t] < graph.w_lower[t])
					ic3 = myRound(ic3 + graph.w_lower[t] - clst3_attrs[t]);		// 下溢的部分
			}
			ifsb_change = myRound(ic1 + ic2 + ic3 - sol.cluster_infeasible[npid1] - sol.cluster_infeasible[npid2] - sol.cluster_infeasible[npid3]); 						// 如果是infeasible下降阶段就不是-0了

			// 计算cost变化
			dist_change = myRound(sol.dist_vtx_to_cluster[nvtx1][npid2] - sol.dist_vtx_to_cluster[nvtx1][npid1]
						+ sol.dist_vtx_to_cluster[nvtx2][npid3] - sol.dist_vtx_to_cluster[nvtx2][npid2]
						+ sol.dist_vtx_to_cluster[nvtx3][npid1] - sol.dist_vtx_to_cluster[nvtx3][npid3]
						- graph.edge[nvtx1][nvtx2] - graph.edge[nvtx1][nvtx3] - graph.edge[nvtx2][nvtx3]);

			// 更新move gain结点
			if (tabu_list[nvtx1][npid2] > ts_iter && tabu_list[nvtx2][npid3] > ts_iter && tabu_list[nvtx3][npid1] > ts_iter)  // 禁忌点
			{
				if (myRound(dist_change + param_beta*ifsb_change) < myRound(tabu_min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					tabu_min_gain = myRound(dist_change + param_beta*ifsb_change);
					tabu_candidate_gain.resize(1);
					tabu_candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change + param_beta*ifsb_change) == myRound(tabu_min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					tabu_candidate_gain.push_back(new_node);
				}
			}
			else  // 非禁忌点
			{
				if (myRound(dist_change + param_beta*ifsb_change) < myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					min_gain = myRound(dist_change + param_beta*ifsb_change);
					candidate_gain.resize(1);
					candidate_gain[0] = new_node;
				}
				else if (myRound(dist_change + param_beta*ifsb_change) == myRound(min_gain))
				{
					Gain_node new_node = Gain_node(nvtx1, nvtx2, nvtx3, ifsb_change, dist_change, 3);
					candidate_gain.push_back(new_node);
				}
			}
		}
#endif
	}

	// 释放辅助数组的内存
	delete[] idx_list;
	delete[] clst1_attrs;
	delete[] clst2_attrs;
	delete[] clst3_attrs;
	clst1_attrs = NULL;
	clst2_attrs = NULL;
	clst3_attrs = NULL;
	idx_list = NULL;

	/* 4.选择邻域移动 */
	// 从备选的邻域移动中随机选择一个
	int size_cl = candidate_gain.size();
	int ts_size_cl = tabu_candidate_gain.size();

	// 特赦准则：①非禁忌点中没有能改进的，禁忌点中有②禁忌点中的最大改进比best_cost要好（换成ts_best_cost也行）
	if (ts_size_cl > 0 && min_gain > 0 && tabu_min_gain < 0)
	{
		int rnd = rand() % ts_size_cl;
		min_c_gain = tabu_candidate_gain[rnd];
	}
	else if (size_cl > 0)
	{
		int rnd = rand() % size_cl;
		min_c_gain = candidate_gain[rnd];
	}
	else
		min_c_gain.clear();

//	sol.verify_clst_attrs(graph, "禁忌邻域移动");
//	sol.verify_ifsb_change(graph, "禁忌邻域移动");
//	sol.verify_infeasible(graph, "禁忌邻域移动");											// DEBUG
//	sol.verify_dist(graph, "禁忌邻域移动");
}

/*
 * 将解决方案解压缩为更原始一级的对等sol
 * graph1 → graph2（是层次更低的图）
 */
void uncoarsen_sol(Solution &sol, Graph &graph1, Graph &graph2)
{
	Solution sol2 = Solution(graph2);
	sol2.level = graph2.cur_level;
	sol2.ifsb = myRound(sol.ifsb);
	sol2.dist = myRound(sol.dist);									//

	for (int i = 0; i < graph1.n; i++)										// 每个折叠点
	{
		int cur_ptn = sol.ptn[i];
		sol2.cluster_size[cur_ptn] += graph1.vertices[i].pre_vts.size();
		int vts_size = graph1.vertices[i].pre_vts.size();
		for (int j = 0; j < vts_size; j++)			// 每个折叠点中的上一级点
		{
			int cur_vtx = graph1.vertices[i].pre_vts[j];
			sol2.ptn[cur_vtx] = cur_ptn;
		}
	}

	for (int i = 0; i < graph1.k; i++)
		memcpy(sol2.cluster_attrs[i], sol.cluster_attrs[i], sizeof(double)*graph1.tau);

	// 改变大小
	sol.free_memory(graph1);
	sol = Solution(sol2, graph2);
//	sol.cpy(sol2, graph2);
	sol2.free_memory(graph2);

	// 验证
//	sol.verify_cluster_attrs(graph2, "解压缩之后");
//	sol.verify_ifsb_change(graph2, "解压缩之后");
//	sol.verify_infeasible(graph2, "解压缩之后");											// DEBUG
//	if(abs(sol.ifsb) < epsilon)
//		sol.verify_dist(graph2, "解压缩之后");
}


/*
 * 生成一个infeasible为0的解
 */
void infeasible_tabu(Solution &cur_sol, Graph &graph)
{
	// 初始化
	for(int i = 0; i < graph.n; i++)
		memset(tabu_list[i], 0, sizeof(short int)*graph.k);

	// 计算目标函数，计算不可行程度
	init_i_move_gain(cur_sol, graph, min_i_gain);

	ts_best_ifsb = myRound(cur_sol.ifsb);
	ts_best_sol = Solution(cur_sol, graph);

	int non_improve = 0;
	int ts_iter = 0;
	while (1)
	{
		non_improve++;
		ts_iter++;

		// 选点和分区
		if (min_i_gain.elem1 == -1)
		{
			break;
		}
#if(DEBUG_VRF)
	printf("\n before——type=%d ifsb_change=%f cur_ifsb=%f", min_i_gain.type, min_i_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后2");
#endif
#if(DEBUG_VRF)
		if(min_i_gain.type == 0)
			printf("\n%d(%d)→(%d)\n", min_i_gain.elem1, cur_sol.ptn[min_i_gain.elem1], min_i_gain.elem2);
		else if(min_i_gain.type == 1)
			printf("\n%d(%d)→(%d)， %d(%d)→(%d)\n", min_i_gain.elem1, cur_sol.ptn[min_i_gain.elem1], cur_sol.ptn[min_i_gain.elem2], min_i_gain.elem2, cur_sol.ptn[min_i_gain.elem2], cur_sol.ptn[min_i_gain.elem1]);
		fflush(stdout);
#endif
		update_i_move_gain_ts(cur_sol, graph, min_i_gain);
#if(DEBUG_VRF)
		if(min_i_gain.type == 0)
			printf("\n%d(%d)→(%d)\n", min_i_gain.elem1, cur_sol.ptn[min_i_gain.elem1], min_i_gain.elem2);
		else if(min_i_gain.type == 1)
			printf("\n%d(%d)→(%d)， %d(%d)→(%d)\n", min_i_gain.elem1, cur_sol.ptn[min_i_gain.elem1], cur_sol.ptn[min_i_gain.elem2], min_i_gain.elem2, cur_sol.ptn[min_i_gain.elem2], cur_sol.ptn[min_i_gain.elem1]);
		fflush(stdout);
#endif
#if(DEBUG_VRF)
	printf("\n after——type=%d ifsb_change=%f cur_ifsb=%f", min_i_gain.type, min_i_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后1");
#endif
	// 必更新
		if (myRound(cur_sol.ifsb) < myRound(ts_best_ifsb))
		{
			non_improve = 0;
			ts_best_sol.cpy(cur_sol, graph);
			ts_best_ifsb = myRound(cur_sol.ifsb);					 
		}
#if(DEBUG_PRT_ifsb_ts)
		printf("iter=%d, cur_ifsb=%f\n", ts_iter, cur_sol.ifsb);
		fflush(stdout);
#endif

		// 判断是否继续（去掉这里程序仍然可以停止，但是加上可以减少搜索时间）
		if (abs(ts_best_ifsb) < epsilon || non_improve >= tabu_depth)
		{
			break;
		}
	}

	cur_sol.cpy(ts_best_sol, graph);

	// 作用完成，释放
	ts_best_sol.free_memory(graph);
	// 验证
#if(DEBUG_VRF)
	cur_sol.verify_ifsb_change(graph, "infeasible_tabu之后");
	cur_sol.verify_cluster_attrs(graph, "infeasible_tabu之后");
	cur_sol.verify_infeasible(graph, "infeasible_tabu之后");											// DEBUG
#endif
//	if(abs(sol.ifsb) < epsilon)
//		sol.verify_dist(graph2, "解压缩之后");
}


/*
 * 禁忌搜索改进cost
 * 进入该函数前需要保证：空的dist_vtx_to_cluster
 */
void cost_descent(Solution &cur_sol, Graph &graph)
{
	// 初始化
//	printf("\npre=%.9f", cur_sol.dist);
	init_c_move_gain(cur_sol, graph, min_c_gain);													// 初始化cost相关
//	printf(", init=%.9f", cur_sol.dist);
	// 改进
	bool improved = true;
	while (improved)
	{
		improved = false;

		// 选点和分区
		if (myRound(min_c_gain.ifsb_mg) < 0)
		{
			update_c_move_gain2(cur_sol, graph, min_c_gain);
			improved = true;
		}
	}
//	printf(", end=%.9f\n", cur_sol.dist);

	// 验证
//	cur_sol.verify_cluster_attrs(graph, "cost descent之后");
//	cur_sol.verify_infeasible(graph,  "cost descent之后");							// DEBUG
//	cur_sol.verify_dist(graph,  "cost descent之后");
}


/*
 * 禁忌搜索改进cost（根据学习矩阵扰动）
 * 进入该函数前需要保证：空的dist_vtx_to_cluster
 */
void cost_tabu(Solution &cur_sol, Graph &graph, int *move_flag)
{
#if(DEBUGx)
	printf(" <tabu开始");
	fflush(stdout);
#endif
	// 初始化
	for(int i = 0; i < graph.n; i++)
		memset(tabu_list[i], 0, sizeof(short int)*graph.k);
	init_c_move_gain(cur_sol, graph, min_c_gain);													// 初始化cost相关

	if(abs(cur_sol.ifsb) < epsilon)
	{
		ts_best_cost = cur_sol.dist;
		ts_best_sol = Solution(cur_sol, graph);
	}
	else
	{
		ts_best_cost = MAX_VALUE;
		ts_best_sol = Solution(graph);
	}

	if(myRound(ts_best_cost) < myRound(rbest_cost) && abs(ts_best_sol.ifsb) < epsilon)
	{
		rbest_cost = myRound(ts_best_cost);
		rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
		rcd_cost.push_back(rbest_cost);
		rcd_time.push_back(rbest_time);
	}

	// 改进
	ts_iter = 0;
	int non_improve = 0;
	while (1)
	{
		non_improve++;

		// 更新，选下一个点和分区
		if (min_c_gain.elem1 == -1)
			break;
#if(DEBUG_cost_ts)
	printf("\n before——type=%d ifsb_change=%f cur_ifsb=%f", min_c_gain.type, min_c_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后00");
#endif
		update_tabu_move_gain2(cur_sol, graph, min_c_gain, move_flag);
#if(DEBUG_cost_ts)
	printf("\n after——type=%d ifsb_change=%f cur_ifsb=%f", min_c_gain.type, min_c_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后11");
#endif
//		update_tabu_move_gain(cur_sol, graph, min_c_gain);

		// 验证
//		cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
//		cur_sol.verify_ifsb_change(graph, "tabu更新之后");
//		cur_sol.verify_infeasible(graph, "tabu更新之后");											// DEBUG
//		cur_sol.verify_dist(graph, "tabu更新之后");

		if (myRound(cur_sol.dist) < myRound(ts_best_cost))
		{
			non_improve = 0;
			Solution tsol = Solution(cur_sol, graph);			// 复制cur_sol
			cost_descent(tsol, graph);
			ts_best_sol.cpy(tsol, graph);
			ts_best_cost = myRound(tsol.dist);					 
			tsol.free_memory(graph);

			if(myRound(ts_best_cost) < myRound(rbest_cost))
			{
				rbest_cost = myRound(ts_best_cost);
				rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
				rcd_cost.push_back(rbest_cost);
				rcd_time.push_back(rbest_time);
			}
		}

#if(DEBUG_PRT1)
		printf("cost tabu --- iter=%d, cur cost=%f, ts best cost=%f, ils best cost=%.6f\n", ts_iter, cur_sol.dist, ts_best_cost, ils_best_cost);
		fflush(stdout);
#endif
		if(non_improve >= tabu_depth)
			break;

		ts_iter++;
	}

	cur_sol.cpy(ts_best_sol, graph);

	// 作用完成，释放
	ts_best_sol.free_memory(graph);

	// 验证
//	cur_sol.verify_cluster_attrs(graph, "tabu中while循环之后");
//	cur_sol.verify_infeasible(graph,  "tabu中while循环之后");							// DEBUG
//	cur_sol.verify_dist(graph,  "tabu中while循环之后");
}


/*
 * 禁忌搜索改进cost
 */
void cost_tabu_SO(Solution &cur_sol, Graph &graph, int *move_flag)
{
#if(DEBUGx)
	printf(" 【SO开始");
	fflush(stdout);
#endif
#if(DEBUG_VRF)
	// 验证
	cur_sol.verify_ifsb_change(graph, "cost tabu之前");
	cur_sol.verify_cluster_attrs(graph, "cost tabu之前");
	cur_sol.verify_infeasible(graph, "cost tabu之前");
#endif
	int rcd[5];
	memset(rcd, 0, sizeof(int)*5);
	int sum_rcd = 0;
	// 初始化
	param_beta = 1.0;
	for(int i = 0; i < graph.n; i++)
		memset(tabu_list[i], 0, sizeof(short int)*graph.k);
	init_i_move_gain(cur_sol, graph, min_i_gain);
	init_c_move_gain_so(cur_sol, graph, min_c_gain);													// 初始化cost相关

	if(abs(cur_sol.ifsb) < epsilon)
	{
		ts_best_cost = cur_sol.dist;
		ts_best_sol = Solution(cur_sol, graph);
	}
	else
	{
		ts_best_cost = MAX_VALUE;
		ts_best_sol = Solution(graph);
	}
	if(myRound(ts_best_cost) < myRound(rbest_cost) && abs(ts_best_sol.ifsb) < epsilon)
	{
		rbest_cost = myRound(ts_best_cost);
		rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
		rcd_cost.push_back(rbest_cost);
		rcd_time.push_back(rbest_time);
	}

	// 改进
	ts_iter = 0;
	int non_improve = 0;
	while (1)
	{
		non_improve++;

		// 更新，选下一个点和分区
		if (min_c_gain.elem1 == -1)
			break;
#if(DEBUG_cost_tsso)
	printf("\n before——type=%d ifsb_change=%f cur_ifsb=%f", min_c_gain.type, min_c_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后11");
#endif
		update_tabu_move_gain3_1(cur_sol, graph, min_c_gain, move_flag);
#if(DEBUG_cost_tsso)
	printf("\n after——type=%d ifsb_change=%f cur_ifsb=%f", min_c_gain.type, min_c_gain.ifsb_mg, cur_sol.ifsb);
	fflush(stdout);
	cur_sol.verify_cluster_attrs(graph, "tabu更新之后");
	cur_sol.verify_ifsb_change(graph, "tabu更新之后");
	cur_sol.verify_infeasible(graph, "更新之后22");
#endif

		// 自适应更新param_beta
		sum_rcd += abs(cur_sol.ifsb) < epsilon ? 1:0 - rcd[ts_iter % 5];
		if(sum_rcd == 5 && ts_iter > 5)
			param_beta /= param_brate;
		else if(sum_rcd == 0 && ts_iter > 5)
			param_beta *= param_brate;

		if(param_beta < 1)
			param_beta = 1;

		rcd[ts_iter % 5] = abs(cur_sol.ifsb) < epsilon ? 1:0;

		// 更新
		if (myRound(cur_sol.dist) < myRound(ts_best_cost) && abs(cur_sol.ifsb) < epsilon)
		{

			non_improve = 0;
			Solution tsol = Solution(cur_sol, graph);			// 复制cur_sol
			cost_descent(tsol, graph);
			ts_best_sol.cpy(tsol, graph);
			ts_best_cost = myRound(tsol.dist);					 
			tsol.free_memory(graph);
			if(myRound(ts_best_cost) < myRound(rbest_cost))
			{
				rbest_cost = myRound(ts_best_cost);
				rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
				rcd_cost.push_back(rbest_cost);
				rcd_time.push_back(rbest_time);
			}
		}

#if(DEBUG_PRT1)
		printf("cost tabu SO --- iter=%d, cur cost=%.2f, ts best cost=%.2f, ils best cost=%.2f\n", ts_iter, cur_sol.dist, ts_best_cost, ils_best_cost);
		fflush(stdout);
#endif
		if(non_improve >= tabu_depth)
			break;

		ts_iter++;
	}

	cur_sol.cpy(ts_best_sol, graph);

	// 作用完成，释放
	ts_best_sol.free_memory(graph);

#if(DEBUGx)
	printf(" tabu结束>");
	fflush(stdout);
#endif
#if(DEBUGx)
	printf("，SO结束】");
	fflush(stdout);
#endif
	// 验证
#if(DEBUG_cost_tsso)
	cur_sol.verify_cluster_attrs(graph, "SO之后");
	cur_sol.verify_infeasible(graph,  "SO之后");							// DEBUG
	cur_sol.verify_dist(graph,  "SO之后");
#endif
}


/*
 * 计算常用指标并输出：平均cost，平均时间，hit次数，cost的标准差
 */
void cal_indicators()
{
	double sum_cost = 0, sum_time = 0;
	hit = 0;
	std_dev = 0;
	for (int i = 0; i < runs; i++)
	{
		sum_cost += each_run_rlt[i];
		sum_time += each_run_time[i];
//		sum_htime += each_hit_time[i];

		if (abs(myRound(each_run_rlt[i], 2) - myRound(glb_best_cost, 2)) < 0.01)
			hit++;
	}

	avg_cost = sum_cost / runs;
	avg_time = sum_time / runs;
//	avg_htime = sum_htime / runs;

	sum_avg_cost += avg_cost;
	sum_avg_time += avg_time;

	for (int i = 0; i < runs; i++)
		std_dev += pow(each_run_rlt[i] - avg_cost, 2) / runs;

	std_dev = sqrt(std_dev);

	// 输出
//	cout << "best cost: " << glb_best_cost << "\nhit:" << hit << "\navg_cost: " << avg_cost << "\navg_time: " << avg_time << "\nstd_dev: " << std_dev << endl;
}


/*
 * 要保证cluster_attrs、cluster_size、cur_sol是对的
 */
void perturb_update(Solution &cur_sol, int vtx, int clst_move_to, Graph &graph)
{
	int clst_move_from = cur_sol.ptn[vtx];

	// 更新权重,两个分区的不可解程度
	for (int t = 0; t < graph.tau; t++)
	{
		cur_sol.cluster_attrs[clst_move_from][t] = myRound(cur_sol.cluster_attrs[clst_move_from][t] - graph.vertices[vtx].attrs[t]);
		cur_sol.cluster_attrs[clst_move_to][t] = myRound(cur_sol.cluster_attrs[clst_move_to][t] + graph.vertices[vtx].attrs[t]);
	}

	// 移动
	cur_sol.cluster_size[clst_move_from]--;
	cur_sol.ptn[vtx] = clst_move_to;
	cur_sol.cluster_size[clst_move_to]++;
}


/*
 * 施加随机扰动，只更新sol和部分相关变量
 */
bool perturb(double np_pct, Solution &sol, Graph &graph, int *move_flag)
{
	int num_ptb = ceil(np_pct * graph.n);
	bool cur_flag = false;							// cur_flag=true说明本次扰动成功了
	bool once_flag = false;							// once_flag=true说明num_ptb次扰动中有一次成功了

	int *idxlist = new int[graph.n];
	for(int i = 0; i < graph.n; i++)
		idxlist[i] = i;
	Quick_Sort_up(idxlist, move_flag, 0, graph.n-1);
	for (int i = 0; i < num_ptb; i++)
	{
		cur_flag = false;
		int rnd_vtx = idxlist[i];
		for (int j = 0; j < graph.n; j++)
		{
			if (sol.cluster_size[sol.ptn[(rnd_vtx + j) % graph.n]] > 1)	// 保证每个cluster中至少有一个点
			{
				rnd_vtx = (rnd_vtx + j) % graph.n;
				cur_flag = true;
				once_flag = true;
				break;
			}
		}

		if (cur_flag)
		{
			int rnd_clst = rand() % graph.k;
			if (sol.ptn[rnd_vtx] == rnd_clst)						// 如果随机cluster是随机点所在的cluster，简单处理一下
			{
				int rnd = rand() % (graph.k - 2) + 1;					// 选择加一个随机数字（不可以是0，也不可以是k）
				rnd_clst = (rnd_clst + rnd) % graph.k;
			}

			// 移动（只更新了部分）
			perturb_update(sol, rnd_vtx, rnd_clst, graph);
		}
		else
			break;													// 不可能再有扰动的组合了，不用继续扰动下去了
	}

	delete[] idxlist;
	return once_flag;
}


/*
 * 迭代地构造初始解并进行局部搜索改进（这个函数内的每一个解决方案都是同一level的）
 */
void itered_local_search(Solution &sol, Graph &graph)
{
#if(DEBUGx)
	printf("{ILS开始：");
	fflush(stdout);
#endif
	// 重置ils_best相关
	if(abs(sol.ifsb) < epsilon)			// 如果初始解sol是一个可行解
	{
		ils_best_cost = myRound(sol.dist);
		ils_best_sol = Solution(sol, graph);
	}
	else								// 如果初始解是一个不可行解
	{
		ils_best_sol = Solution(sol, graph);
		ils_best_cost = MAX_VALUE;
		ils_best_sol.ifsb = MAX_VALUE;
	}

	//迭代改进
//	int iter = 0;
	int count_frozen = 0;
	bool is_perturb = true;
	int *move_flag = new int[graph.n];														// 记录每个点的移动次数
	memset(move_flag, 0, sizeof(int)*graph.n);

	bool improve = true;
	while (improve && (clock() - start) / (double)CLOCKS_PER_SEC < tmax)
	{
		if (graph.k >= graph.n)			// 如果分区数和点数相同就没有改进的必要了
			break;

		if(1)
		{
			improve = false;

			count_frozen = -1;

			// cost_tabu(sol, graph, move_flag);						// 改进目标函数
			cost_tabu_SO(sol, graph, move_flag);						// 改进目标函数
			if(abs(sol.ifsb) < epsilon)
				cost_tabu(sol, graph, move_flag);	


			if (myRound(sol.dist) < myRound(ils_best_cost) && abs(sol.ifsb) < epsilon)
			{
				ils_best_cost = myRound(sol.dist);
				ils_best_sol.cpy(sol, graph);			// sol中保存ts搜索过程中最好的解

				improve = true;

				if(myRound(sol.dist) < myRound(rbest_cost))
				{
					rbest_cost = myRound(sol.dist);
					rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
					rcd_cost.push_back(rbest_cost);
					rcd_time.push_back(rbest_time);
				}

//				printf("ils_best_cost = %.4f, rbest_cost=")
			}
		}

#if(DEBUG_PRT2)
		printf("itered_local_search --- iter=%d, cur cost=%.2f, ils best cost=%.2f\n", iter, sol.dist, ils_best_cost);
		fflush(stdout);
#endif

		// 记录infeasible最小的那个解
		if(count_frozen != -1 && sol.ifsb < ils_best_sol.ifsb)
		{
			ils_best_sol.cpy(sol, graph);
		}
		count_frozen++;

		// 补充扰动逻辑，如果improve = false，不用扰动
		if(improve)
		{
			is_perturb = perturb(np_pct, sol, graph, move_flag);			// 施加扰动
			if (is_perturb == false)						// 如果不能扰动，提前停止迭代
				break;
		}

		// 如果连续max_frozen次没有改进，就退出（说明当前压缩图大概率是没有可行解的）
		if(count_frozen == max_frozen)
		{
			nexit++;
			break;
		}
	}

	// 每次搜索完都要把sol更新为该搜索过程中最好的那个解
	sol.cpy(ils_best_sol, graph);

	// 释放ils_best_sol的内存（ils_best_sol的作用域只在这一个函数中）
	ils_best_sol.free_memory(graph);
	delete[] move_flag;

	// 验证
#if(DEBUG_ils)
	sol.verify_cluster_attrs(graph, "ils之后");
	sol.verify_ifsb_change(graph, "ils之后");
	sol.verify_infeasible(graph, "ils之后");											// DEBUG
	if(abs(sol.ifsb) < epsilon)
		sol.verify_dist(graph, "ils之后");
#endif
#if(DEBUGx)
	printf(" ILS结束}\n");
	fflush(stdout);
#endif
}


/*
 * 多级优化
 */
Solution multi_level(Graph &graph)
{
	/* 1.coarsening phase */
	ml_best_cost = MAX_VALUE;
	coarsen_process(graph);

	/* 2.Initial partition and its refinement */
	Solution sol = Solution(coarsen_graph[cur_lvl]);
//	sol.init_sol_greedy(coarsen_graph[cur_lvl]);
	sol.init_sol(coarsen_graph[cur_lvl]);

	itered_local_search(sol, coarsen_graph[cur_lvl]);
	if (cur_lvl > 0)
	{
		uncoarsen_sol(sol, coarsen_graph[cur_lvl], coarsen_graph[cur_lvl - 1]);
		cur_lvl--;
	}

	/* 3.uncoarsening Phase */
	while (cur_lvl > 0)
	{
		// cout << "Multi-level:" << cur_lvl;
		if(cur_lvl > 1)
			itered_local_search(sol, coarsen_graph[cur_lvl]);
		if (cur_lvl > 0)
			uncoarsen_sol(sol, coarsen_graph[cur_lvl], coarsen_graph[cur_lvl - 1]);
		cur_lvl--;
	}
	if(cur_lvl == 0)
		itered_local_search(sol, coarsen_graph[cur_lvl]);
	//itered_local_search(sol, coarsen_graph[cur_lvl]);

	// 每一层都跑完之后，复制ils的最好结果（sol）
	if (myRound(sol.dist) < myRound(ml_best_cost) && abs(sol.ifsb) < epsilon)		// 只接受可行解
	{
		ml_best_cost = myRound(sol.dist);
		ml_best_sol.cpy(sol, graph);

		if(myRound(sol.dist) < myRound(rbest_cost))
		{
			rbest_cost = myRound(sol.dist);
			rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
			rcd_cost.push_back(rbest_cost);
			rcd_time.push_back(rbest_time);
		}
	}

	// 如果multi-level没能找到一个可行解，直接对原始图进行ils
	if(abs(myRound(sol.ifsb)) > epsilon)
	{
		printf("由于未找到可行解，进入备选搜索阶段\n");
		fflush(stdout);
		Solution rst_sol = Solution(coarsen_graph[0]);
//		rst_sol.init_sol_greedy(coarsen_graph[0]);
		rst_sol.init_sol(coarsen_graph[0]);
		itered_local_search(rst_sol, coarsen_graph[0]);
		ml_best_sol.cpy(rst_sol, graph);
		ml_best_cost = myRound(rst_sol.dist);
		if(myRound(ml_best_cost) < myRound(rbest_cost) && abs(ml_best_sol.ifsb) < epsilon)
		{
			rbest_cost = myRound(ml_best_cost);
			rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
			rcd_cost.push_back(rbest_cost);
			rcd_time.push_back(rbest_time);
		}
		rst_sol.free_memory(coarsen_graph[0]);
	}


	// 剩余的时间用来跑原始图
	while((clock() - start) / static_cast<double>(CLOCKS_PER_SEC) < tmax)
	{
		itered_local_search(sol, coarsen_graph[0]);
		if (myRound(sol.dist) < myRound(ml_best_cost)  && abs(sol.ifsb) < epsilon)		// 只接受可行解
		{
			ml_best_cost = myRound(sol.dist);
			ml_best_sol.cpy(sol, graph);

			if(myRound(sol.dist) < myRound(rbest_cost))
			{
				rbest_cost = myRound(sol.dist);
				rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
				rcd_cost.push_back(rbest_cost);
				rcd_time.push_back(rbest_time);
			}
		}
	}
	// 验证
	//ml_best_sol.verify_cluster_attrs(graph, "multi-level之后");
	//ml_best_sol.verify_ifsb_change(graph, "multi-level之后");
//	ml_best_sol.verify_infeasible(graph, "multi-level之后");											// DEBUG
//	if(myRound(ml_best_sol.ifsb) == 0)
//		ml_best_sol.verify_dist(graph, "multi-level之后");


//	test_print("multi-level结尾");
	return sol;
}


/*
 * 输出文件中的标题信息 MDFD
 */
void output_header(string dir, string rltfile)
{
	FILE *writefile;

	// 记录每次运行的结果

	// 记录结果概要
//	string outfile_overall = dir + rltfile + "_overall" + ".txt";
//	writefile = fopen(outfile_overall.c_str(), "a+");
	char outfile_overall[1000];
	sprintf(outfile_overall, "./results/overall/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_overall, "a+");
 
	fprintf(writefile, "\n\n------------------------- test_MLWT_v2.2.2.1.2.h ---------------------------------------------");
	fprintf(writefile, "\nILS_iter=%d \tnp_pct=%f \tmax_cv_count=%d \tmax_lvl=%d \ttt_pct=%.2f \ttd_pct=%.2f \tnp_pct=%.2f \tparam_beta=%.2f \tparam_brate=%.2f", max_iter, np_pct, max_cv_count, max_lvl, tt_pct, td_pct, np_pct, param_beta, param_brate);
	fprintf(writefile, "\ninstance \t\t\t\thit compare  best_cost \t trg_cost \t avg_cost \t avg_time \t hit_time \t std_dev \t std_dev%%");
	fprintf(writefile, "\n------------------------------------------------------------------------------------------\n");

	fclose(writefile);
}


void output_header(char *filename, string rltfile)
{
	FILE *writefile;

	// 记录每次运行的结果
	// 记录结果概要
//	string outfile_overall = rltfile + "_overall" + ".txt";
//	writefile = fopen(outfile_overall.c_str(), "a+");
	char outfile_overall[1000];
	sprintf(outfile_overall, "./results/overall/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_overall, "a+");
 
	fprintf(writefile, "\n\n------------------------- test_MLWT_v2.2.2.1.2.h ---------------------------------------------");
	fprintf(writefile, "\nILS_iter=%d \tnp_pct=%f \tmax_cv_count=%d \tmax_lvl=%d \ttt_pct=%.2f \ttd_pct=%.2f \tnp_pct=%.2f \tparam_beta=%.2f \tparam_brate=%.2f", max_iter, np_pct, max_cv_count, max_lvl, tt_pct, td_pct, np_pct, param_beta, param_brate);
	fprintf(writefile, "\ninstance \t\t\t\thit compare \t best_cost \t trg_cost \t avg_cost \t avg_time \t hit_time \t std_dev \t std_dev%%");
	fprintf(writefile, "\n------------------------------------------------------------------------------------------\n");

	fclose(writefile);
}


void output_tail(string dir, string rltfile)
{
	FILE *writefile;



	// 记录结果概要
	char outfile_overall[1000];
	sprintf(outfile_overall, "./results/overall/%s_%d.txt", graph_name.c_str(), seed);
//	string outfile_overall = dir + rltfile + "_overall" + ".txt";
	writefile = fopen(outfile_overall, "a+");
 
	fprintf(writefile, "-------------------------------------------------------------------------------------------");
	fprintf(writefile, "\nsum_avg_cost=%.2f\nsum_avg_time=%.2f\n", sum_avg_cost, sum_avg_time);

	fclose(writefile);
}

/*
 * 对partition的编号重新排序
 */
void sort_pid(short int *ptn, Graph &graph)
{
	// 获取编号对应
	int size = 0;
	int *projection = new int[graph.k];
	memset(projection, -1, sizeof(int)*graph.k);
	for(int i = 0; i < graph.n; i++)
	{
		int id = ptn[i];
		if(projection[id] == -1)
		{
			projection[id] = size;
			size++;
			if(size == graph.k)
				break;
		}
	}

	// 根据新编号与旧编号的映射重置ptn
	for(int i = 0; i < graph.n; i++)
		ptn[i] = projection[ptn[i]];

	delete[] projection;
}


/*
 * 对某个文件filename执行ILS算法
 */
void execute_instance(string dir, string filename, string rltfile)
{
	string path = dir + filename + ".txt";

	// 读图
	Graph graph;
	graph.read_file2(path);
	char *ftemp = basename(const_cast<char*>(filename.c_str()));
	graph_name = string(ftemp);
//	delete[] ftemp;
	//	print_graph();

	// 初始化
	allocate_memery(graph);															// 分配内存

	glb_best_cost = MAX_VALUE;
	for (int i = 0; i < runs; i++)
	{
		is_hit = false;
		rbest_time = 99999999;
		rbest_cost = 99999999;

		start = clock();

		Solution best_sol = multi_level(graph);										// best_sol没有用，要删掉

		if (myRound(ml_best_cost) < myRound(glb_best_cost) && abs(ml_best_sol.ifsb) < epsilon)
		{
			glb_best_cost = myRound(ml_best_cost);
			glb_best_sol.cpy(ml_best_sol, graph);

			if(myRound(ml_best_cost) < myRound(rbest_cost))
			{
				rbest_cost = myRound(ml_best_cost);
				rbest_time = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
				rcd_cost.push_back(rbest_cost);
				rcd_time.push_back(rbest_time);
			}
		}

		// 更新记录
//		each_run_time[i] = (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
		each_run_time[i] = rbest_time;
		each_run_rlt[i] = ml_best_cost;
//		each_hit_time[i] = hit_time;
		each_best_sol[i] = Solution(ml_best_sol, graph);

		printf("Round %d: best=%.4f, rbest=%.4f, all time=%.2f||", i+1, ml_best_cost, rbest_cost, each_run_time[i]);
		ml_best_sol.verify_infeasible(graph, "一轮run之后");											// DEBUG
		ml_best_sol.verify_dist(graph, "一轮run之后");

		fflush(stdout);
		best_sol.free_memory(graph);
	}
//	print_record();
	// 计算各种指标
	cal_indicators();

	FILE *writefile;

	// 记录每次运行的结果
//	string outfile_detail = dir + rltfile + "_detail" + ".txt";
//	writefile = fopen(outfile_detail.c_str(), "a+");

	// 记录结果概要
	char outfile_detail[1000];
	sprintf(outfile_detail, "./results/detail/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_detail, "a+");

	fprintf(writefile, "===========%s\n", filename.c_str());
	for (int i = 0; i < runs; i++)
	{
		fprintf(writefile, "\t%d--", i);
		fprintf(writefile, "\t%.4f \t%.4f \t%.4f", each_run_rlt[i], each_best_sol[i].ifsb, each_run_time[i]);

		fprintf(writefile, "\tpartitioning[] = ");
		sort_pid(each_best_sol[i].ptn, graph);
		for(int j = 0; j < graph.n; j++)
		{
			fprintf(writefile, "%d ", each_best_sol[i].ptn[j]);
		}
		fprintf(writefile, "\n");
	}
	fprintf(writefile, "\n");
	fclose(writefile);

	// 记录结果概要
//	string outfile_overall = dir + rltfile + "_overall" + ".txt";
//	writefile = fopen(outfile_overall.c_str(), "a+");

	char outfile_overall[1000];
	sprintf(outfile_overall, "./results/overall/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_overall, "a+");
 

	// 与target最好对比
	string cmp = " ";

	fprintf(writefile, "%s\t", filename.c_str());
	fprintf(writefile, "%d/%d\t", hit, runs);
	fprintf(writefile, "graph_name:%s \nbest_cost:%.4f \navg_cost:%.4f \navg_time:%.6f\n",
			cmp.c_str(),
			myRound(glb_best_cost, 4),
			myRound(avg_cost, 4),
			rbest_time);
	fclose(writefile);

//	printf("Best %f\n", glb_best_cost);
	// 释放内存
	free_memory(graph);
}


/*
 * 读算例列表
 */
void execute_each_instance(string dir, string listfile, string rltfile)
{
	string path = dir + listfile;
	string *instance_list;
	double *target_cost_list;
	int num_instances;

	ifstream fin(path);		// 打开文件

	// 判断是否成功打开文件
	if (!fin.is_open())
	{
		cerr << "Can not open the file!" << path << endl;
		exit(-14);
	}

	// 检查文件流是否处于错误状态
	if (fin.fail())
	{
		cerr << "Error occurred during file operation." << path << endl;
		exit(-15);
	}

	// 判断文件是否为空
	if (fin.eof())
	{
		cerr << "Empty file" << path << endl;
		exit(-16);
	}

	// 开始读内容
	fin >> num_instances;
	instance_list = new string[num_instances];  // 分配内存，并调整大小
	target_cost_list = new double[num_instances];

	int temp = 0;
	while (!fin.eof())
	{
		fin >> instance_list[temp];
		fin >> target_cost_list[temp++];
	}

	fin.close();

	// 遍历测试所有算例
	for (int i = 0; i < num_instances; i++)
	{
		target_cost = target_cost_list[i];
		execute_instance(dir, instance_list[i], rltfile);
	}


	output_tail(dir, rltfile);

	// 释放内存
	delete[] target_cost_list;
	delete[] instance_list;
}


void execute_instance(char * filename, string rltfile)
{
	// 读图
	Graph graph;
	graph.read_file2(filename);
	char *ftemp = basename(filename);
	graph_name = string(ftemp);
	//	print_graph();

	// 初始化
	allocate_memery(graph);															// 分配内存

	glb_best_cost = MAX_VALUE;
	for (int i = 0; i < runs; i++)
	{
		start = clock();
		rbest_time = 99999999;
		rbest_cost = 99999999;
		Solution best_sol = multi_level(graph);										// best_sol没有用，要删掉

		if (myRound(ml_best_cost) < myRound(glb_best_cost) && abs(ml_best_sol.ifsb) < epsilon)
		{
			glb_best_cost = myRound(ml_best_cost);
			glb_best_sol.cpy(ml_best_sol, graph);

			if(myRound(ml_best_cost) < myRound(rbest_cost))
			{
				rbest_cost = myRound(ml_best_cost);
				rbest_time =  (clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
				rcd_cost.push_back(rbest_cost);
				rcd_time.push_back(rbest_time);
			}
		}

		// 更新记录
		each_run_time[i] = rbest_time;
		each_run_rlt[i] = ml_best_cost;
//		each_hit_time[i] = hit_time;
		each_best_sol[i] = Solution(ml_best_sol, graph);

		printf("Round %d: best=%.4f, all time=%.2f||", i+1, ml_best_cost, each_run_time[i]);
		ml_best_sol.verify_infeasible(graph, "一轮run之后");											// DEBUG
		ml_best_sol.verify_dist(graph, "一轮run之后");
		fflush(stdout);
		best_sol.free_memory(graph);
	}
//	print_record();
	// 计算各种指标
	cal_indicators();

	FILE *writefile;

	// 记录每次运行的结果
//	string outfile_detail = rltfile + "_detail" + ".txt";
//	writefile = fopen(outfile_detail.c_str(), "a+");
	char outfile_detail[1000];
	sprintf(outfile_detail, "./results/detail/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_detail, "a+");
 

	fprintf(writefile, "===========%s\n", filename);
	for (int i = 0; i < runs; i++)
	{
		fprintf(writefile, "\t%d--", i);
		fprintf(writefile, "\t%.4f \t%.4f \t%.4f", each_run_rlt[i], each_best_sol[i].ifsb, each_run_time[i]);

		fprintf(writefile, "\tpartitioning[] = ");
		sort_pid(each_best_sol[i].ptn, graph);
		for(int j = 0; j < graph.n; j++)
		{
			fprintf(writefile, "%d ", each_best_sol[i].ptn[j]);
		}
		fprintf(writefile, "\n");
	}
	fprintf(writefile, "\n");
	fclose(writefile);

	// 记录结果概要
//	string outfile_overall = rltfile + "_overall" + ".txt";
//	writefile = fopen(outfile_overall.c_str(), "a+");
	char outfile_overall[1000];
	sprintf(outfile_overall, "./results/overall/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(outfile_overall, "a+");
 

	// 与target最好对比
	string cmp = " ";

	fprintf(writefile, "%s\t", filename);
	fprintf(writefile, "%d/%d\t", hit, runs);
//	fprintf(writefile, "%s \t%.4f \t%.4f \t%.4f \t\t%.6f \t%.4f \t%.4f \t%.4f\n", cmp.c_str(), glb_best_cost, target_cost, avg_cost, avg_time, avg_htime, std_dev, std_dev/glb_best_cost);
	fprintf(writefile, "graph_name:%s \nbest_cost:%.4f \navg_cost:%.4f \navg_time:%.6f\n",
				cmp.c_str(),
				myRound(glb_best_cost, 4),
				myRound(avg_cost, 4),
				rbest_time);
	fclose(writefile);

	char cvg_file[1000];
	sprintf(cvg_file, "./results/convergence/%s_%d.txt", graph_name.c_str(), seed);
	writefile = fopen(cvg_file, "a+");

	int size = int(rcd_cost.size());
	for(int i = 0; i < size; i++)
	{
		fprintf(writefile, "%.4f, %.4f\n", rcd_time[i], rcd_cost[i]);
	}

	fclose(writefile);


	// 释放内存
	free_memory(graph);
}

/* for tuning*/
void execute_each_instance(char * filename, string rltfile)
{
	// 输出标题信息
//	output_header(filename, rltfile);
	execute_instance(filename, rltfile);

	printf("\nnexit:%d\n", nexit);fflush(stdout);
}
#endif /* TEST_MLWT_V2_0_0_H_ */
