#pragma once
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "TERRAUtility.h"

using std::vector;
using std::deque;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::all;
using Eigen::last;
using TERRAConfig::uavData;

class Vertex {
public:
    const Vertex* parent;
    Point_2 p;
    double t;
    double g;
    double h;
    double f;
    int r;
public:
    Vertex();
    ~Vertex();
    Vertex(const Vertex* parent, Point_2& _p, double t, int g, double h, double f, int r);
    double x();
    double y();
};

int SearchUAVOperations();
double ComputeH_LKH2(const vector<Point_2>& path, const Vertex& vertexNode);
Vertex FindLowest(vector<Vertex>& OPEN);
void DeleteNode(vector<Vertex>& OPEN);
bool isSamePoint(const Point_2& p1, const Point_2& p2);
vector<Vertex> ExpandGraph(const Vertex& currentNode, const vector<Point_2>& path, const Vertex& vertexNode);
int isVisited(const vector<Vertex>& list, const Vertex& node);
void remove(vector<Vertex>& list, const Vertex& node);
void UpdateVertex(vector<Vertex>& OPEN, const Vertex& currentNode, const Vertex& suc, const Vertex& vertexNode);
void ReconstructPath(const Vertex& currentNode, const vector<Point_2>& path, deque<int>& rte, double& dist);
bool isInPath(const Vertex* cur, const Point_2& suc);
bool isSamePoint(const Point_2& p1, const Point_2& p2);
double getDistance(const Point_2& p1, const Point_2& p2);
double getDistance(double x1, double y1, double x2, double y2);
double getTime(double dis);
VectorXi LKH_TSP(const MatrixXd& costMatrix, double CostMatrixMulFactor, const std::string& fileName);