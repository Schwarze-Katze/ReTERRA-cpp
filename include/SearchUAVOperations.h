#pragma once
#include "TERRAUtility.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::all;
using Eigen::last;

class Vertex {
public:
    Vertex* parent;
    Point_2 p;
    int t;
    int g;
    double h;
    double f;
    int r;
public:
    Vertex();
    ~Vertex();
    Vertex(Vertex* parent, Point_2& _p, int t, int g, double h, double f, int r);
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
void UpdateVertex(vector<Vertex>& OPEN, const Vertex& currentNode, const Vertex& suc, const Vertex& vertexNode);
void ReconstructPath(const Vertex& currentNode, const vector<Point_2>& path, VectorXi& rte, double& dist);