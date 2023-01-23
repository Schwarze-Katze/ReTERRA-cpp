#include "SearchUAVOperations.h"

int SearchUAVOperations(vector<Point_2>& path, VectorXi& rte, double& dist, int& time, int& stop) {
    Point_2 home = path.at(0);
    int stops = 0;
    auto N = path.size();
    auto h = N - 1;
    Vertex vertexNode(nullptr, home, 0, 0, 0, 0, N - 1);
    vertexNode.h = ComputeH_LKH2(path, vertexNode);
    vertexNode.f = vertexNode.h;
    bool solFound = false;
    int expandedNodes = 0;
    dist = 0;
    time = 0;
    vector<Vertex> OPEN(1,vertexNode);
    vector<Vertex> CLOSE;
    Vertex currentNode;
    while (!OPEN.empty()) {
        currentNode = FindLowest(OPEN);
        OPEN[0].f = INFINITY;
        DeleteNode(OPEN);
        if (isSamePoint(currentNode.p, home) and currentNode.r == 0) {
            solFound = true;
            break;
        }
        vector<Vertex> successor = ExpandGraph(currentNode, path, vertexNode);
        for (auto& tmp : successor) {
            Vertex suc;
            int inClosed = isVisited(CLOSE, tmp);
            if (inClosed == -1) {
                int inOpen = isVisited(OPEN, tmp);
                if (inOpen == -1) {
                    tmp.t = INFINITY;
                    tmp.parent = nullptr;
                    suc = tmp;//need predefine
                }
                else {
                    suc = OPEN[inOpen];
                }
                UpdateVertex(OPEN, currentNode, suc, vertexNode);
            }
        }
        expandedNodes += successor.size();
        CLOSE.push_back(currentNode);
    }
    if (solFound) {
        ReconstructPath(currentNode, path, rte, dist);
        for (int i = 1;i < rte.size() - 1;++i) {
            if (rte(i) == 1) {
                stops++;
            }
        }
        time = currentNode.t;
    }
    else {
        std::cout << "UAV Subpath solution NOT found!" << std::endl;
    }
    return 0;
}

Vertex::Vertex(Vertex* _parent, Point_2& _p, int _t, int _g, double _h, double _f, int _r):parent(_parent), p(_p), t(_t), g(_g), h(_h), f(_f), r(_r) { }

double Vertex::x() {
    return p.x();
}
double Vertex::y() {
    return p.y();
}

double ComputeH_LKH2(const vector<Point_2>& path, const Vertex& vertexNode) {
    return 0.0;
}
