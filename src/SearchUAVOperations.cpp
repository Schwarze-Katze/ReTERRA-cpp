#include "SearchUAVOperations.h"

int SearchUAVOperations(vector<Point_2>& path, VectorXi& rte, double& dist, double& time, int& stop) {
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
        expandedNodes += successor. size();  
        CLOSE.push_back(currentNode);
    }
    if (solFound) {
        deque<int> tmpRte(rte.begin(), rte.end());
        ReconstructPath(currentNode, path, tmpRte, dist);
        rte = Eigen::Map<VectorXi>(&tmpRte[0], tmpRte.size());//CANNOT DO THIS
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

Vertex::Vertex(const Vertex* _parent, Point_2& _p, double _t, int _g, double _h, double _f, int _r):parent(_parent), p(_p), t(_t), g(_g), h(_h), f(_f), r(_r) { }

double Vertex::x() {
    return p.x();
}
double Vertex::y() {
    return p.y();
}

double ComputeH_LKH2(const vector<Point_2>& path, const Vertex& suc) {
    double h = 0.0;
    vector<Point_2> rnodes(1, suc.p);
    double costToHome = 0;
    VectorXi TSPSolution;
    for (int i = 1;i < path.size();++i) {
        if (!isInPath(&suc, path[i])) {
            rnodes.push_back(path[i]);
            costToHome += getTime(getDistance(path[0], path[i]));
        }
    }
    int So = std::ceil(costToHome / TERRAConfig::uavData.TimeBudget) + 0.1;
    if (So > 0) {
        const int n = rnodes.size();
        MatrixXd costMatrix(n, n);
        for (int i = 0;i < n;++i) {
            for (int j = 0;j < n;++j) {
                costMatrix(i, j) = getTime(getDistance(rnodes[i], rnodes[j]));
            }
        }
        if (n > 2) {
            TSPSolution = LKH_TSP(costMatrix, 1, "tsp_solution");
        }
        else {
            TSPSolution << 1, 2, 1;
        }
        vector<Point_2> orderedNodes(n);
        for (int i = 0;i < n;++i) {
            orderedNodes[i] = rnodes[TSPSolution[i]];
        }
        cv::Mat cvOrderedNodes(2, n, CV_64F);
        for(int i=0;i<n;++i){
            cvOrderedNodes.at<double>(0, 0) = orderedNodes[i].x();
            cvOrderedNodes.at<double>(1, 0) = orderedNodes[i].y();
        }
        cv::Mat labels;
        
        //NOTE: apply 3 attempts then pick the minimum. Set attempts to 1 for better performance.
        cv::kmeans(cvOrderedNodes, So, labels, cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 100, 1e-6), 3, cv::KMEANS_PP_CENTERS);
        MatrixXd clusterMatrix(3, So);
        clusterMatrix.row(0).setConstant(0);
        clusterMatrix.row(1).setConstant(orderedNodes[0].x());
        clusterMatrix.row(2).setConstant(orderedNodes[0].y());
        for (int i = 1;i < n;++i) {
            int idx = labels.at<int>(i - 1);
            clusterMatrix(0, idx) += getTime(getDistance(clusterMatrix(1, idx), clusterMatrix(2, idx), orderedNodes[idx].x(), orderedNodes[idx].y()));
        }
        double Tch = costToHome * 2;
        for (int i = 0;i < So;++i) {
            clusterMatrix(0, i) += getTime(getDistance(clusterMatrix(1, i), clusterMatrix(2, i), path[0].x(), path[0].y()));
            double Se = clusterMatrix(0, i) / TERRAConfig::uavData.TimeBudget;
            h += Se * Tch + 2 * clusterMatrix(0, i);
        }
    }
    return h;
}

Vertex FindLowest(vector<Vertex>& OPEN) {
    int co = OPEN.size();
    vector<Vertex> L;
    int cl = 0;
    int lowest = 0;
    Vertex node;
    bool isNodeAvailable = false;
    while (cl < co)
    {
        for (int i = 0;i < co;++i) {
            if (OPEN[i].f != INFINITY) {
                node = OPEN[i];
                isNodeAvailable = true;
                lowest = i;
                break;
            }
        }
        if (isNodeAvailable) {
            for (int i = 0;i < co;++i) {
                if (OPEN[i].f < node.f) {
                    node = OPEN[i];
                    lowest = i;
                }
            }
            OPEN[lowest].f = INFINITY;
            L.push_back(node);
            isNodeAvailable = false;
            cl = L.size();
        }
    }
    node = L[0];
    return node;
}

void DeleteNode(vector<Vertex>& OPEN) {
    vector<Vertex> L;
    for (int i = 0;i < OPEN.size();++i) {
        if (OPEN[i].f != INFINITY) {
            L.push_back(OPEN[i]);
        }
    }
    OPEN = L;
    return;
}

inline bool isSamePoint(const Point_2& p1, const Point_2& p2) {
    return (p1 - p2).squared_length() < eps;
}

vector<Vertex> ExpandGraph(const Vertex& currentNode, const vector<Point_2>& path, const Vertex& vertexNode) {
    vector<Vertex> L;
    bool isRel = true;
    for (int i = 0;i < path.size();++i) {
        auto suc = path[i];
        if (!isSamePoint(currentNode.p, suc)) {
            if (isSamePoint(vertexNode.p, suc)) {
                isRel = false;
            }
            else {
                isRel = isInPath(&currentNode, suc);
            }
            double dToSuc, tToSuc, g;
            if (!isRel) {
                dToSuc = getDistance(currentNode.p, suc);
                tToSuc = getTime(dToSuc);
                if (isSamePoint(currentNode.p, vertexNode.p)) {
                    g = uavData.TimeTakeoff + tToSuc;
                    tToSuc += currentNode.t + uavData.TimeTakeoff;
                }
                else {
                    g = currentNode.g + tToSuc;
                    if (isSamePoint(suc, vertexNode.p)) {
                        auto tc = uavData.TimeBudget - g;
                        tToSuc += currentNode.t + uavData.TimeLanding + tc;
                    }
                    else {
                        tToSuc += currentNode.t;
                    }
                }
            }
            auto dToHome = getDistance(vertexNode.p, suc);
            auto tToHome = getTime(dToHome);
            if (uavData.TimeBudget >= std::floor(g + tToHome + uavData.TimeLanding)) {
                Vertex successor(&(currentNode), suc, tToSuc, g, 0, 0, 0);
                if (isSamePoint(suc, vertexNode.p)) {
                    successor.r = currentNode.r;
                }
                else {
                    successor.r = currentNode.r - 1;
                }
                successor.h = ComputeH_LKH2(path, successor);
                successor.f = tToSuc + successor.h;
                L.push_back(successor);
            }
        }
    }
    return L;
}

int isVisited(const vector<Vertex>& _list, const Vertex& node) {
    int nodeIdx = 0;
    auto list = _list;
    if (!list.empty()) {
        for (int i = 0;i < list.size();++i) {
            if (isSamePoint(list[i].p, node.p) and list[i].h - node.h < eps and list[i].r == node.r) {
                nodeIdx = i;
                break;
            }
        }
    }
    return nodeIdx;
}

void remove(vector<Vertex>& list, const Vertex& node) {
    int idx = 0;
    if (!list.empty()) {
        for (int i = 0;i < list.size();++i) {
            if (isSamePoint(list[i].p, node.p) and list[i].h - node.h < eps and list[i].r == node.r) {
                idx = i;
                break;
            }
        }
        list[idx].f = INFINITY;
        DeleteNode(list);
    }
}

void UpdateVertex(vector<Vertex>& OPEN, const Vertex& currentNode, const Vertex& _suc, const Vertex& vertexNode) {
    auto suc = _suc;
    auto d = getDistance(suc.p, currentNode.p);
    auto t = getTime(d);
    if (isSamePoint(currentNode.p, vertexNode.p)) {
        t += uavData.TimeTakeoff;
    }
    else{
        if (isSamePoint(suc.p, vertexNode.p)) {
            t = uavData.TimeLanding + uavData.TimeBudget;
        }
    }
    if (currentNode.t + t < suc.t) {
        suc.t = currentNode.t + t;
        suc.parent = &currentNode;
        if (isVisited(OPEN, suc)) {
            remove(OPEN, suc);
        }
        suc.f = suc.t + suc.h;
        OPEN.push_back(suc);
    }
}

void ReconstructPath(const Vertex& currentNode, const vector<Point_2>& path, deque<int>& rte, double& dist) {
    
    for (int i = 0;i < path.size();++i) {
        if (isSamePoint(path[i], currentNode.p)) {
            
        }
    }
}
