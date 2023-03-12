#include "Voronoi.h"

using std::vector;

int VoronoiCoveringTimeOptimize(vector<Point_2>& V1, vector<Point_2>& coveredTarget) {
    vector<Point_2> uncoveredTarget = TERRAConfig::problemParam.Target;
    vector<Point_2> trustedVertices;
    vector<Point_2> usedVertices;
    vector<Point_2> artificialVertices{ TERRAConfig::problemParam.Home };

    bool duplicated = false;
    while (!uncoveredTarget.empty()) {
        // std::cout << "--uncoveredTarget:"<<uncoveredTarget.size() << std::endl;
        // for (auto& tmp : uncoveredTarget) {
        //     std::cout << tmp << std::endl;
        // }
        vector<Point_2> nearestVertices;
        if (!trustedVertices.empty()) {
            for (int i = 0;i < uncoveredTarget.size();++i) {
                int idx = 0;
                double tmin = sqrt((uncoveredTarget[i] - trustedVertices[idx]).squared_length()) / TERRAConfig::uavData.Vuav;
                for (int j = 1;j < trustedVertices.size();++j) {
                    double tcur = sqrt((uncoveredTarget[i] - trustedVertices[j]).squared_length()) / TERRAConfig::uavData.Vuav;
                    if (tcur < tmin) {
                        tmin = tcur;
                        idx = j;
                    }
                }
                nearestVertices.push_back(trustedVertices[idx]);
            }
            std::sort(nearestVertices.begin(), nearestVertices.end());
            auto endPos = std::unique(nearestVertices.begin(), nearestVertices.end(), isUnique);
            nearestVertices.resize(nearestVertices.end() - endPos);
            duplicated = isDuplicated(nearestVertices, usedVertices);
            usedVertices.insert(usedVertices.end(), nearestVertices.begin(), nearestVertices.end());
        }
        vector<Point_2> voronoiVertices = uncoveredTarget;
        voronoiVertices.insert(voronoiVertices.end(), nearestVertices.begin(), nearestVertices.end());
        if (voronoiVertices.size() == 2 or duplicated) {
            double midX = (voronoiVertices[0].x() + voronoiVertices[1].x()) / 2;
            double midY = (voronoiVertices[0].y() + voronoiVertices[1].y()) / 2;
            double dist = sqrt((voronoiVertices[0] - Point_2(midX, midY)).squared_length());
            while (dist > TERRAConfig::problemParam.Radius) {
                midX = (voronoiVertices[0].x() + midX) / 2;
                midY = (voronoiVertices[0].y() + midY) / 2;
                dist = sqrt((voronoiVertices[0] - Point_2(midX, midY)).squared_length());
            }
            artificialVertices.push_back(Point_2(midX, midY));
        }
        else {
            //voronoi
            DT delaunayTriangulation;
            delaunayTriangulation.insert(voronoiVertices.begin(), voronoiVertices.end());
            for (auto tmp = delaunayTriangulation.finite_faces_begin();tmp < delaunayTriangulation.finite_faces_end();++tmp) {
                auto obj = delaunayTriangulation.dual(tmp);
                trustedVertices.push_back(obj);
            }
            
        }
        trustedVertices.insert(trustedVertices.end(), artificialVertices.begin(), artificialVertices.end());
        std::sort(trustedVertices.begin(), trustedVertices.end());
        auto endPos = std::unique(trustedVertices.begin(), trustedVertices.end(), isUnique);
        trustedVertices.resize(endPos - trustedVertices.begin());
        std::cout << "----trustedVertices:" << trustedVertices.size() << std::endl;
        for (auto& tmp : trustedVertices) {
            std::cout << tmp << std::endl;
        }
        std::cout << "----uncoveredTarget:" << uncoveredTarget.size() << std::endl;
        for (auto& tmp : uncoveredTarget) {
            std::cout << tmp << std::endl;
        }
        for (int i = 0;i < uncoveredTarget.size();++i) {
            int bestVertice = -1;
            int targetCovered = -1;
            double Radius = TERRAConfig::problemParam.Radius;
            for (int j = 0;j < trustedVertices.size();++j) {
                if (uncoveredTarget[i].x() != INFINITY) {
                    double d = sqrt((uncoveredTarget[i] - trustedVertices[j]).squared_length());
                    if (d < Radius) {
                        Radius = d;
                        bestVertice = j;
                        targetCovered = i;
                    }
                }
            }
            if (bestVertice != -1) {
                coveredTarget.push_back(uncoveredTarget[targetCovered]);
                uncoveredTarget[targetCovered] = Point_2(INFINITY, INFINITY);
                V1.push_back(trustedVertices[bestVertice]);
            }
        }
        int delCnt = 0;
        for (auto &tmp : uncoveredTarget) {
            if (isinf(tmp.x())) {
                ++delCnt;
            }
        }
        if (delCnt) {
            auto tmpWaypointNotCovered = uncoveredTarget;
            uncoveredTarget.clear();
            for (auto &tmp : tmpWaypointNotCovered) {
                if (!isinf(tmp.x())) {
                    uncoveredTarget.push_back(tmp);
                }
            }
        }
        else {
            trustedVertices.insert(trustedVertices.end(), uncoveredTarget.begin(), uncoveredTarget.end());
        }
    }
    std::sort(V1.begin(), V1.end());
    return 0;
}

bool isDuplicated(vector<Point_2> const v1, vector<Point_2> const v2){
    for (auto &tmp1 : v1) {
        for (auto &tmp2 : v2) {
            if (tmp1 == tmp2) {
                return true;
            }
        }
    }
    return false;
}

bool isUnique(const Point_2& p1, const Point_2& p2) {
    return(p1 - p2).squared_length() < eps;
}
