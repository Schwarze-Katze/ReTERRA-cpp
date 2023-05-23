#include "UAVComputePath.h"

int UAVComputePath(const vector<Point_2>& coveredTarget, const MatrixXi& setCoverTable, const VectorXi& solutionSetsLabelsV, const vector<Point_2>& V1, const vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& distance, double& time, int& stops) {
    vector<bool> isCoveredTargetReached(coveredTarget.size(), false);
    time = 0;
    distance = 0;
    stops = 0;
    //此处修改
    VD voronoiDiagram;
    for (auto k : solutionSetsLabelsV) {
        voronoiDiagram.insert(V1[k]);
    }
    std::vector<std::pair<Point_2, Point_2>> targetToVertice;
    for (auto& p : coveredTarget) {
        VD::Locate_result loc = voronoiDiagram.locate(p);
        if (VD::Face_handle* fh = boost::get<VD::Face_handle>(&loc)) {
            targetToVertice.emplace_back((*fh)->dual()->point(), p);
            // std::cout << "--pair: \n" << targetToVertice.back().first << '\n' << targetToVertice.back().second << std::endl;
        }
        else if (VD::Vertex_handle* vh = boost::get<VD::Vertex_handle>(&loc)) {
            Point_2 p1((*vh)->point().x() - eps, (*vh)->point().y() - eps);//avoid vertex
            VD::Locate_result loc1 = voronoiDiagram.locate(p1);
            if (VD::Face_handle* fh = boost::get<VD::Face_handle>(&loc1)) {
                targetToVertice.emplace_back((*fh)->dual()->point(), p);
                // std::cout << "--pair: \n" << targetToVertice.back().first << '\n' << targetToVertice.back().second << std::endl;
            }
        }
    }
    for (auto k : solutionSetsLabelsV) {
        vector<Point_2> subPath = { V1[k] };
        for (auto& pair : targetToVertice) {
            if (isSamePoint(pair.first, V1[k])) {
                subPath.push_back(pair.second);
                // std::cout << "--pair: \n" << pair.first << '\n' << pair.second << std::endl;
            }
        }
        // std::cout << "--subPath: " << subPath.size() << std::endl;
        // for (auto& tmp : subPath) {
        //     std::cout << tmp << std::endl;
        // }
        //TIME BASED SEARCHING
        VectorXi rteT;
        double dis, t;
        int stop;
        SearchUAVOperations(subPath, rteT, dis, t, stop);
        std::cout << "--rteT: " << rteT.size() << std::endl;
        std::cout << rteT << std::endl;
        std::cout << "--UGVPath: " << UGVPath.size() << std::endl;
        for (auto& tmp : UGVPath) {
            std::cout << tmp << std::endl;
        }

        time += t;
        distance += dis;
        stops += stop;
        //UAV Path for Distance-Based Searching or Time-Based Searching
        for (int z = 0;z < UGVPath.size();++z) {
            if (isSamePoint(UGVPath[z], V1[k])) {
                vector<Point_2> singlePathPack;
                for (auto tmp : rteT) {
                    singlePathPack.push_back(subPath[tmp]);
                }
                uavPath2.push_back(singlePathPack);
            }
        }
    }
    // original uav subpath generating
    // vector<Point_2> subPath = { V1[i] };
    // for (int j = 0;j < setCoverTable.rows();++j) {
    //     if (setCoverTable(j, i) and !isCoveredTargetReached[j]) {
    //         subPath.push_back(coveredTarget[j]);
    //         isCoveredTargetReached[j] = true;
    //     }
    // }

    return 0;
}

inline bool isSamePoint(const Point_2& p1, const Point_2& p2) {
    return (p1 - p2).squared_length() < eps;
}