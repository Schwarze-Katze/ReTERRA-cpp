#include "UAVComputePath.h"

int UAVComputePath(const vector<Point_2>& coveredTarget, const MatrixXi& setCoverTable, const VectorXi& solutionSetsLabelsV, const vector<Point_2>& V1, const vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& distance, double& time, int& stops) {
    vector<bool> isCoveredTargetReached(coveredTarget.size(), false);
    time = 0;
    distance = 0;
    stops = 0;
    for (int i = 0;i < setCoverTable.cols();++i) {
        for (auto k : solutionSetsLabelsV) {
            if (i == k) {
                vector<Point_2> subPath = { V1[i] };
                for (int j = 0;j < setCoverTable.rows();++j) {
                    if (setCoverTable(j, i) and !isCoveredTargetReached[j]) {
                        subPath.push_back(coveredTarget[j]);
                        isCoveredTargetReached[j] = true;
                    }
                }
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
                    if (isSamePoint(UGVPath[z], V1[i])) {
                        vector<Point_2> singlePathPack;
                        for (auto tmp : rteT) {
                            singlePathPack.push_back(subPath[tmp]);
                        }
                        uavPath2.push_back(singlePathPack);
                    }
                }
            }
        }
    }
    return 0;
}

inline bool isSamePoint(const Point_2& p1, const Point_2& p2) {
    return (p1 - p2).squared_length() < eps;
}