#include "UAVComputePath.h"

int UAVComputePath(vector<Point_2>& coveredTarget, MatrixXi& setCoverTable, VectorXi& solutionSetsLabelsV, vector<Point_2>& V1, std::vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& distance, int& time, int& stops) {
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
                double dis;
                int stop, t;
                //SearchUAVOperations(subPath, rteT, dis, t, stop);
                time += t;
                distance += dis;
                stops += stop;
                //UAV Path for Distance-Based Searching or Time-Based Searching
                for (int z = 0;z < UGVPath.size();++z) {
                    if (UGVPath[z] == V1[i]) {
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
