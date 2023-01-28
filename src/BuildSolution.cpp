#include "BulidSolution.h"

int BuildSolution(const vector<Point_2>& ugvPath, const std::vector<std::vector<Point_2>>& uavPath1, const std::vector<std::vector<Point_2>>& uavPath2) {
    const int n = ugvPath.size();
    pathSol.resize(n);
    for (int i = 0;i < n;++i) {
        pathSol[i].ugv = ugvPath[i];
    }
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < uavPath1.size();++j) {
            if (!uavPath1[j].empty() and isSamePoint(pathSol[i].ugv, uavPath1[j][0])) {
                pathSol[i].uav1 = uavPath1[j];
            }
        }
        for (int j = 0;j < uavPath2.size();++j) {
            if (!uavPath2[j].empty() and isSamePoint(pathSol[i].ugv, uavPath2[j][0])) {
                pathSol[i].uav2 = uavPath2[j];
            }
        }
    }
    return 0;
}