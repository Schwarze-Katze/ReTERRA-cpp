#pragma once
#include <iostream>
#include "TERRAUtility.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using std::vector;

int UAVComputePath(const vector<Point_2>& coveredTarget, const MatrixXi& setCoverTable, const VectorXi& solutionSetsLabelsV, const vector<Point_2>& V1, const vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& distance, double& time, int& stops);
int SearchUAVOperations(const vector<Point_2>& path, VectorXi& rte, double& dist, double& time, int& stop);
bool isSamePoint(const Point_2& p1, const Point_2& p2);