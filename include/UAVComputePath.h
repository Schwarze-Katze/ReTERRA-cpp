#pragma once
#include <iostream>
#include "TERRAUtility.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using std::vector;

int UAVComputePath(vector<Point_2>& coveredTarget, MatrixXi& setCoverTable, VectorXi& solutionSetsLabelsV, vector<Point_2>& V1, vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& uav_distance, int& uav_time, int& stops);