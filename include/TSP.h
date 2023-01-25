#pragma once
#include "TERRAUtility.h"

using std::vector;
using std::string;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXi LKH_TSP(const MatrixXd& costMatrix, double CostMatrixMulFactor, const string& fileName);

extern "C" {
    int LKHmain(char* ParameterFileName);
}