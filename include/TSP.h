#pragma once
#include "TERRAUtility.h"

using std::vector;
using std::string;
using std::fstream;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXi LKH_TSP(const MatrixXd& costMatrix, double CostMatrixMulFactor, const string& fileName);
inline void WriteTSPLibFile(const string& fileName, const MatrixXd& costMatrix);
inline VectorXi ReadSolution(const string& fileName);

extern "C" {
    int LKHmain(int argc, char* argv[]);
}