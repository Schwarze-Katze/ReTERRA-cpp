#pragma once
#include <cmath>
#include "TERRAUtility.h"


using std::vector;
using TERRAConfig::ugvData;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::all;
using Eigen::last;
using Eigen::seq;
using Eigen::lastN;
using Eigen::NoChange;
using Eigen::Logical;

int tspGaUgv(//1 input
    const vector<Point_2>& V3,
    //3 output
    double& minDist, vector<Point_2>& UGVPath, VectorXi& rte);

int OX(MatrixXi& child, const VectorXi& p1, const VectorXi& p2);

int CX(MatrixXi& child, VectorXi p1, VectorXi p2);

int OBX(MatrixXi& child, const VectorXi& p1, const VectorXi& p2);