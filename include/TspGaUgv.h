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
    vector<Point_2>& V3,
    //3 output
    double& minDist, vector<Point_2>& UGVPath, VectorXi& rte);

MatrixXi& OX(//1 output
    MatrixXi& child, 
    //2 input
    VectorXi& firstParent, VectorXi& secondParent);

MatrixXi& CX(MatrixXi& child, VectorXi& firstParent, VectorXi& secondParent);

MatrixXi& OBX(MatrixXi& child, VectorXi& firstParent, VectorXi& secondParent);