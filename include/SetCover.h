#pragma once
#include <iostream>
#include "TERRAUtility.h"


using std::vector;
using Eigen::MatrixXd;
using Eigen::Dynamic;
using Eigen::all;
using Eigen::last;
using Eigen::NoChange;
using Eigen::MatrixXi;
using Eigen::VectorXi;
using Eigen::seq;


int GreedySetCovering(//2 input
    std::vector<Point_2>& V1, std::vector<Point_2>& coveredTarget, 
    //3 output
    std::vector<Point_2>& V2, MatrixXi& setCoverTable, VectorXi& solutionSetsLabelsV);

inline int SetCoveringProblem(//4 input
    MatrixXi& A, VectorXi& setsLabelsV, VectorXi& setsCardinalitiesV, VectorXi& setsCardinalitiesL,
    //2 output
    MatrixXi& solutionA, VectorXi& solutionSetsLabelsV);
