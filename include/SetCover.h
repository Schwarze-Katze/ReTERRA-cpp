#pragma once
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <eigen3/Eigen/Eigen>
#include "TERRA.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::all;
using Eigen::last;
using Eigen::NoChange;
using Eigen::MatrixXi;
using Eigen::VectorXi;

int GreedySetCovering(std::vector<Point_2>& V1, std::vector<Point_2>& coveredTarget, std::vector<Point_2>& V2);

inline VectorXi SetCoveringProblem(MatrixXi& A, VectorXi& setsLabelsV, VectorXi& setsCardinalitiesV, VectorXi& setsCardinalitiesL);
