#pragma once
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "TERRA.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
using std::vector;
using TERRAConfig::ugvData;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::last;
using Eigen::seq;

int tspGaUgv(vector<Point_2>& V3, double& minDist, vector<Point_2>& UGVPath);