#pragma once
#include <iostream>
#include "TERRAUtility.h"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using std::vector;
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP> VD;

int UAVComputePath(const vector<Point_2>& coveredTarget, const MatrixXi& setCoverTable, const VectorXi& solutionSetsLabelsV, const vector<Point_2>& V1, const vector<Point_2>& UGVPath, vector<vector<Point_2>>& uavPath1, vector<vector<Point_2>>& uavPath2, double& distance, double& time, int& stops);
int SearchUAVOperations(const vector<Point_2>& path, VectorXi& rte, double& dist, double& time, int& stop);
bool isSamePoint(const Point_2& p1, const Point_2& p2);