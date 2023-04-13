#pragma once
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Circular_kernel_intersections.h>
#include "TERRAUtility.h"


using std::vector;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;

typedef CGAL::Exact_circular_kernel_2 CK;
typedef CK::Line_2 Line_2;
typedef CK::Circle_2 Circle_2;
typedef CK::Intersect_2 Intersect_2;
typedef typename CGAL::CK2_Intersection_traits<CK, Line_2, Circle_2>::type IntersectionResult;


int GravitationalOptimization(const vector<Point_2>& V1, const Point_2& avg, const vector<Point_2>& coveredTarget, VectorXi solutionSetsLabelsV, MatrixXi setCoverTable, vector<Point_2>& VOpt, vector<Point_2>& VRes);
bool isSamePoint(const Point_2& p1, const Point_2& p2);
Point_2 FOptimus(const vector<Point_2>& pa, const Point_2& pm, const Point_2& avg);