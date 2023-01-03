#pragma once
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "TERRA.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

using std::vector;

int GreedySetCovering(std::vector<Point_2>& V1, std::vector<Point_2>& coveredTarget, std::vector<Point_2>& V2);

inline vector<int> SetCoveringProblem(vector<vector<bool>>& A, vector<int>& setsLabelsV, vector<int>& setsCardinalitiesV, vector<int>& setsCardinalitiesL);
