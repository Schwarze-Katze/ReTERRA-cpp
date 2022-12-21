#pragma once
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "TERRA.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

int GreedySetCovering(std::vector<Point_2> &V1,std::vector<Point_2> &coveredTarget);