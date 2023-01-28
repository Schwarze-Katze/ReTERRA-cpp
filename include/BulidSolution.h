#pragma once

#include "TERRAUtility.h"

using std::vector;
using namespace TERRAResult;

int BuildSolution(const vector<Point_2>& ugvPath, const std::vector<std::vector<Point_2>>& uavPath1, const std::vector<std::vector<Point_2>>& uavPath2);
bool isSamePoint(const Point_2& p1, const Point_2& p2);
