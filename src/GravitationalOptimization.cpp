#include "GravitationalOptimization.h"

int GravitationalOptimization(const vector<Point_2>& V1, const Point_2& avg, const vector<Point_2>& coveredTarget, VectorXi solutionSetsLabelsV,  MatrixXi setCoverTable, vector<Point_2>& VOpt, vector<Point_2>& VRes) {
    std::vector<bool> coveredTmp(coveredTarget.size(), true);
    // std::cout << "--coveredTarget: " << coveredTarget.size() << std::endl;
    // for (auto& tmp : coveredTarget) {
    //     std::cout << tmp << std::endl;
    // }
    VOpt.clear();
    VOpt.push_back(TERRAConfig::problemParam.Home);
    VRes.resize(setCoverTable.cols());
    for (int i = 0;i < setCoverTable.cols();++i) {
        for (int k = 0;k < solutionSetsLabelsV.size();++k) {
            vector<Point_2> pa;
            if (i == solutionSetsLabelsV(k) and !isSamePoint(V1[i], TERRAConfig::problemParam.Home)) {
                Point_2 pm = V1[i], pSol;
                if (!isSamePoint(pm, avg)) {
                    for (int j = 0;j < setCoverTable.rows();++j) {
                        if (setCoverTable(j, i) and coveredTmp[j]) {
                            pa.push_back(coveredTarget[j]);
                            coveredTmp[j] = false;
                        }
                    }
                    // std::cout << "--pa: " << pa.size() << std::endl;
                    // for (auto& tmp : pa) {
                    //     std::cout << tmp << std::endl;
                    // }
                    pSol = FOptimus(pa, pm, avg);
                }
                else {
                    pSol = V1[i];
                }
                // std::cout << "--pSol: " << pSol << std::endl;
                if (!isinf(pSol.x()) and !isinf(pSol.y())) {
                    VOpt.push_back(pSol);
                VRes[i] = pSol;
                }
            }
            else if (i == solutionSetsLabelsV(k) and isSamePoint(V1[i], TERRAConfig::problemParam.Home)) {
                VRes[i] = V1[i];
            }
        }
    }
    return 0;
}

Point_2 FOptimus(const vector<Point_2>& pa, const Point_2& pm, const Point_2& avg) {
    vector<CK::Point_2> cpa;
    CK::Point_2 cpm(pm.x(), pm.y());
    CK::Point_2 cavg(avg.x(), avg.y());
    for (int i = 0;i < pa.size();++i) {
        cpa.push_back(CK::Point_2(pa[i].x(), pa[i].y()));
    }
    Line_2 line(cpm, cavg);
    // std::cout << "--Line: " << std::endl;
    // std::cout << "k=" << -CGAL::to_double(line.a() / line.b()) << ", b=" << -CGAL::to_double(line.c() / line.b()) << std::endl;
    vector<Point_2> intersect;
    for (auto& centre : cpa) {
        Circle_2 circle(centre, CGAL::square(TERRAConfig::problemParam.Radius));
        // std::cout << "--Circle: " << std::endl;
        // std::cout << '(' << CGAL::to_double(circle.center().x()) << ' ' << CGAL::to_double(circle.center().y()) << "), R=" << CGAL::sqrt(CGAL::to_double(circle.squared_radius())) << std::endl;
        std::vector<IntersectionResult> ans;
        CGAL::intersection(line, circle, std::back_inserter(ans));
        for (const auto& elem : ans) {
            const auto& tmp = boost::get<std::pair<CGAL::Circular_arc_point_2<CGAL::Exact_circular_kernel_2>, unsigned int>>(elem).first;
            intersect.emplace_back(CGAL::to_double(tmp.x()), CGAL::to_double(tmp.y()));
        }
    }
    // std::cout << "--intersect: " << intersect.size() << std::endl;
    for (auto& tmp : intersect) {
        std::cout << tmp << std::endl;
    }
    vector<Point_2> validPoint;
    for (int j = 0;j < intersect.size();++j) {
        int cnt = 0;
        for (int k = 0;k < pa.size();++k) {
            if (sqrt((intersect[j] - pa[k]).squared_length()) <= TERRAConfig::problemParam.Radius) {
                ++cnt;
            }
        }
        if (cnt == pa.size()) {
            validPoint.push_back(intersect[j]);
        }
    }
    // std::cout << "--validPoint: " << validPoint.size() << std::endl;
    // for (auto& tmp : validPoint) {
    //     std::cout << tmp << std::endl;
    // }
    double minDist = INFINITY, minX = INFINITY, minY = INFINITY;
    for (const auto& tmp : validPoint) {
        if (sqrt((tmp - avg).squared_length()) <= minDist) {
            minX = tmp.x();
            minY = tmp.y();
            minDist = sqrt((tmp - avg).squared_length());
        }
    }
    return Point_2(minX, minY);
}

inline bool isSamePoint(const Point_2& p1, const Point_2& p2) {
    return (p1 - p2).squared_length() < eps;
}