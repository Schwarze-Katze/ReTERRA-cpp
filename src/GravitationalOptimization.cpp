#include "GravitationalOptimization.h"

int GravitationalOptimization(const vector<Point_2>& V1, const Point_2& avg, const vector<Point_2>& coveredTarget, const Eigen::VectorXi& solutionSetsLabelsV, const  MatrixXi& setCoverTable, vector<Point_2>& VOpt, vector<Point_2>& VRes) {
    auto coveredTmp = coveredTarget;
    VOpt.clear();
    VRes.resize(setCoverTable.cols());
    for (int i = 0;i < setCoverTable.cols();++i) {
        for (int k = 0;k < solutionSetsLabelsV.size();++k) {
            vector<Point_2> pa;
            if (i == solutionSetsLabelsV(k) and !isSamePoint(V1[i], TERRAConfig::problemParam.Home)) {
                Point_2 pm = V1[i], pSol;
                if (!isSamePoint(pm, avg)) {
                    for (int j = 0;j < setCoverTable.rows();++j) {
                        if (setCoverTable(j, i) and coveredTmp[j].x() >-0.1) {
                            pa.push_back(coveredTarget[j]);
                            coveredTmp[j] = Point_2(-1, -1);
                        }
                    }
                    //pSol = FOptimus(pa, pm, avg);
                }
                else {
                    pSol = V1[i];
                }
                VOpt.push_back(pSol);
                VRes[i] = pSol;
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
    for (auto &centre : cpa) {
        Circle_2 circle(centre, CGAL::square(TERRAConfig::problemParam.Radius));
        typedef typename CGAL::CK2_Intersection_traits<CK, Line_2, Circle_2>::type IntersectionResult;
        std::vector<IntersectionResult> ans;
        CGAL::intersection(line, circle, std::back_inserter(ans));
        
    }
    
}