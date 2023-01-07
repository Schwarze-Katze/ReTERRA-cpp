#include "SetCover.h"

// vx_sol is the 'x' coordinate of the vertexes of the solution and vy is the
// 'y' coordinate of the vertexes of the solution, and wp_c the list of all
// the covered waypoints. 
// The OUTPUT is the same of the previous version (see below).
//
// GREEDYSCP Greedy SCP algorithm:: function [solution_setsCell, solution_setsLabelsV] = greedyscp(setsCell, setsLabelsV) .
// 	[SolC,SolL] = GREEDYSCP(C, L) if C is an array, creates a cell array SolC that is a solution of Set Cover Problem defined by C, where C{i} = S_i, an input set made by some of the elements we want to cover; 
//    SolC is made by the cells of C selected by the algorithm. The elements that we want to cover are indicates by numbers from 1 to n, where n is the number of elements we want to cover; 
//    therefore, C{i} is a vector of integers between 1 and n.
//
// 	If C is a logical or numerical array of n rows, where C(j,i) > 0 iff element j is contained in set S_i, the output SolC will be a logical array made by the column of log(C) corresponding to the solution
//
// 	If a vector L of integer labels of the elements of C is provided, SolL contains the labels corresponding to SolC. Otherwise SolL contains the positions of elements of SolC in C. SolC and SolL elements are sorted in ascending order of SolL.
//
// 	This is an implementation of the well-known greedy algorithm (Chvátal, 1979), with two small modifications:
// 	* In case of more than one possible choice at one step, the biggest set is chosen.
// 	* Once the solution is found, we check the selected sets to find a better cover solution, removing a set if is a subset of the union of the other set.
//
// 	If you use this code, please cite:
// 	F. Gori, G. Folino, M.S.M. Jetten, E. Marchiori
// 	"MTR: Taxonomic annotation of short metagenomic reads using clustering at multiple taxonomic ranks", Bioinformatics 2010.
// 	doi = 10.1093/bioinformatics/btq649 
int GreedySetCovering(std::vector<Point_2>& V1, std::vector<Point_2>& coveredTarget, std::vector<Point_2>& V2) {
    const int tlen = coveredTarget.size(), vlen = V1.size();
    //vector<vector<bool>> setCoverTable(tlen);
    
    MatrixXi setCoverTable(tlen, vlen);
    for (int i = 0;i < tlen;++i) {
        for (int j = 0;j < vlen;++j) {
            double d = sqrt((coveredTarget[i] - V1[j]).squared_length());
            if (d <= TERRAConfig::problemParam.Radius) {
                setCoverTable(i, j) = 1;
            }
            else {
                setCoverTable(i, j) = 0;
            }
        }
    }
//  Building scp_table, without empty columns. Deleting the vertices which
//  don't surround any waypoint
    std::vector<int> deletedVertices;
    MatrixXi setCoverTable_s(tlen, 0);
    for (int i = 0;i < vlen;++i) {
        
        if (setCoverTable(all, i).any()) {
            //setCoverTable_s(all, Eigen::lastp1) = setCoverTable(all, i);
            setCoverTable_s.conservativeResize(NoChange, setCoverTable_s.cols() + 1);
            setCoverTable_s(all, last) = setCoverTable(all, i);
        }
        else {
            deletedVertices.push_back(i);
        }
    }
    
    // Generate an adecuate labeling of the vertices, if it has been some vertice
    // deleted so Lp is the complement of deleted_vertices
    std::vector<int>Lp;
    for (int i = 1;i <= vlen;++i) {
        Lp.push_back(i);
        if (std::find(deletedVertices.begin(), deletedVertices.end(), i) != deletedVertices.end()) {//find i in deleted vertices
            Lp.back() = 0;
        }
    }

    std::sort(Lp.rbegin(), Lp.rend());
    while (Lp.back() == 0) {
        Lp.pop_back();
    }
    std::reverse(Lp.begin(), Lp.end());
    auto tmpv = setCoverTable_s.colwise().sum().array();
    VectorXi setsCardinalitiesV = tmpv;
    auto tmpl = setCoverTable_s.rowwise().sum().array();
    VectorXi setsCardinalitiesL = tmpl;
    auto& A = setCoverTable_s;
    VectorXi setsLabelsV(Lp.data(), Lp.size());
    auto solution_setsLabelsV = SetCoveringProblem(A, setsLabelsV, setsCardinalitiesV, setsCardinalitiesL);

    return 0;
}

VectorXi SetCoveringProblem(MatrixXi& A, VectorXi& setsLabelsV, VectorXi& setsCardinalitiesV, VectorXi& setsCardinalitiesL) {
    //Identify elements covered by just one set and select them
    vector<int> UniquelyCoveredL;
    for (int i = 0;i < setsCardinalitiesL.size();++i) {
        if (setsCardinalitiesL[i] == 1) {
            UniquelyCoveredL.push_back(i);
        }
    }
    //construct solutions
    MatrixXi solutionA;
    VectorXi solutionSetsLabelsV;
    VectorXi remainingElementsL;
    if (!UniquelyCoveredL.empty()) {
        //toBeSelectL = logical(sum(A(UniquelyCoveredL,:), 1)) ;
        vector<int> toBeSelectedL, notToBeSelectedL;
        auto tmp = A(UniquelyCoveredL, all).colwise().sum().array();
        for (int i = 0;i < tmp.size();++i) {
            if (tmp(i)>0) {
                toBeSelectedL.push_back(i);
            }
            else {
                notToBeSelectedL.push_back(i);
            }
        }
        //solution_A = A(:, toBeSelectL) ;
        solutionA = A(all, toBeSelectedL);
        solutionSetsLabelsV = setsLabelsV(toBeSelectedL);
        A = A(all, notToBeSelectedL).eval();
        setsLabelsV = setsLabelsV(notToBeSelectedL).eval();
        setsCardinalitiesV = setsCardinalitiesV(notToBeSelectedL).eval();
        remainingElementsL = solutionA.rowwise().sum().array() == 0;
    }
    else {
        remainingElementsL.resize(A.rows());
        remainingElementsL.fill(1);
    }
    if (remainingElementsL.size()) {
        std::sort(setsCardinalitiesV.data(), setsCardinalitiesV.data() + setsCardinalitiesV.size(), std::greater<int>());
        VectorXi sortedIndexVector(setsCardinalitiesV.size());
        std::iota(sortedIndexVector.data(), sortedIndexVector.data() + sortedIndexVector.size(), 0);
        std::sort(sortedIndexVector.data(), sortedIndexVector.data() + sortedIndexVector.size(), [&](int i, int j) {return setsCardinalitiesV(i) > setsCardinalitiesV(j);});
        setsLabelsV = setsLabelsV(sortedIndexVector).eval();
        A = A(all, setsLabelsV).eval();
    }
    while (remainingElementsL.size()) {
        int thresholdN = 0;
        for (int i = 0;i < remainingElementsL.size();++i) {
            thresholdN += A(i, 1) and remainingElementsL(i);
        }
        VectorXi indexFocusedSetsLogical = setsCardinalitiesV.array() >= thresholdN;
        MatrixXi focusedSetsA = A(all, indexFocusedSetsLogical).eval();
        VectorXi focusedLabelsV = setsLabelsV(indexFocusedSetsLogical).eval();
        //TODO:greedy_scp.m:376
    }
}