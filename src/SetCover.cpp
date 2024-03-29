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
int GreedySetCovering(const std::vector<Point_2>& V1, const std::vector<Point_2>& coveredTarget, std::vector<Point_2>& V2, MatrixXi& setCoverTable, VectorXi& solutionSetsLabelsV) {
    const int tlen = coveredTarget.size(), vlen = V1.size();
    //vector<vector<bool>> setCoverTable(tlen);
    
    setCoverTable.resize(tlen, vlen);
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
    // std::cout << "--setCoverTable: " << setCoverTable.rows() << ',' << setCoverTable.cols() << std::endl;
    // std::cout << setCoverTable << std::endl;
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
    for (int i = 0;i < vlen;++i) {
        Lp.push_back(i);
        if (std::find(deletedVertices.begin(), deletedVertices.end(), i) != deletedVertices.end()) {//find i in deleted vertices
            Lp.back() = -1;
        }
    }

    std::sort(Lp.rbegin(), Lp.rend());
    while (Lp.back() == -1) {
        Lp.pop_back();
    }
    std::reverse(Lp.begin(), Lp.end());
    auto tmpv = setCoverTable_s.colwise().sum().array();
    VectorXi setsCardinalitiesV = tmpv;
    auto tmpl = setCoverTable_s.rowwise().sum().array();
    VectorXi setsCardinalitiesL = tmpl;
    auto& A = setCoverTable_s;
    VectorXi setsLabelsV = Eigen::Map<VectorXi>(Lp.data(), Lp.size());
    MatrixXi solutionA;
    // std::cout << "--A1: " << A.rows() << ',' << A.cols() << std::endl;
    // std::cout << A << std::endl;
    // std::cout << "--setsCardinalitiesV: " << setsCardinalitiesV.rows() << ',' << setsCardinalitiesV.cols() << std::endl;
    // std::cout << setsCardinalitiesV << std::endl;
    // std::cout << "--setsCardinalitiesL: " << setsCardinalitiesL.rows() << ',' << setsCardinalitiesL.cols() << std::endl;
    // std::cout << setsCardinalitiesL << std::endl;
    SetCoveringProblem(A, setsLabelsV, setsCardinalitiesV, setsCardinalitiesL, solutionA, solutionSetsLabelsV);
    for (int i = 0;i < vlen;++i) {
        for (int j = 0;j < solutionSetsLabelsV.size();++j) {
            if (i == solutionSetsLabelsV(j)){
                V2.push_back(V1[i]);
            }
        }
    }
    return 0;
}

int SetCoveringProblem(MatrixXi A, VectorXi setsLabelsV, VectorXi setsCardinalitiesV, VectorXi setsCardinalitiesL, MatrixXi& solutionA, VectorXi& solutionSetsLabelsV) {
    //Identify elements covered by just one set and select them
    vector<int> UniquelyCoveredL;
    for (int i = 0;i < setsCardinalitiesL.size();++i) {
        if (setsCardinalitiesL[i] == 1) {
            UniquelyCoveredL.push_back(i);
        }
    }
    
    //construct solution
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
        solutionA = A(all, toBeSelectedL);
        solutionSetsLabelsV = setsLabelsV(toBeSelectedL);
        A = A(all, notToBeSelectedL).eval();
        // std::cout << "--A2: " << A.rows() << ',' << A.cols() << std::endl;
        // std::cout << A << std::endl;
        // std::cout << "--solutionSetsLabelsV: " << solutionSetsLabelsV.rows() << ',' << solutionSetsLabelsV.cols() << std::endl;
        // std::cout << solutionSetsLabelsV << std::endl;
        setsLabelsV = setsLabelsV(notToBeSelectedL).eval();
        setsCardinalitiesV = setsCardinalitiesV(notToBeSelectedL).eval();
        remainingElementsL = (solutionA.rowwise().sum().array() == 0).cast<int>();
        // std::cout << "--remainingElementsL: " << remainingElementsL.rows() << ',' << remainingElementsL.cols() << std::endl;
        // std::cout << remainingElementsL << std::endl;
        // std::cout << "--setsLabelsV: " << setsLabelsV.rows() << ',' << setsLabelsV.cols() << std::endl;
        // std::cout << setsLabelsV << std::endl;

    }
    else {
        remainingElementsL.resize(A.rows());
        remainingElementsL.fill(1);
    }
    if (remainingElementsL.size()) {
        VectorXi sortedIndexVector(setsCardinalitiesV.size());
        std::iota(sortedIndexVector.data(), sortedIndexVector.data() + sortedIndexVector.size(), 0);
        std::stable_sort(sortedIndexVector.data(), sortedIndexVector.data() + sortedIndexVector.size(), [&](int i, int j) {return setsCardinalitiesV(i) > setsCardinalitiesV(j);});
        // std::cout << "--sortedIndexVector: " << sortedIndexVector.rows() << ',' << sortedIndexVector.cols() << std::endl;
        // std::cout << sortedIndexVector << std::endl;
        std::sort(setsCardinalitiesV.data(), setsCardinalitiesV.data() + setsCardinalitiesV.size(), std::greater<int>());
        // std::cout << "--setsCardinalitiesV: " << setsCardinalitiesV.rows() << ',' << setsCardinalitiesV.cols() << std::endl;
        // std::cout << setsCardinalitiesV << std::endl;
        setsLabelsV = setsLabelsV(sortedIndexVector).eval();
        // std::cout << "--setsLabelsV: " << setsLabelsV.rows() << ',' << setsLabelsV.cols() << std::endl;
        // std::cout << setsLabelsV << std::endl;
        A = A(all, sortedIndexVector).eval();
        // std::cout << "--A3: " << A.rows() << ',' << A.cols() << std::endl;
        // std::cout << A << std::endl;

        while (remainingElementsL.count()) {
            int thresholdN = 0;
            for (int i = 0;i < remainingElementsL.size();++i) {
                thresholdN += A(i, 0) and remainingElementsL(i);
            }
            Eigen::Logical indexFocusedSetsLogical(setsCardinalitiesV.array() >= thresholdN);
            MatrixXi focusedSetsA = A(all, indexFocusedSetsLogical).eval();
            VectorXi focusedLabelsV = setsLabelsV(indexFocusedSetsLogical).eval();
            // std::cout << "--setsCardinalitiesV: " << setsCardinalitiesV.rows() << ',' << setsCardinalitiesV.cols() << std::endl;
            // std::cout << setsCardinalitiesV << std::endl;
            // std::cout << "--A: " << A.rows() << ',' << A.cols() << std::endl;
            // std::cout << A << std::endl;
            VectorXi intV = remainingElementsL.transpose() * focusedSetsA;
            // std::cout << "--intV: " << intV.rows() << ',' << intV.cols() << std::endl << intV << std::endl;
            int selectedSetPosN;
            intV.maxCoeff(&selectedSetPosN);
            // std::cout << "--focusedSetsA: " << focusedSetsA.rows() << ',' << focusedSetsA.cols() << std::endl;
            // std::cout << focusedSetsA << std::endl;
            // std::cout << "--focusedLabelsV: " << focusedLabelsV.rows() << ',' << focusedLabelsV.cols() << std::endl;
            // std::cout << focusedLabelsV << std::endl;
            solutionA.conservativeResize(focusedSetsA.rows(), solutionA.cols() + 1);
            solutionA(all, last) = focusedSetsA(all, selectedSetPosN);
            // std::cout << "--solutionA: " << solutionA.rows() << ',' << solutionA.cols() << std::endl;
            // std::cout << solutionA << std::endl;
            solutionSetsLabelsV.conservativeResize(solutionSetsLabelsV.size() + 1);
            solutionSetsLabelsV(solutionSetsLabelsV.size() - 1) = focusedLabelsV(selectedSetPosN);
            remainingElementsL = (remainingElementsL.eval().array().cast<bool>() && !(focusedSetsA(all, selectedSetPosN).array().cast<bool>())).cast<int>();
            
            vector<int> vectorRelativePosV;
            for (int i = 0;i < setsCardinalitiesV.size();++i) {
                if (setsCardinalitiesV(i) >= thresholdN) {
                    vectorRelativePosV.push_back(i);
                }
            }
            int target = vectorRelativePosV[selectedSetPosN];
            vector<int> vectorNotRelativePosV;
            for (int i = 0;i < setsCardinalitiesV.size();++i) {
                if (target != i) {
                    vectorNotRelativePosV.push_back(i);
                }
            }
            VectorXi notRelativePosV = Eigen::Map<VectorXi>(vectorNotRelativePosV.data(), vectorNotRelativePosV.size());
            A = A(all, notRelativePosV).eval();
            setsLabelsV = setsLabelsV(notRelativePosV).eval();
            setsCardinalitiesV = setsCardinalitiesV(notRelativePosV).eval();
            // std::cout << "--remainingElementsL: " << remainingElementsL.rows() << ',' << remainingElementsL.cols() << std::endl;
            // std::cout << remainingElementsL << std::endl;
            // std::cout << "--A4: " << A.rows() << ',' << A.cols() << std::endl;
            // std::cout << A << std::endl;
            // std::cout << "--setsLabelsV: " << setsLabelsV.rows() << ',' << setsLabelsV.cols() << std::endl;
            // std::cout << setsLabelsV << std::endl;
            // std::cout << "--setsCardinalitiesV: " << setsCardinalitiesV.rows() << ',' << setsCardinalitiesV.cols() << std::endl;
            // std::cout << setsCardinalitiesV << std::endl;
        }
        vector<int> reserveIdx;
        MatrixXi tmpA = solutionA;
        for (int setIdx = solutionSetsLabelsV.size() - 1;setIdx >= 0;--setIdx) {
            tmpA(all, setIdx).fill(0);
            if (!(tmpA.rowwise().sum().all())) {
                reserveIdx.push_back(setIdx);
                tmpA(all, setIdx) = solutionA(all, setIdx);
            }
        }
        solutionA = solutionA(all, reserveIdx).eval();
        solutionSetsLabelsV = solutionSetsLabelsV(reserveIdx).eval();
    }
    VectorXi sortedSolutionIndexVector(solutionSetsLabelsV.size());
    std::iota(sortedSolutionIndexVector.data(), sortedSolutionIndexVector.data() + sortedSolutionIndexVector.size(), 0);
    std::sort(sortedSolutionIndexVector.data(), sortedSolutionIndexVector.data() + sortedSolutionIndexVector.size(), [&](int i, int j) {return solutionSetsLabelsV(i) < solutionSetsLabelsV(j);});
    std::sort(solutionSetsLabelsV.data(), solutionSetsLabelsV.data() + solutionSetsLabelsV.size());
    solutionA = solutionA(all, sortedSolutionIndexVector).eval();
    return 0;
}