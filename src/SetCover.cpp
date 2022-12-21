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
int GreedySetCovering(std::vector<Point_2> &V1,std::vector<Point_2> &coveredTarget) {
    
    return 0;
}