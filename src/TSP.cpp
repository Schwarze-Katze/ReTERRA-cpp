#include "TSP.h"
#include "SearchUAVOperations.h"

//   Syntax:
//   TSPsolution = LKH_TSP(CostMatrix,pars_struct,fname_tsp,LKHdir,TSPLIBdir)
//
//   This functions solves TSP problems using the Lin-Kernighan-Helsgaun
//   solver. It assumes that a compiled executable of the LKH solver as
//   found at: http://www.akira.ruc.dk/~keld/research/LKH/ is available at
//   the LKHdir directory. Furthermore a TSPLIB directory is assumed.
//   For the definition of the TSPLIB and the compilation of the LKH code
//   check the aforementioned url. 
//
//   Inputs:
//   CostMatrix      : the Cost Matrix of the (asymmetric) TSP. [e.g. it can be an NxN matrix of distances]
//   pars_struct     : parameters structure with
//                   : -> CostMatrixMulFactor (value that makes Cost Matrix
//                        almost integer. [eg. pars_struct.CostMatrixMulFactor = 1000; ]
//                     -> user_comment (a user comment for the problem) [optional]
//   fname_tsp       : the filename to save the tsp problem
//   LKHdir          : the directory of the LKH executable
//   TSPLIBdir       : the directory of the TSPLIB files
//
//   Outputs:
//   TSPsolution     : the TSP solution
//
//   Authors:
//   Kostas Alexis (kalexis@unr.edu)
VectorXi LKH_TSP(const MatrixXd& costMatrix, double costMatrixMulFactor, const string& fileName) {
    MatrixXd costMatrix4Tsp = costMatrixMulFactor * costMatrix;
    costMatrix4Tsp.unaryExpr<double(*)(double)>(&floor);
    WriteTSPLibFile(fileName, costMatrix4Tsp);
    std::cout << (std::string(_LKH_DIR) + "lkhExec " + fileName + ".par") << std::endl;
    system((std::string(_LKH_DIR) + "lkhExec " + fileName + ".par").c_str());
    return ReadSolution((fileName + ".txt").c_str());
}

void WriteTSPLibFile(const string& fileName, const MatrixXd& costMatrix) {
    int dims = costMatrix.cols();
    string name = "NAME : " + fileName + '\n';
    string type = "TYPE: TSP\n";
    string comment = "COMMENT : \n";
    string dimension = "DIMENSION : " + std::to_string(dims) + '\n';
    string edgeWeightTime = "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    string edgeWeightFormat = "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n";
    string edgeWeightSection = "EDGE_WEIGHT_SECTION\n";
    string eof_ = "EOF\n";
    string costMatrixStr;
    for (int i = 0;i < dims;++i) {
        for (int j = 0;j < dims;++j) {
            costMatrixStr += std::to_string(int(costMatrix(i, j) + 0.5)) + ' ';
        }
        costMatrixStr += '\n';
    }

    fstream fParam, fTsp;
    fParam.open((fileName + ".tsp").c_str(), std::ios::out | std::ios::trunc);
    fParam << name << comment << type << dimension << edgeWeightTime << edgeWeightFormat << edgeWeightSection << costMatrixStr << eof_;
    fParam.close();

    string problemFile = "PROBLEM_FILE = " + fileName + ".tsp\n";
    string optimum = "OPTIMUM = 378032\n";
    string moveType = "MOVE_TYPE = 5\n";
    string patchingC = "PATCHING_C = 3\n";
    string patchingA = "PATCHING_A = 2\n";
    string runs = "RUNS = 1\n";
    string tourFile = "TOUR_FILE = " + fileName + ".txt\n";
    fTsp.open((fileName + ".par").c_str(), std::ios::out | std::ios::trunc);
    fTsp << problemFile << optimum << moveType << patchingC << patchingA << runs << tourFile;
    fTsp.close();
    return;
}

inline VectorXi ReadSolution(const string& fileName) {
    vector<int> tspSol;
    double tspcost = 0;
    fstream fSol;
    fSol.open(fileName, std::ios::in);
    string tmpStr;
    while (tmpStr != "TOUR_SECTION") {
        fSol >> tmpStr;
    }
    int tmpIdx;
    fSol >> tmpIdx;
    while (tmpIdx != -1) {
        tspSol.push_back(tmpIdx);
        fSol >> tmpIdx;
    }
    return Eigen::Map<VectorXi>(tspSol.data(), tspSol.size());
}
