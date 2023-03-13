#include "TERRA.h"
#include "TERRAUtility.h"



namespace TERRAConfig {
    
    ConfigParam::ConfigParam() { }

    ConfigParam::ConfigParam(int iter, bool printAns, bool saveAns, bool vrep, const std::string& _fullname):iterations(iter), printResults(printAns), saveResults(saveAns), Vrep(vrep), fullName(_fullname) {
        time_t now = time(0);
        tm* ltm = localtime(&now);
        saveDir = "Results/Test-" + std::to_string(ltm->tm_year + 1900) + "." + std::to_string(ltm->tm_mon + 1) + "." + std::to_string(ltm->tm_mday) + "_" + std::to_string(ltm->tm_hour) + "." + std::to_string(ltm->tm_min) + "." + std::to_string(ltm->tm_sec) + "/";

        if (saveResults) {
            std::cout << "mkdir -p ./" + saveDir << std::endl;
            system(("mkdir -p ./" + saveDir).c_str());//Create directory if saveresults = true
        }
    }

    ConfigParam::~ConfigParam() { }


    ProblemParam::ProblemParam(int targetCnt, int radius, int delta, double homeX, double homeY, int area, const std::string& _Gp):Radius(radius), TargetCnt(targetCnt), Delta(delta), AreaSize(area), Gp(_Gp) {
        Home = Point_2(homeX, homeY);
    }
    ProblemParam::ProblemParam() { }
    ProblemParam::~ProblemParam() { }

    // This method creates a map ensuring that 'n' target points are distributed
    // in 'delta' groups inside an specific radius 2x'r', around a map
    // with an specific area.
    // Mandatory: 'n' must be multiple of 'delta'
    //
    // INPUTS
    // problem_params
    //
    // OUTPUTS
    // xy = coordinates of the target points
    //
    // Example :
    // n = 12
    // r = 3
    // area = 50 km2
    // delta = 2
    // Solution : This algorithm will create 6 groups of 2 target points inside a
    // r = 3 km of radius in a 50km2 area
    int ProblemParam::SceneGenerator(const std::string& fileName, const std::string& iterDir) {
        std::random_device realRng;
        std::mt19937 rng(realRng());

        std::vector<Point_2> groups;
        std::vector<Point_2> targetRaw;

        if (TargetCnt > 0 and AreaSize > 0 and Delta > 0) {
            int groupCnt = TargetCnt / Delta;
            for (int i = 0;i < Delta;++i) {
                std::uniform_real_distribution<double> uniformRng(0, 1);
                Point_2 g(AreaSize * uniformRng(rng), AreaSize * uniformRng(rng));
                groups.push_back(g);

                double sigma = 1.5 * uniformRng(rng) + 0.5;//0.5-2
                sigma *= Radius;
                std::normal_distribution<double> normalRng(0, sigma);
                for (int j = 0;j < groupCnt;++j) {
                    Point_2 targetPoint(g.x() + normalRng(rng), g.y() + normalRng(rng));
                    targetRaw.push_back(targetPoint);
                }
            }
            Target = targetRaw;
#ifdef _RES_DIR
            std::string filePath = _RES_DIR + fileName;
            std::cout << "Scene file at:" << filePath << std::endl;
            std::fstream fOut(filePath, std::ios::out | std::ios::trunc);
            for (auto& tmp : targetRaw) {
                fOut << tmp.x() << ", " << tmp.y() << "\n";
            }
            fOut.close();
            std::cout << "cp " + filePath + " ./" + iterDir + fileName<< std::endl;
            system(("cp " + filePath + " ./" + iterDir + fileName).c_str());
#endif
        }
        else {
            std::cerr << "Invalid Parameter" << std::endl;
        }
        return 0;
    }

    int ProblemParam::ReadScene(const std::string& fileName) {
        try {
#ifdef _RES_DIR
            std::fstream fIn(_RES_DIR + fileName, std::ios::in);
            std::string curLine;
            double tmpx, tmpy;
            while (std::getline(fIn, curLine)) {
                std::istringstream readStr(curLine);
                std::string tmp;
                getline(readStr, tmp, ',');
                tmpx = std::stod(tmp);
                getline(readStr, tmp, ',');
                tmpy = std::stod(tmp);
                Target.push_back(Point_2(tmpx, tmpy));
            }
            fIn.close();
#endif
            std::cout << "--Read Scene: " << Target.size() << std::endl;
            for (auto& tmp : Target) {
                std::cout << tmp << std::endl;
            }
            if (Target.size() != TargetCnt) {
                throw "Scene size not match";
            }
        }
        catch (const char* e) {
            std::cerr << e << '\n';
            exit(1);
        }
        return 0;
    }


    UGVData::UGVData(double vugv, int popSize, int tour, int mutOpr, double mutRate, int crossOpr, double eliteP):Vugv(vugv), PopulationSize(popSize), Tournaments(tour), MutationOperator(mutOpr), MutationRate(mutRate), CrossoverOperator(crossOpr), EliteProbability(eliteP) { }
    UGVData::UGVData() { }
    UGVData::~UGVData() { }

    UAVData::UAVData(double timeBudget, double vuav, double timeLanding, double timeTakeoff):TimeBudget(timeBudget), Vuav(vuav), TimeLanding(timeLanding), TimeTakeoff(timeTakeoff) {
        UAVRadius = vuav * (timeBudget - timeTakeoff - timeLanding) / 2;
        problemParam.Radius = UAVRadius;
    }
    UAVData::UAVData() { }
    UAVData::~UAVData() { }

}

namespace TERRAResult {
    DataSolution::DataSolution(/* args */) { }

    DataSolution::~DataSolution() { }

    DataSolution::DataSolution(double ugvDist, double ugvTime, double uavDist, double uavTime, double totalDist, double totalTime, int stop):ugvDist(ugvDist), ugvTime(ugvTime), uavDist(uavDist), uavTime(uavTime), totalDist(totalDist), totalTime(totalTime), stop(stop) { }
    
    PathSolution::PathSolution(/* args */) { }

    PathSolution::~PathSolution() { }
}

int TERRA() {
    vector<Point_2> V1, coveredTarget, V2;
    MatrixXi setCoverTable;
    VectorXi solutionSetsLabelsV;
    // (1) Compute Voronoi Diagram
    VoronoiCoveringTimeOptimize(V1, coveredTarget);
    std::cout << "--V1: " << V1.size() << std::endl;
    for (auto& tmp : V1) {
        std::cout << tmp << std::endl;
    }
    std::cout << "--coveredTarget: " << coveredTarget.size() << std::endl;
    for (auto& tmp : coveredTarget) {
        std::cout << tmp << std::endl;
    }
    // (2) Compute Set Covering Problem
    GreedySetCovering(V1, coveredTarget, V2, setCoverTable, solutionSetsLabelsV);
    std::cout << "--V2: " << V2.size() << std::endl;
    for (auto& tmp : V2) {
        std::cout << tmp << std::endl;
    }
    std::cout << "--setCoverTable: " << setCoverTable.rows() << ", " << setCoverTable.cols() << std::endl;
    std::cout << setCoverTable << std::endl;
    std::cout << "--solutionSetsLabelsV: " << solutionSetsLabelsV.size() << std::endl;
    std::cout << solutionSetsLabelsV << std::endl;
    // Check If Home is a vertex of the solution
    vector<Point_2> V3 = CheckHome(V2);
    // std::cout << "--V3: " << V3.size() << std::endl;
    // for (auto& tmp : V3) {
    //     std::cout << tmp << std::endl;
    // }
    double ugvDist, ugvTime, uavDist, uavTime;
    int uavStop;
    std::vector<Point_2> ugvPath;
    std::vector<std::vector<Point_2>> uavPath1, uavPath2;
    Eigen::VectorXi rte;
    if (TERRAConfig::problemParam.Gp.empty()) {
        //UGV's Path without Gp
        tspGaUgv(V3, ugvDist, ugvPath, rte);
        // std::cout << "--ugvPath: " << ugvPath.size() << std::endl;
        // for (auto& tmp : ugvPath) {
        //     std::cout << tmp << std::endl;
        // }
        ugvTime = ugvDist / TERRAConfig::ugvData.Vugv;
        //UAV's Path without Gp
        UAVComputePath(coveredTarget, setCoverTable, solutionSetsLabelsV, V1, ugvPath, uavPath1, uavPath2, uavDist, uavTime, uavStop);
    }
    else {
        Point_2 avg;
        MatrixXd V3Mat(2, V3.size());
        for (int i=0;i<V3.size();++i) {
            V3Mat(0, i) = V3[i].x();
            V3Mat(1, i) = V3[i].y();
        }
        if (TERRAConfig::problemParam.Gp == "GravityCenter") {
            avg = { V3Mat.row(0).mean(), V3Mat.row(1).mean() };
        }
        else if (TERRAConfig::problemParam.Gp == "HomeCenter"){
            avg = TERRAConfig::problemParam.Home;
        }
        else if (TERRAConfig::problemParam.Gp == "MedianCenter"){
            avg = { Eigen::median(V3Mat.row(0)), Eigen::median(V3Mat.row(1)) };
        }
        std::vector<Point_2> VOpt, VRes;
        GravitationalOptimization(V1,avg,coveredTarget,solutionSetsLabelsV,setCoverTable,VOpt,VRes);
        //Genetic Algorithm to UGV Path
        tspGaUgv(VOpt, ugvDist, ugvPath, rte);
        ugvTime = ugvDist / TERRAConfig::ugvData.Vugv;
        //Search Algorithm to UAV Path
        UAVComputePath(coveredTarget, setCoverTable, solutionSetsLabelsV, VRes, ugvPath, uavPath1, uavPath2, uavDist, uavTime, uavStop);
    }
    auto totalTime = uavTime + ugvTime;
    auto totalDist = uavDist + ugvDist;
    TERRAResult::dataSol = TERRAResult::DataSolution(ugvDist, ugvTime, uavDist, uavTime, totalDist, totalTime, uavStop);
    // std::cout << "--uavPath1: " << uavPath1.size() << std::endl;
    // for (auto& tmpv : uavPath1) {
    //     for (auto& tmp1 : tmpv) {
    //         std::cout << tmp1 << " -> ";
    //     }
    //     std::cout << std::endl;
    // }
    std::cout << "--uavPath2: " << uavPath2.size() << std::endl;
    for (auto& tmpv : uavPath2) {
        for (auto& tmp1 : tmpv) {
            std::cout << tmp1 << " -> ";
        }
        std::cout << "END OF PATH" << std::endl;
    }
    BuildSolution(ugvPath, uavPath1, uavPath2);
    //TODO:visualisation
    return 0;
}

vector<Point_2> CheckHome(vector<Point_2> &V2) {
    vector<Point_2> ugvPath;
    int pos = 0;
    bool enc = false;
    for (int i = 0;i < V2.size();++i) {
        if ((V2[i] - TERRAConfig::problemParam.Home).squared_length() <= eps) {
            enc = true;
            pos = i;
        }
    }
    if (enc) {
        ugvPath = V2;
        std::swap(ugvPath[pos], ugvPath[0]);//ugvPath[pos]==Home
    }
    else {
        ugvPath.push_back(TERRAConfig::problemParam.Home);
        ugvPath.insert(ugvPath.end(), V2.begin(), V2.end());
    }
    return ugvPath;
}