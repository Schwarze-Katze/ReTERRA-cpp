#include "TERRA.h"

namespace TERRAConfig {
    
    ConfigParam::ConfigParam() { }

    ConfigParam::ConfigParam(int iter, bool printAns, bool saveAns, bool vrep, char* fullname):iterations(iter), printResults(printAns), saveResults(saveAns), Vrep(vrep), fullName(fullname) {
        time_t now = time(0);
        tm* ltm = localtime(&now);
        saveDir = "Results\\Test-" + std::to_string(ltm->tm_year + 1900) + "." + std::to_string(ltm->tm_mon + 1) + "." + std::to_string(ltm->tm_mday) + "." + std::to_string(ltm->tm_hour) + "." + std::to_string(ltm->tm_min) + "\\";

        if (saveResults) {
            std::cout << ("mkdir .\\" + saveDir).c_str() << std::endl;
            system(("mkdir .\\" + saveDir).c_str());//Create directory if saveresults = true
            system("pwd");
        }
    }

    ConfigParam::~ConfigParam() { }


    ProblemParam::ProblemParam(int targetCnt, int radius, int delta, double homeX, double homeY, int area):Radius(radius), TargetCnt(targetCnt), Delta(delta), AreaSize(area) {
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
    int ProblemParam::SceneGenerator() {
        //CGAL::Random rng;

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

                double sigma = 1.5 * uniformRng(rng) + 0.5;
                sigma *= Radius;
                std::normal_distribution<double> normalRng(0, sigma);
                for (int j = 0;j < groupCnt;++j) {
                    Point_2 targetPoint(g.x() + normalRng(rng), g.y() + normalRng(rng));
                    targetRaw.push_back(targetPoint);
                }
            }
            Target = targetRaw;
        }
        else {
            std::cerr << "Invalid Parameter" << std::endl;
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

    PathSolution::PathSolution(/* args */) { }

    PathSolution::~PathSolution() { }
}

int TERRA() {
    std::vector<Point_2> V1, coveredTarget, V2;
    MatrixXi setCoverTable;
    VectorXi solutionSetsLabelsV;
    // (1) Compute Voronoi Diagram
    VoronoiCoveringTimeOptimize(V1, coveredTarget);
    // (2) Compute Set Covering Problem
    GreedySetCovering(V1, coveredTarget, V2, setCoverTable, solutionSetsLabelsV);
    // Check If Home is a vertex of the solution
    
    return 0;
}