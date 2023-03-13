#pragma once
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <yaml-cpp/yaml.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;


const double eps = 1e-8;

namespace TERRAConfig {

    // Configuration of the testing parameters
    //     iterations = n ? of executions with different scenarios
    //     Vrep = to launch the V - REP simulation
    //     saveDir = o.O
    class ConfigParam {
    public:
        int iterations;//n ? of executions with different scenarios
        bool printResults;
        bool saveResults;
        bool Vrep;//to launch the V - REP simulation
        std::string saveDir;
        std::string fullName;
    private:
        /* data */
    public:
        ConfigParam();
        ConfigParam(int iter, bool printAns, bool saveAns, bool vrep, const std::string& fullname);
        ~ConfigParam();
    };

    // ECU_CSURP Parameters
    //     Target = Target Points of the scenario
    //     TargetCnt = Number of target points
    //     Radius = Farthest distance the UAV can travel[m]
    //     Delta = Number of targets per group for the random map generator
    //     Home = Home location
    //     AreaSize = map area
    //     Gp = gravitational point(if null, then no GOA is applied) - Gps = 'GravityCenter', 'HomeCenter', 'MedianCenter'
    class ProblemParam {
    public:
        std::vector<Point_2> Target;//Target Points of the scenario
        double Radius;//Farthest distance the UAV can travel[m]
        int TargetCnt;//Number of target points
        int Delta;//Number of targets per group for the random map generator
        Point_2 Home;//Home location
        double AreaSize;//map area
        std::string Gp;
    private:
        /* data */
    public:
        ProblemParam();
        ProblemParam(int targetCnt, int radius, int delta, double homeX, double homeY, int area, const std::string& Gp);
        ~ProblemParam();
        int SceneGenerator(const std::string& fileName, const std::string& iterDir);
        int ReadScene(const std::string& fileName);
    };
    
    // UGV Path Planning (GA) Parameters
    //     popSize = size of the population of chromosomes.
    //     tournaments = number of chromosomes to select for the tournament selection
    //     mutOper = mutation operator (1 = Flip; 2 = Swap; 3 = Slide)
    //     mutRate = mutation rate from 0.0 to 1.0 (0 - 100 %)
    //     crossOper = crossover operator (1 = OX; 2 = CX; 3 = CBX)
    //     eliteP = elitism selection from 0 to 100 %
    class UGVData {
    public:
        double Vugv;//UGV's speed [m/s]
        int PopulationSize;//size of the population of chromosomes.
        int Tournaments;//number of chromosomes to select for the tournament selection.
        int MutationOperator;//mutation operator (1 = Flip; 2 = Swap; 3 = Slide)
        double MutationRate;//mutation rate from 0 to 100 %
        int CrossoverOperator;//crossover operator (1 = OX; 2 = CX; 3 = CBX)
        double EliteProbability;//elitism selection from 0 to 100 %
    private:
        /* data */
    public:
        UGVData();
        UGVData(double vugv, int popSize, int tour, int mutOpr, double mutRate, int crossOpr, double eliteP);
        ~UGVData();
    };

    // UAV Path Planning Parameters
    //     Tt : Total time budget[budget] budget = flight seconds
    //     Vuav : UAV's speed [m/s]
    //     Tl : Landing time[s]
    //     To : Taking off time[s]
    //     R : UAV working radius[m]
    class UAVData {
    public:
        double TimeBudget;//Total time budget[budget] budget = flight seconds
        double Vuav;//UAV's speed [m/s]
        double TimeLanding;//Tl: Landing time[s]
        double TimeTakeoff;//To : Taking off time[s]
        double UAVRadius;//UAV working radius[m]
    private:
        /* data */
    public:
        UAVData();
        UAVData(double timeBudget, double vuav, double timeLanding, double timeTakeoff);
        ~UAVData();
    };

} // namespace TERRAConfig

namespace TERRAConfig {
    extern ConfigParam configParam;
    extern ProblemParam problemParam;
    extern UGVData ugvData;
    extern UAVData uavData;
}

namespace TERRAResult {
    class DataSolution
    {
    public:
        double ugvDist;
        double ugvTime;
        double uavDist;
        double uavTime;
        double totalDist;
        double totalTime;
        int stop;
    private:
        /* data */
    public:
        DataSolution(/* args */);
        DataSolution(double ugvDist, double ugvTime, double uavDist, double uavTime, double totalDist, double totalTime, int stop);
        ~DataSolution();
    };
    

    class PathSolution
    {
    public:
        Point_2 ugv;
        std::vector<Point_2> uav1;
        std::vector<Point_2> uav2;
    private:
        /* data */
    public:
        PathSolution(/* args */);
        ~PathSolution();
    };
    
}

namespace TERRAResult{
    extern DataSolution dataSol;
    extern std::vector<PathSolution> pathSol;
}
namespace Eigen {
    class Logical {
    private:
        Index new_size;
        Array<Index, Dynamic, 1> old_inds;

    public:
        Logical(const Array<bool, Dynamic, 1>& keep): new_size(keep.count()), old_inds(new_size) {
            Index j = 0;
            for (Index i = 0; i < keep.size(); i++)
                if (keep(i))
                    old_inds(j++) = i;
            old_inds.conservativeResize(j);
            new_size = j;
        }
        Index size() const { return new_size; }
        Index operator[](Index new_ind) const { return old_inds(new_ind); }
        bool empty() const { return new_size == 0; }
    };

    template<typename Derived>
    typename Derived::Scalar median(Eigen::DenseBase<Derived>& d) {
        auto r{ d.reshaped() };
        std::sort(r.begin(), r.end());
        return r.size() % 2 == 0 ?
            r.segment((r.size() - 2) / 2, 2).mean() :
            r(r.size() / 2);
    }

    template<typename Derived>
    typename Derived::Scalar median(const Eigen::DenseBase<Derived>& d) {
        typename Derived::PlainObject m{ d.replicate(1,1) };
        return median(m);
    }
}


//  This function executes the whole solution designed and developed to solve
//  the Energy Constrained UAV and Charging Station UGV Routing Problem (ECU-CSURP).

//  The goal is to visit all the waypoints with the UAV, supporting it with
//  UGV as a moving recharging station.

//  - INPUTS -
//  cfg_params     :configuration parameters
//  problem_params :ECU_CSURP Parameters
//  map_data       :xy coordinates
//  ugv_data       :ugv parameters
//  uav_data       :uav parameters

//  - OUTPUTS -
//  data_sol :list of structures with the cost of the paths of each solution.
//   struct {
//      f_ugv_d       :Cost of the UGV as distance travelled in meters.
//      f_ugv_t       :Cost of the UGV as time in seconds.
//      f_uav_d       :Cost of the UAV as distance travelled in meters. (SEARCH_UAV_PATH)
//      f_uav_t       :Cost of the UAV as time in seconds. (SEARCH_UAV_OPERATIONS)
//      ftotal_d      :Total cost in meters.
//   }
//  path_sol :paths of both UGV-UAV to reach the solutions computed in the algorithm.
//  figures  :list of figures with the solutions of each algorithm stage
int TERRA();
