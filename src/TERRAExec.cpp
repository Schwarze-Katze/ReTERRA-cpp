#include "TERRAUtility.h"
#include <set>
inline int initTERRAParam();
inline int TERRALaunch();
inline int ResultOutput(const std::string& iterDir);

namespace TERRAConfig {
    ConfigParam configParam;
    ProblemParam problemParam;
    UGVData ugvData;
    UAVData uavData;
}

namespace TERRAResult {
    DataSolution dataSol;
    std::vector<PathSolution> pathSol;
}

int main() {
    std::cout << "hello world" << std::endl;
    initTERRAParam();
    TERRALaunch();
    return 0;
}

inline int initTERRAParam() {
    using namespace TERRAConfig;
    configParam = ConfigParam(1, true, true, false, "");
    problemParam = ProblemParam(20, 0, 1, 0.5, 0.5, 200, "GravityCenter");
    ugvData = UGVData(0.4, 430, 9, 2, 0.06, 2, 2.7);
    uavData = UAVData(308, 1, 4, 4);
    return 0;
}

inline int TERRALaunch() {
    using namespace TERRAConfig;
    std::string iterDir;
    for (int iter = 0;iter < configParam.iterations;++iter) {
        std::cout << "Start computing scenario " << iter << std::endl;
        if (configParam.iterations > 1) {
            iterDir = configParam.saveDir + "Iteration" + std::to_string(iter) + "/";
            system(("mkdir -p ./" + iterDir).c_str());
        }
        else {
            iterDir = configParam.saveDir;
        }
        problemParam.SceneGenerator("TestScene.in", iterDir);
        // problemParam.ReadScene("TestScene.in", iterDir);
        TERRA();
        std::cout << "Finished computing scenario " << iter << std::endl;
        ResultOutput(iterDir);
    }
    return 0;
}

inline int ResultOutput(const std::string& iterDir) {
    std::unordered_set < std::pair<int, int>, boost::hash<std::pair<int, int>>> scene, target;
    for (auto& tmp : TERRAConfig::problemParam.Target) {
        scene.insert(std::make_pair(int(tmp.x()), int(tmp.y())));
    }
    YAML::Node config, map, agents;
    map["dimensions"].push_back(TERRAConfig::problemParam.AreaSize);
    map["dimensions"].push_back(TERRAConfig::problemParam.AreaSize);
    for (int i = 0; i < 200;++i) {
        std::vector<int> obs(2);
        do {
            obs = { rand() % int(TERRAConfig::problemParam.AreaSize),rand() % int(TERRAConfig::problemParam.AreaSize) };
        } while (scene.count(std::make_pair(obs[0], obs[1])) or target.count(std::make_pair(obs[0], obs[1])));
        map["obstacles"].push_back(obs);
    }
    for (int i = 0;i < TERRAResult::pathSol.size();++i) {
        auto& path = TERRAResult::pathSol[i];
        YAML::Node agent;
        agent["name"] = "agent" + std::to_string(i);
        std::vector<int> point = { int(path.ugv.x()),int(path.ugv.y()) };
        agent["start"].push_back(point);
        for (auto& tmp : path.uav2) {
            point = { int(tmp.x()),int(tmp.y()) };
            agent["potentialGoals"].push_back(point);
        }
        agents.push_back(agent);
    }
    config["map"] = map;
    config["agents"] = agents;
    std::fstream fOut(iterDir + "output.yaml", std::ios::out | std::ios::trunc);
    fOut << config;
    // std::fstream fOut(iterDir + "ugv.out", std::ios::out | std::ios::trunc);
    // for (auto& tmp : TERRAResult::pathSol) {
    //     fOut << tmp.ugv.x() << (&tmp.ugv == &TERRAResult::pathSol.back().ugv ? "\n" : ", ");
    // }
    // for (auto& tmp : TERRAResult::pathSol) {
    //     fOut << tmp.ugv.y() << (&tmp.ugv == &TERRAResult::pathSol.back().ugv ? "\n" : ", ");
    // }
    // fOut.close();
    // fOut.open(iterDir + "uav.out", std::ios::out | std::ios::trunc);
    // for (auto& tmpUGV : TERRAResult::pathSol) {
    //     for (auto& tmpUAV : tmpUGV.uav2) {
    //         fOut << tmpUAV.x() << (&tmpUAV == &tmpUGV.uav2.back() ? "\n" : ", ");
    //     }
    //     for (auto& tmpUAV : tmpUGV.uav2) {
    //         fOut << tmpUAV.y() << (&tmpUAV == &tmpUGV.uav2.back() ? "\n" : ", ");
    //     }
    // }
    fOut.close();
    return 0;
}