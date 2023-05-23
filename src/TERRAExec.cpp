#include "TERRAUtility.h"
#include <set>
inline int initTERRAParam();
inline int TERRALaunch();
inline int CBSResultOutput(const std::string& iterDir);
inline int VisualizeResultOutput(const std::string& iterDir);
inline int CompareUAVDist();

namespace TERRAConfig {
    ConfigParam configParam;
    ProblemParam problemParam;
    UGVData ugvData;
    UAVData uavData;
}

namespace TERRAResult {
    DataSolution dataSol;
    std::vector<PathSolution> pathSol;
    std::vector<double> newdist, origindist;
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
    problemParam = ProblemParam(100, 0, 1, 1, 1, 200, "");
    ugvData = UGVData(0.4, 430, 9, 2, 0.06, 2, 2.7);
    uavData = UAVData(308, 0.5, 4, 4);
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
        CBSResultOutput(iterDir);
        VisualizeResultOutput(iterDir);
        CompareUAVDist();
    }
    for (int i = 0;i < configParam.iterations;++i) {
        std::cout << TERRAConfig::problemParam.TargetCnt << ", " << TERRAResult::origindist[i] * 2 << ", " << TERRAResult::newdist[i] << std::endl;
    }
    return 0;
}

inline int CBSResultOutput(const std::string& iterDir) {
    std::unordered_set < std::pair<int, int>, boost::hash<std::pair<int, int>>> scene, target;
    auto pathSol = TERRAResult::pathSol;
    auto GetMinimumPoint = [](const std::vector<TERRAResult::PathSolution>& paths) {
        double minx = 0, miny = 0;
        for (auto& path : paths) {
            minx = std::min(minx, path.ugv.x());
            miny = std::min(miny, path.ugv.y());
            for (auto& p : path.uav2) {
                minx = std::min(minx, p.x());
                miny = std::min(miny, p.y());
            }
        }
        return Point_2(minx, miny);
    };
    auto GetMaximumPoint = [](const std::vector<TERRAResult::PathSolution>& paths) {
        double maxx = 0, maxy = 0;
        for (auto& path : paths) {
            maxx = std::max(maxx, path.ugv.x());
            maxy = std::max(maxy, path.ugv.y());
            for (auto& p : path.uav2) {
                maxx = std::max(maxx, p.x());
                maxy = std::max(maxy, p.y());
            }
        }
        return Point_2(maxx, maxy);
    };
    Point_2 minimumPoint = GetMinimumPoint(pathSol), maximunPoint = GetMaximumPoint(pathSol);
    for (auto& tmp : TERRAConfig::problemParam.Target) {
        scene.insert(std::make_pair(int(tmp.x() - minimumPoint.x()), int(tmp.y() - minimumPoint.y())));
    }
    for (auto& path : pathSol) {
        auto tmpP = path.ugv - minimumPoint;
        path.ugv = Point_2(tmpP.x(), tmpP.y());
        for (auto &p : path.uav2) {
            auto tmpP = p - minimumPoint;
            p = Point_2(tmpP.x(), tmpP.y());
        }
    }
    for (auto& tmp : pathSol) {
        target.insert(std::make_pair(int(tmp.ugv.x()), int(tmp.ugv.y())));
    }
    YAML::Node config, map, agents;
    map["dimensions"].push_back(minimumPoint.x() < 0 ? int(std::ceil(maximunPoint.x() - minimumPoint.x())) : int(std::ceil(maximunPoint.x())));
    map["dimensions"].push_back(minimumPoint.y() < 0 ? int(std::ceil(maximunPoint.y() - minimumPoint.y())) : int(std::ceil(maximunPoint.y())));;
    for (int i = 0; i < 20;++i) {
        std::vector<int> obs(2);
        do {
            obs = { rand() % int(TERRAConfig::problemParam.AreaSize),rand() % int(TERRAConfig::problemParam.AreaSize) };
        } while (scene.count(std::make_pair(obs[0], obs[1])) or target.count(std::make_pair(obs[0], obs[1])));
        map["obstacles"].push_back(obs);
    }
    bool isdone = false;
    int idx = 0;
    while (!isdone) {
        isdone = true;
        for (int i = 0;i < pathSol.size() - 1;++i) {
            auto& path = pathSol[i];
            YAML::Node agent;
            agent["name"] = "agent" + std::to_string(i);
            if (idx < int(path.uav2.size()) - 1) {
                agent["start"].push_back(int(path.uav2[idx].x()));
                agent["start"].push_back(int(path.uav2[idx].y()));
                isdone = false;
            }
            else if (path.uav2.size()) {
                agent["start"].push_back(int(path.uav2.back().x()));
                agent["start"].push_back(int(path.uav2.back().y()));
            }
            std::vector<int> point(2);
            if (idx >= path.uav2.size() - 1) {
                agent["goal"] = YAML::Null;
            }
            else if (path.uav2.size()) {
                point = { int(path.uav2[idx + 1].x()),int(path.uav2[idx + 1].y()) };
                // agent["potentialGoals"].push_back(point);
                agent["goal"].push_back(int(path.uav2[idx + 1].x()));
                agent["goal"].push_back(int(path.uav2[idx + 1].y()));
            }
            if (agent.size()>1) {
                agents.push_back(agent);
            }
        }
        config["map"] = map;
        config["agents"] = agents;
        std::fstream fOut(iterDir + "output" + std::to_string(idx) + ".yaml", std::ios::out | std::ios::trunc);
        fOut << config;
        fOut.close();/* code */
        idx++;
        agents = YAML::Null;
    }
    
    return 0;
}

inline int VisualizeResultOutput(const std::string& iterDir) {
    YAML::Node visualize, schedule, baseline;
    int cnt = 0;
    for (int i = 0; i < TERRAResult::pathSol.size(); ++i) {
        auto& tmpUAV = TERRAResult::pathSol[i].uav2;
        auto agentName = "agent" + std::to_string(cnt);
        if (tmpUAV.empty()) {
            continue;
        }
        else {
            for (int j = 0; j < tmpUAV.size(); ++j) {

                schedule[agentName][j]["x"] = tmpUAV[j].x();
                schedule[agentName][j]["y"] = tmpUAV[j].y();
                schedule[agentName][j]["z"] = 0;
                schedule[agentName][j]["t"] = 0;
            }
            cnt++;
        }
    }
    visualize["schedule"] = schedule;
    std::fstream fOut(iterDir + "Visualize.yaml", std::ios::out | std::ios::trunc);
    fOut << visualize;
    fOut.close();
    return 0;
}

inline int CompareUAVDist() {
    // double newDist = TERRAResult::dataSol.uavDist;
    double newDist = 0;
    double originDist = 0;
    auto GetDist = [](const Point_2& a, const Point_2& b) {
        return sqrt((a - b).squared_length());
    };
    for (const auto& paths : TERRAResult::pathSol) {
        for (int i = 1;i < paths.uav2.size();++i) {
            originDist += GetDist(paths.uav2[i], paths.ugv);
            newDist += GetDist(paths.uav2[i], paths.uav2[i - 1]);
        }
    }
    TERRAResult::origindist.push_back(originDist);
    TERRAResult::newdist.push_back(newDist);
    return 0;
}
