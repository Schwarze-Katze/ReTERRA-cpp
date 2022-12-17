#include <iostream>
#include "TERRA.h"

using namespace TERRAConfig;

int initTERRAParam();
int TERRALaunch();

int main() {
    std::cout << "hello world" << std::endl;
    initTERRAParam();
    TERRALaunch();
    return 0;
}

int initTERRAParam() {
    TERRAConfig::configParam = TERRAConfig::ConfigParam(1, true, true, false, "placeholder");
    TERRAConfig::problemParam = TERRAConfig::ProblemParam(20, 0, 1, 0.5, 0.5, 200);
    TERRAConfig::ugvData = TERRAConfig::UGVData(0.4, 430, 9, 2, 0.06, 1, 2.7);
    TERRAConfig::uavData = TERRAConfig::UAVData(308, 0.5, 4, 4);
    return 0;
}

int TERRALaunch() {
    std::string iterDir;
    for (int iter = 0;iter < configParam.iterations;++iter) {
        std::cout << "Start computing scenario " << iter << std::endl;
        if (configParam.iterations > 1) {
            iterDir = configParam.saveDir + "Iteration" + std::to_string(iter) + "\\";
            system(("mkdir .\\" + iterDir).c_str());
        }
        problemParam.SceneGenerator();
        TERRA();
        std::cout << "Finished computing scenario " << iter << std::endl;
    }
    return 0;
}