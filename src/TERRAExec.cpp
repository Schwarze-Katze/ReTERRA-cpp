#include "TERRAUtility.h"

inline int initTERRAParam();
inline int TERRALaunch();

namespace TERRAConfig{
    ConfigParam configParam = ConfigParam(1, true, true, false, "");
    ProblemParam problemParam = ProblemParam(20, 0, 1, 0.5, 0.5, 200,"");
    UGVData ugvData = UGVData(0.4, 430, 9, 2, 0.06, 2, 2.7);
    UAVData uavData = UAVData(308, 0.5, 4, 4);
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
    // using namespace TERRAConfig;
    
    return 0;
}

inline int TERRALaunch() {
    using namespace TERRAConfig;
    std::string iterDir;
    for (int iter = 0;iter < configParam.iterations;++iter) {
        std::cout << "Start computing scenario " << iter << std::endl;
        if (configParam.iterations > 1) {
            iterDir = configParam.saveDir + "Iteration" + std::to_string(iter) + "\\";
            system(("mkdir .\\" + iterDir).c_str());
        }
        //problemParam.SceneGenerator("TestScene.in");
        problemParam.ReadScene("TestScene.in");
        TERRA();
        std::cout << "Finished computing scenario " << iter << std::endl;
    }
    return 0;
}