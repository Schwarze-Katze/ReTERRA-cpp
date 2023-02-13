#include "TspGaUgv.h"

int tspGaUgv(vector<Point_2>& V3, double& minDist, vector<Point_2>& UGVPath, VectorXi& rte) {
    auto popSize = ugvData.PopulationSize;
    auto eliteP = ugvData.EliteProbability;
    auto muteRate = ugvData.MutationRate;
    auto tournaments = ugvData.Tournaments;
    auto mutOp = ugvData.MutationOperator;
    auto crossOp = ugvData.CrossoverOperator;
    int lastDepot = 1;
    double seval = 30000;
    int numIter = seval / popSize;
    double elite = popSize * eliteP / 100;

    double elitesMut = muteRate * popSize;
    const int n = V3.size();
    MatrixXd dmat(n, n);
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < n;++j) {
            dmat(i, j) = (V3[i] - V3[j]).squared_length();
        }
    }
    dmat = dmat.eval().cwiseSqrt();
    // std::cout << "--dmat: " << dmat.rows() << ", " << dmat.cols() << std::endl;
    // std::cout << dmat << std::endl;
    popSize = std::lround(4 * std::ceil(double(popSize) / 4));
    numIter = std::max(1, numIter);
    
    MatrixXi pop(popSize, n);
    Eigen::VectorXi permTmp(n);
    std::iota(permTmp.begin(), permTmp.end(), 0);
    for (int i = 0;i < popSize;++i) {
        std::random_shuffle(permTmp.begin(), permTmp.end());
        pop(i, all) = permTmp;
    }
    double globalMin = INFINITY;
    VectorXd totalDist(popSize);
    MatrixXd distHistory(3, numIter);
    MatrixXi newPop(int(elitesMut + elite + popSize), n);
    VectorXi optRoute;
    int newPopCnt = 0;
    for (int iter = 0;iter < numIter;++iter) {
        newPopCnt = 0;
        //Evaluate Each Population Member (Calculate Total Distance)
        for (int p = 0;p < popSize;++p) {
            double d = 0;
            if (lastDepot == 1) {
                d = dmat(pop(p, last), pop(p, 0));
            }
            else {
                d = 0;
            }
            for (int k = 1;k < n;++k) {
                d += dmat(pop(p, k - 1), pop(p, k));
            }
            totalDist(p) = d;
        }
        //Find the Best Route in the Population
        int minIndex;
        minDist = totalDist.minCoeff(&minIndex);
        
        distHistory(0, iter) = minDist;
        if (minDist < globalMin) {
            globalMin = minDist;
            optRoute = pop(minIndex, all);
        }
        distHistory(1, iter) = totalDist.mean();
        auto maxDist = totalDist.maxCoeff();
        distHistory(2, iter) = maxDist;
        //Only for demonstrating to a reviewer an issue
        if (iter > 40) {
            elitesMut = 0;
        }
        //Mutation Operator for the Genetic Algorithm
        for (int p = 0;p < elitesMut;++p) {
            VectorXi tmpPop(n);
            Eigen::ArrayXd dRandIdx(tournaments);
            dRandIdx.setRandom();
            dRandIdx = dRandIdx.eval().abs() * (popSize - 1);
            VectorXi randIdx = dRandIdx.cast<int>();
            // randIdx = randIdx.eval().unaryExpr<double(*)(double)>(&std::ceil);
            MatrixXi rtes = pop(randIdx, all);
            VectorXd dists = totalDist(randIdx);
            int minDistIdx;
            dists.minCoeff(&minDistIdx);
            VectorXi best = rtes(minDistIdx, all);
            int I = rand() % n, J = rand() % n;
            while (I == J) {
                J = rand() % n;
            }
            if (I > J) {
                std::swap(I, J);
            }
            tmpPop = best;
            switch (mutOp)
            {
            case 1://flip
                tmpPop(seq(I, J)) = tmpPop(seq(J, I, -1)).eval();
                break;
            case 2://swap
                std::swap(tmpPop(I), tmpPop(J));
                break;
            case 3://slide
            {
                auto tmpCell = tmpPop(I);
                tmpPop(seq(I, J - 1)) = tmpPop(seq(I + 1, J));
                tmpPop(J) = tmpCell;
                break;
            }
            default:
                break;
            }
            newPop(newPopCnt, all) = tmpPop;
            ++newPopCnt;
        }
        int restPop = 0;
        if (elitesMut > 0) {
            restPop = popSize - newPopCnt;
        }
        else {//elitesMut==0
            restPop = popSize;
        }
        //Elitism Selection Mechanism
        if (elite > 0) {
            restPop -= elite;
            auto tmpTotalDist = totalDist;
            for (int i = 0;i < elite;++i) {
                int minIdx = 0;
                tmpTotalDist.minCoeff(&minIdx);
                tmpTotalDist(minIdx) = INFINITY;
                newPop(newPopCnt, all) = pop(minIdx, all);
                newPopCnt++;
            }
        }
        //
        while (restPop > 0) {
            MatrixXi tmpPop;
            VectorXi randOrder(popSize);
            std::iota(randOrder.begin(), randOrder.end(), 0);
            std::random_shuffle(randOrder.begin(), randOrder.end());
            VectorXi firstParent(n), secondParent(n);
            int randIdx = tournaments + rand() % (popSize - tournaments);//interval (a,b) with the formula r = a + (b-a).*rand(N,1).
            MatrixXi rtes = pop(randOrder(seq(randIdx - tournaments, randIdx)).array(), all);
            VectorXd dists = totalDist(randOrder(seq(randIdx - tournaments, randIdx)));
            int idx = 0;
            dists.minCoeff(&idx);
            firstParent = rtes(idx, all);
            dists(idx) = INFINITY;
            dists.minCoeff(&idx);
            secondParent = rtes(idx, all);

            switch (crossOp)
            {
            case 1://Order Crossover
                OX(tmpPop, firstParent, secondParent);
                --restPop;
                break;
            case 2://Cycle Crossover
                CX(tmpPop, firstParent, secondParent);
                restPop -= 2;
                if (restPop < 0) {
                    tmpPop.conservativeResize(1, NoChange);
                }
                break;
            case 3://Order Base Crossover
                OBX(tmpPop, firstParent, secondParent);
                --restPop;
                break;
            default:
                break;
            }
            for (int i = 0;i < tmpPop.rows();++i) {
                newPop(newPopCnt, all) = tmpPop(i, all);
                ++newPopCnt;
            }
        }
        pop = newPop;
    }
    //Translate the path to a Home to Home path i.e., 1,..,..,..,1
    VectorXi idxOptRoute = optRoute.cwiseEqual(1).cast<int>();
    int idx = std::find(idxOptRoute.begin(), idxOptRoute.end(), 1) - idxOptRoute.begin();
    int len = optRoute.size();
    int cycles = len - idx + 1;
    rte.resize(optRoute.size());
    rte(seq(cycles, last)) = optRoute(seq(0, last - cycles));
    if (cycles) {
        rte(seq(0, cycles - 1)) = optRoute(lastN(cycles));
    }//rte = circshift(optRoute,cycles);
    if (lastDepot == 1) {
        rte.conservativeResize(optRoute.size() + 1);
        rte(optRoute.size()) = 0;
    }
    for (auto tmpRte : rte) {
        UGVPath.push_back(V3[tmpRte]);
    }
    return 0;
}

int OX(MatrixXi& child, const VectorXi& p1, const VectorXi& p2) {
    const int c1 = p1.size();
    const int c2 = p2.size();
    child.resize(Eigen::fix<1>, c1);
    child.fill(-1);
    int randIdx1 = rand() % c1;
    int randIdx2 = rand() % c1;
    while (randIdx1 == randIdx2) {
        randIdx2 = rand() % c1;
    }
    if (randIdx1 >= randIdx2) {
        std::swap(randIdx1, randIdx2);
    }
    child(0, seq(randIdx1, randIdx2)) = p1(seq(randIdx1, randIdx2));
    for (int i = 0;i < c2;++i) {
        Eigen::Logical k(child(0, all).array() == p2(i));
        if (k.empty()) {
            int pos = 0;
            for (int j = 0;pos == 0;++j) {
                if (child(0, j) == -1) {
                    pos = j;
                }
            }
            child(0, pos) = p2(i);
        }
    }
    return 0;
}

int CX(MatrixXi& child, VectorXi p1, VectorXi p2) {
    const int c1 = p1.size();
    const int c2 = p2.size();
    child.resize(Eigen::fix<2>, c1);
    MatrixXi totalCycles(0, c1);
    bool isCompleted = false;
    while (!isCompleted) {
        MatrixXi cycle(2, c1);
        cycle.fill(-1);
        //Find first element of the new parent
        int r = 0;
        while (p1(r) == -1) {
            ++r;
        }
        auto findCycle = [&](auto&& self, int curVal, int initVal)->MatrixXi {
            Logical idx(p1.array() == curVal);
            cycle(0, idx[0]) = curVal;
            curVal = p2(idx[0]);
            cycle(1, idx[0]) = curVal;
            int t = cycle.cols();
            if (t == 1) {
                initVal = cycle(0, idx[0]);
            }
            else if (t > 1) {
                if (initVal != p2(idx[0])) {
                    cycle = self(self, curVal, initVal);
                }
            }
            return cycle;
        };
        cycle = findCycle(findCycle, p1(r), p1(r));
        std::vector<int> newP1, newP2;
        for (int i = 0;i < c1;++i) {
            auto idx(cycle.row(0).array() == p1(i));
            newP1.push_back(idx.any() ? -1 : p1(i));
        }
        for (int i = 0;i < c2;++i) {
            auto idx(cycle.row(1).array() == p2(i));
            newP2.push_back(idx.any() ? -1 : p2(i));
        }
        p1 = Eigen::Map<VectorXi>(newP1.data(), newP1.size());
        p2 = Eigen::Map<VectorXi>(newP2.data(), newP2.size());
        totalCycles.conservativeResize(totalCycles.rows() + 2, NoChange);
        totalCycles(lastN(2), all) = cycle;
        if ((p1.array() == -1).all()) {
            isCompleted = true;
        }
    }
    //Filling the new chromosomes
    int c = 0;
    for (int i = 0;i < totalCycles.rows() / 2;++i) {
        if (i % 2) {
            for (int j = 0;j < totalCycles.cols();++j) {
                if (totalCycles(i + c, j) != -1) {
                    child(0, j) = totalCycles(i + c, j);
                }
                if (totalCycles(i + c + 1, j) != -1) {
                    child(1, j) = totalCycles(i + c + 1, j);
                }
            }
        }
        else {
            for (int j = 0;j < totalCycles.cols();++j) {
                if (totalCycles(i + c, j) != -1) {
                    child(1, j) = totalCycles(i + c, j);
                }
                if (totalCycles(i + c + 1, j) != -1) {
                    child(0, j) = totalCycles(i + c + 1, j);
                }
            }
        }
        ++c;
    }
    return 0;
}

int OBX(MatrixXi& child, const VectorXi& p1, const VectorXi& p2) {
    const int c1 = p1.size();
    const int c2 = p2.size();
    child.resize(Eigen::fix<1>, c1);
    child.fill(-1);
    //Step 1: Random indexs
    int set = rand() % c1;
    VectorXi selectedIdx(set);
    std::iota(selectedIdx.begin(), selectedIdx.end(), 0);
    std::random_shuffle(selectedIdx.begin(), selectedIdx.end());
    VectorXi p1Selected = p1(selectedIdx);
    //Step 2: Copy parent2 to child
    for (int i = 0;i < c1;++i) {
        if (!(p1Selected.array() == p2(i)).any()) {
            child(0, i) = p2(i);
        }
    }
    //Step 3: Copy parent1 to child
    for (int i = 0;i < c1;++i) {
        if ((p1Selected.array() == p1(i)).any()) {
            for (int j = 0;j < c1;++j) {
                if (child(0, j) == -1) {
                    child(0, j) = p1(i);
                    break;
                }
            }
        }
    }
    return 0;
}
