#include "TspGaUgv.h"

int tspGaUgv(vector<Point_2>& V3, double& minDist, vector<Point_2>& UGVPath) {
    auto popSize = ugvData.PopulationSize;
    auto eliteP = ugvData.EliteProbability;
    auto muteRate = ugvData.MutationRate;
    auto tournaments = ugvData.Tournaments;
    auto mutOp = ugvData.MutationOperator;
    int lastDepot = 1;
    double seval = 30000;
    int numIter = seval / popSize;
    double elite = popSize * eliteP / 100;

    double minDist = 0;
    double elitesMut = muteRate * popSize;
    const int n = V3.size();
    MatrixXd dmat(n, n);
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < n;++j) {
            dmat(i, j) = (V3[i] - V3[j]).squared_length();
        }
    }
    dmat = dmat.sqrt().eval();

    popSize = std::lround(4 * std::ceil(double(popSize) / 4));
    numIter = std::max(1, numIter);
    
    MatrixXi pop(popSize, n);
    Eigen::VectorXi permTmp(n);
    std::iota(permTmp.begin(), permTmp.end(), 0);
    for (int i = 0;i < popSize;++i) {
        std::random_shuffle(permTmp.begin(), permTmp.end());
        pop(i, all) = permTmp;
    }
    double globalMin = 1e18;
    VectorXd totalDist(popSize);
    MatrixXd distHistory(3, numIter);
    VectorXd tmpPop(n);
    for (int iter = 1;iter < numIter;++iter) {
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
        double minDist = totalDist.minCoeff(&minIndex);
        VectorXi optRoute;
        distHistory(1, iter) = minDist;
        if (minDist < globalMin) {
            globalMin = minDist;
            optRoute = pop(minIndex, all);
        }
        distHistory(2, iter) = totalDist.mean();
        auto maxDist = totalDist.maxCoeff();
        distHistory(3, iter) = maxDist;
        //Only for demonstrating to a reviewer an issue
        if (iter > 40) {
            elitesMut = 0;
        }
        //Mutation Operator for the Genetic Algorithm
        for (int p = 0;p < elitesMut;++p) {
            VectorXd randIdx(tournaments);
            randIdx.Random(tournaments);
            randIdx = 1 + double(popSize - 1) / 2 * (randIdx.array().eval() + 0.5);
            randIdx = randIdx.eval().unaryExpr<double(*)(double)>(&std::ceil);
            VectorXi rtes = pop(randIdx(seq(0, tournaments)), all);
            VectorXd dists = totalDist(randIdx(seq(0, tournaments)));
            int minDistIdx;
            dists.minCoeff(&minDistIdx);
            VectorXd best = rtes(minDistIdx, all);
            int I = rand() % n, J = rand() % n;
            while (I == J) {
                J = rand() % n;
            }
            if (I > J) {
                std::swap(I, J);
            }
            tmpPop(0, all) = best;
            switch (mutOp)
            {
            case 1://flip
                tmpPop(0, seq(I, J)) = tmpPop(0, seq(J, I, -1)).eval();
                break;
            case 2://swap
                std::swap(tmpPop(0, I), tmpPop(0, J));
                break;
            case 3://slide
                auto tmpCell = tmpPop(0, I);
                tmpPop(0, seq(I, J - 1)) = tmpPop(0, seq(I + 1, J));
                tmpPop(0, J) = tmpCell;
                break;
            default:
                break;
            }
        }
    }

}