//
// Created by chkwon on 3/23/22.
//

#include "AlgorithmParameters.h"
#include <iostream>

extern "C"
struct AlgorithmParameters default_algorithm_parameters() {
	struct AlgorithmParameters ap{};

	ap.ficEvol="evolution/";

	ap.nbGranular = 20;
	ap.mu = 25;
	ap.lambda = 40;
	ap.nbElite = 4;
	ap.nbClose = 5;

	ap.nbIterPenaltyManagement = 100;
	ap.targetFeasible = 0.2;
	ap.penaltyDecrease = 0.85;
	ap.penaltyIncrease = 1.2;

	ap.seed = 0;
	ap.nbIter = 20000;
	ap.nbIterTraces = 500;
	ap.timeLimit = 0;
	ap.useSwapStar = 1;

	ap.testmode = 0;

	ap.nombreTry=1;
    
    ap.crossoverType = "OX";
	ap.useBestMove=false;

	ap.chromTOrder = 0; // 0 = original, 1 = shuffle routes, 2 = shuffle + reverse
	
	ap.nbCut = 2;
	ap.segmentDivisor = 2;
    
    ap.eqSeg = 1;
    ap.useCostBenefit = 1;
    ap.randSelect = 1;
    ap.insertSeg = 1;
    

    

	return ap;
}

void print_algorithm_parameters(const AlgorithmParameters & ap)
{
	std::cout << "=========== Algorithm Parameters =================" << std::endl;
	std::cout << "---- nbGranular              is set to " << ap.nbGranular << std::endl;
	std::cout << "---- mu                      is set to " << ap.mu << std::endl;
	std::cout << "---- lambda                  is set to " << ap.lambda << std::endl;
	std::cout << "---- nbElite                 is set to " << ap.nbElite << std::endl;
	std::cout << "---- nbClose                 is set to " << ap.nbClose << std::endl;
	std::cout << "---- nbIterPenaltyManagement is set to " << ap.nbIterPenaltyManagement << std::endl;
	std::cout << "---- targetFeasible          is set to " << ap.targetFeasible << std::endl;
	std::cout << "---- penaltyDecrease         is set to " << ap.penaltyDecrease << std::endl;
	std::cout << "---- penaltyIncrease         is set to " << ap.penaltyIncrease << std::endl;
	std::cout << "---- seed                    is set to " << ap.seed << std::endl;
	std::cout << "---- nbIter                  is set to " << ap.nbIter << std::endl;
	std::cout << "---- nbIterTraces            is set to " << ap.nbIterTraces << std::endl;
	std::cout << "---- timeLimit               is set to " << ap.timeLimit << std::endl;
	std::cout << "---- useSwapStar             is set to " << ap.useSwapStar << std::endl;
	std::cout << "---- crossoverType           is set to " << ap.crossoverType << std::endl;
	std::cout << "---- nbCut				   is set to " << ap.nbCut << std::endl;
	std::cout << "---- segmentDivisor				   is set to " << ap.segmentDivisor << std::endl;
	std::cout << "---- eqSeg				   is set to " << ap.eqSeg << std::endl;
	std::cout << "---- useCostBenefit          is set to " << ap.useCostBenefit << std::endl;
	std::cout << "---- randSelect              is set to " << ap.randSelect << std::endl;
	std::cout << "---- insertSeg               is set to " << ap.insertSeg << std::endl;
	std::cout << "---- useBestMove               is set to " << ap.useBestMove << std::endl;
	std::cout << "---- testmode               is set to " << ap.testmode << std::endl;
	std::cout << "---- chromTOrder               is set to " << ap.chromTOrder << std::endl;
    std::cout << "==================================================" << std::endl;
}
