//
// Created by chkwon on 3/23/22.
//

// This header file must be readable in C.

#ifndef ALGORITHMPARAMETERS_H
#define ALGORITHMPARAMETERS_H
#include <string>

struct AlgorithmParameters {

	// GLOBAL PARAMETERS
	int nbGranular;			// Granular search parameter, limits the number of moves in the RI local search
	int mu;					// Minimum population size
	int lambda;				// Number of solutions created before reaching the maximum population size (i.e., generation size)
	int nbElite;			// Number of elite individuals
	int nbClose;			// Number of closest solutions/individuals considered when calculating diversity contribution

	int nbIterPenaltyManagement;  // Number of iterations between penalty updates
	double targetFeasible;	      // Reference proportion for the number of feasible individuals, used for the adaptation of the penalty parameters
	double penaltyDecrease;	      // Multiplier used to decrease penalty parameters if there are sufficient feasible individuals
	double penaltyIncrease;	      // Multiplier used to increase penalty parameters if there are insufficient feasible individuals

	int seed;				// Random seed. Default value: 0
	int nbIter;				// Nb iterations without improvement until termination (or restart if a time limit is specified). Default value: 20,000 iterations
	int nbIterTraces;       // Number of iterations between traces display during HGS execution
	double timeLimit;		// CPU time limit until termination in seconds. Default value: 0 (i.e., inactive)
	int useSwapStar;		// Use SWAP* local search or not. Default value: 1. Only available when coordinates are provided.

	bool testmode;			// If true, enables collection and export of detailed statistics (move counts, distances, LS traces). Disable for production runs.

	long nombreTry;			// Experiment index, used to differentiate output files across repeated runs


	std::string ficEvol;	// Output directory path where statistics files are written

	// CROSSOVER PARAMETERS

	std::string crossoverType;		// Crossover operator to use.

	bool useBestMove;				// If true, selects the best possible move during crossover (e.g. best E-set in EAX, best relocate in PL) instead of a random one

	int chromTOrder ;				// Ordering mode for the giant tour when exporting an individual:
									//   0 = sorted by route polar angle barycenter
									//   1 = shuffled randomly
									//   2 = shuffled with random per-route reversal


	// GOX/ GLOX PARAMETERS
	
	int nbCut;						// Number of cut points used to split the giant tour into segments in GOX
	int segmentDivisor;				//controls the segment size as segSize = nbClients / segmentDivisor.
									// A higher value produces shorter segments, a lower value produces longer ones for GLOXss.
    
    int eqSeg;						// If 1, cut points are drawn randomly (equal-chance segments); if 0, segments are equally spaced
	int useCostBenefit;				// If 1, segments are scored by benefit/cost ratio (demand-weighted distance from depot / route cost); if 0, scored by raw route cost
	int randSelect;					// If 1, segments are selected randomly instead of greedily by score
	int insertSeg;					// Insertion mode for selected segments into the offspring:
									//   0 = in selection order
									//   1 = in random order
									//   2 = sorted by proximity to the last inserted client

    
    
};

#ifdef __cplusplus
extern "C"
#endif
struct AlgorithmParameters default_algorithm_parameters();

#ifdef __cplusplus
void print_algorithm_parameters(const AlgorithmParameters & ap);
#endif

#endif //ALGORITHMPARAMETERS_H
