/*MIT License

Copyright(c) 2020 Thibaut Vidal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <iostream>
#include <string>
#include <climits>
#include "AlgorithmParameters.h"

class CommandLine
{
public:
	AlgorithmParameters ap = default_algorithm_parameters();

	int nbVeh		 = INT_MAX;		// Number of vehicles. Default value: infinity
	std::string pathInstance;		// Instance path
	
	
	std::string pathSolution;		// Solution path
	bool verbose     = true;
	bool isRoundingInteger = true;

	// Reads the line of command and extracts possible options
	CommandLine(int argc, char* argv[])
	{
		if (argc % 2 != 1|| argc > 35 || argc < 3)
		{
			std::cout << "----- NUMBER OF COMMANDLINE ARGUMENTS IS INCORRECT: " << argc << std::endl;
			display_help(); throw std::string("Incorrect line of command");
		}
		else
		{
			pathInstance = std::string(argv[1]);
            
			pathSolution = "results/"+std::string(argv[2]);
            
			for (int i = 3; i < argc; i += 2)
			{
				if (std::string(argv[i]) == "-t")
					ap.timeLimit = atof(argv[i+1]);
                else if (std::string(argv[i]) == "-crossoverType")
					ap.crossoverType  = argv[i+1];
                else if (std::string(argv[i]) == "-nbCut")
					ap.nbCut  = atoi(argv[i+1]);
                else if (std::string(argv[i]) == "-segmentDivisor")
                	ap.segmentDivisor  = atoi(argv[i+1]);
                else if (std::string(argv[i]) == "-eqSeg")
					ap.eqSeg  = atoi(argv[i+1]);
                else if (std::string(argv[i]) == "-useCostBenefit")
					ap.useCostBenefit  = atoi(argv[i+1]);
                else if (std::string(argv[i]) == "-randSelect")
					ap.randSelect  = atoi(argv[i+1]);
                else if (std::string(argv[i]) == "-insertSeg")
					ap.insertSeg  = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-it")
					ap.nbIter  = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-seed")
					ap.seed    = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-veh")
					nbVeh = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-round")
					isRoundingInteger = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-log")
					verbose = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-nbGranular")
					ap.nbGranular = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-mu")
					ap.mu = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-lambda")
					ap.lambda = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-nbElite")
					ap.nbElite = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-nbClose")
					ap.nbClose = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-nbIterPenaltyManagement")
					ap.nbIterPenaltyManagement = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-nbIterTraces")
					ap.nbIterTraces = atoi(argv[i + 1]);
				else if (std::string(argv[i]) == "-targetFeasible")
					ap.targetFeasible = atof(argv[i+1]);
				else if (std::string(argv[i]) == "-penaltyIncrease")
					ap.penaltyIncrease = atof(argv[i+1]);
				else if (std::string(argv[i]) == "-penaltyDecrease")
					ap.penaltyDecrease = atof(argv[i+1]);
				else if (std::string(argv[i]) == "-nbtry")
					ap.nombreTry=  atol( argv[i+1]);
				else if (std::string(argv[i]) == "-ficEvol")
					ap.ficEvol=   argv[i+1];
				else if (std::string(argv[i]) == "-useBestMove")
					ap.useBestMove=   (std::string(argv[i+1]) == "true" || std::string(argv[i+1]) == "1");
				else if (std::string(argv[i]) == "-chromTOrder")
					ap.chromTOrder = atoi(argv[i+1]);
				else if (std::string(argv[i]) == "-testmode")
					ap.testmode=   (std::string(argv[i+1]) == "true" || std::string(argv[i+1]) == "1");
				else
				{
					std::cout << "----- ARGUMENT NOT RECOGNIZED: " << std::string(argv[i]) << std::endl;
					display_help(); throw std::string("Incorrect line of command");
				}
			}
		}
	}

	// Printing information about how to use the code
	void display_help()
		{
	    std::cout << std::endl;
	    std::cout << "-------------------------------------------------- HGS-CVRP algorithm (2020) ----------------------------------------------------------------------" << std::endl;
	    std::cout << "Call with: ./hgs instancePath solPath [-it nbIter] [-t myCPUtime] [-seed mySeed] [-veh nbVehicles] [-log verbose] [-crossoverType type]            " << std::endl;
	    std::cout << "[-it <int>] sets a maximum number of iterations without improvement. Defaults to 20,000                                                            " << std::endl;
	    std::cout << "[-t <double>] sets a time limit in seconds. If this parameter is set the code will be run iteratively until the time limit                         " << std::endl;
	    std::cout << "[-seed <int>] sets a fixed seed. Defaults to 0                                                                                                     " << std::endl;
	    std::cout << "[-veh <int>] sets a prescribed fleet size. Otherwise a reasonable UB on the the fleet size is calculated                                           " << std::endl;
	    std::cout << "[-round <bool>] rounding the distance to the nearest integer or not. It can be 0 (not rounding) or 1 (rounding). Defaults to 1.                    " << std::endl;
	    std::cout << "[-log <bool>] sets the verbose level of the algorithm log. It can be 0 or 1. Defaults to 1.                                                        " << std::endl;
	    std::cout << "[-crossoverType <string>] crossover operator to use. Options: OX, LOX, AOX, UOX, PMX, EAX, PL, GrPX, GrPX2, GOX, GLOX, RANDOM.                   " << std::endl;
	    std::cout << std::endl;
	    std::cout << "Additional Arguments:                                                                                                                               " << std::endl;
	    std::cout << "[-nbIterTraces <int>] number of iterations between traces display during HGS execution. Defaults to 500                                            " << std::endl;
	    std::cout << "[-nbGranular <int>] granular search parameter, limits the number of neighbors evaluated in the RI local search. Defaults to 20                     " << std::endl;
	    std::cout << "[-mu <int>] minimum population size. Defaults to 25                                                                                                " << std::endl;
	    std::cout << "[-lambda <int>] number of offspring generated per generation before culling the population back to mu + nbElite. Defaults to 40                    " << std::endl;
	    std::cout << "[-nbElite <int>] number of elite individuals protected from elimination. Defaults to 5                                                              " << std::endl;
	    std::cout << "[-nbClose <int>] number of closest individuals considered when computing the diversity contribution of a solution. Defaults to 4                   " << std::endl;
	    std::cout << "[-nbIterPenaltyManagement <int>] number of iterations between two penalty parameter updates. Defaults to 100                                       " << std::endl;
	    std::cout << "[-targetFeasible <double>] target proportion of feasible individuals in the population, used to adapt penalty parameters. Defaults to 0.2          " << std::endl;
	    std::cout << "[-penaltyIncrease <double>] multiplier applied to penalty parameters when the feasibility rate falls below the target. Defaults to 1.2             " << std::endl;
	    std::cout << "[-penaltyDecrease <double>] multiplier applied to penalty parameters when the feasibility rate exceeds the target. Defaults to 0.85                " << std::endl;
	    std::cout << "[-useBestMove <bool>] if true, selects the best possible move during crossover (e.g. best E-set in EAX, best relocate in PL). Defaults to 0       " << std::endl;
	    std::cout << "[-chromTOrder <int>] ordering mode for the giant tour when exporting an individual:                                                                 " << std::endl;
	    std::cout << "                     0 = sorted by route polar angle barycenter, 1 = shuffled, 2 = shuffled with random per-route reversal. Defaults to 0          " << std::endl;
	    std::cout << "[-nbtry <long>] experiment index, used to differentiate output files across repeated runs. Defaults to 0                                           " << std::endl;
	    std::cout << "[-ficEvol <string>] output directory path where statistics files are written. Defaults to an 'evolution' directory			                      " << std::endl;
	    std::cout << "[-testmode <bool>] if true, enables collection and export of detailed statistics (move counts, distances, LS traces). Defaults to 0                " << std::endl;
	    std::cout << std::endl;
	    std::cout << "GOX Parameters:                                                                                                                                     " << std::endl;
	    std::cout << "[-nbCut <int>] number of cut points used to split the giant tour into segments in GOX. Defaults to 2                                               " << std::endl;
	    std::cout << "[-eqSeg <int>] if 1, cut points are drawn randomly (equal-chance segments); if 0, segments are equally spaced. Defaults to 1                      " << std::endl;
	    std::cout << "[-useCostBenefit <int>] if 1, segments are scored by benefit/cost ratio; if 0, scored by raw route cost. Defaults to 1                            " << std::endl;
	    std::cout << "[-randSelect <int>] if 1, segments are selected randomly instead of greedily by score. Defaults to 1                                              " << std::endl;
	    std::cout << "[-insertSeg <int>] insertion mode for selected segments: 0 = in selection order, 1 = random order, 2 = sorted by proximity. Defaults to 1         " << std::endl;
	    std::cout << std::endl;
	    std::cout << "GLOX Parameters:                                                                                                                                    " << std::endl;
	    std::cout << "[-segmentDivisor <int>] controls segment size as segSize = nbClients / segmentDivisor. Higher = shorter segments. Defaults to 2                    " << std::endl;
	    std::cout << "[-useCostBenefit <int>] if 1, the best segment is selected by benefit/cost ratio; if 0, by raw route cost. Defaults to 1                          " << std::endl;
	    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	    std::cout << std::endl;
	};
};
#endif
