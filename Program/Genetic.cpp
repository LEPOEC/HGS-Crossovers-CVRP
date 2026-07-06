#include "Genetic.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <unordered_map>
#include <chrono>
#include <map>
#include <queue>

std::vector<std::vector<int>> evalDistRoute(std::vector<int>  v1, std::vector<int>  v2){

    int c;

    int n = v1.size();
    int m = v2.size();

    //D matrix (n+1) x (m+1)
    std::vector<std::vector<int> > D (n+1, std::vector<int>(m+1));

    // Init D matrix
    for (int i = 0; i < n+1; i++){
        D[i][0] = i;
    }


    for (int j = 0; j < m+1; j++){
        D[0][j] = j;
    }


    for (int i = 1; i < n+1; i++){
        for (int j = 1; j < m+1; j++){
            if (v1[i-1]  == v2[j-1]){
                c = 0;
            } else {
                c = 99999;
            }
            D[i][j] = std::min(std::min(D[i - 1][j] + 1, D[i][j - 1] + 1), D[i - 1][j - 1] + c);
        }
    }
    return D;
}


void Genetic::run()
{
        /* FONCTION CALCUL DISTANCE */
        auto edgeDistance = [&](const Individual& indiv, const Individual& parent) -> double {
            // Construire le set d'arêtes du parent
            std::set<std::pair<int,int>> parentEdges;
            for (auto& route : parent.chromR) {
                if (route.empty()) continue;
                parentEdges.insert({0, route.front()});
                for (int i = 0; i + 1 < (int)route.size(); i++) {
                    int u = route[i], v = route[i+1];
                    parentEdges.insert({std::min(u,v), std::max(u,v)});
                }
                parentEdges.insert({0, route.back()});
            }
            // Compter les arêtes de l'enfant absentes du parent
            int diff = 0, total = 0;
            for (auto& route : indiv.chromR) {
                if (route.empty()) continue;
                auto check = [&](int u, int v) {
                    total++;
                    if (!parentEdges.count({std::min(u,v), std::max(u,v)})) diff++;
                };
                check(0, route.front());
                for (int i = 0; i + 1 < (int)route.size(); i++) check(route[i], route[i+1]);
                check(0, route.back());
            }
            return (total > 0) ? (double)diff / total : 0.0;
        };

        auto plDistance = [&](const Individual& indiv, const Individual& parent) -> double {
            // Matching greedy entre routes comme dans crossoverPathRelinking
            std::vector<std::vector<int>> routesA = indiv.chromR;
            std::vector<std::vector<int>> routesB = parent.chromR;
            int nV = params.nbVehicles;

            // Matrice de distances
            std::vector<std::vector<int>> A(nV, std::vector<int>(nV, 0));
            for (int i = 0; i < nV; i++)
                for (int j = 0; j < nV; j++) {
                    auto v1 = routesA[i], v2 = routesB[j];
                    int d1 = evalDistRoute(v1, v2)[v1.size()][v2.size()];
                    std::reverse(v2.begin(), v2.end());
                    int d2 = evalDistRoute(v1, v2)[v1.size()][v2.size()];
                    A[i][j] = std::min(d1, d2);
                }

            // Matching greedy
            double totalDist = 0;
            int matched = 0;
            std::vector<bool> usedI(nV, false), usedJ(nV, false);
            for (int c = 0; c < nV; c++) {
                int minVal = INT_MAX, minI = -1, minJ = -1;
                for (int i = 0; i < nV; i++)
                    if (!usedI[i])
                        for (int j = 0; j < nV; j++)
                            if (!usedJ[j] && A[i][j] < minVal) {
                                minVal = A[i][j]; minI = i; minJ = j;
                            }
                if (minI == -1) break;
                usedI[minI] = usedJ[minJ] = true;
                totalDist += minVal;
                matched++;
            }
            return (matched > 0) ? totalDist / matched : 0.0;
        };



        /* INITIAL POPULATION */
        population.generatePopulation();


        // Periodically store the best cost found
        std::ofstream Evolutionfile,CSVfile,DistFile;
        if (params.ap.testmode) {
            Evolutionfile.open(params.ap.ficEvol+"Instance_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + params.ap.crossoverType + "_seed_" +  std::to_string(params.ap.seed)+"_"+std::to_string(params.ap.nombreTry) + ".txt");

            CSVfile.open(params.ap.ficEvol + "LS_stats_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + params.ap.crossoverType + "_seed_" + std::to_string(params.ap.seed) + "_" + std::to_string(params.ap.nombreTry) + ".csv");

            CSVfile << "iteration,crossover_type,nbMoves,nbLoopsLS," << "m1,m2,m3,m4,m5,m6,m7,m8,m9,swapStar,"<<"nbMovesRepair,nbLoopsRepair,"<<"m1r,m2r,m3r,m4r,m5r,m6r,m7r,m8r,m9r,swapStarR,"<<"time_crossover,time_LS,ratio_crossover_LS\n";

            DistFile.open(params.ap.ficEvol + "dist_stats_" + std::to_string(params.nbClients + 1) + "_"+ "crossover_type_" + params.ap.crossoverType + "_seed_" + std::to_string(params.ap.seed) + "_" + std::to_string(params.ap.nombreTry)+ "_chromTOrder_" + std::to_string(params.ap.chromTOrder) + ".csv");
            DistFile << "iteration,crossover_type,"
             << "feasible_before,feasible_after,"
             << "dist_p1p2_edge,dist_p1p2_pl,"
             << "edge_p1_before,edge_p2_before,"
             << "edge_p1_after,edge_p2_after,"
             << "pl_p1_before,pl_p2_before,"
             << "pl_p1_after,pl_p2_after\n";
        }


        const double interval = 0.2; // Intervals in s
        double lastOutputTime = (double)clock() / (double)CLOCKS_PER_SEC; // Last save in s


        // For Random crossover
        const std::vector<std::string> allCrossovers = {"OX", "LOX", "AOX", "UOX", "PMX", "EAX", "PL", "GrPX", "GOX"};
        std::uniform_int_distribution<> crossoverDistr(0, allCrossovers.size() - 1);

        std::map<std::string, int> crossoverUsageCount;
        for (const auto& c : allCrossovers)
            crossoverUsageCount[c] = 0;

        int nbIter;
        int nbIterNonProd = 1;
        if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;

        std::cout << " Nombreb Iteration " <<   params.ap.nbIter << std::endl;
        std::cout << " Type Crossover " <<   params.ap.crossoverType << std::endl;


        Individual Parent1 ;
        Individual Parent2 ;

        for (nbIter = 0 ;  (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit ; nbIter++) {
            /* SELECTION AND CROSSOVER */
            std::string activeCrossover = params.ap.crossoverType;
            Parent1=population.getBinaryTournament();
            Parent2=population.getBinaryTournament();

            auto startCrossover = std::chrono::high_resolution_clock::now();

            if (params.ap.crossoverType == "RANDOM") {
                activeCrossover = allCrossovers[crossoverDistr(params.ran)];
            }



            if (activeCrossover == "OX") {

                //Original crossover of Vidal.

                crossoverOX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "LOX") {

                // Linear order crossover
                // https://www.researchgate.net/publication/304400311_A_Study_of_Crossover_Operators_for_Genetic_Algorithm_and_Proposal_of_a_New_Crossover_Operator_to_Solve_Open_Shop_Scheduling_Problem

                crossoverLOX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "AOX") {

                // Adjust version of the OX crossover from the paper [Hvattum 2022] "Adjusting the order crossover operator for capacitated vehicle routing problems".
                //https://www.sciencedirect.com/science/article/pii/S0305054822002246

                crossoverAOX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "UOX") {

                // version of the OX crossover using bit mask
                // https://www.mecs-press.org/ijisa/ijisa-v7-n11/IJISA-V7-N11-3.pdf

                crossoverUOX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "PMX") {

                // adapted version of Alleles, loci, and the traveling salesman problem.
                // https://www.mecs-press.org/ijisa/ijisa-v7-n11/IJISA-V7-N11-3.pdf
                // https://www.researchgate.net/publication/304400311_A_Study_of_Crossover_Operators_for_Genetic_Algorithm_and_Proposal_of_a_New_Crossover_Operator_to_Solve_Open_Shop_Scheduling_Problem

                crossoverPMX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "EAX") {

                // Edge assembly crossover
                // [Nagata 2009] Edge assembly-based memetic algorithm for the capacitated vehicle routing problem
                crossoverEAX(offspring, Parent1, Parent2);


            }  else if (activeCrossover == "PL") {

                // [Sorensen 2013] Pathrelinking for CVRP : https://www.sciencedirect.com/science/article/pii/S0305054813000373

                crossoverPathRelinking(offspring, Parent1, Parent2);

            } else if (activeCrossover == "GrPX") {

                // Graph-based crossover inspired by the classical graph coloring crossover.
                // Routes are selected greedily based on a cost/benefit ratio,
                // alternating between the two parents at each step.
                crossoverGrPX(offspring, Parent1, Parent2);

            } else if (activeCrossover == "GrPX2") {

                // Variant of GrPX where Parent1 is used twice before switching to Parent2:
                // Parent1 -> Parent1 -> Parent2 -> Parent1 -> ...

                crossoverGrPX2(offspring, Parent1, Parent2);

            } else if (activeCrossover == "GOX") {

                // splits both parents' giant tours into segments,
                // evaluates each segment by cost or cost/benefit ratio, then greedily selects
                // the best segments alternating between parents. Remaining clients are fille from Parent2.

                crossoverGOX(offspring, Parent1, Parent2    );

            } else if (activeCrossover == "GLOX") {

                // finds the single best segment in Parent1
                // according to a cost/benefit score, copies it in place (like LOX),
                // then fills the remaining positions linearly from Parent2.

                crossoverGLOX(offspring, Parent1, Parent2    );
            }
            else {
                throw std::runtime_error("Erreur : crossover invalide");
            }


            auto endCrossover = std::chrono::high_resolution_clock::now();
            double timeCrossover = std::chrono::duration<double>(endCrossover - startCrossover).count();

            double edgeP1Before=0, edgeP2Before=0, plP1Before=0, plP2Before=0;
            double edgeP1After=0,  edgeP2After=0,  plP1After=0,  plP2After=0;
            double distP1P2Edge = 0;
            double distP1P2PL   =0;
            bool feasibleBefore ;

            if (params.ap.testmode) {


                double currentTime = (double)clock() / (double)CLOCKS_PER_SEC;
                if (currentTime - lastOutputTime >= interval)
                {
                    Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << currentTime<< " - "<<nbIter << std::endl;
                    lastOutputTime = currentTime;
                }

                /* Distances before LS */
                if (nbIter % 50 ==0) {
                    edgeP1Before = edgeDistance(offspring, Parent1);
                    edgeP2Before = edgeDistance(offspring, Parent2);
                    plP1Before   = plDistance(offspring, Parent1);
                    plP2Before   = plDistance(offspring, Parent2);

                    distP1P2Edge = edgeDistance(Parent1, Parent2);
                    distP1P2PL   = plDistance(Parent1, Parent2);
                    feasibleBefore = offspring.eval.isFeasible;
                }

            }

            /* LOCAL SEARCH */
            auto startLS = std::chrono::high_resolution_clock::now();
            auto [nbMoves,nbloopsLS,statsLS]=localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
            auto endLS = std::chrono::high_resolution_clock::now();
            double timeLS = std::chrono::duration<double>(endLS - startLS).count();
            if (params.ap.testmode) {


                /* Distances  after LS */
                if (nbIter % 50 ==0) {
                    edgeP1After = edgeDistance(offspring, Parent1);
                    edgeP2After = edgeDistance(offspring, Parent2);
                    plP1After   = plDistance(offspring, Parent1);
                    plP2After   = plDistance(offspring, Parent2);

                    bool feasibleAfter = offspring.eval.isFeasible;

                    DistFile << nbIter << ","
                     << activeCrossover << ","
                     << feasibleBefore << "," << feasibleAfter << ","
                     << distP1P2Edge << "," << distP1P2PL << ","
                     << edgeP1Before << "," << edgeP2Before << ","
                     << edgeP1After  << "," << edgeP2After  << ","
                     << plP1Before   << "," << plP2Before   << ","
                     << plP1After    << "," << plP2After    << "\n";
                }
            }




            bool isNewBest = population.addIndividual(offspring,true);


            int nbMovesRepair = 0, nbLoopsRepair = 0;
            LSStats statsRepair;
            double timeRepair = 0.0;


            if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
            {
                auto startRepair = std::chrono::high_resolution_clock::now();
                auto [nbMovesR, nbLoopsR, sR] =localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
                auto endRepair = std::chrono::high_resolution_clock::now();
                timeRepair = std::chrono::duration<double>(endRepair - startRepair).count();
                nbMovesRepair  = nbMovesR;
                nbLoopsRepair  = nbLoopsR;
                statsRepair    = sR;
                if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
            }

            if (params.ap.testmode) {
                double totalLS = timeLS + timeRepair;
                double ratio = (totalLS > 1e-9) ? timeCrossover / totalLS : 0.0;

                if (nbIter % 10 == 0)
                    CSVfile << nbIter << "," << activeCrossover << ","
                            << nbMoves << "," << nbloopsLS << ","
                            << statsLS.move1 << "," << statsLS.move2 << "," << statsLS.move3 << ","
                            << statsLS.move4 << "," << statsLS.move5 << "," << statsLS.move6 << ","
                            << statsLS.move7 << "," << statsLS.move8 << "," << statsLS.move9 << ","
                            << statsLS.swapStar << ","
                            << nbMovesRepair << "," << nbLoopsRepair << ","
                            << statsRepair.move1 << "," << statsRepair.move2 << "," << statsRepair.move3 << ","
                            << statsRepair.move4 << "," << statsRepair.move5 << "," << statsRepair.move6 << ","
                            << statsRepair.move7 << "," << statsRepair.move8 << "," << statsRepair.move9 << ","
                            << statsRepair.swapStar << ","
                            << timeCrossover << "," << totalLS << "," << ratio << "\n";
            }



            /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
            if (isNewBest) nbIterNonProd = 1;
            else nbIterNonProd ++ ;

            /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
            if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
            if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

            /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/

            if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
            {
            	population.restart();
            	nbIterNonProd = 1;
                if (params.ap.testmode) Evolutionfile << "Restart" << std::endl;

            }
            crossoverUsageCount[activeCrossover]++;
        }


        if (params.ap.testmode) Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        if (params.ap.testmode) {
            std::ofstream UsageFile(params.ap.ficEvol + "crossover_usage_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + params.ap.crossoverType + "_seed_" + std::to_string(params.ap.seed) + "_" + std::to_string(params.ap.nombreTry) + ".csv");
            UsageFile << "crossover_type,count\n";
            for (const auto& [name, count] : crossoverUsageCount)
                UsageFile << name << "," << count << "\n";
            UsageFile.close();
        }

        if (params.ap.testmode) {
            Evolutionfile.close();
            CSVfile.close();
            DistFile.close();
        }

}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
        // Frequency table to track the customers which have been already inserted
        std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

        // Picking the beginning and end of the crossover zone
        std::uniform_int_distribution<> distr(0, params.nbClients-1);
        int start = distr(params.ran);
        int end = distr(params.ran);

        // Avoid that start and end coincide by accident
        while (end == start) end = distr(params.ran);

        // Copy from start to end
        int j = start;
        while (j % params.nbClients != (end + 1) % params.nbClients)
        {
                result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
                freqClient[result.chromT[j % params.nbClients]] = true;
                j++;
        }

        // Fill the remaining elements in the order given by the second parent
        for (int i = 1; i <= params.nbClients; i++)
        {
                int temp = parent2.chromT[(end + i) % params.nbClients];
                if (freqClient[temp] == false)
                {
                        result.chromT[j % params.nbClients] = temp;
                        j++;
                }
        }
        // Complete the individual with the Split algorithm
        split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverLOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Ensure start is smaller than end
    if (start > end)
        std::swap(start, end);

    // Copy the crossover segment from parent1 to the offspring
    for (int i = start; i <= end; i++) {
        result.chromT[i] = parent1.chromT[i];
        freqClient[result.chromT[i]] = true;
    }

    // Fill the remaining elements using parent2
    int j = 0 ; // Start from the left
    for (int i = 0; i < params.nbClients; i++) {
        if (!freqClient[parent2.chromT[i]]) {
            if (j==start) j=end+1;                  //skip the segment of parent1
            result.chromT[j] = parent2.chromT[i];
            freqClient[parent2.chromT[i]] = true;
            j ++;
        }
    }
    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverAOX(Individual & result, const Individual & p1, const Individual & p2)
{
    // Frequency table to track customers already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start1 = distr(params.ran);
    int end1 = distr(params.ran);
    while (end1 == start1) end1 = distr(params.ran);

    // Shift zone in p2 to match final customer of zone in p1
    int start2 = start1, end2 = end1;
    while (p2.chromT[end2 % params.nbClients] != p1.chromT[end1 % params.nbClients])
        start2++, end2++;

    // Test if zone in p1 is different to zone in p2
    bool same = true;
    int size = (start1 < end1 ? end1 - start1 : params.nbClients - start1 + end1);
    for (int j = 0; j < size && same; j++) {
        if (p1.chromT[(start1 + j) % params.nbClients] != p2.chromT[(start2 + j) % params.nbClients])
            same = false;
    }

    // If same, randomize point in p2
    if (same) end2 = end2 + rand() % (params.nbClients - size);

    // Copy in place the elements from start to end
    int j = start1;
    while (j % params.nbClients != (end1 + 1) % params.nbClients) {
        result.chromT[j % params.nbClients] = p1.chromT[j % params.nbClients];
        freqClient[result.chromT[j % params.nbClients]] = true;
        j++;
    }

    // Fill the remaining elements in the order given by p2
    for (int i = 1; i <= params.nbClients; i++) {
        int temp = p2.chromT[(end2 + i) % params.nbClients];
        if (freqClient[temp] == false) {
            result.chromT[j % params.nbClients] = temp;
            j++;
        }
    }

    // Completing the individual with the Split algorithm
    split.generalSplit(result, p1.eval.nbRoutes);
}

void Genetic::crossoverUOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);


    // Generate the bits mask that will be use
    double probability=0.5;  // probability to have 1
    std::bernoulli_distribution dist(probability);

    std::vector<bool> mask(params.nbClients);
    for(int i = 0; i < params.nbClients; ++i)
    {
        mask[i] = dist(params.ran);
        // copy the element from parent1 to result
        if (mask[i]) {
            result.chromT[i] = parent1.chromT[i];
            freqClient[result.chromT[i]] = true;
        }
    }

    // Fill remaining elements from parent2
    int j = 0;
    while (j < params.nbClients && mask[j]) j++; // find the first free j
    for (int e : parent2.chromT) {
        if (!freqClient[e]) {
            result.chromT[j] = e;
            freqClient[e] = true;
            j++;
            while (j < params.nbClients && mask[j]) j++;
        }
    }
    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverPMX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Avoid that start and end coincide by accident
    while (end == start)
        end = distr(params.ran);

    if (start > end)
        std::swap(start, end);

    // Copy the segment from parent2 to the result and create mapping table
    std::unordered_map<int, int> mapping;
    for (int i = start; i <= end; ++i) {
        result.chromT[i] = parent2.chromT[i];
        mapping[parent2.chromT[i]] = parent1.chromT[i];
    }
    // Map elements from parent1 to result, avoiding elements already in the segment
    for (int i = 0; i < params.nbClients; ++i) {
        if (i < start || i > end) {
            int value = parent1.chromT[i];
            while (mapping.find(value) != mapping.end()) {
                value = mapping[value];
            }
            result.chromT[i] = value;
        }
    }
    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverEAX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    auto getEdges = [](const Individual& parent) {
        std::set<std::pair<int,int>> edges;
        for (auto& route : parent.chromR) {
            if (route.empty()) continue;
            edges.insert({0, route[0]});
            for (int i = 0; i + 1 < route.size(); i++) {
                int u = route[i], v = route[i+1];
                if (u > v) std::swap(u, v);
                edges.insert({u, v});
            }
            edges.insert({0, route.back()});
        }
        return edges;
    };



    // step 1 :  GAB = edgesA ∪ edgesB\  (edgesA ∩ edgesB)

    std::set<std::pair<int,int>> edgesA = getEdges(parent1);
    std::set<std::pair<int,int>> edgesB = getEdges(parent2);

    // separate to keep the origine for step 2
    std::set<std::pair<int,int>> GAB_A;
    std::set<std::pair<int,int>> GAB_B;
    for (auto& e : edgesA)
        if (edgesB.find(e) == edgesB.end())
            GAB_A.insert(e);
    for (auto& e : edgesB)
        if (edgesA.find(e) == edgesA.end())
            GAB_B.insert(e);


    //step 2
    // Structure pour représenter un AB-cycle
    struct ABCycle {
        std::vector<std::pair<int,int>> edges; // edge of the cycle
        std::vector<bool> fromA;              // edge comes from A
    };

    auto buildAdjacency = [](const std::set<std::pair<int,int>>& edges) {
        std::unordered_map<int, std::vector<int>> adj;
        for (auto& [u, v] : edges) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        return adj;
    };

    auto adjA = buildAdjacency(GAB_A);
    auto adjB = buildAdjacency(GAB_B);


    std::vector<ABCycle> abCycles;

    while (!GAB_A.empty()) {
        ABCycle cycle;

        // select a Random starting point for the Cycle
        auto it = GAB_A.begin();
        std::advance(it, params.ran() % GAB_A.size());
        int startNode = it->first;
        int currentNode = startNode;
        bool useA = true;

        while (true) {
            if (useA) {
                // looking for a edge from A
                auto& neighbors = adjA[currentNode];
                if (neighbors.empty()) break;

                // random selection from the available neighbors
                int idx = params.ran() % neighbors.size();
                int nextNode = neighbors[idx];


                std::pair<int,int> edge = std::make_pair(
                    std::min(currentNode, nextNode),
                    std::max(currentNode, nextNode)
                );
                cycle.edges.push_back(edge);
                cycle.fromA.push_back(true);


                GAB_A.erase(edge);
                neighbors.erase(neighbors.begin() + idx);
                auto& backNeighbors = adjA[nextNode];
                backNeighbors.erase(
                    std::find(backNeighbors.begin(), backNeighbors.end(), currentNode)
                );

                currentNode = nextNode;
            } else {
                // looking for a edge from B
                auto& neighbors = adjB[currentNode];
                if (neighbors.empty()) break;

                int idx = params.ran() % neighbors.size();
                int nextNode = neighbors[idx];

                auto edge = std::make_pair(
                    std::min(currentNode, nextNode),
                    std::max(currentNode, nextNode)
                );
                cycle.edges.push_back(edge);
                cycle.fromA.push_back(false);

                GAB_B.erase(edge);
                neighbors.erase(neighbors.begin() + idx);
                auto& backNeighbors = adjB[nextNode];
                backNeighbors.erase(
                    std::find(backNeighbors.begin(), backNeighbors.end(), currentNode)
                );

                currentNode = nextNode;
            }

            useA = !useA;

            // if start node is the same as current node then we have a cycle
            if (currentNode == startNode && cycle.edges.size() >= 2)
                break;
        }

        if (!cycle.edges.empty())
            abCycles.push_back(cycle);
    }



    // step 3
    // Si pas d'AB-cycles, fallback vers OX
    if (abCycles.empty()) {
        crossoverOX(result, parent1, parent2);
        return;
    }

    int esetIdx;

    if (params.ap.useBestMove) {
        // Choisir l'AB-cycle qui minimise le nombre de clients dans des subtours
        int bestSubtourCount = std::numeric_limits<int>::max();
        esetIdx = 0;

        for (int idx = 0; idx < (int)abCycles.size(); idx++) {
            // Simuler l'application de cet E-set sur edgesA
            std::set<std::pair<int,int>> candidateEdges = edgesA;
            for (int e = 0; e < (int)abCycles[idx].edges.size(); e++) {
                if (abCycles[idx].fromA[e])
                    candidateEdges.erase(abCycles[idx].edges[e]);
                else
                    candidateEdges.insert(abCycles[idx].edges[e]);
            }

            // Construire l'adjacence candidate
            std::unordered_map<int, std::vector<int>> adjCandidate;
            for (auto& [u, v] : candidateEdges) {
                adjCandidate[u].push_back(v);
                adjCandidate[v].push_back(u);
            }

            // BFS depuis le dépôt pour trouver les clients atteignables
            std::vector<bool> vis(params.nbClients + 1, false);
            std::queue<int> q;
            q.push(0);
            vis[0] = true;
            while (!q.empty()) {
                int cur = q.front(); q.pop();
                for (int nb : adjCandidate[cur])
                    if (!vis[nb]) { vis[nb] = true; q.push(nb); }
            }

            // Compter les clients dans des subtours
            int subtourCount = 0;
            for (int c = 1; c <= params.nbClients; c++)
                if (!vis[c]) subtourCount++;

            if (subtourCount < bestSubtourCount) {
                bestSubtourCount = subtourCount;
                esetIdx = idx;
            }
        }
    } else {
        esetIdx = params.ran() % abCycles.size(); // original
    }
    ABCycle& eset = abCycles[esetIdx];


    //step 4

    // Partir des arêtes de PA
    std::set<std::pair<int,int>> childEdges = edgesA;

    // Pour chaque arête dans l'E-set
    for (auto& e : eset.edges) {
        if (eset.fromA[&e - &eset.edges[0]]) {
            // Supprimer l'arête de PA
            childEdges.erase(e);
        } else {
            // Ajouter l'arête de PB
            childEdges.insert(e);
        }
    }


    // step 5
    // Construire la liste d'adjacence de l'enfant
    std::unordered_map<int, std::vector<int>> childAdj;
    for (auto& [u, v] : childEdges) {
        childAdj[u].push_back(v);
        childAdj[v].push_back(u);
    }




    // Identifier les routes (passent par le dépôt) et les subtours
    std::vector<std::vector<int>> routes;    // passent par dépôt (0)
    std::vector<std::vector<int>> subtours;  // ne passent pas par dépôt

    std::vector<bool> visited(params.nbClients + 1, false);

    // Reconstruire les chemins depuis la liste d'adjacence
    auto buildPath = [&](int start, int prev) {
        std::vector<int> path;
        int curr = start;
        int p = prev;
        while (true) {
            path.push_back(curr);
            visited[curr] = true;
            auto& neighbors = childAdj[curr];
            int next = -1;
            for (int n : neighbors) {
                if (n != p && !visited[n]) {
                    next = n;
                    break;
                }
            }
            if (next == -1 || next == start) break;
            p = curr;
            curr = next;
        }
        return path;
    };



    // Les routes partent du dépôt
    for (int neighbor : childAdj[0]) {
        if (!visited[neighbor]) {
            auto path = buildPath(neighbor, 0);
            routes.push_back(path);
            visited[0]=false;
        }

    }

    // Les subtours sont les cycles restants
    for (int c = 1; c <= params.nbClients; c++) {
        if (!visited[c]) {
            auto path = buildPath(c, -1);
            subtours.push_back(path);
        }
    }

    // Mélanger les subtours aléatoirement (comme dit le papier)
    std::shuffle(subtours.begin(), subtours.end(), params.ran);

    // Fusionner chaque subtour avec la meilleure route via 2-opt*
    for (auto& subtour : subtours) {
        double bestGain = std::numeric_limits<double>::infinity();
        int bestRoute = -1;
        int bestPosRoute = -1;   // position dans la route
        int bestPosSubtour = -1; // position dans le subtour

        // tester toutes les combinaisons route x subtour

        for (int r = 0; r < routes.size(); r++) {
            auto& route = routes[r];

            for (int i = 0; i < route.size(); i++) {
                // 1 arête supprimée dans la route
                int r1 = (i == 0) ? 0 : route[i-1];
                int r2 = route[i];

                for (int j = 0; j < subtour.size(); j++) {
                    // 1 arête supprimée dans le subtour
                    int s1 = subtour[j];
                    int s2 = subtour[(j+1) % subtour.size()];

                    // On supprime (r1→r2) et (s1→s2)
                    // On ajoute (r1→s1) et (s2→r2)  -- variante 1
                    // Ou        (r1→s2) et (s1→r2)  -- variante 2 (subtour inversé)
                    double oldCost = params.timeCost[r1][r2]
                                   + params.timeCost[s1][s2];

                    double newCost1 = params.timeCost[r1][s1]  // r1 connecté à début subtour
                                    + params.timeCost[s2][r2]; // fin subtour connecté à r2

                    double newCost2 = params.timeCost[r1][s2]  // r1 connecté à fin subtour
                                    + params.timeCost[s1][r2]; // début subtour connecté à r2

                    double gain = std::min(newCost1, newCost2) - oldCost;

                    if (gain < bestGain) {
                        bestGain = gain;
                        bestRoute = r;
                        bestPosRoute = i;
                        bestPosSubtour = j;
                    }
                }
            }
        }

        // Appliquer le meilleur move trouvé
        if (bestRoute != -1) {
            auto& route = routes[bestRoute];

            // Insérer le subtour dans la route à la meilleure position
            std::vector<int> newRoute;
            for (int i = 0; i <= bestPosRoute; i++)
                newRoute.push_back(route[i]);

            // Ajouter le subtour (potentiellement inversé)
            for (int i = bestPosSubtour; i < subtour.size(); i++)
                newRoute.push_back(subtour[i]);
            for (int i = 0; i < bestPosSubtour; i++)
                newRoute.push_back(subtour[i]);

            for (int i = bestPosRoute + 1; i < route.size(); i++)
                newRoute.push_back(route[i]);

            routes[bestRoute] = newRoute;
        }
    }

    for (int k = 0; k < routes.size() && k < params.nbVehicles; k++){
        routes[k].erase(
            std::remove(routes[k].begin(), routes[k].end(), 0),
            routes[k].end()
        );
        result.chromR[k] = routes[k];
    }
    // Vider les routes restantes
    for (int k = routes.size(); k < params.nbVehicles; k++){
        result.chromR[k].clear();
    }
    result.evaluateCompleteCost(params);



}

void Genetic::crossoverGrPX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;

                currentClient = customer;

            }

            routeCost[k][l] += params.timeCost[0][currentClient];


        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
        double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];
        valMax=-1;
        vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }

            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }
        }

        for (int customer : currentParent.chromR[vehMax]) {

            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;

                result.chromR[i].push_back(customer);

            }
        }


        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;


        int indice2 = (indice==0) ? 1 : 0;


        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;

            double totalDemande = 0;


            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer] ;
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;

                    //routeBenefit[indice2][l] += 1;
                    //routeBenefit[indice2][l] += params.cli[customer].demand;

                    currentClient = customer;

                }


            }

            routeCost[indice2][l] += params.timeCost[0][currentClient];

        }

        if(indice == 0){
            indice = 1;

        } else {
            indice = 0;
        }
    }
    for (int i=0; i<params.nbClients; i++){

                if (tabAffectedClients[i]<0) {
                        int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
            result.chromR[randVeh].push_back(i+1);
                }
        }

    result.evaluateCompleteCost(params);

}

class Dmatrix {
    public:
        std::vector<std::vector<int> > tab;


};

void Genetic::crossoverPathRelinking(Individual & result, const Individual & parent1, const Individual & parent2)
{

        std::vector<std::vector<int> > A (params.nbVehicles, std::vector<int>(params.nbVehicles));

        std::vector <int> affectation = std::vector <int> (params.nbVehicles);
        std::vector <int> affectationBis = std::vector <int> (params.nbVehicles);

        std::vector<std::vector<int> > ifRoadReversed (params.nbVehicles, std::vector<int>(params.nbVehicles));

        std::vector<std::vector<Dmatrix>> storeDistMatrix(params.nbVehicles, std::vector<Dmatrix>(params.nbVehicles));


        // Calculate Dmatrix and construct of the matrix A

        for (int r1 = 0; r1 < params.nbVehicles; r1++){

            for (int r2 = 0; r2 < params.nbVehicles; r2++){


               std::vector<int> v1 = parent1.chromR[r1];
               std::vector<int> v2 = parent2.chromR[r2];

               unsigned long n = v1.size();
               unsigned long m = v2.size();


               std::vector<std::vector<int> > D1 = evalDistRoute(v1, v2);

               int d1 =  D1[n][m];

               std::reverse(v2.begin(),v2.end());

               std::vector<std::vector<int> > D2 = evalDistRoute(v1, v2);

               int d2 =  D2[n][m];

                // Aij = min (d1,d2)

                if(d1 < d2){
                    A[r1][r2] = d1;
                    Dmatrix d;
                    d.tab = D1;

                    storeDistMatrix[r1][r2] = d;


                } else{
                    A[r1][r2] = d2;
                    ifRoadReversed[r1][r2] = 1;
                    Dmatrix d;
                    d.tab = D2;

                    storeDistMatrix[r1][r2] = d;
                }


            }

        }




        for (int c = 0; c < params.nbVehicles; c++){

            int minVal = 999999;
            int minI = -1;
            int minJ = -1;

            for (int i = 0; i < params.nbVehicles; i++){

                for (int j = 0; j < params.nbVehicles; j++){

                    if (A[i][j] < minVal){

                        minVal = A[i][j];
                        minI = i;
                        minJ = j;

                    }

                }

            }

            for (int k = 0; k < params.nbVehicles; k++){

                A[minI][k] = 99999;
                A[k][minJ] = 99999;

            }

            affectation[minI] = minJ;

            affectationBis[minJ] = minI;


        }


        std::vector<std::set <int>> toRemove (params.nbVehicles);

        std::vector<std::set <int>> toAdd (params.nbVehicles);

        std::vector<std::set <int>> F (params.nbVehicles);

        std::set <int> M;


        for (int k = 0; k < params.nbVehicles; k++){

            int r1 = k;
            int r2 = int(affectation[k]);

            std::vector<int> v1 = parent1.chromR[r1];
            std::vector<int> v2 = parent2.chromR[r2];

            int n = v1.size();
            int m = v2.size();



            if(ifRoadReversed[r1][r2] == 1){
                std::reverse(v2.begin(),v2.end());
            }
            
            std::vector<std::vector<int> > D = storeDistMatrix[r1][r2].tab;

            int pos1 = n;
            int pos2 = m;


            // backtracking on the Dmatrix to find the path
            while (pos1 != 0 || pos2!= 0){



                if(pos1 == 0){
                    toAdd[r2].insert(v2[pos2-1]);
                    M.insert(v2[pos2-1]);
                    pos2 = pos2 - 1;

                } else if(pos2 == 0) {
                    //M.push_back(v1[pos1-1]);
                    toRemove[r1].insert(v1[pos1-1]);
                    M.insert(v1[pos1-1]);
                   // M.push_back(0);
                    pos1 = pos1 - 1;
                } else if(D[pos1-1][pos2 - 1] == D[pos1][pos2]  && D[pos1-1][pos2 - 1] < D[pos1][pos2 - 1] && D[pos1-1][pos2-1] < D[pos1-1][pos2]){

                    F[r2].insert(v2[pos2-1]);
                    //F.push_back(0);
                    pos1 = pos1 - 1;
                    pos2 = pos2 - 1;

                } else if( D[pos1-1][pos2] < D[pos1-1][pos2 - 1]  && D[pos1-1][pos2] < D[pos1][pos2 - 1]){

                    //M.push_back(v1[pos1-1]);
                   // M.push_back(0);
                    toRemove[r1].insert(v1[pos1-1]);

                    M.insert(v1[pos1-1]);

                    pos1 = pos1 - 1;

                } else {

                    //M.push_back(v2[pos2-1]);

                    toAdd[r2].insert(v2[pos2-1]);

                    M.insert(v2[pos2-1]);


                    //M.push_back(0);
                    pos2 = pos2 - 1;

                }



            }



        }


        std::set< std::tuple<int,int,int>> relocateMoves;


        for (auto m : M) {

            for (int i = 0; i < params.nbVehicles; i++){

                if (toRemove[i].count(m)) {

                    for (int j = 0; j < params.nbVehicles; j++){

                        if (toAdd[j].count(m)) {

                            std::tuple <int, int, int> move = std::make_tuple(m, i, j);

                            relocateMoves.insert(move);
                        }

                    }

                }

            }

        }

        for (int k = 0; k < params.nbVehicles; k++){

            result.chromR[k] = parent1.chromR[k];

        }



        if (params.ap.useBestMove) {
            std::vector<std::tuple<int,int,int>> remainingMoves(relocateMoves.begin(), relocateMoves.end());
            int nbToApply = std::max(1, (int)remainingMoves.size() / 2);

            for (int k = 0; k < nbToApply && !remainingMoves.empty(); k++) {
                // Réévaluer tous les moves restants sur l'état courant
                double bestCost = std::numeric_limits<double>::max();
                int bestIdx = 0;

                for (int m = 0; m < (int)remainingMoves.size(); m++) {
                    int c  = std::get<0>(remainingMoves[m]);
                    int r1 = std::get<1>(remainingMoves[m]);
                    int r2 = std::get<2>(remainingMoves[m]);

                    auto& route1 = result.chromR[r1];
                    auto& route2 = result.chromR[r2];

                    // c est peut-être déjà parti de r1 si un move précédent l'a déplacé
                    auto itC = std::find(route1.begin(), route1.end(), c);
                    if (itC == route1.end()) continue; // move invalide, on skip

                    int idxC = itC - route1.begin();
                    int prev1 = (idxC > 0)                    ? route1[idxC-1] : 0;
                    int next1 = (idxC < (int)route1.size()-1) ? route1[idxC+1] : 0;

                    double costRemove = params.timeCost[prev1][next1]
                                      - params.timeCost[prev1][c]
                                      - params.timeCost[c][next1];

                    double bestInsert = std::numeric_limits<double>::max();
                    int prevR2 = 0;
                    for (int idx = 0; idx <= (int)route2.size(); idx++) {
                        int nextR2 = (idx < (int)route2.size()) ? route2[idx] : 0;
                        double insertCost = params.timeCost[prevR2][c]
                                          + params.timeCost[c][nextR2]
                                          - params.timeCost[prevR2][nextR2];
                        if (insertCost < bestInsert) bestInsert = insertCost;
                        prevR2 = nextR2;
                    }

                    double totalDelta = costRemove + bestInsert;
                    if (totalDelta < bestCost) {
                        bestCost = totalDelta;
                        bestIdx = m;
                    }
                }

                // Appliquer le meilleur move trouvé
                auto bestMove = remainingMoves[bestIdx];
                int c  = std::get<0>(bestMove);
                int r1 = std::get<1>(bestMove);
                int r2 = std::get<2>(bestMove);

                result.chromR[r1].erase(
                    std::remove(result.chromR[r1].begin(), result.chromR[r1].end(), c),
                    result.chromR[r1].end());

                double bestInsert = std::numeric_limits<double>::max();
                int bestPos = 0, prevR2 = 0;
                for (int idx = 0; idx <= (int)result.chromR[r2].size(); idx++) {
                    int nextR2 = (idx < (int)result.chromR[r2].size()) ? result.chromR[r2][idx] : 0;
                    double insertCost = params.timeCost[prevR2][c]
                                      + params.timeCost[c][nextR2]
                                      - params.timeCost[prevR2][nextR2];
                    if (insertCost < bestInsert) { bestInsert = insertCost; bestPos = idx; }
                    prevR2 = nextR2;
                }
                result.chromR[r2].insert(result.chromR[r2].begin() + bestPos, c);

                // Retirer le move appliqué de la liste
                remainingMoves.erase(remainingMoves.begin() + bestIdx);
            }
        } else {
            //keep 50% of the change
            unsigned long nbRelocateMovesTODO =   relocateMoves.size()/2;



            std::vector< std::tuple<int,int,int>> relocateMovesToPerform;

            std::sample(relocateMoves.begin(), relocateMoves.end(), std::back_inserter(relocateMovesToPerform), nbRelocateMovesTODO, std::mt19937{std::random_device{}()});


            std::random_shuffle(relocateMovesToPerform.begin(), relocateMovesToPerform.end());



            for (auto move : relocateMovesToPerform){

                int c = std::get<0>(move);
                int r1 = std::get<1>(move);
                int r2 = std::get<2>(move);

                result.chromR[r1].erase(std::remove(result.chromR[r1].begin(), result.chromR[r1].end(), c), result.chromR[r1].end());;

                std::vector<int> v2 = parent2.chromR[r2];

                if(ifRoadReversed[affectationBis[r2]][r2] == 1){
                    std::reverse(v2.begin(),v2.end());
                }
                int pivot = -1;

                bool notFound = true;

                int idx = 0;

                while(notFound){

                    if (F[r2].count(v2[idx])) {

                        pivot = v2[idx];
                    }

                    if(v2[idx] == c){

                        notFound = false;
                    }

                    idx = idx+1;

                }
                F[r2].insert(c);
                notFound = true;
                idx = 0;
                if(pivot != -1){
                    while(notFound){
                        if(result.chromR[affectationBis[r2]][idx] == pivot){
                            notFound = false;
                        }
                        idx = idx+1;
                    }
                }
                result.chromR[affectationBis[r2]].insert(result.chromR[affectationBis[r2]].begin() + idx, c);
            }

        }



        // Build up the rest of the Individual structure
        result.evaluateCompleteCost(params);

}

void Genetic::crossoverGrPX2(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;
                //routeBenefit[k][l] += params.timeCost[0][customer];

                currentClient = customer;




            }

            routeCost[k][l] += params.timeCost[0][currentClient];

        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    int cpt = 0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
        double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];

        valMax=-1;
        vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }


            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }

        }



        for (int customer : currentParent.chromR[vehMax]) {


            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;
                result.chromR[i].push_back(customer);
            }
        }

        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;

        int indice2 = (indice==0) ? 1: 0;

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }

        for (int l=0; l<params.nbVehicles; l++){

            int currentClient = 0;

            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer];
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;
                    //routeBenefit[indice2][l] += params.timeCost[0][customer];
                    currentClient = customer;

                }
            }
            routeCost[indice2][l] += params.timeCost[0][currentClient];
        }

        if(indice == 0 and cpt == 1){
            indice = 1;

        } else if(indice == 0 and cpt == 0){
            cpt = 1;

        } else if(indice == 1){
            indice = 0;
            cpt = 0;
        }
    }


    for (int i=0; i<params.nbClients; i++){

            if (tabAffectedClients[i]<0) {
                    int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
                result.chromR[randVeh].push_back(i+1);
            }
        }
    result.evaluateCompleteCost(params);
}

void Genetic::crossoverGOX(Individual &result, const Individual &parent1, const Individual &parent2)
{

    for (auto i =0 ; i<result.chromT.size(); i++) {
        result.chromT[i]=-1;
    }
    int nbClients = params.nbClients;
    
    // Point de coupe
    int pointCut = params.ap.nbCut;
    
//     std::cout << " pointCut " <<   pointCut << std::endl;

    // Variables booléennes pour les options de crossover
    bool segmentsEqual = false;
    if(params.ap.eqSeg == 1){
        segmentsEqual = true;
    }
    
    bool useCostBenefit = false;
    if(params.ap.useCostBenefit == 1){
        useCostBenefit = true;
    }

    bool randomSegmentSelection = false;
    if(params.ap.randSelect == 1){
        randomSegmentSelection = true;
    }
    
    int insertionSegment = params.ap.insertSeg;

    // Vérification que pointCut est valide
    if (pointCut < 2 || pointCut > nbClients) {
        throw std::invalid_argument("pointCut doit être entre 2 et nbClients.");
    }

    std::vector<int> cutPoints(pointCut);

    if (segmentsEqual) {
        std::uniform_int_distribution<> distr(0, nbClients - 1);
        std::set<int> cuts;
        while (cuts.size() < pointCut) {
            cuts.insert(distr(params.ran));
        }
        cutPoints.assign(cuts.begin(), cuts.end());
    } else {
        int step = nbClients / pointCut;
        for (int i = 0; i < pointCut; ++i) {
            cutPoints[i] = i * step;
        }
    }

    std::sort(cutPoints.begin(), cutPoints.end());

    // Initialiser un tableau pour suivre les éléments déjà ajoutés
    std::vector<bool> freqClient(nbClients + 1, false);

    // Lambda pour évaluer un segment
    auto evaluateSegment = [&](const std::vector<int>& segment) -> double {
        double totalCost = 0.0;
        double totalBenefit = 0.0;
        std::vector<bool> tempFreqClient(freqClient);
        int currentClient = 0;
        for (int customer : segment) {
            if (!tempFreqClient[customer]) {
                double cost = params.timeCost[currentClient][customer];
                double benefit = params.timeCost[0][customer] * params.cli[customer].demand;
                totalCost += cost;
                totalBenefit += benefit;
                currentClient = customer;
                tempFreqClient[customer] = true;
            }
        }
        totalCost += params.timeCost[currentClient][0];
        if (useCostBenefit) {
            if (totalCost < 1e-9) return 0.0;
            return -(totalBenefit / totalCost);  // in the following code the score is minimized. So we return the negative
        }
        return totalCost;
    };

    // Création et évaluation des segments pour les deux parents
    std::vector<std::vector<int>> segments1(pointCut);
    std::vector<std::vector<int>> segments2(pointCut);


    for (int k = 0; k < pointCut; ++k) {
        int start = cutPoints[k];
        int end = (k + 1 < pointCut) ? cutPoints[k + 1] : nbClients;
        for (int i = start; i < end; ++i) {
            segments1[k].push_back(parent1.chromT[i % nbClients]);
            segments2[k].push_back(parent2.chromT[i % nbClients]);
        }
    }

    std::vector<double> scores1(pointCut), scores2(pointCut);
    for (int i = 0; i < pointCut; ++i) {
        scores1[i] = evaluateSegment(segments1[i]);
        scores2[i] = evaluateSegment(segments2[i]);
    }

    // Sélection des segments
    std::vector<std::vector<int>> selectedSegmentsList;
    std::vector<bool> segmentSelected(pointCut, false);
    std::unordered_set<int> newlySelectedClients;

    for (int i = 0; i < pointCut; ++i) {
        int parentIndex = (i % 2 == 0) ? 0 : 1;
        int minScoreSegmentIndex = -1;
        double minScore = std::numeric_limits<double>::max();

        if (randomSegmentSelection) {
            do {
                minScoreSegmentIndex = std::rand() % pointCut;
            } while (segmentSelected[minScoreSegmentIndex]);
        } else {
            for (int j = 0; j < pointCut; ++j) {
                double score = (parentIndex == 0 ? scores1[j] : scores2[j]);
                if (!segmentSelected[j] && score < minScore) {
                    minScoreSegmentIndex = j;
                    minScore = score;
                }
            }
        }

        if (minScoreSegmentIndex == -1) break;

        auto& chosenSegment = (parentIndex == 0) ? segments1[minScoreSegmentIndex] : segments2[minScoreSegmentIndex];


        chosenSegment.erase(
         std::remove_if(
                chosenSegment.begin(),
                chosenSegment.end(),
                [&](int client) {
                    if (!freqClient[client]) {
                        freqClient[client] = true;
                        newlySelectedClients.insert(client);
                        return false; // on garde
                    }
                    return true; // on supprime
                }
            ),
            chosenSegment.end()
         );
        selectedSegmentsList.push_back(chosenSegment);
        segmentSelected[minScoreSegmentIndex] = true;

        for (int j = 0; j < pointCut; ++j) {
            if (!segmentSelected[j]) {
                bool affected = false;
                for (int client : segments1[j]) {
                    if (newlySelectedClients.find(client) != newlySelectedClients.end()) {
                        affected = true;
                        break;
                    }
                }
                if (!affected) {
                    for (int client : segments2[j]) {
                        if (newlySelectedClients.find(client) != newlySelectedClients.end()) {
                            affected = true;
                            break;
                        }
                    }
                }
                if (affected) {
                    scores1[j] = evaluateSegment(segments1[j]);
                    scores2[j] = evaluateSegment(segments2[j]);
                }
            }
        }
    }

    // Gestion de l'insertion des segments dans le résultat
    if (insertionSegment == 0) {
        // Insertion simple des segments
        int j = 0;
        for (const auto& segment : selectedSegmentsList) {
            for (int client : segment) {
                result.chromT[j++] = client;
            }
        }
        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!freqClient[client]) {
                result.chromT[j++] = client;
                freqClient[client] = true;
            }
        }
    } else if (insertionSegment == 1) {
        // Insertion aléatoire des segments
        std::shuffle(selectedSegmentsList.begin(), selectedSegmentsList.end(), std::mt19937{std::random_device{}()});
        int j = 0;
        for (const auto& segment : selectedSegmentsList) {
            for (int client : segment) {
                result.chromT[j++] = client;
            }
        }

        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!freqClient[client]) {
                result.chromT[j++] = client;
            }
        }
    } else if (insertionSegment == 2) {
        
        // Méthode 2: Insertion avancée optimisée
        // Pré-calculer les distances entre les clients
        std::vector<std::vector<double>> distances(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));
        for (int i = 0; i <= nbClients; ++i) {
            for (int j = 0; j <= nbClients; ++j) {
                distances[i][j] = params.timeCost[i][j];
            }
        }

        // Créer un tableau pour stocker les clients déjà ajoutés
        std::vector<bool> added(nbClients + 1, false);

        // Initialiser le dernier client inséré
        int lastInsertedClient = 0;
        result.chromT.clear();  // Réinitialiser le chromosome résultant

        // Priorité des segments à insérer en fonction des distances minimales
        std::vector<size_t> segmentOrder(selectedSegmentsList.size());
        std::iota(segmentOrder.begin(), segmentOrder.end(), 0); // Remplir avec 0, 1, 2, ..., numCuts - 1
        std::sort(segmentOrder.begin(), segmentOrder.end(), [&](size_t a, size_t b) {
            double distA = distances[lastInsertedClient][selectedSegmentsList[a].front()];
            double distB = distances[lastInsertedClient][selectedSegmentsList[b].front()];
            return distA < distB;
        });

        // Insérer les segments sélectionnés
        for (size_t i : segmentOrder) {
            const auto& segment = selectedSegmentsList[i];
            for (int client : segment) {
                if (!added[client]) {
                    result.chromT.push_back(client);
                    added[client] = true;
                    lastInsertedClient = client;
                }
            }
        }


        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!added[client]) {
                result.chromT.push_back(client);
                added[client] = true;
            }
        }
    }


    // Compléter l'individu avec l'algorithme Split
    split.generalSplit(result, parent1.eval.nbRoutes);

}


void Genetic::crossoverGLOX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    std::vector<bool> freqClient(params.nbClients + 1, false);
    int nbClients = params.nbClients;

    // Évaluer tous les segments possibles de parent1 par score coût/bénéfice (comme GOX)
    int bestStart = 0, bestEnd = nbClients / 2;
    double bestScore = std::numeric_limits<double>::max();

    // Taille de segment fixe = nbClients / nbCut comme dans GOX
    int segSize = nbClients / params.ap.segmentDivisor;

    for (int start = 0; start < nbClients; start++) {
        int end = start + segSize;
        if (end > nbClients) break;

        // Calcul du score coût/bénéfice du segment
        double totalCost = 0.0, totalBenefit = 0.0;
        int currentClient = 0;
        for (int i = start; i < end; i++) {
            int customer = parent1.chromT[i];
            totalCost    += params.timeCost[currentClient][customer];
            totalBenefit += params.timeCost[0][customer] * params.cli[customer].demand;
            currentClient = customer;
        }
        totalCost += params.timeCost[currentClient][0];

        double score = (params.ap.useCostBenefit && totalCost > 1e-9)
                       ? -(totalBenefit / totalCost)
                       : totalCost;

        if (score < bestScore) {
            bestScore = score;
            bestStart = start;
            bestEnd   = end;
        }
    }

    // Copier le meilleur segment de P1 (comme LOX)
    for (int i = bestStart; i < bestEnd; i++) {
        result.chromT[i] = parent1.chromT[i];
        freqClient[result.chromT[i]] = true;
    }

    // Remplissage linéaire depuis P2 (comme LOX)
    int j = 0;
    for (int i = 0; i < nbClients; i++) {
        if (!freqClient[parent2.chromT[i]]) {
            if (j == bestStart) j = bestEnd; // sauter le segment
            result.chromT[j] = parent2.chromT[i];
            freqClient[parent2.chromT[i]] = true;
            j++;
        }
    }

    split.generalSplit(result, parent1.eval.nbRoutes);
}

Genetic::Genetic(Params & params) :
        params(params),
        split(params),
        localSearch(params),
        population(params,this->split,this->localSearch),
        offspring(params){}
