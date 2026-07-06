# HGS-Crossovers-CVRP: A study of crossovers for the CVRP



## Compiling the executable 

You need [`CMake`](https://cmake.org) to compile.

Build with:
```console
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles"
make bin
```
This will generate the executable file `hgs` in the `build` directory.

Test with:
```console
ctest -R bin --verbose
```

## Running the algorithm

After building the executable, try an example: 
```console
./hgs ../Instances/CVRP/X-n157-k13.vrp mySolution.sol -seed 1 -t 30
```

The following options are supported:
```
Call with: ./hgs instancePath solPath [-it nbIter] [-t myCPUtime] [-seed mySeed] [-veh nbVehicles] [-log verbose] [-crossoverType type]
[-it <int>] sets a maximum number of iterations without improvement. Defaults to 20,000
[-t <double>] sets a time limit in seconds. If this parameter is set, the code will be run iteratively until the time limit
[-seed <int>] sets a fixed seed. Defaults to 0
[-veh <int>] sets a prescribed fleet size. Otherwise a reasonable UB on the fleet size is calculated
[-round <bool>] rounding the distance to the nearest integer or not. It can be 0 (not rounding) or 1 (rounding). Defaults to 1.
[-log <bool>] sets the verbose level of the algorithm log. It can be 0 or 1. Defaults to 1.
[-crossoverType <string>] crossover operator to use. Options: OX, LOX, AOX, UOX, PMX, EAX, PL, GrPX, GrPX2, GOX, GLOX, RANDOM. Defaults to OX.

Additional Arguments:
[-nbIterTraces <int>] number of iterations between traces display during HGS execution. Defaults to 500
[-nbGranular <int>] granular search parameter, limits the number of neighbors evaluated in the RI local search. Defaults to 20
[-mu <int>] minimum population size. Defaults to 25
[-lambda <int>] number of offspring generated per generation before culling the population back to mu + nbElite. Defaults to 40
[-nbElite <int>] number of elite individuals protected from elimination. Defaults to 5
[-nbClose <int>] number of closest individuals considered when computing the diversity contribution of a solution. Defaults to 4
[-nbIterPenaltyManagement <int>] number of iterations between two penalty parameter updates. Defaults to 100
[-targetFeasible <double>] target proportion of feasible individuals in the population. Defaults to 0.2
[-penaltyIncrease <double>] multiplier applied to penalty parameters when feasibility rate falls below target. Defaults to 1.2
[-penaltyDecrease <double>] multiplier applied to penalty parameters when feasibility rate exceeds target. Defaults to 0.85
[-useBestMove <bool>] if true, selects the best possible move during crossover (e.g. best E-set in EAX, best relocate in PL). Defaults to 0
[-chromTOrder <int>] ordering mode for the giant tour: 0 = polar angle, 1 = shuffled, 2 = shuffled with random reversal. Defaults to 0
[-nbtry <long>] experiment index, used to differentiate output files across repeated runs. Defaults to 0
[-ficEvol <string>] output directory path where result and statistics files are written. Defaults to current directory
[-testmode <bool>] if true, enables collection and export of detailed statistics (move counts, distances, LS traces). Defaults to 0

GOX Parameters:
[-nbCut <int>] number of cut points used to split the giant tour into segments. Defaults to 3
[-eqSeg <int>] if 1, cut points are drawn randomly; if 0, segments are equally spaced. Defaults to 1
[-useCostBenefit <int>] if 1, segments scored by benefit/cost ratio; if 0, by raw route cost. Defaults to 1
[-randSelect <int>] if 1, segments selected randomly instead of greedily. Defaults to 0
[-insertSeg <int>] insertion mode: 0 = selection order, 1 = random order, 2 = sorted by proximity. Defaults to 0

GLOX Parameters:
[-segmentDivisor <int>] controls segment size as segSize = nbClients / segmentDivisor. Defaults to 2
[-useCostBenefit <int>] if 1, best segment selected by benefit/cost ratio; if 0, by raw route cost. Defaults to 1
```

There exist different conventions regarding distance calculations in the academic literature.
The default code behavior is to apply integer rounding, as it should be done on the X instances of Uchoa et al. (2017).
To change this behavior (e.g., when testing on the CMT or Golden instances), give a flag `-round 0`, when you run the executable.

The progress of the algorithm in the standard output will be displayed as:

``
It [N1] [N2] | T(s) [T] | Feas [NF] [BestF] [AvgF] | Inf [NI] [BestI] [AvgI] | Div [DivF] [DivI] | Feas [FeasC] [FeasD] | Pen [PenC] [PenD]
``
```
[N1] and [N2]: Total number of iterations and iterations without improvement
[T]: CPU time spent until now
[NF] and [NI]: Number of feasible and infeasible solutions in the subpopulations 
[BestF] and [BestI]: Value of the best feasible and infeasible solution in the subpopulations 
[AvgF] and [AvgI]: Average value of the solutions in the feasible and infeasible subpopulations 
[DivF] and [DivI]: Diversity of the feasible and infeasible subpopulations
[FC] and [FD]: Percentage of naturally feasible solutions in relation to the capacity and duration constraints
[PC] and [PD]: Current penalty level per unit of excess capacity and duration
```

## Code structure

The main classes containing the logic of the algorithm are the following:
* **Params**: Stores the main data structures for the method
* **Individual**: Represents an individual solution in the genetic algorithm, also provides I/O functions to read and write individual solutions in CVRPLib format.
* **Population**: Stores the solutions of the genetic algorithm into two different groups according to their feasibility. Also includes the functions in charge of diversity management.
* **Genetic**: Contains the main procedures of the genetic algorithm as well as the crossover
* **LocalSearch**: Includes the local search functions, including the SWAP* neighborhood
* **Split**: Algorithms designed to decode solutions represented as giant tours into complete CVRP solutions
* **CircleSector**: Small code used to represent and manage arc sectors (to efficiently restrict the SWAP* neighborhood)

In addition, additional classes have been created to facilitate interfacing:
* **AlgorithmParameters**: Stores the parameters of the algorithm
* **CVRPLIB** Contains the instance data and functions designed to read input data as text files according to the CVRPLIB conventions
* **commandline**: Reads the line of command
* **main**: Main code to start the algorithm
* **C_Interface**: Provides a C interface for the method

## Compiling the shared library

You can also build a shared library to call the HGS-CVRP algorithm from your code.

```console
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles"
make lib
```
This will generate the library file, `libhgscvrp.so` (Linux), `libhgscvrp.dylib` (macOS), or `hgscvrp.dll` (Windows),
in the `build` directory.

To test calling the shared library from a C code:
```console
make lib_test_c
ctest -R lib --verbose
```

## Contributing

Thank you very much for your interest in this code.
This code is still actively maintained and evolving. Pull requests and contributions seeking to improve the code in terms of readability, usability, and performance are welcome. Development is conducted in the `dev` branch. I recommend to contact me beforehand at <thibaut.vidal@polymtl.ca> before any major rework.

As a general guideline, the goal of this code is to stay **simple**, **stand-alone**, and **specialized** to the CVRP. 
Contributions that aim to extend this approach to different variants of the vehicle routing problem should usually remain in a separate repository.
Similarly, additional libraries or significant increases of conceptual complexity will be avoided. Indeed, when developing (meta-)heuristics, it seems always possible to do a bit better at the cost of extra conceptual complexity. The overarching goal of this code is to find a good trade-off between algorithm simplicity and performance.

There are two main types of contributions:
* Changes that do not impact the sequence of solutions found by the HGS algorithm when running `ctest` and testing other instances with a fixed seed.
This is visible by comparing the average solution value in the population and diversity through a test run. Such contributions include refactoring, simplification, and code optimization. Pull requests of this type are likely to be integrated more quickly.
* Changes that impact the sequence of solutions found by the algorithm.
In this case, I recommend to contact me beforehand with (i) a detailed description of the changes, (ii) detailed results on 10 runs of the algorithm for each of the 100 instances of Uchoa et al. (2017) before and after the changes, using the same termination criterion as used in [2](https://arxiv.org/abs/2012.10384).

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
- Copyright(c) 2020 Thibaut Vidal

