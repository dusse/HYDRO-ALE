//
//  SimulationManager.hpp

#ifndef SimulationManager_hpp
#define SimulationManager_hpp

#include <stdio.h>
#include <chrono>
#include <mpi.h>
#include <Python.h>
#include "../grid/GridManager.hpp"
#include "../physics/laser/LaserManager.hpp"
#include "../physics/heat/HeatManager.hpp"
#include "../misc/Logger.hpp"
#include "../output/Writer.hpp"
#include "../input/Loader.hpp"
#include "../solvers/LagrangianSolver.hpp"

#include "../common/variables/ScalarVar.hpp"
#include "../common/variables/VectorVar.hpp"

class SimulationManager {
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<LagrangianSolver> solver;
    std::shared_ptr<GridManager> gridMng;
    std::shared_ptr<LaserManager> laserMng;
    std::shared_ptr<HeatManager> heatMng;
    std::shared_ptr<Loader> loader;
    std::unique_ptr<Writer> writer;

	double computeTimeStep(int);
public:
	SimulationManager(int ac, char **av);
	void initialize();
	void initSubdomains();
	void runSimulation(int ac, char **av);
	void finilize();
	~SimulationManager();
};

#endif
