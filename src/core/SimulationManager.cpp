#include "SimulationManager.hpp"
#include <iostream>
#include <string>
#include <thread>
using namespace std;
using namespace chrono;


SimulationManager::SimulationManager(int ac, char **av) {
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (ac > 1){
	    setenv("PYTHONPATH", av[1] , 1);
	}else{
        if(rank == 0){
            string msg ="[SimulationManager] Use default input file path = src/input/";
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
	    const char  *PRJ_PATH = "src/input/";
	    setenv("PYTHONPATH", PRJ_PATH , 1);
	}
	initialize();
}

void SimulationManager::initialize() {
    logger.reset(new Logger());
    loader.reset(new Loader());
    loader->load();
    gridMng.reset(new GridManager(loader));
    laserMng.reset(new LaserManager(loader, gridMng));
    heatMng.reset(new HeatManager(loader, gridMng));
    solver.reset(new LagrangianSolver(loader, gridMng, laserMng, heatMng));
    writer.reset(new Writer(loader, gridMng, laserMng, heatMng));
    
    logger->writeMsg("[SimulationManager] init...OK");
}

void SimulationManager::runSimulation(int ac, char **av) {
	logger->writeMsg("[SimulationManager] run Simulation...OK");
    auto start_time_tot = high_resolution_clock::now();
    int xRes = loader->resolution[0], yRes = loader->resolution[1];
    int maxTimeStep = loader->getMaxTimestepsNum();
    int maxTimeStep2Write = loader->getTimestepsNum2Write();
    int fileNumCount = 0;
    int i_time;
    int STOP_SIMULATION = 1;

    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for ( i_time=0; i_time < maxTimeStep; i_time++ ){
        auto start_time = high_resolution_clock::now();
        logger->writeMsg("*****************************************************");
        logger->writeMsg(string("[SimulationManager] step = "+to_string(i_time)+"; time = "+to_string(i_time*loader->getTimeStep())).c_str());

        if( i_time % maxTimeStep2Write == 0 ){
            writer->write(fileNumCount);
            fileNumCount++;
        }


        if(STOP_SIMULATION == 0){
             break;
        }
        if( solver->solve(i_time*loader->getTimeStep()) == SOLVE_FAIL ){
            string msg ="[SimulationManager] STOP SIMULATION!!!";
            logger->writeMsg(msg.c_str(), CRITICAL);
            STOP_SIMULATION = 0;
        }
        
        if(rank == 0){
            auto end_time = high_resolution_clock::now();
            string msg = "[SimulationManager] Step duration = "+to_string(duration_cast<seconds>(end_time - start_time).count())+" s";
            logger->writeMsg(msg.c_str(), INFO);
        }
    }
    if(rank == 0){
    	auto end_time_tot = high_resolution_clock::now();
        string msg = "[SimulationManager] Total duration for "+to_string(i_time)+" steps = "+to_string(duration_cast<minutes>(end_time_tot - start_time_tot).count())+" min";
        logger->writeMsg(msg.c_str(), INFO);
    }



}

void SimulationManager::initSubdomains(){

}

double SimulationManager::computeTimeStep(int tmp) {
	return tmp + 1;
}

void SimulationManager::finilize() {
	logger->writeMsg("[SimulationManager] finalize...OK");
}

SimulationManager::~SimulationManager() {
	finilize();
}

