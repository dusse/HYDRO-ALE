//
//  LagrangianSolver.hpp

#ifndef LagrangianSolver_hpp
#define LagrangianSolver_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>
#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"
#include "../physics/laser/LaserManager.hpp"
#include "../physics/heat/HeatManager.hpp"
#include "Remaper.hpp"
#include "VariableNames.hpp"

#include "../grid/Node.hpp"
#include "../grid/Zone.hpp"
#include "../grid/Subzone.hpp"

const static int  SOLVE_OK   = 0;
const static int  SOLVE_FAIL = 1;

class LagrangianSolver
{
    
private:
    
    std::unique_ptr<Logger> logger;
    std::unique_ptr<Remaper> remaper;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<LaserManager> laserMgr;
    std::shared_ptr<HeatManager> heatMgr;
    
    int timeStepNum = 0;
    
    void initNodes();
    void initZones();
    void calculateForces();
    void updateNodeVelocities();
    void updateNodeCoordinates();
    void calculateInternalEnergy();
    void calculateDensity();
    void calculatePressure();
    void applyPressureBC();
    void updateBoundaryZoneVariable();
    void checkConservativeLaws();

    
public:
    LagrangianSolver(std::shared_ptr<Loader>, std::shared_ptr<GridManager>, std::shared_ptr<LaserManager>, std::shared_ptr<HeatManager> );
    void initialize();
    void initNodalMass();
    void initZonalMass();
    int solve(double);
    void finilize();
    ~LagrangianSolver();
};
#endif
