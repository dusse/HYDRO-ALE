//
//  Remaper.hpp

#ifndef Remaper_hpp
#define Remaper_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>
#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"
#include "../misc/Logger.hpp"

#include "VariableNames.hpp"

#include "../grid/Node.hpp"
#include "../grid/Zone.hpp"
#include "../grid/Subzone.hpp"



const  static std::string  DENSITY_RECONSTRUCTED_COEF ="dens_coef";
const  static std::string  MOMENTUM_X_DENSITY_RECONSTRUCTED_COEF ="mom_x_coef";
const  static std::string  MOMENTUM_Y_DENSITY_RECONSTRUCTED_COEF ="mom_y_coef";
const  static std::string  TOT_ERG_RECONSTRUCTED_COEF ="tot_coef";

const  static std::string  DENSITY_MIN_MAX ="dens_min_max";
const  static std::string  MOMENTUM_X_DENSITY_MIN_MAX ="mom_x_min_max";
const  static std::string  MOMENTUM_Y_DENSITY_MIN_MAX ="mom_y_min_max";
const  static std::string  TOT_ERG_MIN_MAX ="tot_min_max";

const  static std::string  DENSITY_INT_REMAP = "densintrem";
const  static std::string  MOMENTUM_X_DENSITY_INT_REMAP = "momXintrem";
const  static std::string  MOMENTUM_Y_DENSITY_INT_REMAP = "momYintrem";
const  static std::string  TOT_ERG_DENSITY_INT_REMAP ="totergintrem";


class Remaper
{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::vector<std::vector<VectorVar>> reconstructedCoefs;
    std::vector<std::vector<VectorVar>> minMaxDensities;
    std::vector<std::vector<ScalarVar>> densitiesIntegrated;
    
    void reconstructInsideOldZones();
    void integrateInsideNewZones(std::vector<VectorVar>, std::vector<VectorVar>);
    void calculateSpatialIntegrals( std::vector<std::vector<std::vector<double>>>*, std::vector<VectorVar> );
    void calculateMomentsForAllQuantities(std::vector<std::vector<double>>* , std::vector<std::vector<double>>* );
    void reconstructSolverVariables();
    void calculateAvgDerivativeInNeighborhood(std::vector<std::vector<double>>, std::vector<VectorVar>);
    void minimizeAvgDerivativeInNeighborhood(std::vector<std::vector<double>>, std::vector<VectorVar>);
    
public:
    
    Remaper(std::shared_ptr<Loader>, std::shared_ptr<GridManager>);
    void runRemapping();
    void initialize();
    void finilize();
    ~Remaper();
    
};
#endif
