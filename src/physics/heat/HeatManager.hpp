#ifndef HeatManager_hpp
#define HeatManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include "../../grid/GridManager.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../../common/variables/ScalarVar.hpp"

const double INITIAL_TEMPERATURE = 0.08499;// eV
const std::string  TEMPERATURE   = "temp";
const std::string  FLUX          = "flux";
const std::string  SINUS         = "sin" ;
const std::string  VOLUME        = "vol" ;
const std::string  BASIS         = "bas" ;
const std::string  DIV_Q         = "divq" ;


class HeatManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::vector<ScalarVar> currentTemperature;
    std::vector<VectorVar> currentFlux;
    std::vector<ScalarVar> currentDivFlux;
    void initialize();
    void initCurrentTemperature();
    void calculateDivQAndTemperature(std::vector<VectorVar>, std::vector<ScalarVar> );
    std::vector<VectorVar> getKsiCoeficients(std::vector<VectorVar> , std::vector<VectorVar> , std::vector<ScalarVar>  );
    std::vector<VectorVar> getEttaCoeficients(std::vector<VectorVar> , std::vector<VectorVar> , std::vector<ScalarVar>  );
    double getConductivity(double temp);
    std::vector<double> solveTridiagonalMatrix(std::vector<VectorVar> );
    
    
public:
    HeatManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>);
    void solve();
    std::vector<ScalarVar> getCurrentTemperature();
    std::vector<VectorVar> getCurrentFlux();
    std::vector<ScalarVar> getFluxDivergence();

};


#endif /* HeatManager_hpp */
