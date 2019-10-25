#ifndef LaserManager_hpp
#define LaserManager_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include "../../grid/GridManager.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../../common/variables/ScalarVar.hpp"
#include "../../solvers/VariableNames.hpp"


const std::string  DIV_I         = "divi" ;

class LaserManager{


    std::string SURFACE_POSITION = "sface";//TODO add to output
    std::vector<int> surfacePosition;
//
//    //I(t_FWHM/2) = Imax/2 define t_FWHM and give τ = t_FWHM /(2 sqrt(ln(2))).
    double SQRT_LN_2 = sqrt(log(2));
    double t_FWHM = 1.0E-9;
    double laserWaveL = 1.315;//1st harmonic 1.315 μm
    double tau = 0.5*t_FWHM/SQRT_LN_2;//τ
    double totalLaserPulseErg = 5.0e-3;//10 J
    double laserSpotRadius = 150.0E-4;//300 μm for Gauss
    const double laserSpotRadiusTot = 400.0E-4;//300 μm
    const int initialFrontPosIn = 3;
    
    //material Aluminum
    double A = 27.0;//atomic mass
    double Z =  3.0;//plasma mean ion charge
    
    double Ca_1 = 0.5; // absorption coefficient for Al[the first harmonic]
//    double Ca_3 = 0.75;// absorption coefficient for Al[the third harmonic frequency]
    double Ca = Ca_1;
    
    double xf = 0.906;//root of equation 5 erf(xf) = 4
    double x0 = laserSpotRadius/xf;

    double I_max = Ca*xf*SQRT_LN_2*(totalLaserPulseErg/(t_FWHM*PI*laserSpotRadius*laserSpotRadius));
    /**
    ** 2D                               EL
    ** Imax = Ca*xf*sqrt(ln(2)) [  -------------   ]
    **                           t_FWHM * π * rf^2
    **
    **  xf ≈ 0.906 rf - laser spot radius
    **  EL - total laser pulse energy
    **  x0 = rf / ln(5)
    **/
    
//  The laser beam penetrates the material till the critical density
    double rho_critical =1.86E-3*A/(Z*laserWaveL*laserWaveL);//laser wavelength in μm

    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::vector<ScalarVar> currentDivI;
    
    void initialize();
    void updateSurfacePosition();
    //I(x,y,t) = Imax exp(−(t/τ)^2 −(x/x0)^2)
    double getIntensity2D(double, double);
public:
    LaserManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>);
    void solve(double);
    std::vector<ScalarVar> getIntensityDivergence();
};


#endif /* LaserManager_hpp */
