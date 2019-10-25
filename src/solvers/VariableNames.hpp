//
//  VariableNames.hpp
#ifndef VariableNames_hpp
#define VariableNames_hpp

#include <stdio.h>
#include <iostream>
#include <string>

//material Aluminum
const double A = 27.0;//atomic mass
const double Z =  3.0;//plasma mean ion charge
const double Cpt = 0.9648e12;// [erg/eV/g] pressure/temperature conversion factor
const double Tout = 0.08499;// eV
const double Tboil = 0.25;// eV for Aluminum
const double Tmelt = 0.085;// melting eV for Aluminum

const  static double  GAMMA = 1.4;
const  static double  BOUNDARY_PRESSURE = 1.0e-5;

const  std::string  NODE_TOT_PRES_FORCE  = "Fpl" ;
const  std::string  SUBZONAL_PRES_FORCE  = "Fdpl";
const  std::string  VISCOSITY_PRES_FORCE = "Fql" ;

//scalar
const  static std::string  DENSITY = "dens";
const  static std::string  PRESURE = "pres";
const  static std::string  INTERNAL_ERG = "specIntErg";


const  static std::string  DELTA_PRES  = "delta_pres";
/** SCALAR **/
const  static std::string  NODAL_MASS = "nodem";
const  static std::string  ZONAL_MASS = "zonem";
const  static std::string  SUBZONAL_MASS = "szonem";

//vector
const  static std::string  VELOCITY ="vel";

const static std::string  VELOCITY_CTS ="vel_cts";//velocity current time step
const static std::string  VELOCITY_DELTA ="vel_delta";//velocity delta

const static std::string  VELOCITY_HTS ="vel_hts";//velocity half    time step
const static std::string  VELOCITY_NTS ="vel_nts";//velocity next    time step
const static std::string  VELOCITY_OUT ="vel_out";//velocity current time step output



#endif /* VariableNames_hpp */
