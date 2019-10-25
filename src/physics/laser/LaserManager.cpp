#include "LaserManager.hpp"

using namespace std;

LaserManager::LaserManager(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr):loader(move(load)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());
    initialize();
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        string msg ="[LaserManager] rho_critical = "+to_string(rho_critical);
        logger->writeMsg(msg.c_str(), INFO);
    }
    logger->writeMsg("[LaserManager] create...OK");
}


void LaserManager::initialize(){
    for (int idx=0; idx < gridMgr->getTotalZoneNumber(); idx++){
        currentDivI.push_back(ScalarVar(DIV_I, 0.0));
    }

    for (int i = 0; i < loader->resolution[0]+1; i++){
        surfacePosition.push_back(initialFrontPosIn);
    }
}

double LaserManager::getIntensity2D(double x, double t){
    return I_max*exp(-pow((t/tau),2)-pow(x/x0,2));
}
/***The components of laser intensity I(x,y) are projected on the normals to the edge in the center of each edge
 as either Iξ[i,j+1/2] or I[ηi+1/2,j]

 special treatment is used to get values of laser intensity at each edge.
 The density at nodes ρ[i,j] is obtained by interpolation (weighted by subzonal volumes)
 from two neighboring zones on right ρ[i+1/2,j+1/2] and ρ[i+1/2,j−1/2].
 If both these densities are subcritical the laser penetrates 
 to the node i, j and thus density ρ[i,j] should be also subcritical.
 **/



/***
 div(I) = 0 if density in all four nodes of the cell is either subcritical or supercritical.
 Otherwise div(I) = 1/Vz (Iξ SξL −Iξ SξL +Iη SηL −Iη SηL)
 **/
/*
 The divergence of the laser intensity in cell c is zero, 
 if all four nodal densities are either subcritical or supercritical.
 If the values are mixed (some sub- and supercritical),
 we set the cell divergence to the value of the integral of the divergence over cell c divided by its volume,
 and after applying the Green formula, to the integral over the cell boundary ∂c. 
 It is evaluated as the sum (over all cell edges) of the edge intensities Ie 
 projected to the direction of the edge outer normal, 
 multiplied by the subcritical edge length Ls(e), and divided by the cell volume.
 */

void LaserManager::solve(double time){
    if(time > t_FWHM){
        for (int idx=0; idx < gridMgr->getTotalZoneNumber(); idx++){
            currentDivI[idx] = ScalarVar(DIV_I, 0.0);
        }
        return;
    }
    
    vector<ScalarVar> divI;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int xResExt = xRes+2, yResExt = yRes+2, zResExt = 1;
    int i,j,k=0, ijZone;
    
    vector<VectorVar> side;
    vector<ScalarVar> vol;
    vector<VectorVar> sin;// with horizontal line for ksi and etta
    
    int i1jNode, ij1Node;
    int ijNode, i1j1Node;
    double Xij, Xi1j, Xij1, Xi1j1;
    double Yij, Yi1j, Yij1, Yi1j1;
    
    for(int tmp = 0; tmp < xResExt*yResExt; tmp++){
        vol.push_back(ScalarVar("vol", 0.0));
        sin.push_back(VectorVar("sin", {0.0, 0.0}));
        side.push_back(VectorVar("bas", {0.0, 0.0}));
    }
    double laserSpotRadiusOld = laserSpotRadius;
    vector<VectorVar> coords = gridMgr->getExtendedGridCoordinates();
    double ksi_ij, etta_ij, etta_ij1, ksi_i1j, volij;
    double sinijKsi, sinijEtta, sini1jKsi, sinij1Etta;
    for ( i=0; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            ijZone   = IDX(i  , j  , k, xRes+1 , yRes+1 , zRes);
            ijNode   = IDX(i  , j  , k, xResExt, yResExt, zResExt);
            i1jNode  = IDX(i+1, j  , k, xResExt, yResExt, zResExt);
            ij1Node  = IDX(i  , j+1, k, xResExt, yResExt, zResExt);
            i1j1Node = IDX(i+1, j+1, k, xResExt, yResExt, zResExt);
            
            Xij   = coords[ijNode].getValue()[0];
            Yij   = coords[ijNode].getValue()[1];
            Xi1j  = coords[i1jNode].getValue()[0];
            Yi1j  = coords[i1jNode].getValue()[1];
            Xij1  = coords[ij1Node].getValue()[0];
            Yij1  = coords[ij1Node].getValue()[1];
            Xi1j1 = coords[i1j1Node].getValue()[0];
            Yi1j1 = coords[i1j1Node].getValue()[1];
            
            ksi_ij  = sqrt(pow((Xi1j - Xij),2) + pow((Yi1j - Yij),2));// sξij = √{[x(i+1,j)-x(i,j)]^2+[y(i+1,j)-y(i,j)]^2}
            sinijKsi  = (Xi1j - Xij)/ksi_ij;
            etta_ij = sqrt(pow((Xij1 - Xij),2) + pow((Yij1 - Yij),2));// sηij = √{[x(i,j+1)-x(i,j)]^2+[y(i,j+1)-y(i,j)]^2}
            sinijEtta  = (Xij1 - Xij)/etta_ij;
            side[ijNode]  = VectorVar("bas", { ksi_ij, etta_ij });
            volij = 0.5 * ((Xi1j1 - Xij)*(Yij1 - Yi1j) - (Xij1 - Xi1j)*(Yi1j1 - Yij));//ij
            vol[ijNode]  = ScalarVar("vol", volij);
            sin[ijNode]  = VectorVar("sin", { sinijKsi, sinijEtta });
            if( volij < 0.0 ){
                string msg ="[LaserManager] volume problem volij = "+to_string(volij);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }
        }
    }
    
    updateSurfacePosition();
    
    double fluxi1jKsi, fluxijKsi;
    double fluxij1Etta, fluxijEtta;
    double fluxDiv;
    for ( i = 0; i < xRes+1; i++){
        for( j = 0; j < surfacePosition[i]+2; j++ ){
            ijZone = IDX(i, j, k, xRes+1, yRes+1, zRes);
            
            ijNode   = IDX(i  , j  , k, xResExt, yResExt, zResExt);
            i1jNode  = IDX(i+1, j  , k, xResExt, yResExt, zResExt);
            ij1Node  = IDX(i  , j+1, k, xResExt, yResExt, zResExt);
            
            Xij   = coords[ijNode].getValue()[0];
            Xi1j  = coords[i1jNode].getValue()[0];
            Xij1  = coords[ij1Node].getValue()[0];
            
            if(abs(Xij) > laserSpotRadiusTot){
                fluxDiv = 0.0;
                currentDivI[ijZone] = ScalarVar(DIV_I, Ca*fluxDiv);
                continue;
            }
            ksi_ij  = side[ijNode].getValue()[0];
            etta_ij = side[ijNode].getValue()[1];
            ksi_i1j = side[i1jNode].getValue()[0];
            etta_ij1  = side[ij1Node].getValue()[1];
            volij  = vol[ijNode].getValue();
            sinijKsi  = sin[ijNode].getValue()[0];
            sinijEtta = sin[ijNode].getValue()[1];
            sini1jKsi  = sin[ijNode].getValue()[0];
            sinij1Etta = sin[ijNode].getValue()[1];
            
            laserSpotRadius = laserSpotRadius*(1-j/6);
            
            fluxijKsi   = getIntensity2D( 0.5*(Xij+Xi1j) , time)*sinijKsi;
            fluxijEtta  = getIntensity2D( 0.5*(Xij+Xij1) , time)*sinijEtta;

//            if( j == surfacePosition[i]+1){
                fluxi1jKsi  = getIntensity2D(  0.5*(Xij+Xi1j), time)*sini1jKsi;
                fluxij1Etta = getIntensity2D(  0.5*(Xij+Xij1), time)*sinij1Etta;
//            }else{
//                fluxi1jKsi  = 0.0;
//                fluxij1Etta = 0.0;
//            }
            /*
             *    Ω(i,j) x  div(I)  = - ( Iξ(i+1,j) x sξ(i+1,j) - Iξ(i,j) x sξ(i,j) + Iη(i,j+1) x sη(i,j+1) - Iη(i,j) x sη(i,j) )
             *
             */
//            fluxDiv = (fluxi1jKsi*ksi_i1j - fluxijKsi*ksi_ij + fluxij1Etta*etta_ij1 - fluxijEtta*etta_ij)/volij;
            
            fluxDiv = fluxi1jKsi*ksi_i1j/volij;
           
            currentDivI[ijZone] = ScalarVar(DIV_I, Ca*fluxDiv);
        }
    }
    laserSpotRadius = laserSpotRadiusOld;
}

void LaserManager::updateSurfacePosition(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k=0, ij1Zone;
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    double dens;
    for ( i = 0; i < xRes+1; i++ ){
        j = surfacePosition[i];
        ij1Zone = IDX(i, j+1, k, xRes+1, yRes+1, zRes);
        dens = densities[ij1Zone].getValue();
        if( dens < rho_critical ){
            surfacePosition[i] = surfacePosition[i] + 1;
        }
    }
}

vector<ScalarVar> LaserManager::getIntensityDivergence(){
    return currentDivI;
}

