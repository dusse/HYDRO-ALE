#include "Remaper.hpp"

#include <iostream>
#include <string>
using namespace std;
using namespace chrono;


Remaper::Remaper(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr):loader(move(load)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());
    initialize();
    for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        vector<VectorVar> reconstructedCoefsForZone;
        reconstructedCoefsForZone.push_back(VectorVar(DENSITY_RECONSTRUCTED_COEF, {0.0, 0.0}));
        reconstructedCoefsForZone.push_back(VectorVar(MOMENTUM_X_DENSITY_RECONSTRUCTED_COEF, {0.0, 0.0}));
        reconstructedCoefsForZone.push_back(VectorVar(MOMENTUM_Y_DENSITY_RECONSTRUCTED_COEF, {0.0, 0.0}));
        reconstructedCoefsForZone.push_back(VectorVar(TOT_ERG_RECONSTRUCTED_COEF, {0.0, 0.0}));
        reconstructedCoefs.push_back(reconstructedCoefsForZone);
        
        vector<VectorVar> minMaxDensitiesForZone;
        minMaxDensitiesForZone.push_back(VectorVar(DENSITY_MIN_MAX, {0.0, 0.0}));
        minMaxDensitiesForZone.push_back(VectorVar(MOMENTUM_X_DENSITY_MIN_MAX, {0.0, 0.0}));
        minMaxDensitiesForZone.push_back(VectorVar(MOMENTUM_Y_DENSITY_MIN_MAX, {0.0, 0.0}));
        minMaxDensitiesForZone.push_back(VectorVar(TOT_ERG_MIN_MAX, {0.0, 0.0}));
        minMaxDensities.push_back(minMaxDensitiesForZone);
        
        vector<ScalarVar> densitiesIntegratedForZone;
        densitiesIntegratedForZone.push_back(ScalarVar(DENSITY_INT_REMAP, 0.0));
        densitiesIntegratedForZone.push_back(ScalarVar(MOMENTUM_X_DENSITY_INT_REMAP, 0.0));
        densitiesIntegratedForZone.push_back(ScalarVar(MOMENTUM_Y_DENSITY_INT_REMAP, 0.0));
        densitiesIntegratedForZone.push_back(ScalarVar(TOT_ERG_DENSITY_INT_REMAP, 0.0));
        densitiesIntegrated.push_back(densitiesIntegratedForZone);
    }
}

void Remaper::initialize(){
    logger->writeMsg("[Remaper] initialize...OK");
}


void Remaper::runRemapping(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[Remaper] remapping...");
//    reconstructInsideOldZones();
    vector<VectorVar> coordsOld        = gridMgr->getVectorVariableForAllNodes(COORDINATE);
    vector<VectorVar> zonalCentroidsOld = gridMgr->getVectorVariableForAllZones(ZONAL_CENTROID);
    gridMgr->smoothMesh();
    gridMgr->updateGrid();
//    integrateInsideNewZones(coordsOld, zonalCentroidsOld);
//    reconstructSolverVariables();
    //todo do in parallel???
    logger->writeMsg("[Remaper] remapping...OK");
    auto end_time = high_resolution_clock::now();
    string msg ="[Remaper] remapping duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
    
}

//g(x,y)= g(ij) + (∂g/∂x)(x−xc)+(∂g/∂y)(y−yc).
void Remaper::reconstructInsideOldZones(){
    auto start_time = high_resolution_clock::now();
    int ijZone, ijNode;
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<VectorVar> velocities  = gridMgr->getVectorVariableForAllNodes(VELOCITY_CTS);
    vector<ScalarVar> energies = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    vector<VectorVar> zonalCentroids = gridMgr->getVectorVariableForAllZones(ZONAL_CENTROID);
    vector<vector<int>> neighborhoodIndices = gridMgr->getNeighborhoodZoneIndicesForAllZone();
    int neighborZoneIDX;
    vector<vector<double>> coefs;
    vector<double> densCoefs;
    vector<double> momXCoefs;
    vector<double> momYCoefs;
    vector<double> ergCoefs;
    
    vector<vector<double>> vars;
    vector<double> dens;
    vector<double> ergTot;
    vector<double> momentumX;
    vector<double> momentumY;
    
    double densMin, momentumXMin, momentumYMin, ergTotMin;
    double densMax, momentumXMax, momentumYMax, ergTotMax;
    
    double velX, velY;
    
    for ( ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        ijNode = gridMgr->getNodeIDXForZone(ijZone);
        velX = velocities[ijNode].getValue()[0];
        velY = velocities[ijNode].getValue()[1];
        dens.push_back( densities[ijZone].getValue() );
        momentumX.push_back(dens[ijZone]*velX);
        momentumY.push_back(dens[ijZone]*velY);
        ergTot.push_back(dens[ijZone]*( energies[ijZone].getValue() + 0.5*( velX*velX + velY*velY ) ));
    }
    vector<int> zonesInBoundaryIndexes = gridMgr->getZonesInBoundaryIndexes();
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexes.begin(); it != zonesInBoundaryIndexes.end(); it++) {
        	ijZone = *it;
            //find min/max value
            densMin = dens[ijZone], momentumXMin = momentumX[ijZone], momentumYMin = momentumY[ijZone], ergTotMin = ergTot[ijZone];
            densMax = dens[ijZone], momentumXMax = momentumX[ijZone], momentumYMax = momentumY[ijZone], ergTotMax = ergTot[ijZone];
            vector<int>::iterator iterator;
            for( iterator = neighborhoodIndices[ijZone].begin(); iterator != neighborhoodIndices[ijZone].end(); iterator++) {
                neighborZoneIDX = (*iterator);
                densMin      = dens[neighborZoneIDX]     <densMin     ?dens[neighborZoneIDX]     :densMin;
                momentumXMin = momentumX[neighborZoneIDX]<momentumXMin?momentumX[neighborZoneIDX]:momentumXMin;;
                momentumYMin = momentumY[neighborZoneIDX]<momentumYMin?momentumY[neighborZoneIDX]:momentumYMin;;
                ergTotMin    = ergTot[neighborZoneIDX]   <ergTotMin   ?ergTot[neighborZoneIDX]   :ergTotMin;
                
                densMax      = dens[neighborZoneIDX]     >densMax     ?dens[neighborZoneIDX]     :densMax;
                momentumXMax = momentumX[neighborZoneIDX]>momentumXMax?momentumX[neighborZoneIDX]:momentumXMax;;
                momentumYMax = momentumY[neighborZoneIDX]>momentumYMax?momentumY[neighborZoneIDX]:momentumYMax;;
                ergTotMax    = ergTot[neighborZoneIDX]   >ergTotMax   ?ergTot[neighborZoneIDX]   :ergTotMax;
                
            }
            minMaxDensities[ijZone][0] = VectorVar(DENSITY_MIN_MAX, {densMin, densMax});
            minMaxDensities[ijZone][1] = VectorVar(MOMENTUM_X_DENSITY_MIN_MAX, {momentumXMin, momentumXMax});
            minMaxDensities[ijZone][2] = VectorVar(MOMENTUM_Y_DENSITY_MIN_MAX, {momentumYMin, momentumYMax});
            minMaxDensities[ijZone][3] = VectorVar(TOT_ERG_MIN_MAX, {ergTotMin, ergTotMax});

    }
    
    minimizeAvgDerivativeInNeighborhood({ dens, momentumX, momentumY, ergTot}, zonalCentroids);
//    calculateAvgDerivativeInNeighborhood({ dens, momentumX, momentumY, ergTot}, zonalCentroids);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[Remaper] reconstructInsideOldZones() duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
}


/**
 
 ∂G = ∫g(x,y)dxdy =  (g(ij) - (∂g/∂x)xc - (∂g/∂y)yc)∫1dxdy + (∂g/∂x)∫xdxdy + (∂g/∂y)∫ydxdy

 **/
void Remaper::integrateInsideNewZones(vector<VectorVar> coordsOld, vector<VectorVar> zonalCentroidsOld){
    auto start_time = high_resolution_clock::now();
    int ijZone, regionNum, ijZoneDonor, varNum;
    vector<ScalarVar> zonalVolumes   = gridMgr->getScalarVariableForAllZones(ZONAL_VOLUME);// after smoothing
    
    vector<vector<double>> deltas = { {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0} };
    vector<double> deltasTotal = {0.0, 0.0, 0.0, 0.0};    
    vector<vector<double>> zoneDensities;
    vector<vector<double>> zoneMasses;
    calculateMomentsForAllQuantities( &zoneDensities, &zoneMasses );
    for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        for( varNum = 0; varNum < 4; varNum++){
            densitiesIntegrated[ijZone][varNum] = ScalarVar(densitiesIntegrated[ijZone][varNum].getName(),  zoneMasses[ijZone][varNum]);
        }
    }
    double sweptRegionMass, sweptRegionMassTot, minDens, maxDens, minDonorDens, maxDonorDens;
    double volIntegral, xIntegral, yIntegral;
    double idealDens, idealMass, currentMass, deltaMass, availDeltaMass, currentDonorDens, currentDonorMass;
    double xc, yc, xCoef, yCoef;
    double smoothedZoneVolume, smoothedDonorZoneVolume;
    double ratio;
    vector<vector<int>> neighborhoodIndicesAll = gridMgr->getNeighborhoodZoneIndicesForAllZone();
    vector<int> neighborNodeIDXs, neighborZoneIDXs, donorZoneIDXs;
    vector<vector<vector<double>>> spatialIntegrals;
    calculateSpatialIntegrals( &spatialIntegrals, coordsOld );
    vector<int> zonesInBoundaryIndexesSpiral = gridMgr->getZonesInBoundaryIndexesSpiral();
    //calculate volume and decide which zone will be used for integration
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexesSpiral.begin(); it != zonesInBoundaryIndexesSpiral.end(); it++) {
        ijZone = *it;
        neighborZoneIDXs = neighborhoodIndicesAll[ijZone];
        donorZoneIDXs = {neighborZoneIDXs[1], neighborZoneIDXs[3], neighborZoneIDXs[5], neighborZoneIDXs[7]};
        // calculate smoothed volume
        smoothedZoneVolume = zonalVolumes[ijZone].getValue();
        deltasTotal = {0.0, 0.0, 0.0, 0.0};
        for ( regionNum = 0; regionNum < 4; regionNum++ ){// 4 quadrilateral swept regions
                
                volIntegral = spatialIntegrals[ijZone][regionNum][0];
                xIntegral   = spatialIntegrals[ijZone][regionNum][1];
                yIntegral   = spatialIntegrals[ijZone][regionNum][2];
                
                ijZoneDonor = (volIntegral > 0.0) ? donorZoneIDXs[regionNum] : ijZone;
                xc = zonalCentroidsOld[ijZoneDonor].getValue()[0];
                yc = zonalCentroidsOld[ijZoneDonor].getValue()[1];
                
                for( varNum = 0; varNum < 4; varNum++){
                    xCoef =  reconstructedCoefs[ijZoneDonor][varNum].getValue()[0];
                    yCoef =  reconstructedCoefs[ijZoneDonor][varNum].getValue()[1];
                    //todo back to normal calcul
//                    sweptRegionMass = ((zoneDensities[ijZoneDonor][varNum] - xCoef*xc - yCoef*yc)*volIntegral + xCoef*xIntegral + yCoef*yIntegral);
                    sweptRegionMass = zoneDensities[ijZoneDonor][varNum]*volIntegral;
                    deltas[regionNum][varNum] = sweptRegionMass;
                    deltasTotal[varNum]      += sweptRegionMass;
                }
        }
            
            //todo check neibor before action
            for( varNum = 0; varNum < 4; varNum++){
                 
                minDens = minMaxDensities[ijZone][varNum].getValue()[0];
                maxDens = minMaxDensities[ijZone][varNum].getValue()[1];
                
                // we want to save old density
                idealDens   = zoneDensities[ijZone][varNum];
                idealMass   = idealDens*smoothedZoneVolume;
//                zoneMasses[ijZone][varNum] = idealMass;
                currentMass = zoneMasses[ijZone][varNum];//old value or changed by previous step
                double currentDens = currentMass/smoothedZoneVolume;
                deltaMass   = idealMass - currentMass;// how much need to add/extract to save the same density
                
                if( deltaMass < EPSILON || areSame( idealDens, currentDens )){
                    continue;
                }
                sweptRegionMassTot = deltasTotal[varNum];
                if( fabs(sweptRegionMassTot) >= fabs(deltaMass) ){
                    ratio = fabs(deltaMass/sweptRegionMassTot);
                }else{
                    ratio = 1.0;
                }
                
                // extract from donor zones with smart algorithm
                availDeltaMass = 0.0;
                for ( regionNum = 0; regionNum < 4; regionNum++ ){
                    ijZoneDonor = donorZoneIDXs[regionNum];
                    smoothedDonorZoneVolume = zonalVolumes[ijZoneDonor].getValue();
                    minDonorDens = minMaxDensities[ijZoneDonor][varNum].getValue()[0];
                    maxDonorDens = minMaxDensities[ijZoneDonor][varNum].getValue()[1];
                    currentDonorMass = zoneMasses[ijZoneDonor][varNum];
                    sweptRegionMass = deltas[regionNum][varNum]*ratio;
                    currentDonorMass -= sweptRegionMass;
                    
                    currentDonorDens = currentDonorMass/smoothedDonorZoneVolume;
                    if( currentDonorDens > maxDonorDens || currentDonorDens < minDonorDens){
                        sweptRegionMassTot -= sweptRegionMass;
                        ratio = ratio == 1.0 ? 1.0 : fabs(deltaMass/sweptRegionMassTot);
                    }else{
                        zoneMasses[ijZoneDonor][varNum] -= sweptRegionMass;
                        availDeltaMass += sweptRegionMass;
                        
                    }
                }
                zoneMasses[ijZone][varNum] += availDeltaMass;
                
                densitiesIntegrated[ijZone][varNum] = ScalarVar(densitiesIntegrated[ijZone][varNum].getName(),  zoneMasses[ijZone][varNum]);
        }
    }
    auto end_time = high_resolution_clock::now();
    string msg ="[Remaper] Integration inside new zones duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
}

void Remaper::calculateMomentsForAllQuantities(vector<vector<double>> *zoneDensities, vector<vector<double>> *zoneDensitiesIntegrated){
    int ijZone, ijNode;
    double velX, velY;
    double dens, mass, momxDens, momyDens, momx, momy, erg, tote, toteDens;
    vector<ScalarVar> densities   = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<VectorVar> velocities  = gridMgr->getVectorVariableForAllNodes(VELOCITY_CTS);
    vector<ScalarVar> energies    = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    vector<ScalarVar> zonalMasses = gridMgr->getScalarVariableForAllZones(ZONAL_MASS);
    
    for ( ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        ijNode = gridMgr->getNodeIDXForZone(ijZone);
        velX = velocities[ijNode].getValue()[0];
        velY = velocities[ijNode].getValue()[1];
        dens = densities[ijZone].getValue();
        if( dens < 0.0 ){
            string msg ="[Remaper] negative density problem ijZone = "+to_string(ijZone)+" dens = "+to_string(dens);
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        mass = zonalMasses[ijZone].getValue();
        if( mass <= 0.0 ){
            string msg ="[Remaper] negative mass problem ijZone = "+to_string(ijZone)+" mass = "+to_string(mass);
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        erg  = energies[ijZone].getValue();
        if( erg < 0.0 ){
            string msg ="[Remaper] negative energy problem ijZone = "+to_string(ijZone)+" erg = "+to_string(erg);
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        momxDens = dens*velX;
        momyDens = dens*velY;
        toteDens = dens*( erg + 0.5*( velX*velX + velY*velY ));
        (*zoneDensities).push_back( {dens, momxDens, momyDens, toteDens} );
        
        momx = mass*velX;
        momy = mass*velY;
        tote = mass*( erg + 0.5*( velX*velX  +velY*velY ));
        (*zoneDensitiesIntegrated).push_back( {mass, momx, momy, tote} );
    }
}



// ∫1dxdy =  ∮xdy    =  Σ1/2(x1 + x2)(y2 - y1)
// ∫xdxdy =  ∮x^2dy  =  Σ1/6(x1^2 + x1x2 + x2^2)(y2-y1)
// ∫ydxdy = -∮y^2dx  = -Σ1/6(y1^2 + y1y2 + y2^2)(x2-x1)
void Remaper::calculateSpatialIntegrals(vector<vector<vector<double>>> *spatialIntegrals, vector<VectorVar> coordsOld){
    double xleft, xright, yleft, yright, dx, dy;
    double volIntegral, xIntegral, yIntegral;
    int  cornerNumber, regionNum;
    vector<vector<int>> indicesOfNeighborNodes = gridMgr->getIndicesOfNeighborNodesForAllZones();
    vector<VectorVar> coordsSmoothed         = gridMgr->getVectorVariableForAllNodes(COORDINATE);// already smoothed
    
    vector<vector<double>> integrals = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
    double sweptRegionVolTest;

    //calculate volume
    vector<vector<VectorVar>> right;
    vector<vector<VectorVar>> left;
    
    /**              
     *                                                    x Node'(i,j)
     *                          swept region 3
     *           .Node(i-1,j)___________________.Node(i,j)
     *           |               |               |
     *     x Node'(i-1,j)        .               |  s
     *     |     |               |               |  w
     *     |     |                                  e
     *     |  s  |               |               |  p
     *     |  w  |               .               |  t
     *     |  e  |               |               |
     *     |  p  |               |               |  r
     *     |  t  |- - - - - - - -x Zone(i,j) - - |  e
     *     |     |               |               |  g
     *     |  r  |               .               |  i
     *     |  e  |               |               |  o
     *     |  g  |               .               |  n
     *     |  i  |
     *     |  o  |               .               |  2
     *     |  n  |               |               |
     *     |    .Node(i-1,j-1)__________________.Node(i,j-1)
     *     |  4/                                \
     *      \ /                swept region 1    \
     *       x _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _x
     *      Node'(i-1,j-1)                      Node'(i,j-1)
     **/
    int indexR, indexL;
    vector<int> neighborNodeIDX;
    

    int ijZone;
    for ( ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        (*spatialIntegrals).push_back(integrals);
    }
    vector<int> zonesInBoundaryIndexes = gridMgr->getZonesInBoundaryIndexes();
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexes.begin(); it != zonesInBoundaryIndexes.end(); it++) {
        	ijZone = *it;
            neighborNodeIDX = indicesOfNeighborNodes[ijZone];
            integrals = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
            for ( regionNum = 0; regionNum < 4; regionNum++ ){// 4 quadrilateral swept regions
                indexL = regionNum;
                indexR = (regionNum == 3) ? 0 : regionNum + 1;
                left.push_back({coordsSmoothed[neighborNodeIDX[indexL]], coordsSmoothed[neighborNodeIDX[indexR]], coordsOld[neighborNodeIDX[indexR]],coordsOld[neighborNodeIDX[indexL]]});
                right.push_back({coordsSmoothed[neighborNodeIDX[indexR]], coordsOld[neighborNodeIDX[indexR]], coordsOld[neighborNodeIDX[indexL]],coordsSmoothed[neighborNodeIDX[indexL]]});
                volIntegral = 0.0; // ∫1dxdy =  ∮xdy    =  Σ1/2(x1 + x2)(y2 - y1)
                xIntegral = 0.0;   // ∫xdxdy =  ∮x^2dy  =  Σ1/6(x1^2 + x1x2 + x2^2)(y2-y1)
                yIntegral = 0.0;   // ∫ydxdy = -∮y^2dx  = -Σ1/6(y1^2 + y1y2 + y2^2)(x2-x1)
                
                sweptRegionVolTest  = 0.0;
                for ( cornerNumber=0; cornerNumber<4; cornerNumber++ ){// 4 vertices
                    xleft  = left [regionNum][cornerNumber].getValue()[0];
                    xright = right[regionNum][cornerNumber].getValue()[0];
                    yleft  = left [regionNum][cornerNumber].getValue()[1];
                    yright = right[regionNum][cornerNumber].getValue()[1];
                    
                    dy = yright-yleft;
                    dx = xright-xleft;
                    volIntegral += 0.5*(xright+xleft)*dy;
                    xIntegral += ( xright*xright + xright*xleft + xleft*xleft )*dy/6;
                    yIntegral -= ( yright*yright + yright*yleft + yleft*yleft )*dx/6;
            
                    sweptRegionVolTest -= 0.5*(yright+yleft)*dx;
                }
                if( !areSame(volIntegral, sweptRegionVolTest)){
                    logger->writeMsg("[Remaper] GRID is broken...", CRITICAL);
                }else{
                    integrals[regionNum][0] = volIntegral;
                    integrals[regionNum][1] = xIntegral;
                    integrals[regionNum][2] = yIntegral;
                }
            }
            (*spatialIntegrals)[ijZone] = integrals;
    }

}

void Remaper::reconstructSolverVariables(){
    auto start_time = high_resolution_clock::now();
    int ijZone, ijNode, subzoneNum;
    double subzonalMass;
    double zonalMass;
    double density = 0.0;
    double subzonalVolume;
    double velX, velY;
    double intErg;
    double zoneVol;
    
    vector<ScalarVar> zonalVolumes = gridMgr->getScalarVariableForAllZones(ZONAL_VOLUME);
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();

    for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
            ijNode = gridMgr->getNodeIDXForZone(ijZone);
            zonalMass = 0.0;
            zoneVol = zonalVolumes[ijZone].getValue();
            density = densitiesIntegrated[ijZone][0].getValue()/zoneVol;
            if( density < 0.0 ){
                string msg ="[Remaper] set abs for negative density ijZone = "+to_string(ijZone)+" density = "+to_string(density);
                logger->writeMsg(msg.c_str(), CRITICAL);
                density = fabs(density);
            }

            for ( subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++){
                subzonalVolume = subzones[ijZone][subzoneNum].getVariable(SUBZONAL_VOLUME).getValue();
                subzonalMass = density*subzonalVolume;
                zonalMass += subzonalMass;
                gridMgr->setScalarVariableForSubzone(ijZone, subzoneNum, ScalarVar(SUBZONAL_MASS, subzonalMass));
                gridMgr->setScalarVariableForSubzone(ijZone, subzoneNum, ScalarVar(DENSITY, density));
            }
            if( zonalMass != zonalMass || zonalMass <= 0.0 || !areSame(density*zoneVol, zonalMass)){
                string msg ="[Remaper]  negative zonalMass ijZone = "+to_string(ijZone)+" zonalMass = "+to_string(zonalMass);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(ZONAL_MASS, zonalMass));
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(DENSITY, density));
            
            velX = densitiesIntegrated[ijZone][1].getValue()/zonalMass;
            velY = densitiesIntegrated[ijZone][2].getValue()/zonalMass;
            
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_CTS, {velX, velY, 0.0}));
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_OUT, {velX, velY, 0.0}));
            
            intErg = (densitiesIntegrated[ijZone][3].getValue()/zonalMass - 0.5*( velX*velX + velY*velY ));
            if( intErg != intErg || intErg < 0.0 ){
                string msg ="[Remaper] set negative ENERGY to zero for ijZone = "+to_string(ijZone)+" intErg = "+to_string(intErg);
                logger->writeMsg(msg.c_str(), CRITICAL);
                intErg = 0.0;
            }
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(INTERNAL_ERG, intErg));
            

    }
    
    //fill border set mass to zero
    map<int, int> zonesBoundaryIndexes = gridMgr->getZonesBoundaryIndexes();
    map<int, int>::iterator it;
    for( it = zonesBoundaryIndexes.begin(); it != zonesBoundaryIndexes.end(); it++) {
            ijZone = it->first;
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(ZONAL_MASS, 0.0));
    }

    auto end_time = high_resolution_clock::now();
    string msg ="[Remaper] Reconstruct variables duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
}




void Remaper::calculateAvgDerivativeInNeighborhood(vector<vector<double>> variables, vector<VectorVar> zonalCentroids){
    int  cornerNumber, dim;
    int zoneIDX;
    int const pairNumberTot = 8;
    double xleft, xright, yleft, yright, dx, dy;
    double neighborhoodVolume = 0.0, neighborhoodVolume1 = 0.0;
    vector<int> neighborhoodIndices;
    vector<vector<double>> derivative;
    for(int var = 0; var < variables.size(); var++){
        derivative.push_back({0.0, 0.0});
    }
    int firstTermIDX, secondTermIDX;
    double firstCoord, secondCoord;
    int deltaDim;
    vector<vector<int>> neighborhoodIndicesAll = gridMgr->getNeighborhoodZoneIndicesForAllZone();
    //calculate volume
    vector<vector<int>> right;
    vector<vector<int>> left;
    int regionNum, pairNumber;
    vector<int> zonesInBoundaryIndexes = gridMgr->getZonesInBoundaryIndexes();
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexes.begin(); it != zonesInBoundaryIndexes.end(); it++) {
    	zoneIDX = *it;
            neighborhoodIndices = neighborhoodIndicesAll[zoneIDX];
            right.push_back({neighborhoodIndices[0], neighborhoodIndices[1], zoneIDX, neighborhoodIndices[7]});
            left.push_back ({neighborhoodIndices[7], neighborhoodIndices[0], neighborhoodIndices[1], zoneIDX});

            right.push_back({neighborhoodIndices[2], neighborhoodIndices[3], zoneIDX, neighborhoodIndices[1]});
            left.push_back ({neighborhoodIndices[1], neighborhoodIndices[2], neighborhoodIndices[3], zoneIDX});

            right.push_back({neighborhoodIndices[4], neighborhoodIndices[5], zoneIDX, neighborhoodIndices[3]});
            left.push_back ({neighborhoodIndices[3], neighborhoodIndices[4], neighborhoodIndices[5], zoneIDX});

            right.push_back({neighborhoodIndices[6], neighborhoodIndices[7], zoneIDX, neighborhoodIndices[5]});
            left.push_back ({neighborhoodIndices[5], neighborhoodIndices[6], neighborhoodIndices[7], zoneIDX});


            for ( regionNum = 0; regionNum < 4; regionNum++ ){// 4 quadrilateral regions
                for ( cornerNumber=0; cornerNumber<4; cornerNumber++ ){// 4 points
                    xleft  = zonalCentroids[left [regionNum][cornerNumber]].getValue()[0];
                    xright = zonalCentroids[right[regionNum][cornerNumber]].getValue()[0];
                    yleft  = zonalCentroids[left [regionNum][cornerNumber]].getValue()[1];
                    yright = zonalCentroids[right[regionNum][cornerNumber]].getValue()[1];

                    dy = yright-yleft;
                    dx = xright-xleft;
                    neighborhoodVolume  += 0.5*(xright+xleft)*dy;
                    neighborhoodVolume1 += -0.5*(yright+yleft)*dx;
                }
            }
            if( !areSame(neighborhoodVolume, neighborhoodVolume1)){
                string msg ="[Remaper] PROBLEM neighborhoodVolume = "+to_string(neighborhoodVolume)+"  ; neighborhoodVolume1 = "+to_string(neighborhoodVolume1);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }

            for (dim = 0; dim<2; dim++){
                deltaDim = (dim == 0) ? 1 : 0 ;
                neighborhoodVolume = (dim == 0) ? fabs(neighborhoodVolume) : -fabs(neighborhoodVolume) ;

                for ( pairNumber = 0; pairNumber < pairNumberTot; pairNumber++ ){
                    firstTermIDX  = neighborhoodIndices[pairNumber];
                    secondTermIDX = (pairNumber == (pairNumberTot-1))?neighborhoodIndices[0]:neighborhoodIndices[pairNumber+1];
                    firstCoord  = zonalCentroids[firstTermIDX].getValue()[deltaDim];
                    secondCoord = zonalCentroids[secondTermIDX].getValue()[deltaDim];
                    for(int var = 0; var < variables.size(); var++){
                        derivative[var][dim] += 0.5*(variables[var][firstTermIDX]+variables[var][secondTermIDX])*(secondCoord-firstCoord)/neighborhoodVolume;
                    }
                }
            }
            for(int var = 0; var < variables.size(); var++){
                reconstructedCoefs[zoneIDX][var] = VectorVar(reconstructedCoefs[zoneIDX][var].getName(), derivative[var]);
            }
    }
}

/**
 
 | Axx Axy | . | x' | = | Bx |
 | Axy Ayy |   | y' |   | By |
 

 | x' | =  _1_  |  Ayy -Axy | . | Bx |
 | y' |   detA  | -Axy  Axx |   | By |
 

 | x' | = |  AyyBx - AxyBy | .  _1_
 | y' |   | -AxyBx + AxxBy |    detA

 
 **/
void Remaper::minimizeAvgDerivativeInNeighborhood(vector<vector<double>> variables, vector<VectorVar> zonalCentroids){
    vector<vector<double>> derivativeMin;
    for(int var = 0; var < variables.size(); var++){
        derivativeMin.push_back({0.0, 0.0});
    }
    int zoneIDX;
    int const neighborsNumber = 8;
    double Axx = 0.0, Axy = 0.0, Ayy = 0.0, Bx = 0.0, By = 0.0;
    double detA = 0.0;
    int neighborNum, neighborIDX;
    double xZone, xNeighbor,yZone, yNeighbor;
    double varZone, varNeighbor;
    double dx, dy;
    vector<vector<int>> neighborhoodIndicesAll = gridMgr->getNeighborhoodZoneIndicesForAllZone();
    vector<int> zonesInBoundaryIndexes = gridMgr->getZonesInBoundaryIndexes();
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexes.begin(); it != zonesInBoundaryIndexes.end(); it++) {
    	zoneIDX = *it;
            vector<int> neighborhoodIndices = neighborhoodIndicesAll[zoneIDX];
            xZone = zonalCentroids[zoneIDX].getValue()[0];
            yZone = zonalCentroids[zoneIDX].getValue()[1];

            for(int var = 0; var < variables.size(); var++){
                varZone = variables[var][zoneIDX];
                for ( neighborNum=0; neighborNum<neighborsNumber; neighborNum++ ){
                    neighborIDX  = neighborhoodIndices[neighborNum];
                    xNeighbor = zonalCentroids[neighborIDX].getValue()[0];
                    yNeighbor = zonalCentroids[neighborIDX].getValue()[1];
                    varNeighbor = variables[var][neighborIDX];
                    dx = (xNeighbor-xZone);
                    dy = (yNeighbor-yZone);
                    Axx += 2*dx*dx;
                    Axy += 2*dx*dy;
                    Ayy += 2*dy*dy;
                    Bx  += 2*dx*(varNeighbor-varZone);
                    By  += 2*dy*(varNeighbor-varZone);
                }
                detA = Axx*Ayy - Axy*Axy;
                if( detA<EPSILON ){
                    string msg ="[Remaper] DET(A) less then EPSILON; detA = "+to_string(detA);
                    logger->writeMsg(msg.c_str(), CRITICAL);
                    derivativeMin[var] = { 0.0, 0.0 };
                }else{
                    derivativeMin[var] = {(Bx*Ayy-By*Axy)/detA, (By*Axx-Bx*Axy)/detA};
                }
                reconstructedCoefs[zoneIDX][var] = VectorVar(reconstructedCoefs[zoneIDX][var].getName(), derivativeMin[var]);
            }
    }
}

void Remaper::finilize()
{
    logger->writeMsg("[Remaper] finilize...OK");
}

Remaper::~Remaper()
{
	finilize();
}


