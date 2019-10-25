#include "LagrangianSolver.hpp"



using namespace std;
using namespace chrono;

/**             1 Dρ  = - ∇w
 *              ρ Dt
 *
 *
 *              ρ Dw  = - ∇p
 *                Dt
 *
 *              ρ De  = - p∇w  − Ca div(I) + div(k∇T) //Ca = 0.5 - 1st,  Ca = 0.75 3d laser frequency
 *                Dt
 *
 *            |  D   =   ∂  + W x ∇ |
 *            |  Dt      ∂t         |
 *
 *
 *
 **/


LagrangianSolver::LagrangianSolver(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr, shared_ptr<LaserManager> laserMnr, shared_ptr<HeatManager> heatMnr):loader(move(load)), gridMgr(move(gridMnr)),laserMgr(move(laserMnr)),heatMgr(move(heatMnr))
{
    logger.reset(new Logger());
    
    remaper.reset(new Remaper(loader, gridMgr));
    
    initialize();
    logger->writeMsg("[LagrangianSolver] create...OK");
}

void LagrangianSolver::initialize()
{
    initNodes();
    initZones();
    initZonalMass();
    initNodalMass();
}

int LagrangianSolver::solve( double time)
{
    int ijNode, dimNum;
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    if( time > 0.2E-9 && time < 20.0E-9 ){ // todo hardcoded times
        remaper->runRemapping();
        initNodalMass();
        calculatePressure();
        applyPressureBC();
    }

//    checkConservativeLaws();
    calculateForces();
    updateNodeVelocities();
    updateNodeCoordinates();
    
    // Update the geometry according to the new nodal position
    if ( gridMgr->updateGrid() == GRID_FAIL ){
        //try again
        //rollback old coordinates
        vector<VectorVar> coordsPrev = gridMgr->getVectorVariableForAllNodes(COORDINATE_OLD);
        for(ijNode=0; ijNode<gridMgr->getTotalNodeNumber();ijNode++){
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(COORDINATE, coordsPrev[ijNode].getValue()));

        }
        remaper->runRemapping();
        initNodalMass();
        calculatePressure();
        applyPressureBC();
        
        calculateForces();
        updateNodeVelocities();
        updateNodeCoordinates();
        if( gridMgr->updateGrid() == GRID_FAIL ){
            logger->writeMsg("[LagrangianSolver] remapping has not helped...(", CRITICAL);
            return SOLVE_FAIL;
        }
    }
    
//    heatMgr->solve(); todo heat manager switch on
    laserMgr->solve(time);

    calculateInternalEnergy();
    calculateDensity();

    calculatePressure();
    applyPressureBC();
    

//    checkConservativeLaws();
    
    
    vector<VectorVar> velocitiesNTS = gridMgr->getVectorVariableForAllNodes(VELOCITY_NTS);
    for (ijNode=0; ijNode<gridMgr->getTotalNodeNumber(); ijNode++){
        vector<double> vel;
        for (dimNum = 0; dimNum < loader->dim; dimNum++){
            vel.push_back(velocitiesNTS[ijNode].getValue()[dimNum]);
        }
        gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_OUT, vel));
        gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_CTS, vel));
       
    }

    if (rank == 0){
        string msg ="[LagrangianSolver] SOLVER time ="+to_string(time)+"; timeStepNum = "+to_string(timeStepNum);
        logger->writeMsg(msg.c_str(), INFO);
    }
    timeStepNum++;
    return SOLVE_OK;
}


void LagrangianSolver::initNodes(){
    double x,y,z;
    vector<double> vel;
    vector<VectorVar> coords = gridMgr->getVectorVariableForAllNodes(COORDINATE);
    for (int ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
         x = coords[ijNode].getValue()[0];
         y = coords[ijNode].getValue()[1];
         z = coords[ijNode].getValue()[2];
         vel = loader->getVelocity(x, y, z);

         gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_CTS, vel));
         gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_HTS, vel));
         gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_OUT, vel));
    }
}

void LagrangianSolver::initZones(){
    
    double x,y,z;
    int ijZone, ijNode;
    double density = 0.0;
    double presure = 0.0;
    double intErg = 0.0;
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<VectorVar> coords = gridMgr->getVectorVariableForAllNodes(COORDINATE);
    vector<int> zonesInBoundaryIndexes = gridMgr->getZonesInBoundaryIndexes();
    vector<int>::iterator it;
    for( it = zonesInBoundaryIndexes.begin(); it != zonesInBoundaryIndexes.end(); it++) {
            ijZone = *it;
            ijNode = gridMgr->getNodeIDXForZone(ijZone);
            x = coords[ijNode].getValue()[0];
            y = coords[ijNode].getValue()[1];
            z = coords[ijNode].getValue()[2];
            
            density = loader->getDensity(x, y, z);
            //      Z + 1
            // p =  ----- x C(p/T) x T x ρ 
            //        A
            presure = 1/A*Cpt*Tout*density;
            intErg  = presure/( ( GAMMA - 1 )*density);
            
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(DENSITY, density));
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(PRESURE, presure));
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(INTERNAL_ERG, intErg));

    }
    updateBoundaryZoneVariable();
}

void LagrangianSolver::updateBoundaryZoneVariable(){
    map<int, int> zonesBoundaryIndexes = gridMgr->getZonesBoundaryIndexes();
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<ScalarVar> presures  = gridMgr->getScalarVariableForAllZones(PRESURE);
    double dens, pres, erg;
    int zoneIDX, donorIDX;

    map<int, int>::iterator it;
    for( it = zonesBoundaryIndexes.begin(); it != zonesBoundaryIndexes.end(); it++) {
        zoneIDX  = it->first;
        donorIDX = it->second;
        dens = densities[donorIDX].getValue();
        pres = presures [donorIDX].getValue();
        erg  = pres/( ( GAMMA - 1 )*dens );
        gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(DENSITY, dens));
        gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(PRESURE, pres));
        gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(INTERNAL_ERG, erg));
    }
}

void LagrangianSolver::initZonalMass(){
    int subzoneNum, ijZone;
    double subzonalMass;
    double zonalMass;
    double density;
    double subzonalVolume;
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    for ( ijZone=0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        subzonalMass = 0.0;
        zonalMass = 0.0;
        density = densities[ijZone].getValue();
        for ( subzoneNum=0; subzoneNum < loader->subzoneNum; subzoneNum++ ){
            subzonalVolume = subzones[ijZone][subzoneNum].getVariable(SUBZONAL_VOLUME).getValue();
            subzonalMass = density*subzonalVolume;
            zonalMass += subzonalMass;
                gridMgr->setScalarVariableForSubzone(ijZone, subzoneNum, ScalarVar(SUBZONAL_MASS, subzonalMass));
                gridMgr->setScalarVariableForSubzone(ijZone, subzoneNum, ScalarVar(DENSITY, density));
        }
        if(zonalMass <= 0){
            string msg ="[LagrangianSolver] init zonal mass problem!!! zonalMass = "+to_string(zonalMass);
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
        gridMgr->setScalarVariableForZone(ijZone, ScalarVar(ZONAL_MASS, zonalMass));
    }
}

void LagrangianSolver::initNodalMass(){
    int subIDX, subzoneNum;
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<vector<int>> indicesOfNeighbors  = gridMgr->getIndicesOfNeighborZonesForAllNodes();
    double nodalMass, subzonalMass;
    for (int ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
         nodalMass = 0.0;
         for ( subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++ ){
              subIDX = indicesOfNeighbors[ijNode][subzoneNum];
              subzonalMass = subzones[subIDX][subzoneNum].getVariable(SUBZONAL_MASS).getValue();
              nodalMass += subzonalMass;
          }
          gridMgr->setScalarVariableForNode(ijNode, ScalarVar(NODAL_MASS, nodalMass));
    }
    
}



void LagrangianSolver::calculateForces(){
    //1. forces =  - ∇p.
    //1. For each subzone compute  zonal, subzonal and viscosity pressure forces
    double nodeTotPresForce_X, nodeTotPresForce_Y;
    double pres;//Fplc
    
    double xEdgeCenter1, xEdgeCenter2, xZoneCenter;
    double yEdgeCenter1, yEdgeCenter2, yZoneCenter;
    
    double zonalDens ,zonalDensL, zonalDensR;
    double subzonalDens, subzonalDensL, subzonalDensR;
    double deltaPres, deltaPresL, deltaPresR;
    double integralForceX, integralForceX1, integralForceX2;
    double integralForceY, integralForceY1, integralForceY2;
    double subzonePresForce_X, subzonePresForce_Y;
    double erg, ergL, ergR;
    int ijNode, subIDX;
    int neighborZoneIDX, subzoneNum;
    vector<ScalarVar> presures = gridMgr->getScalarVariableForAllZones(PRESURE);
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<vector<int>> indicesOfNeighborZones  = gridMgr->getIndicesOfNeighborZonesForAllNodes();
    vector<vector<double>> verticesCoordinates;
    vector<VectorVar> zonalCentroids = gridMgr->getVectorVariableForAllZones(ZONAL_CENTROID);
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<ScalarVar> energies = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    
    for ( ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
        nodeTotPresForce_X = 0.0;
        nodeTotPresForce_Y = 0.0;
        for ( subzoneNum = 0; subzoneNum < loader->subzoneNum; subzoneNum++ ){
                neighborZoneIDX = indicesOfNeighborZones[ijNode][subzoneNum];
                pres = presures[neighborZoneIDX].getValue();
                
                /**
                 *          x Zone(i,j+1) - - - - - - - - - x Zone(i+1,j+1)
                 *          .               |               .
                 *          |               |               |
                 *          .               |               .
                 *          |               |               |
                 *          .  Subzone 2    |   Subzone 1   .
                 *          |               |               |
                 *          .               |               .
                 *          |               |               |
                 *          ________________.Node(i,j)_______
                 *          |               |               |
                 *          .               |               .
                 *          |               |               |
                 *          .    Subzone 3  |   Subzone 4   .
                 *          |               |               |
                 *          .               |               .
                 *          |               |               |
                 *          .               |               .
                 *          x Zone(i,j) - -   - - - - - - - x Zone(i+1,j)
                 **/
                verticesCoordinates = subzones[neighborZoneIDX][subzoneNum].getVerticesCoordinates();

                xEdgeCenter1 = verticesCoordinates[0][1];
                xZoneCenter  = verticesCoordinates[0][2];
                xEdgeCenter2 = verticesCoordinates[0][3];
                
                yEdgeCenter1 = verticesCoordinates[1][1];
                yZoneCenter  = verticesCoordinates[1][2];
                yEdgeCenter2 = verticesCoordinates[1][3];
            
                integralForceX1 = xZoneCenter  - xEdgeCenter1;
                integralForceX2 = xEdgeCenter2 - xZoneCenter;
                integralForceX  = integralForceX1 + integralForceX2;
                
                integralForceY1 = yZoneCenter  - yEdgeCenter1;
                integralForceY2 = yEdgeCenter2 - yZoneCenter;
                integralForceY  = integralForceY1 + integralForceY2;
                
                subzonePresForce_X = -(pres)*integralForceY;
                subzonePresForce_Y =  (pres)*integralForceX;
                
                nodeTotPresForce_X += subzonePresForce_X;
                nodeTotPresForce_Y += subzonePresForce_Y;
                
            /************************************************************************************************/
//                
//            subzonalDens  = subzones [neighborZoneIDX][subzoneNum].getVariable(DENSITY).getValue();
//            zonalDens     = densities[neighborZoneIDX].getValue();
//            erg           = energies [neighborZoneIDX].getValue();
//            
//            int index = subzoneNum == 0 ? loader->subzoneNum-1 : subzoneNum-1;
//            subIDX = indicesOfNeighborZones[ijNode][index];
//            subzonalDensL = subzones[subIDX][index].getVariable(DENSITY).getValue();
//            
//            zonalDensL = densities[subIDX].getValue();
//            ergL       = energies [subIDX].getValue();
//            
//            index = subzoneNum == loader->subzoneNum-1 ? 0 : subzoneNum+1;
//            subIDX = indicesOfNeighborZones[ijNode][index];
//            subzonalDensR = subzones[subIDX][index].getVariable(DENSITY).getValue();
//
//            zonalDensR = densities[subIDX].getValue();
//            ergR       = energies[subIDX].getValue();
//                
//            deltaPres  =  (GAMMA-1)* erg *( subzonalDens  - zonalDens );
//            deltaPresL =   deltaPres - (GAMMA-1)* ergL *( subzonalDensL - zonalDensL );
//            deltaPresR =  (GAMMA-1)* ergR *( subzonalDensR - zonalDensR ) - deltaPres;
//                
//            subzonePresForce_X -= deltaPres*integralForceY + 0.5*( integralForceY2*deltaPresR + integralForceY1*deltaPresL );
//            subzonePresForce_Y += deltaPres*integralForceX + 0.5*( integralForceX2*deltaPresR + integralForceX1*deltaPresL );
//                
//            nodeTotPresForce_X += subzonePresForce_X;
//            nodeTotPresForce_Y += subzonePresForce_Y;
            
            /************************************************************************************************/
                
            gridMgr->setVectorVariableForSubzone(neighborZoneIDX, subzoneNum, VectorVar(SUBZONAL_PRES_FORCE, {subzonePresForce_X, subzonePresForce_Y}));
        }

        gridMgr->setVectorVariableForNode(ijNode, VectorVar(NODE_TOT_PRES_FORCE, {nodeTotPresForce_X, nodeTotPresForce_Y}));
    }
}



void LagrangianSolver::updateNodeVelocities(){
    // forces ->  node velocities.
    // According to the nodal forces compute n+1  and  n+1/2 velocities
    // w[n+1] = w[n] + ∆t/m * F[n]
    // w[n+1/2] = 0.5*( w[n] + w[n+1] )
    int dimNum;
    int ijNode;
    double nodalMass;
    double velNext, velPrev, velDelta;
    double time2massRatio;
    
    vector<VectorVar> zonalPresForces = gridMgr->getVectorVariableForAllNodes(NODE_TOT_PRES_FORCE);
    vector<VectorVar> velocitiesCurrent = gridMgr->getVectorVariableForAllNodes(VELOCITY_CTS);
    vector<ScalarVar> nodalMasses = gridMgr->getScalarVariableForAllNodes(NODAL_MASS);
    
    
    for ( ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
            vector<double> velocitiesNext, velocitiesHalf, velocitiesDelta;
            nodalMass = nodalMasses[ijNode].getValue();
            time2massRatio = loader->getTimeStep()/nodalMass;
            for( dimNum = 0; dimNum < loader->dim; dimNum++ ){
                velPrev   = velocitiesCurrent[ijNode].getValue()[dimNum];
                velDelta = zonalPresForces[ijNode].getValue()[dimNum]*time2massRatio;
                velNext = velPrev + velDelta;
                velocitiesNext.push_back(velNext);
                velocitiesDelta.push_back(velDelta);
                velocitiesHalf.push_back( 0.5*( velNext + velPrev ) );
            }
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_NTS  , velocitiesNext));
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_HTS  , velocitiesHalf));
            gridMgr->setVectorVariableForNode(ijNode, VectorVar(VELOCITY_DELTA, velocitiesDelta));
    }
    
}

void LagrangianSolver::updateNodeCoordinates(){
    // new mesh
    // Move nodes to their new positions according to n+1/2 velocities
    // z[n+1] = z[n] + ∆t * w[n+1/2]
	int dimNum;
    double coordNext, coordPrev, coordDelta ;
    vector<VectorVar> coordsPrev = gridMgr->getVectorVariableForAllNodes(COORDINATE);
    vector<VectorVar> velocitiesHTC = gridMgr->getVectorVariableForAllNodes(VELOCITY_HTS);
    double timeStep = loader->getTimeStep();
    for (int ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
          vector<double> coordinate;
          for( dimNum = 0; dimNum < loader->dim; dimNum++ ){
                coordPrev  = coordsPrev[ijNode].getValue()[dimNum];
                coordDelta = velocitiesHTC[ijNode].getValue()[dimNum]*timeStep;
                coordNext  = coordPrev + coordDelta;
                coordinate.push_back(coordNext);
          }
          gridMgr->setVectorVariableForNode(ijNode, VectorVar(COORDINATE_OLD, coordsPrev[ijNode].getValue()));
          gridMgr->setVectorVariableForNode(ijNode, VectorVar(COORDINATE, coordinate));
    }
}

//E[n] = - Σ[subzones](F[n]*w[n+1/2])
//Eint[n+1] = Eint[n] + ∆t/m * E[n]
void LagrangianSolver::calculateInternalEnergy(){
    // internal energy = - p∇w.
    // Compute the total work done in each zone due to the forces affecting its nodes and compute new internal energy due to this work
    
    int  dimNum, subzoneNum, neighbor;
    int ijZone;
    double zonalMass;
    double time2massRatio;
    vector<VectorVar> nodalTotPresForces = gridMgr->getVectorVariableForAllNodes(NODE_TOT_PRES_FORCE);
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<VectorVar> velocitiesHalf = gridMgr->getVectorVariableForAllNodes(VELOCITY_HTS);
    vector<ScalarVar> energies = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    vector<ScalarVar> zonalMasses = gridMgr->getScalarVariableForAllZones(ZONAL_MASS);
    vector<ScalarVar> divI = laserMgr->getIntensityDivergence();
    vector<ScalarVar> divQ = heatMgr->getFluxDivergence();
    vector<vector<int>> indicesOfNeighborNodes = gridMgr->getIndicesOfNeighborNodesForAllZones();
    VectorVar velHalf;
    VectorVar force;
    double FV;
    double erg, ergPrev;
    double intErg;
    for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
            zonalMass = zonalMasses[ijZone].getValue();
            time2massRatio = loader->getTimeStep()/zonalMass;
            ergPrev = energies[ijZone].getValue();
            erg = 0.0;
            for ( subzoneNum = 0; subzoneNum < loader->subzoneNum; subzoneNum++ ){
                neighbor = indicesOfNeighborNodes[ijZone][subzoneNum];
                velHalf = velocitiesHalf[neighbor];
                force   = subzones[ijZone][subzoneNum].getVectorVariable(SUBZONAL_PRES_FORCE);
                FV = 0.0;
                for( dimNum = 0; dimNum < loader->dim; dimNum++ ){
                    FV += force.getValue()[dimNum]*velHalf.getValue()[dimNum];
                }
                erg -= FV;
                if( ergPrev + erg * time2massRatio < 0.0 ){
                    string msg ="[LagrangianSolver] negative energy problem  ijZone = "+to_string(ijZone);
                    logger->writeMsg(msg.c_str(), CRITICAL);
                }
            }
            

            double laserDelta = fabs(divI[ijZone].getValue());
            erg += laserDelta;//2D
//            erg += divQ[ijZone].getValue(); todo heat manager switch on
            
            intErg = ergPrev + erg * time2massRatio;
            if( intErg < 0.0 ){
                string msg ="[LagrangianSolver] negative energy problem  ijZone = "+to_string(ijZone);
                logger->writeMsg(msg.c_str(), CRITICAL);
                intErg = 0.0;
            }
            gridMgr->setScalarVariableForZone(ijZone, ScalarVar(INTERNAL_ERG, intErg));
    }
}

// ρ=m/v
void LagrangianSolver::calculateDensity(){
    // density = - ∇w
    // Update zone  densities
    int  subzoneNum;
    double zonalMass, subzonalMass;
    double density, subzonalDensity;
    double zonalVolume, subzonalVolume;
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<ScalarVar> volumes = gridMgr->getScalarVariableForAllZones(ZONAL_VOLUME);
    vector<ScalarVar> zonalMasses = gridMgr->getScalarVariableForAllZones(ZONAL_MASS);
     for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        zonalMass   = 0.0;  zonalVolume = 0.0;
        for ( subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++ ){
                subzonalMass   = subzones[ijZone][subzoneNum].getVariable(SUBZONAL_MASS).getValue();
                subzonalVolume = subzones[ijZone][subzoneNum].getVariable(SUBZONAL_VOLUME).getValue();
                subzonalDensity = subzonalMass/subzonalVolume;
                gridMgr->setScalarVariableForSubzone(ijZone, subzoneNum, ScalarVar(DENSITY, subzonalDensity));
                zonalMass   += subzonalMass;
                zonalVolume += subzonalVolume;
        }
        if( !areSame(zonalMasses[ijZone].getValue(), zonalMass) ){
            logger->writeMsg("[LagrangianSolver] zonal mass problem ", CRITICAL);
        }
        if( !areSame(volumes[ijZone].getValue(), zonalVolume) ){
            logger->writeMsg("[LagrangianSolver] zonal volume problem ", CRITICAL);
        }
        density = zonalMass/zonalVolume;
        if( density <= 0.0 ){
            logger->writeMsg("[LagrangianSolver] negative density problem ", CRITICAL);
        }
        gridMgr->setScalarVariableForZone(ijZone, ScalarVar(DENSITY, density));
    }
}

// p = (γ − 1)ρE
void LagrangianSolver::calculatePressure(){
    double pres, dens, erg;
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<ScalarVar> energies  = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    for (int ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
        dens = densities[ijZone].getValue();
        erg  = energies[ijZone].getValue();
        pres = ( GAMMA - 1 )*dens*erg;
        gridMgr->setScalarVariableForZone(ijZone, ScalarVar(PRESURE, pres));
    }
}

void LagrangianSolver::applyPressureBC(){
    map<int, int> zonesBoundaryIndexes = gridMgr->getZonesGlobalBoundaryIndexes();
    vector<ScalarVar> presures  = gridMgr->getScalarVariableForAllZones(PRESURE);
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    double dens, pres, presCur;
    int zoneIDX, donorIDX;
    double mq = A/(Z+1)/Cpt;
    
                //        A       p
    double Tcur;// T =  ----- x ------
                //      Z + 1  C(p/T) ρ
    
    map<int, int>::iterator it;
    for( it = zonesBoundaryIndexes.begin(); it != zonesBoundaryIndexes.end(); it++) {
        zoneIDX  = it->first;
        donorIDX = it->second;
        dens = densities[donorIDX].getValue();
        pres = presures [donorIDX].getValue();
        Tcur = mq*pres/dens;
        if( Tcur > Tboil ){
            gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(PRESURE, BOUNDARY_PRESSURE));
        }else if( Tcur > Tmelt){
            presCur = BOUNDARY_PRESSURE+(Tcur-Tboil)/(Tmelt-Tboil)*(pres-BOUNDARY_PRESSURE);
            gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(PRESURE, presCur));
        }else{
            gridMgr->setScalarVariableForZone(zoneIDX, ScalarVar(PRESURE, pres));
        }
    }

}

void LagrangianSolver::checkConservativeLaws(){
    int ijZone, ijNode, subzoneNum, subIDX, dimNum;
    
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones(DENSITY);
    vector<ScalarVar> volumes = gridMgr->getScalarVariableForAllZones(ZONAL_VOLUME);
    vector<ScalarVar> energies = gridMgr->getScalarVariableForAllZones(INTERNAL_ERG);
    vector<VectorVar> velocities     = gridMgr->getVectorVariableForAllNodes(VELOCITY_CTS);
    vector<ScalarVar> nodalMasses = gridMgr->getScalarVariableForAllNodes(NODAL_MASS);
    vector<ScalarVar> zonalMasses = gridMgr->getScalarVariableForAllZones(ZONAL_MASS);
    
    vector<vector<Subzone>> subzones = gridMgr->getSubzonesForAllZones();
    vector<vector<int>> indicesOfNeighbors  = gridMgr->getIndicesOfNeighborZonesForAllNodes();
    
    double zoneMass, totZoneMass = 0.0, totZoneMass2 = 0.0, totZoneMass3 = 0.0,  subzonalMass, zoneMass1;
    double nodeMass, totNodeMass = 0.0, totNodeMass2 = 0.0 ;
    double totIntErg = 0.0;
    double vel, totKinErg = 0.0;
    double totErg = 0.0;

    for (ijZone = 0; ijZone < gridMgr->getTotalZoneNumber(); ijZone++){
            zoneMass = densities[ijZone].getValue()*volumes[ijZone].getValue();
            totZoneMass += zoneMass;
            totZoneMass2 += zonalMasses[ijZone].getValue();
            zoneMass1 = 0.0;
            for ( subzoneNum = 0; subzoneNum < loader->subzoneNum; subzoneNum++ ){
                subzonalMass = subzones[ijZone][subzoneNum].getVariable(SUBZONAL_MASS).getValue();
                zoneMass1 += subzonalMass;
            }
            if(!areSame(zoneMass1, zonalMasses[ijZone].getValue())){
                string msg ="[LagrangianSolver] negative zonal mass problem "+to_string(ijZone);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }
            totZoneMass3 += zoneMass1;
            if(energies[ijZone].getValue() < 0.0){
                string msg ="[LagrangianSolver] negative ENERGY problem "+to_string(ijZone);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }
            totIntErg += energies[ijZone].getValue()*zoneMass;

    }
    string msg ="[LagrangianSolver] totZoneMass_rhoV = "+to_string(totZoneMass);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totZoneMass_zone = "+to_string(totZoneMass2);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totZoneMass_sub = "+to_string(totZoneMass3);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totIntErg = "+to_string(totIntErg);
    logger->writeMsg(msg.c_str(), INFO);
    

    for (ijNode = 0; ijNode < gridMgr->getTotalNodeNumber(); ijNode++){
            totNodeMass += nodalMasses[ijNode].getValue();
            nodeMass = 0.0;
            for ( subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++ ){
                subIDX = indicesOfNeighbors[ijNode][subzoneNum];
                subzonalMass = subzones[subIDX][subzoneNum].getVariable(SUBZONAL_MASS).getValue();
                nodeMass += subzonalMass;
            }
            if(!areSame(nodeMass, nodalMasses[ijNode].getValue())){
                string msg ="[LagrangianSolver] NODE MASS PROBLEM ijNode = "+to_string(ijNode);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }
            totNodeMass2 += nodeMass;
            
            for( dimNum = 0; dimNum < loader->dim; dimNum++ ){
                vel = velocities[ijNode].getValue()[dimNum];
                totKinErg += 0.5*vel*vel*nodalMasses[ijNode].getValue();
            }

    }
    totErg  = totKinErg + totIntErg;
    
    msg ="[LagrangianSolver] totNodeMass_1 = "+to_string(totNodeMass);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totNodeMass_2 = "+to_string(totNodeMass2);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totKinErg = "+to_string(totKinErg);
    logger->writeMsg(msg.c_str(), INFO);
    msg ="[LagrangianSolver] totErg = "+to_string(totErg);
    logger->writeMsg(msg.c_str(), INFO);
}


LagrangianSolver::~LagrangianSolver(){
    finilize();
    logger->writeMsg("[LagrangianSolver] delete...OK");
}

void LagrangianSolver::finilize(){
    
}
