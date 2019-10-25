#include "GridManager.hpp"

using namespace std;
using namespace chrono;

GridManager::GridManager(shared_ptr<Loader> loader){
    this->loader=loader;
    logger.reset(new Logger());
    initialize();
    string msg ="[GridManager] init...OK {Node number:"+to_string(totalNodeNumber)+"}";
    logger->writeMsg(msg.c_str());
}

int GridManager::getTotalNodeNumber(){
    return totalNodeNumber;
}

int GridManager::getTotalZoneNumber(){
    return totalZoneNumber;
}


void GridManager::initialize(){
    logger->writeMsg("[GridManager] initialize() ...");

    this->totalNodeNumber = loader->resolution[0]*loader->resolution[1]*loader->resolution[2];
    initNodes();
    
    // TODO: change for 3D
    // this->totalZoneNumber = (loader->resolution[0]+1)*(loader->resolution[1]+1)*(loader->resolution[2]+1);
    this->totalZoneNumber = (loader->resolution[0]+1)*(loader->resolution[1]+1);
    for(int idx=0;idx<totalZoneNumber;idx++){
        zones.push_back(unique_ptr<Zone>(new Zone()));
        zones[idx]->setScalarVariable(ScalarVar(ZONAL_VOLUME, 0.0));
        vector<Subzone> subzones;
        for(int subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++ ){
            subzones.push_back(Subzone());
            subzones[subzoneNum].setVariable(ScalarVar(SUBZONAL_VOLUME, 0.0));
        }
        zones[idx]->setSubzones(subzones);
    }
    
    initZones();
    updateGrid();
    
    logger->writeMsg("[GridManager] initialize() ...OK");
}

int GridManager::updateGrid(){
    initGridWithGhostLayer();
    return calculateZoneCentroidsAndVolume();
}
void GridManager::initGridWithGhostLayer(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int xResExt = xRes+2, yResExt = yRes+2, zResExt = 1;
    int i,j,k = 0, ijNode;
    coordNodesAndGhost.clear();
    for(int tmp = 0; tmp < xResExt*yResExt; tmp++){
        coordNodesAndGhost.push_back(VectorVar(COORDINATE, {0.0, 0.0}));
    }
    
    int coordIDX;
    for ( i=1; i<xRes+1; i++){
        for ( j=1; j<yRes+1; j++){
            coordIDX = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            ijNode   = IDX(i-1, j-1, k, xRes  , yRes  , zRes);
            coordNodesAndGhost[coordIDX] = nodes[ijNode]->getVariable(COORDINATE);
        }
    }
    
    int boundaryNodeIDX, afterBoundaryNodeIDX;
    for ( i=1; i<xRes+1; i++){
        coordIDX              = IDX(i  , 0, k, xResExt, yResExt, zResExt);
        boundaryNodeIDX       = IDX(i-1, 0, k, xRes   , yRes   , zRes);
        afterBoundaryNodeIDX  = IDX(i-1, 1, k, xRes   , yRes   , zRes);
        coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
        
        coordIDX              = IDX(i  , yRes+1, k, xResExt, yResExt, zResExt);
        boundaryNodeIDX       = IDX(i-1, yRes-1, k, xRes  , yRes  , zRes);
        afterBoundaryNodeIDX  = IDX(i-1, yRes-2, k, xRes  , yRes  , zRes);
        coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
    }
    
    for ( j=1; j<yRes+1; j++){
        coordIDX              = IDX(0, j  , k, xResExt, yResExt, zResExt);
        boundaryNodeIDX       = IDX(0, j-1, k, xRes  , yRes  , zRes);
        afterBoundaryNodeIDX  = IDX(1, j-1, k, xRes  , yRes  , zRes);
        coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
        
        coordIDX              = IDX(xRes+1, j  , k, xResExt, yResExt, zResExt);
        boundaryNodeIDX       = IDX(xRes-1, j-1, k, xRes  , yRes  , zRes);
        afterBoundaryNodeIDX  = IDX(xRes-2, j-1, k, xRes  , yRes  , zRes);
        coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
    }
    
    coordIDX              = IDX(0, 0, k, xResExt, yResExt, zResExt);
    boundaryNodeIDX       = IDX(0, 0, k, xRes  , yRes  , zRes);
    afterBoundaryNodeIDX  = IDX(1, 1, k, xRes  , yRes  , zRes);
    coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
    
    coordIDX              = IDX(0, yRes+1, k, xResExt, yResExt, zResExt);
    boundaryNodeIDX       = IDX(0, yRes-1, k, xRes  , yRes  , zRes);
    afterBoundaryNodeIDX  = IDX(1, yRes-2, k, xRes  , yRes  , zRes);
    coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
    
    coordIDX              = IDX(xRes+1, 0, k, xResExt, yResExt, zResExt);
    boundaryNodeIDX       = IDX(xRes-1, 0, k, xRes  , yRes  , zRes);
    afterBoundaryNodeIDX  = IDX(xRes-2, 1, k, xRes  , yRes  , zRes);
    coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
    
    coordIDX              = IDX(xRes+1, yRes+1, k, xResExt, yResExt, zResExt);
    boundaryNodeIDX       = IDX(xRes-1, yRes-1, k, xRes  , yRes  , zRes);
    afterBoundaryNodeIDX  = IDX(xRes-2, yRes-2, k, xRes  , yRes  , zRes);
    coordNodesAndGhost[coordIDX] = VectorVar(COORDINATE,  getReflectedCoordinate(boundaryNodeIDX, afterBoundaryNodeIDX ));
}


//z[i,j] = 1/8 (4 z[i,j] + z[i−1,j] + z[i+1,j] + z[i,j−1] + z[i,j+1])
void GridManager::smoothMesh(){
    auto start_time = high_resolution_clock::now();
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int xResExt = xRes+2, yResExt = yRes+2, zResExt = 1;
    int i,j,k = 0, ijNode;
    int i_1jNode,i1jNode, ij_1Node,ij1Node;
    int i1j1Node, i_1j1Node, i1j_1Node, i_1j_1Node;
    double ijNodeCoord, i_1jNodeCoord, i1jNodeCoord, ij_1NodeCoord, ij1NodeCoord;
    double i1j1NodeCoord, i_1j1NodeCoord, i1j_1NodeCoord, i_1j_1NodeCoord;
    vector<double> xDerivative = {0.0, 0.0}, yDerivative = {0.0, 0.0};
    double alpha, beta, gamma;
    double smoothedCoordinate;
    
    
    for ( i=1; i<xRes+1; i++){
        for ( j=1; j<yRes+1; j++){
            ijNode   = IDX(i  , j  , k, xResExt, yResExt, zResExt);
            i_1jNode = IDX(i-1, j  , k, xResExt, yResExt, zResExt);
            i1jNode  = IDX(i+1, j  , k, xResExt, yResExt, zResExt);
            ij_1Node = IDX(i  , j-1, k, xResExt, yResExt, zResExt);
            ij1Node  = IDX(i  , j+1, k, xResExt, yResExt, zResExt);
            
            i_1j_1Node = IDX(i-1, j-1, k, xResExt, yResExt, zResExt);
            i1j_1Node  = IDX(i+1, j-1, k, xResExt, yResExt, zResExt);
            i_1j1Node  = IDX(i-1, j+1, k, xResExt, yResExt, zResExt);
            i1j1Node   = IDX(i+1, j+1, k, xResExt, yResExt, zResExt);

            for(int tmp = 0; tmp < 2; tmp++){//todo correct dimension number
                
                i_1jNodeCoord = coordNodesAndGhost[i_1jNode].getValue()[tmp];
                i1jNodeCoord  = coordNodesAndGhost[i1jNode].getValue()[tmp];
                xDerivative[tmp] = 0.5*(i1jNodeCoord-i_1jNodeCoord);
                
                ij_1NodeCoord = coordNodesAndGhost[ij_1Node].getValue()[tmp];
                ij1NodeCoord  = coordNodesAndGhost[ij1Node].getValue()[tmp];
                yDerivative[tmp] = 0.5*(ij1NodeCoord-ij_1NodeCoord);
            }
            alpha = xDerivative[0]*xDerivative[0]+xDerivative[1]*xDerivative[1];
            beta  = xDerivative[0]*yDerivative[0]+xDerivative[1]*yDerivative[1];
            gamma = yDerivative[0]*yDerivative[0]+yDerivative[1]*yDerivative[1];
            
            vector<double> coordinate;
            for(int tmp = 0; tmp < 2; tmp++){//todo correct dimension number
                ijNodeCoord   = coordNodesAndGhost[ijNode].getValue()[tmp];
                i_1jNodeCoord = coordNodesAndGhost[i_1jNode].getValue()[tmp];
                i1jNodeCoord  = coordNodesAndGhost[i1jNode].getValue()[tmp];
                ij_1NodeCoord = coordNodesAndGhost[ij_1Node].getValue()[tmp];
                ij1NodeCoord  = coordNodesAndGhost[ij1Node].getValue()[tmp];
                
                i_1j_1NodeCoord = coordNodesAndGhost[i_1j_1Node].getValue()[tmp];
                i1j_1NodeCoord  = coordNodesAndGhost[i1j_1Node].getValue()[tmp];
                i_1j1NodeCoord  = coordNodesAndGhost[i_1j1Node].getValue()[tmp];
                i1j1NodeCoord   = coordNodesAndGhost[i1j1Node].getValue()[tmp];
                
//                smoothedCoordinate = 0.125*(4*ijNodeCoord+i_1jNodeCoord+i1jNodeCoord+ij_1NodeCoord+ij1NodeCoord);
                
                smoothedCoordinate = 0.5*(alpha*(ij_1NodeCoord+ij1NodeCoord)+gamma*(i_1jNodeCoord+i1jNodeCoord)-0.5*beta*(i_1j_1NodeCoord-i1j_1NodeCoord-i_1j1NodeCoord+i1j1NodeCoord))/(alpha+gamma);
                
                coordinate.push_back(smoothedCoordinate);
            }
            
//            nodes[ijNode]->setVectorVariable(VectorVar(COORDINATE_SMOOTHED, coordinate));
            ijNode   = IDX(i-1, j-1, k, xRes, yRes, zRes);
            nodes[ijNode]->setVectorVariable(VectorVar(COORDINATE, coordinate));
        }
    }
    auto end_time = high_resolution_clock::now();
    
    string msg ="[GridManager] Smoothing duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());

}

vector<VectorVar> GridManager::getExtendedGridCoordinates(){
    return coordNodesAndGhost;

}
vector<double> GridManager::getReflectedCoordinate(int boundaryNodeIDX, int afterBoundaryNodeIDX){
    vector<double> coord = {0.0, 0.0};
    double boundaruNodeCoord, afterBoundaruNodeCoord;
    for(int tmp = 0; tmp < 2; tmp++){
        boundaruNodeCoord      = nodes[boundaryNodeIDX]->getVariable(COORDINATE).getValue()[tmp];
        afterBoundaruNodeCoord = nodes[afterBoundaryNodeIDX]->getVariable(COORDINATE).getValue()[tmp];
        coord[tmp] = 2*boundaruNodeCoord - afterBoundaruNodeCoord;
    }
    return coord;
}


void GridManager::initNodes(){
	auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node start ");
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResExt = xRes+1, yResExt = yRes+1, zResExt = 1;
    double x,y,z;
    int i,j,k,idx;
    for ( i=0; i<xRes; i++){
        for ( j=0; j<yRes; j++){
        	for ( k=0; k<zRes; k++){
            idx = IDX(i,j,k,xRes,yRes,zRes);
            
            x = i*loader->spatialSteps[0] + loader->boxCoordinates[0][0];
            y = j*loader->spatialSteps[1] + loader->boxCoordinates[1][0];
            z = k*loader->spatialSteps[2] + loader->boxCoordinates[2][0];
            
            map<string,VectorVar> variables;
            variables[COORDINATE] = VectorVar(COORDINATE, {x, y, z});
            
            nodes.push_back(unique_ptr<Node>(new Node(variables)));
            vector<int> indicesOfNeighborZones;
            
            /**
             *     [1]  x Zone(i,j+1) - - - - - - - - - x Zone(i+1,j+1) [0]
             *          .               |               .
             *          |               |               |
             *          .               |               .
             *          |               |               |
             *          .  Subzone 2    |   Subzone 1   .
             *          |               |               |
             *          .               |               .
             *          |               |               |
             *          ________ _______.Node(i,j)_______
             *          |               |        |      |
             *          .               |               .
             *          |               |               |
             *          .    Subzone 3  |   Subzone 4   .
             *          |               |               |
             *          .               |               .
             *          |               |               |
             *          .               |               .
             *     [2]  x Zone(i,j) - -   - - - - - - - x Zone(i+1,j)  [3]
             **/
            //order is important for the pressure force calculation
            indicesOfNeighborZones.push_back(IDX(i+1, j+1 , k, xResExt, yResExt, zResExt));
            indicesOfNeighborZones.push_back(IDX(i  , j+1 , k, xResExt, yResExt, zResExt));
            indicesOfNeighborZones.push_back(IDX(i  , j   , k, xResExt, yResExt, zResExt));
            indicesOfNeighborZones.push_back(IDX(i+1, j   , k, xResExt, yResExt, zResExt));
            
            
            nodes[idx]->setIndicesOfNeighborZones(indicesOfNeighborZones);
        	}
        }
    }
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
}

void GridManager::initZones(){
	auto start_time = high_resolution_clock::now();
    int i,j,k=0,idxZone, idxDonorZone, idxNode;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResExt = xRes+1, yResExt = yRes+1, zResExt = 1;

    for ( i=1; i < xRes; i++){
        for ( j=1; j < yRes; j++){
            
            idxZone = IDX(i,j,k, xResExt, yResExt, zResExt);
            idxNode = IDX(i,j,k, xRes   , yRes   , zRes   );
            zones[idxZone]->setNodeIDX(idxNode);
            zonesInBoundaryIndexes.push_back(idxZone);
            
           
            /**
             *     [3]   .Node(i-1,j)___________________.Node(i,j)   [2]
             *           |               |(x2, y2)       |
             *           |               .               |
             *           |               |               |
             *           |   Subzone 4   .    Subzone 3  |
             *           |               |               |
             *           |               .               |
             *           |               |               |
             *           |               |               |
             *           |- - - - - - - -x Zone(i,j) - - |
             *           |               |               |
             *           |               .               |
             *           |               |               |
             *           |               .               |
             *           |   Subzone 1   |   Subzone 2   |
             *           |               .               |
             *           |               |               |
             *     [0]   .Node(i-1,j-1)__________________.Node(i,j-1) [1]
             *
             *
             *
             **/
            vector<int> indicesOfNeighborNodes;
            indicesOfNeighborNodes.push_back(IDX(i-1, j-1 ,k,xRes,yRes,zRes));
            indicesOfNeighborNodes.push_back(IDX(i  , j-1 ,k,xRes,yRes,zRes));
            indicesOfNeighborNodes.push_back(IDX(i  , j   ,k,xRes,yRes,zRes));
            indicesOfNeighborNodes.push_back(IDX(i-1, j   ,k,xRes,yRes,zRes));
            zones[idxZone]->setIndicesOfNeighborNodes(indicesOfNeighborNodes);
            
            vector<int> neighborhoodIndices;
            neighborhoodIndices.push_back(IDX(i-1, j-1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i  , j-1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i+1, j-1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i+1, j   ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i+1, j+1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i  , j+1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i-1, j+1 ,k, xResExt, yResExt, zResExt));
            neighborhoodIndices.push_back(IDX(i-1, j   ,k, xResExt, yResExt, zResExt));
            zones[idxZone]->setNeighborhoodIndices(neighborhoodIndices);
        }
    }

    int a, b, c;
    for ( i = 1; i < xRes;i++){
        j = 0;
        idxZone       = IDX(i,j  ,k, xResExt, yResExt, zResExt);
        idxDonorZone  = IDX(i,j+1,k, xResExt, yResExt, zResExt);
        idxNode       = IDX(i,j  ,k,xRes  ,yRes  ,zRes);
        zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
          //-1 -> left, bottom, back
          // 0 -> same position
          //+1 -> right, top, front
        a = 0, b = -1, c = 0;
        if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
        	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
        }
        zonesBoundaryIndexes[idxZone] =  idxDonorZone;
        zones[idxZone]->setNodeIDX(idxNode);
        vector<int> neighborhoodIndicesDown;
        neighborhoodIndicesDown.push_back(idxZone);
        neighborhoodIndicesDown.push_back(idxZone);
        neighborhoodIndicesDown.push_back(idxZone);
        neighborhoodIndicesDown.push_back(IDX(i+1, j   ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesDown.push_back(IDX(i+1, j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesDown.push_back(IDX(i  , j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesDown.push_back(IDX(i-1, j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesDown.push_back(IDX(i-1, j   ,k, xResExt, yResExt, zResExt));
        zones[idxZone]->setNeighborhoodIndices(neighborhoodIndicesDown);
        
        j = yRes;
        idxZone       = IDX(i,j  ,k, xResExt, yResExt, zResExt);
        idxDonorZone  = IDX(i,j-1,k, xResExt, yResExt, zResExt);
        idxNode       = IDX(i,j-1,k,xRes  ,yRes  ,zRes);
        zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
        a = 0, b = 1, c = 0;
        if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
        	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
        }
        zonesBoundaryIndexes[idxZone] =  idxDonorZone;
        zones[idxZone]->setNodeIDX(idxNode);
        vector<int> neighborhoodIndicesUp;
        neighborhoodIndicesUp.push_back(IDX(i-1, j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesUp.push_back(IDX(i  , j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesUp.push_back(IDX(i+1, j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesUp.push_back(IDX(i+1, j   ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesUp.push_back(idxZone);
        neighborhoodIndicesUp.push_back(idxZone);
        neighborhoodIndicesUp.push_back(idxZone);
        neighborhoodIndicesUp.push_back(IDX(i-1, j   ,k, xResExt, yResExt, zResExt));
        zones[idxZone]->setNeighborhoodIndices(neighborhoodIndicesUp);
    }
    
    for ( j = 1; j < yRes;j++){
        i = 0;
        idxZone       = IDX(i  ,j,k, xResExt, yResExt, zResExt);
        idxDonorZone  = IDX(i+1,j,k, xResExt, yResExt, zResExt);
        idxNode       = IDX(i  ,j,k,xRes  ,yRes  ,zRes);
        zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
        a = -1, b = 0, c = 0;
        if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
        	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
        }
        zonesBoundaryIndexes[idxZone] =  idxDonorZone;
        zones[idxZone]->setNodeIDX(idxNode);
        vector<int> neighborhoodIndicesLeft;
        neighborhoodIndicesLeft.push_back(idxZone);
        neighborhoodIndicesLeft.push_back(IDX(i  , j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesLeft.push_back(IDX(i+1, j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesLeft.push_back(IDX(i+1, j   ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesLeft.push_back(IDX(i+1, j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesLeft.push_back(IDX(i  , j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesLeft.push_back(idxZone);
        neighborhoodIndicesLeft.push_back(idxZone);
        zones[idxZone]->setNeighborhoodIndices(neighborhoodIndicesLeft);
        
        i = xRes;
        idxZone       = IDX(i  ,j,k, xResExt, yResExt, zResExt);
        idxDonorZone  = IDX(i-1,j,k, xResExt, yResExt, zResExt);
        idxNode       = IDX(i-1,j,k,xRes  ,yRes  ,zRes);
        zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
        a = 1, b = 0, c = 0;
        if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
        	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
        }
        zonesBoundaryIndexes[idxZone] =  idxDonorZone;
        zones[idxZone]->setNodeIDX(idxNode);
        vector<int> neighborhoodIndicesRight;
        neighborhoodIndicesRight.push_back(IDX(i-1, j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesRight.push_back(IDX(i  , j-1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesRight.push_back(idxZone);
        neighborhoodIndicesRight.push_back(idxZone);
        neighborhoodIndicesRight.push_back(idxZone);
        neighborhoodIndicesRight.push_back(IDX(i  , j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesRight.push_back(IDX(i-1, j+1 ,k, xResExt, yResExt, zResExt));
        neighborhoodIndicesRight.push_back(IDX(i-1, j   ,k, xResExt, yResExt, zResExt));
        zones[idxZone]->setNeighborhoodIndices(neighborhoodIndicesRight);
    }
    
    idxZone       = IDX(0, 0, k, xResExt, yResExt, zResExt);
    idxDonorZone  = IDX(1, 1, k, xResExt, yResExt, zResExt);
    idxNode       = IDX(0, 0, k, xRes  ,yRes  ,zRes);
    zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
    a = -1, b = -1, c = 0;
    if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
    	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
    }
    zonesBoundaryIndexes[idxZone] =  idxDonorZone;
    zones[idxZone]->setNodeIDX(idxNode);
    
    idxZone       = IDX(0, yRes  , k, xResExt, yResExt, zResExt);
    idxDonorZone  = IDX(1, yRes-1, k, xResExt, yResExt, zResExt);
    idxNode       = IDX(0, yRes-1, k, xRes  ,yRes  ,zRes);
    zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
    a = -1, b = 1, c = 0;
    if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
    	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
    }
    zonesBoundaryIndexes[idxZone] =  idxDonorZone;
    zones[idxZone]->setNodeIDX(idxNode);
    
    idxZone       = IDX(xRes  , 0, k, xResExt, yResExt, zResExt);
    idxDonorZone  = IDX(xRes-1, 1, k, xResExt, yResExt, zResExt);
    idxNode       = IDX(xRes-1, 0, k, xRes  ,yRes  ,zRes);
    zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
    a = 1, b = -1, c = 0;
    if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
    	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
    }
    zonesBoundaryIndexes[idxZone] =  idxDonorZone;
    zones[idxZone]->setNodeIDX(idxNode);
    
    idxZone       = IDX(xRes  , yRes  , k, xResExt, yResExt, zResExt);
    idxDonorZone  = IDX(xRes-1, yRes-1, k, xResExt, yResExt, zResExt);
    idxNode       = IDX(xRes-1, yRes-1, k, xRes  ,yRes  ,zRes);
    zones[idxZone]->setIndicesOfNeighborNodes({idxNode, idxNode, idxNode, idxNode});
    a = 1, b = 1, c = 0;
    if(loader->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] == MPI_PROC_NULL){
    	zonesGlobalBoundaryIndexes[idxZone] =  idxDonorZone;
    }
    zonesBoundaryIndexes[idxZone] =  idxDonorZone;
    zones[idxZone]->setNodeIDX(idxNode);

    //init spiral indexes
    int rowStart=1;
        int rowLength=xRes-1;
        int colStart=1;
        int colLength = yRes-1;
        while(rowStart <= rowLength && colStart <= colLength){
            for (int i = rowStart; i <= colLength; i++) {
            	idxZone= IDX(rowStart,i,0,xRes+1,yRes+1,1);
            	zonesInBoundaryIndexesSpiral.push_back(idxZone);
            }

            for (int j = rowStart+1; j <= rowLength; j++) {
            	idxZone= IDX(j,colLength,0,xRes+1,yRes+1,1);
            	zonesInBoundaryIndexesSpiral.push_back(idxZone);
            }
            if(rowStart+1 <= rowLength){
                for (int k = colLength-1; k >= colStart; k--) {
                	idxZone= IDX(rowLength,k,0,xRes+1,yRes+1,1);
                	zonesInBoundaryIndexesSpiral.push_back(idxZone);
                }
            }
            if(colStart+1 <= colLength){
                for (int k = rowLength-1; k > rowStart; k--) {
                	idxZone= IDX(k,colStart,0,xRes+1,yRes+1,1);
                	zonesInBoundaryIndexesSpiral.push_back(idxZone);
                }
            }
            rowStart++;
            rowLength--;
            colStart++;
            colLength--;
        }
        auto end_time = high_resolution_clock::now();
        string msg ="[GridManager] init zones duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
        logger->writeMsg(msg.c_str());
}



int GridManager::calculateZoneCentroidsAndVolume(){
    auto start_time = high_resolution_clock::now();
    
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int xResExt = xRes+2, yResExt = yRes+2, zResExt = 1;
    logger->writeMsg("[GridManager] start calculateZoneCentroidsAndVolume()...");
    
    int i,j,k=0, ijZoneIDX;
    int subzoneNum, verticeNum;
    double xleft, xright, yleft, yright, dx, dy;
    double x1,y1,x2,y2,x3,y3;
    double xc ,yc;
    double xc1,yc1;
    double Vc,Vc1,Vc_test;
    double subZonalVolume,subZonalVolume1;
    
    vector<int> idxRight;
    vector<int> idxLeft;
    vector<double> subzoneXcoordRight;
    vector<double> subzoneXcoordLeft;
    vector<double> subzoneYcoordRight;
    vector<double> subzoneYcoordLeft;
    
    vector<Subzone> subzones;
    int ijIDX, i_1jIDX, ij_1IDX, i_1j_1IDX;
    for ( i=0; i < xRes+1; i++){
        for ( j=0; j < yRes+1; j++){
            
            ijZoneIDX = IDX( i , j , k ,xRes+1, yRes+1, zRes);
           
            ijIDX     = IDX( i+1, j+1, k, xResExt, yResExt, zResExt);
            i_1jIDX   = IDX( i  , j+1, k, xResExt, yResExt, zResExt);
            ij_1IDX   = IDX( i+1, j  , k, xResExt, yResExt, zResExt);
            i_1j_1IDX = IDX( i  , j  , k, xResExt, yResExt, zResExt);
            xc = 0.0; yc = 0.0; xc1 = 0.0; yc1 = 0.0; Vc  = 0.0;  Vc1  = 0.0; Vc_test = 0.0;
            
            idxLeft  = {i_1j_1IDX, ij_1IDX, ijIDX, i_1jIDX};
            idxRight = {ij_1IDX, ijIDX, i_1jIDX, i_1j_1IDX};
            

            for ( verticeNum=0; verticeNum<4; verticeNum++ ){
                xleft  = coordNodesAndGhost[idxLeft [verticeNum]].getValue()[0];
                xright = coordNodesAndGhost[idxRight[verticeNum]].getValue()[0];
                yleft  = coordNodesAndGhost[idxLeft [verticeNum]].getValue()[1];
                yright = coordNodesAndGhost[idxRight[verticeNum]].getValue()[1];
                
                dx = xright-xleft;
                dy = yright-yleft;
                
                Vc  +=  0.5*(xright+xleft)*dy;
                Vc1 += -0.5*(yright+yleft)*dx;
                
                xc  +=  (xright*xright+xright*xleft+xleft*xleft)*dy/6;// ∫xdxdy =  ∮x^2dy  =  Σ1/6(x1^2 + x1x2 + x2^2)(y2-y1)
                yc  += -(yright*yright+yright*yleft+yleft*yleft)*dx/6;// ∫ydxdy = -∮y^2dx  = -Σ1/6(y1^2 + y1y2 + y2^2)(x2-x1)
            }
            
            xc = xc/Vc;
            yc = yc/Vc;
            
            if( Vc < EPSILON_VOLUME ){
                string msg ="[GridManager] Zonal volume problem i = "+to_string(i)+"; j = "+to_string(j)+"; Vc = "+to_string(Vc)+"; Vc1 = "+to_string(Vc1);
                logger->writeMsg(msg.c_str(), CRITICAL);
                return GRID_FAIL;
            }
            if( !areSame(Vc, Vc1) ){
                string msg ="[GridManager] Volume problem i = "+to_string(i)+"; j = "+to_string(j)+"; Vc = "+to_string(Vc)+"; Vc1 = "+to_string(Vc1);
                logger->writeMsg(msg.c_str(), CRITICAL);
                return GRID_FAIL;
            }
            zones[ijZoneIDX]->setVectorVariable(VectorVar(ZONAL_CENTROID, { xc, yc }));
            zones[ijZoneIDX]->setScalarVariable(ScalarVar(ZONAL_VOLUME, Vc));
            

            /**
             *            .Node(i-1,j)___________________.Node(i,j) (x1, y1)
             *           |               |(x2, y2)       |
             *           |               .               |
             *           |               |               |
             *           |   Subzone 3   .    Subzone 2  |
             *           |               |               |
             *           |               .               |
             *           |               |               |
             *           |            (xc, yc)           |
             * (x3, y3)  |- - - - - - - -x Zone(i,j) - - |(x3, y3)
             *           |               |               |
             *           |               .               |
             *           |               |               |
             *           |               .               |
             *           |   Subzone 0   |   Subzone 1   |
             *           |               .               |
             *           |               |               |
             *           .Node(i-1,j-1)__________________.Node(i,j-1)
             *        (x1, y1)        (x2, y2)
             *
             *
             **/
            
            subzones = zones[ijZoneIDX]->getSubzones();
            for ( subzoneNum=0; subzoneNum<loader->subzoneNum; subzoneNum++){
                x1 = coordNodesAndGhost[idxLeft[subzoneNum]].getValue()[0];
                y1 = coordNodesAndGhost[idxLeft[subzoneNum]].getValue()[1];
                
                x2 = 0.5*( coordNodesAndGhost[idxRight[subzoneNum]].getValue()[0] + x1 );
                y2 = 0.5*( coordNodesAndGhost[idxRight[subzoneNum]].getValue()[1] + y1 );
                
                int index = (subzoneNum == 0) ? loader->subzoneNum - 1 : subzoneNum - 1;
                
                x3 = 0.5*(coordNodesAndGhost[idxRight[index]].getValue()[0]+coordNodesAndGhost[idxLeft[index]].getValue()[0]);
                y3 = 0.5*(coordNodesAndGhost[idxRight[index]].getValue()[1]+coordNodesAndGhost[idxLeft[index]].getValue()[1]);
                
                subzoneXcoordLeft  = {x1, x2, xc, x3};
                subzoneYcoordLeft  = {y1, y2, yc, y3};
                
                subzoneXcoordRight = {x2, xc, x3, x1};
                subzoneYcoordRight = {y2, yc, y3, y1};
                
                subZonalVolume  = 0.0;
                subZonalVolume1 = 0.0;
                
                subzones[subzoneNum].setVerticesCoordinates({ subzoneXcoordLeft, subzoneYcoordLeft });
                
               
                for ( verticeNum=0; verticeNum<4; verticeNum++ ){
                    xleft  = subzoneXcoordLeft [verticeNum];
                    xright = subzoneXcoordRight[verticeNum];
                    yleft  = subzoneYcoordLeft [verticeNum];
                    yright = subzoneYcoordRight[verticeNum];
                    dy = yright-yleft;
                    dx = xright-xleft;
                    subZonalVolume  +=  0.5*(xright+xleft)*dy;
                    subZonalVolume1 += -0.5*(yright+yleft)*dx;
                }
                if( subZonalVolume < 0.0 ){
                    string msg ="[GridManager] Subzonal volume problem i = "+to_string(i)+"; j = "+to_string(j)+"; subZonalVolume = "+to_string(subZonalVolume);
                    logger->writeMsg(msg.c_str(), CRITICAL);
                    return GRID_FAIL;
                }
                
                if( !areSame(subZonalVolume, subZonalVolume1) ){
                    string msg ="[GridManager] Subzonal volume problem i = "+to_string(i)+"; j = "+to_string(j)+"; subZonalVolume = "+to_string(subZonalVolume)+"; subZonalVolume1 = "+to_string(subZonalVolume1);
                    logger->writeMsg(msg.c_str(), CRITICAL);
                    return GRID_FAIL;
                }
                
                Vc_test += subZonalVolume;

                subzones[subzoneNum].setVariable(ScalarVar(SUBZONAL_VOLUME, subZonalVolume));
                
            }
            zones[ijZoneIDX]->setSubzones(subzones);
            
            if( !areSame(Vc, Vc_test) ){
                string msg ="[GridManager] Volume problem i = "+to_string(i)+"; j = "+to_string(j)+"; Vc = "+to_string(Vc)+"; Vc_test = "+to_string(Vc_test);
                logger->writeMsg(msg.c_str(), CRITICAL);
                return GRID_FAIL;
            }
            
        }
    }
    
    logger->writeMsg("[GridManager] calculateZoneCentroidsAndVolume()...OK");
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] calculateZoneCentroidsAndVolume duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
    return GRID_OK;
}


void GridManager::setSubzones(int zoneIDX, vector<Subzone> subzones){
    zones[zoneIDX]->setSubzones(subzones);
}


vector<vector<ScalarVar>> GridManager::getScalarVariablesForAllNodes(){
    vector<vector<ScalarVar>> result;
    result.reserve(totalNodeNumber);
    int i,j,k = 0,idx;
    for ( i=0; i<loader->resolution[0]; i++){
        for ( j=0; j<loader->resolution[1]; j++){
            idx = IDX(i,j,k,loader->resolution[0]+1,loader->resolution[1]+1,loader->resolution[2]);
            result.push_back(zones[idx]->getAllScalarVariables());
        }
    }
    return result;
}


void GridManager::sendRecvIndecis4MPI(){
	int a, b, c;
	int t;
	int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];

	for(int i=0;i<27;i++){
		sendStartInd_X.push_back(0);
		sendEndInd_X.push_back(0);
		recvStartInd_X.push_back(0);
		recvEndInd_X.push_back(0);
		sendStartInd_Y.push_back(0);
		sendEndInd_Y.push_back(0);
		recvStartInd_Y.push_back(0);
		recvEndInd_Y.push_back(0);
		sendStartInd_Z.push_back(0);
		sendEndInd_Z.push_back(0);
		recvStartInd_Z.push_back(0);
		recvEndInd_Z.push_back(0);
	}

	for (a = -1; a <= 1; a++) {
		for (b = -1; b <= 1; b++) {
			for (c = -1; c <= 1; c++) {
				t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));

				switch (a) {
				case -1:
					sendStartInd_X[t] = 1;
					sendEndInd_X[t] = 1;
					recvStartInd_X[t] = xRes + 1;
					recvEndInd_X[t] = xRes + 1;
					break;

				case 0:
					sendStartInd_X[t] = 1;
					sendEndInd_X[t] = xRes;
					recvStartInd_X[t] = 1;
					recvEndInd_X[t] = xRes;
					break;

				case 1:
					sendStartInd_X[t] = xRes;
					sendEndInd_X[t] = xRes;
					recvStartInd_X[t] = 0;
					recvEndInd_X[t] = 0;
					break;
				} // end switch a

				switch (b) {
				case -1:
					sendStartInd_Y[t] = 1;
					sendEndInd_Y[t] = 1;
					recvStartInd_Y[t] = yRes + 1;
					recvEndInd_Y[t] = yRes + 1;
					break;

				case 0:
					sendStartInd_Y[t] = 1;
					sendEndInd_Y[t] = yRes;
					recvStartInd_Y[t] = 1;
					recvEndInd_Y[t] = yRes;
					break;

				case 1:
					sendStartInd_Y[t] = yRes;
					sendEndInd_Y[t] = yRes;
					recvStartInd_Y[t] = 0;
					recvEndInd_Y[t] = 0;
					break;
				} // end switch b

				switch (c) {
				case -1:
					sendStartInd_Z[t] = 1;
					sendEndInd_Z[t] = 1;
					recvStartInd_Z[t] = zRes + 1;
					recvEndInd_Z[t] = zRes + 1;
					break;

				case 0:
					sendStartInd_Z[t] = 1;
					sendEndInd_Z[t] = zRes;
					recvStartInd_Z[t] = 1;
					recvEndInd_Z[t] = zRes;
					break;

				case 1:
					sendStartInd_Z[t] = zRes;
					sendEndInd_Z[t] = zRes;
					recvStartInd_Z[t] = 0;
					recvEndInd_Z[t] = 0;
					break;
				}
			}
		}
	}


}

void GridManager::sendBoundary2Neighbor(){
	int t;
	int i, j, k;
	int is0, is1, js0, js1, ks0, ks1;
	int ir0, ir1, jr0, jr1, kr0, kr1;
	int cpt = 0;
	MPI_Status st;
	MPI_Datatype mpighost;
	int counter[27];
	for (t = 0; t < 27; t++) {
		is0 = sendStartInd_X[t];
		is1 = sendEndInd_X[t];
		js0 = sendStartInd_Y[t];
		js1 = sendEndInd_Y[t];
		ks0 = sendStartInd_Z[t];
		ks1 = sendEndInd_Z[t];
			counter[t] = (is1 - is0 + 1) * (js1 - js0 + 1) * (ks1 - ks0 + 1);
	}
	int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResExt = xRes+1, yResExt = yRes+1, zResExt = 1;//todo for z
	double *sendBuf[27];
	double *recvBuf[27];

	int count;
	int blocklength[1];
	MPI_Aint displacements[1];
	MPI_Datatype types[1];
	MPI_Aint startaddress, tmpaddress;
	count = 1;

	blocklength[0] = 1;
	displacements[0] = 0;
	types[0] = MPI_DOUBLE;
	MPI_Type_create_struct(count, blocklength, displacements, types, &mpighost);
	MPI_Type_commit(&mpighost);

	for (t = 0; t < 27; t++) {
		sendBuf[t] = new double[counter[t]];
		recvBuf[t] = new double[counter[t]];
	}
    
    int scalarVarsNum = zones[0]->getAllScalarVariables().size();
    scalarVarsNum = 1;//todo check for
    int idxZone;
	for (t = 0; t < 27; t++) {
		cpt = 0;

		is0 = sendStartInd_X[t];
		is1 = sendEndInd_X[t];
		js0 = sendStartInd_Y[t];
		js1 = sendEndInd_Y[t];
		ks0 = sendStartInd_Z[t];
		ks1 = sendEndInd_Z[t];
			// fill the buffer
			for (i = is0; i <= is1; i++) {
				for (j = js0; j <= js1; j++) {
					for (k = ks0; k <= ks1; k++) {
                        double* variablesToSend = new double[scalarVarsNum];
                        idxZone = IDX(i,j,k, xResExt, yResExt, zResExt);
                        for (int idx_var=0;idx_var<scalarVarsNum;idx_var++) {
                            sendBuf[t][cpt] = zones[idxZone]->getAllScalarVariables()[idx_var].getValue();
                        }
                        delete [] variablesToSend;
						cpt++;
					}
				}
			}
		}

		for (t = 0; t < 27; t++) {
			if (t != 13)
					{
				MPI_Sendrecv(sendBuf[t],         // send buffer
						counter[t],          // # of ghost nodes to send
						mpighost,                  // MPI datatype
						loader->neighbors2Send[t],// destination process,
						t,                         // message send tag,
						recvBuf[t],         // receiving buffer
						counter[t],          // # of ghost nodes to recv
						mpighost,                  // recv datatype
						loader->neighbors2Recv[t], // reception process,
						t,                         // message recv tag
						MPI_COMM_WORLD, &st);

			}
		}

		for (t = 0; t < 27; t++) {
			if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
				cpt = 0;

				ir0 = recvStartInd_X[t];
				ir1 = recvEndInd_X[t];
				jr0 = recvStartInd_Y[t];
				jr1 = recvEndInd_Y[t];
				kr0 = recvStartInd_Z[t];
				kr1 = recvEndInd_Z[t];

				for (i = ir0; i <= ir1; i++) {
					for (j = jr0; j <= jr1; j++) {
						for (k = kr0; k <= kr1; k++) {
                            double dataRecv = recvBuf[t][cpt];
                            idxZone = IDX(i,j,k, xResExt, yResExt, zResExt);
                            for (int idx_var=0;idx_var<scalarVarsNum;idx_var++) {
                               string varName = zones[idxZone]->getAllScalarVariables()[idx_var].getName();
                                double var = dataRecv;
                                zones[idxZone]->setScalarVariable(ScalarVar(varName, var));
                            }
							cpt++;
						}
					}
				}
			}
		}

		for (t = 0; t < 27; t++) {
			delete [] sendBuf[t];
			delete [] recvBuf[t];
		}
}

vector<vector<VectorVar>> GridManager::getVectorVariablesForAllNodes(){
    vector<vector<VectorVar>> result;
    result.reserve(totalNodeNumber);
    for (int idx=0; idx < totalNodeNumber; idx++) {
        result.push_back(nodes[idx]->getAllVectorVariables());
    }
    return result;
}


vector<map<string,VectorVar>> GridManager::getVectorVariablesMapForAllNodes(){
    vector<map<string,VectorVar>> result;
    result.reserve(totalNodeNumber);
    for (int idx=0;idx<totalNodeNumber;idx++) {
        result.push_back(nodes[idx]->getVariables());
    }
    return result;
}

map<int, int> GridManager::getZonesBoundaryIndexes(){
    return zonesBoundaryIndexes;
}

map<int,int> GridManager::getZonesGlobalBoundaryIndexes(){
	return zonesGlobalBoundaryIndexes;
}

vector<int> GridManager::getZonesInBoundaryIndexes(){
	return zonesInBoundaryIndexes;
}

vector<int> GridManager::getZonesInBoundaryIndexesSpiral(){
	return zonesInBoundaryIndexesSpiral;
}

vector<vector<int>> GridManager::getIndicesOfNeighborNodesForAllZones(){
    vector<vector<int>> result;
    result.reserve(totalZoneNumber);
    for (int idx=0;idx<totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getIndicesOfNeighborNodes());
    }
    return result;
}

vector<vector<int>> GridManager::getIndicesOfNeighborZonesForAllNodes(){
    vector<vector<int>> result;
    result.reserve(totalNodeNumber);
    for (int idx=0;idx<totalNodeNumber;idx++) {
        result.push_back(nodes[idx]->getIndicesOfNeighborZones());
    }
    return result;
    
}

vector<int> GridManager::getNeighborhoodIndicesForZone(int idx){
    return zones[idx]->getNeighborhoodIndices();
}

vector<vector<int>> GridManager::getNeighborhoodZoneIndicesForAllZone(){
    vector<vector<int>> result;
    result.reserve(totalZoneNumber);
    for (int idx=0;idx<totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getNeighborhoodIndices());
    }
    return result;
    
}

int GridManager::getNodeIDXForZone(int idx){
    return zones[idx]->getNodeIDX();
}
void GridManager::setScalarVariableForZone(int idx, ScalarVar variable){
    zones[idx]->setScalarVariable(variable);
}

void GridManager::setScalarVariableForSubzone(int zoneIDX, int subzoneIDX, ScalarVar variable){
    zones[zoneIDX]->setSubzoneScalarVariable(subzoneIDX, variable);
}

void GridManager::setVectorVariableForSubzone(int zoneIDX, int subzoneIDX, VectorVar variable){
    zones[zoneIDX]->setSubzoneVectorVariable(subzoneIDX, variable);
}

void GridManager::setVectorVariableForZone(int idx, VectorVar variable){
    zones[idx]->setVectorVariable(variable);
}


void GridManager::setScalarVariableForNode(int idx, ScalarVar variable){
    nodes[idx]->setScalarVariable(variable);
}

void GridManager::setVectorVariableForNode(int idx, VectorVar variable){
    nodes[idx]->setVectorVariable(variable);
}



vector<vector<Subzone>> GridManager::getSubzonesForAllZones(){
    vector<vector<Subzone>> result;
    result.reserve(totalZoneNumber);
    for (int idx=0;idx<totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getSubzones());
    }
    
    return result;
    
}

vector<ScalarVar> GridManager::getScalarVariableForAllZones(string varName){
    vector<ScalarVar> result;
    result.reserve(totalZoneNumber);
    for (int idx=0; idx<totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getVariable(varName));
    }
    return result;
}

vector<vector<ScalarVar>> GridManager::getScalarVariablesForAllZones(){
    vector<vector<ScalarVar>> result;
    result.reserve(totalZoneNumber);
    for (int idx=0; idx < totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getAllScalarVariables());
    }
    return result;
}

vector<vector<VectorVar>> GridManager::getVectorVariablesForAllZones(){
    vector<vector<VectorVar>> result;
    result.reserve(totalZoneNumber);
    for (int idx=0; idx < totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getAllVectorVariables());
    }
    return result;
}

vector<VectorVar> GridManager::getVectorVariableForAllZones(string varName){
    vector<VectorVar> result;
    result.reserve(totalZoneNumber);
    for (int idx=0;idx<totalZoneNumber;idx++) {
        result.push_back(zones[idx]->getVectorVariable(varName));
    }
    return result;
}

vector<ScalarVar> GridManager::getScalarVariableForAllNodes(string varName){
    vector<ScalarVar> result;
    result.reserve(totalNodeNumber);
    for (int idx=0;idx<totalNodeNumber;idx++) {
        result.push_back(nodes[idx]->getScalarVariable(varName));
    }
    return result;
}

vector<VectorVar> GridManager::getVectorVariableForAllNodes(string varName){
    vector<VectorVar> result;
    result.reserve(totalNodeNumber);
    for (int idx=0;idx<totalNodeNumber;idx++) {
        result.push_back(nodes[idx]->getVariable(varName));
    }
    return result;
}
