#ifndef GridManager_hpp
#define GridManager_hpp
#include <stdio.h>
#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"
#include "../input/Loader.hpp"
#include "../common/variables/ScalarVar.hpp"
#include "../common/variables/VectorVar.hpp"

#include "Node.hpp"
#include "Zone.hpp"
#include "Subzone.hpp"


//Auxiliary variables:
//ZONE variables
const std::string  ZONAL_CENTROID = "centr_zone";
const std::string  ZONAL_VOLUME   = "vol_zone"  ;

//subzonal variables
const std::string  SUBZONAL_VOLUME   = "vol_subzone"  ;


/** VECTOR **/
const static std::string  COORDINATE ="coord";
const static std::string  COORDINATE_OLD ="coordold";

const static std::string  COORDINATE_SMOOTHED ="coord_smoothed";

const static int  GRID_OK   = 0;
const static int  GRID_FAIL = 1;


class GridManager
{
    
    
private:
    
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    
    
    int totalNodeNumber;
    int totalZoneNumber;
    std::vector<std::unique_ptr<Zone>> zones;
    std::vector<std::unique_ptr<Node>> nodes;
    std::map<int, int> zonesBoundaryIndexes;
    std::map<int, int> zonesGlobalBoundaryIndexes;
    std::vector<int> zonesInBoundaryIndexes;
    std::vector<int> sendStartInd_X;
    std::vector<int> sendEndInd_X;
    std::vector<int> recvStartInd_X;
    std::vector<int> recvEndInd_X;
    std::vector<int> sendStartInd_Y;
    std::vector<int> sendEndInd_Y;
    std::vector<int> recvStartInd_Y;
    std::vector<int> recvEndInd_Y;
    std::vector<int> sendStartInd_Z;
    std::vector<int> sendEndInd_Z;
    std::vector<int> recvStartInd_Z;
    std::vector<int> recvEndInd_Z;
    std::vector<int> zonesInBoundaryIndexesSpiral;
    std::vector<VectorVar> coordNodesAndGhost;
    
    void initialize();
    void initNodes();
    void initZones();
    int  calculateZoneCentroidsAndVolume();
    std::vector<double> getReflectedCoordinate( int , int );
    void initGridWithGhostLayer();
    void sendRecvIndecis4MPI();
    void sendBoundary2Neighbor();
    
public:
    
    GridManager(std::shared_ptr<Loader>);
    
    int getTotalNodeNumber();
    int getTotalZoneNumber();
    
    std::vector<ScalarVar> getScalarVariableForAllZones(std::string);
    std::vector<VectorVar> getVectorVariableForAllZones(std::string);
    std::vector<ScalarVar> getScalarVariableForAllNodes(std::string);
    std::vector<VectorVar> getVectorVariableForAllNodes(std::string);

    
    std::vector<std::vector<int>> getIndicesOfNeighborNodesForAllZones();
    std::vector<std::vector<int>> getIndicesOfNeighborZonesForAllNodes();
    
    std::vector<int> getNeighborhoodIndicesForZone(int);
    std::vector<std::vector<int>> getNeighborhoodZoneIndicesForAllZone();
    
    int getNodeIDXForZone(int);
    
    std::vector<std::vector<ScalarVar>> getScalarVariablesForAllNodes();
    std::vector<std::vector<VectorVar>> getVectorVariablesForAllNodes();
    std::vector<std::vector<ScalarVar>> getScalarVariablesForAllZones();
    std::vector<std::vector<VectorVar>> getVectorVariablesForAllZones();
    
    std::map<int,int> getZonesBoundaryIndexes();
    std::map<int,int> getZonesGlobalBoundaryIndexes();
    std::vector<int> getZonesInBoundaryIndexes();
    std::vector<int> getZonesInBoundaryIndexesSpiral();
    std::vector<std::map<std::string,VectorVar>> getVectorVariablesMapForAllNodes();
    
    std::vector<std::vector<Subzone>> getSubzonesForAllZones();
    void setSubzones(int, std::vector<Subzone>);
    void setScalarVariableForZone(int, ScalarVar);
    void setScalarVariableForSubzone(int, int, ScalarVar);
    void setVectorVariableForSubzone(int, int, VectorVar);
    void setVectorVariableForZone(int, VectorVar);
    void setVectorVariableForNode(int, VectorVar);
    void setScalarVariableForNode(int, ScalarVar);
    int  updateGrid();
    void smoothMesh();
    std::vector<VectorVar> getExtendedGridCoordinates();
    
};
#endif /* GridManager_hpp */
