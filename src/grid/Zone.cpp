#include "Zone.hpp"

using namespace std;

Zone::Zone(std::map<std::string,ScalarVar> scalarVariables, vector<Subzone> subzones, vector<int> indicesOfNeighborNodes, vector<int> neighborhoodIndices){
    this->scalarVariables = scalarVariables;
    this->subzones = subzones;
    this->indicesOfNeighborNodes = indicesOfNeighborNodes;
    this->neighborhoodIndices = neighborhoodIndices;
}

Zone::Zone(){
}

Zone::~Zone(){
}


vector<Subzone> Zone::getSubzones(){
    return this->subzones;
}

vector<int> Zone::getIndicesOfNeighborNodes(){
    return indicesOfNeighborNodes;
}

vector<int> Zone::getNeighborhoodIndices(){
    return neighborhoodIndices;
}

void Zone::setIndicesOfNeighborNodes(vector<int> indicesOfNeighborNodes){
    this->indicesOfNeighborNodes = indicesOfNeighborNodes;
}

void Zone::setNeighborhoodIndices(vector<int> neighborhoodIndices){
 this->neighborhoodIndices = neighborhoodIndices;
}

ScalarVar Zone::getVariable(string varName){
    if ( scalarVariables.find(varName) == scalarVariables.end() ) {
        return ScalarVar(varName, 0.0);
    } else {
        return scalarVariables.find(varName)->second;
    }
}

VectorVar Zone::getVectorVariable(string varName){
    if ( vectorVariables.find(varName) == vectorVariables.end() ) {
        return VectorVar(varName, {0.0,0.0,0.0});
    } else {
        return vectorVariables.find(varName)->second;
    }
}

vector<ScalarVar> Zone::getAllScalarVariables(){
    std::vector<ScalarVar> result;
    result.reserve(scalarVariables.size());
    map<string,ScalarVar>::iterator iterator;
    for( iterator = scalarVariables.begin(); iterator != scalarVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}

vector<VectorVar> Zone::getAllVectorVariables(){
    std::vector<VectorVar> result;
    result.reserve(vectorVariables.size());
    map<string,VectorVar>::iterator iterator;
    for( iterator = vectorVariables.begin(); iterator != vectorVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}

void Zone::setScalarVariable(ScalarVar variable){
    this->scalarVariables[variable.getName()] = variable;
}

void Zone::setVectorVariable(VectorVar variable){
    this->vectorVariables[variable.getName()] = variable;
}

void Zone::setSubzoneScalarVariable(int num, ScalarVar var){
    subzones[num].setVariable(var);
}

void Zone::setSubzoneVectorVariable(int num, VectorVar var){
    subzones[num].setVectorVariable(var);
}

void Zone::setSubzoneVertices(int num, vector<vector<double> > vertices){
    subzones[num].setVerticesCoordinates(vertices );
}

void Zone::setSubzones(vector<Subzone> subzones){
    this->subzones = subzones;
}

void Zone::setNodeIDX(int nodeIDX){
    this->nodeIDX = nodeIDX;
}

int Zone::getNodeIDX(){
    return nodeIDX;
}

