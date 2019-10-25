#include "Node.hpp"

using namespace std;

Node::Node(){
}

Node::Node(map<string,VectorVar> vectorVariables){
    this->vectorVariables = vectorVariables;
}

VectorVar Node::getVariable(string varName){
    if ( vectorVariables.size() == 0 || vectorVariables.find(varName) == vectorVariables.end() ) {
        return VectorVar(varName, {0.0,0.0,0.0});
    } else {
        return vectorVariables.find(varName)->second;
    }
}

ScalarVar Node::getScalarVariable(string varName){
    if ( scalarVariables.find(varName) == scalarVariables.end() ) {
        return ScalarVar(varName, 0.0);
    } else {
        return scalarVariables.find(varName)->second;
    }
}

vector<VectorVar> Node::getAllVectorVariables(){
    vector<VectorVar> result;
    result.reserve(vectorVariables.size());
    map<string,VectorVar>::iterator iterator;
    for( iterator = vectorVariables.begin(); iterator != vectorVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}

vector<ScalarVar> Node::getAllScalarVariables(){
    vector<ScalarVar> result;
    result.reserve(scalarVariables.size());
    map<string,ScalarVar>::iterator iterator;
    for( iterator = scalarVariables.begin(); iterator != scalarVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}

map<string,VectorVar> Node::getVariables(){
    return vectorVariables;
}

void Node::setVectorVariable(VectorVar variable){
    this->vectorVariables[variable.getName()] = variable;
}
    
void Node::setScalarVariable(ScalarVar variable){
    this->scalarVariables[variable.getName()] = variable;
}

vector<int> Node::getIndicesOfNeighborZones(){
    return indicesOfNeighborZones;
}

void Node::setIndicesOfNeighborZones(vector<int> indicesOfNeighborZones){
    this->indicesOfNeighborZones = indicesOfNeighborZones;
}
