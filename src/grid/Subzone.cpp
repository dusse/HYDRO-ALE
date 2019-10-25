#include "Subzone.hpp"

using namespace std;

Subzone::Subzone(){
}


vector<vector<double>> Subzone::getVerticesCoordinates() const{
    return this->verticesCoordinates;
}


void Subzone::setVerticesCoordinates(vector<vector<double>> verticesCoordinates){
    this->verticesCoordinates = verticesCoordinates;
}




ScalarVar Subzone::getVariable(string varName){
    if ( scalarVariables.find(varName) == scalarVariables.end() ) {
        return ScalarVar(varName, 0.0);
    } else {
        return scalarVariables.find(varName)->second;
    }
}

void Subzone::setVariable(ScalarVar var){
    this->scalarVariables[var.getName()] = var;
}


vector<ScalarVar> Subzone::getAllScalarVariables(){
    vector<ScalarVar> result;
    result.reserve(scalarVariables.size());
    map<string,ScalarVar>::iterator iterator;
    for( iterator = scalarVariables.begin(); iterator != scalarVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}



VectorVar Subzone::getVectorVariable(string varName){
    if ( vectorVariables.find(varName) == vectorVariables.end() ) {
        return VectorVar(varName, {0.0,0.0,0.0});
    } else {
        return vectorVariables.find(varName)->second;
    }
}

void Subzone::setVectorVariable(VectorVar var){
    this->vectorVariables[var.getName()] = var;
}

vector<VectorVar> Subzone::getAllVectorVariables(){
    vector<VectorVar> result;
    result.reserve(vectorVariables.size());
    map<string,VectorVar>::iterator iterator;
    for( iterator = vectorVariables.begin(); iterator != vectorVariables.end(); iterator++) {
        result.push_back(iterator->second);
    }
    return result;
}



