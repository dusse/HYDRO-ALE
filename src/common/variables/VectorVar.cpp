#include "VectorVar.hpp"

using namespace std;

VectorVar::VectorVar(){
}

VectorVar::VectorVar(std::string name, std::vector<double> value){
    this->name = name;
    this->value = value;
}

string VectorVar::getName(){
    return name;
}

std::vector<double> VectorVar::getValue(){
    return value;
}

void VectorVar::setValue(std::vector<double> newValue){
    this->value = newValue;
}