#include "ScalarVar.hpp"
using namespace std;

ScalarVar::ScalarVar(){

}

ScalarVar::ScalarVar(std::string name, double value){
    this->name = name;
    this->value = value;
}

string ScalarVar::getName(){
    return name;
}

double ScalarVar::getValue(){    
    return value;
}

void ScalarVar::setValue(double newValue){
    this->value = newValue;
}