#ifndef VectorVar_hpp
#define VectorVar_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

class VectorVar{

private:
    std::string name;
    std::vector<double> value;
    
public:
    VectorVar();
    VectorVar(std::string, std::vector<double>);
    std::string getName();
    std::vector<double> getValue();
    void setValue(std::vector<double>);
};
#endif /* VectorVar_hpp */

