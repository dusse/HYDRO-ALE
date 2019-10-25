#ifndef ScalarVar_hpp
#define ScalarVar_hpp

#include <string>
#include <stdio.h>
#include <iostream>

class ScalarVar{

private:
    std::string name;
    double value;
    
public:
    ScalarVar();
    ScalarVar(std::string, double);
    std::string getName();
    double getValue();
    void setValue(double);
};
#endif /* ScalarVar_hpp */
