#ifndef Node_hpp
#define Node_hpp


#include "../common/variables/VectorVar.hpp"
#include "../common/variables/ScalarVar.hpp"


#include <stdio.h>
#include <iostream>
#include <string>
#include <map>

class Node
{
    /*
     *                           x Zone(i,j+1) - . Node(i,j+1) - - x Zone(i+1,j+1)
     *                           .               |               .
     *                           |               |               |
     *                           .               |               .
     *                           |         a1-<--|-->a2+         |
     *                           .  Subzone 2    |   Subzone 1   .
     *                           |               |               |
     *                           .       a3+     |       a4-     .
     *                           |      /|\      |       /|\     |
     *            .Node(i-1,j)___________|_______.Node(i,j)____________________.Node(i+1,j)
     *                           |       |       |        |      |
     *                           .      \|/      |       \|/     .
     *                           |      a2-      |       a1+     |
     *                           .    Subzone 3  |   Subzone 4   .
     *                      s3<--|         a4+<--|-->a3-         |
     *                           .      s2       |               .
     *                           |      /|\      |               |
     *                           .       |       |               .
     *                           x Zone(i,j) - - . Node(i,j-1) - x Zone(i+1,j)
     */
private:
    std::map<std::string,VectorVar> vectorVariables;
    std::map<std::string,ScalarVar> scalarVariables;
    
    std::vector<int> indicesOfNeighborZones;
    
public:
    Node();
    Node(std::map<std::string,VectorVar> variables);
    VectorVar getVariable(std::string);
    ScalarVar getScalarVariable(std::string);
    std::vector<VectorVar> getAllVectorVariables();
    std::vector<ScalarVar> getAllScalarVariables();
    std::map<std::string,VectorVar> getVariables();
    void setVectorVariable(VectorVar);
    void setScalarVariable(ScalarVar);
    std::vector<int> getIndicesOfNeighborZones();
    void setIndicesOfNeighborZones(std::vector<int>);
    
};
#endif /* Node_hpp */