#ifndef Zone_hpp
#define Zone_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "Subzone.hpp"
#include "../common/variables/ScalarVar.hpp"
#include "../common/variables/VectorVar.hpp"

class Zone
{
    /*
     *                    a4-           a3+
     *                   /|\            /|\
     *            .Node(i-1,j)___________|_______.Node(i,j)
     *           |               |       |       |
     *           |               .      \|/      |
     *           |               |      a2-      |
     *     a4+<--|   Subzone 4   .    Subzone 3  |
     *           |          s3<--|         a4+<--|-->a3-
     *           |               .      s2       |
     *           |               |      /|\      |
     *           |               .       |       |
     *           |- - - - - - - -x Zone(i,j) - - |
     *           |       |       |               |
     *           |      \|/      .               |
     *           |       s4      |               |
     *           |               .-->s1          |
     *     a1-<--|   Subzone 1   |   Subzone 2   |-->a2+
     *           |               .               |
     *           |               |               |
     *           .Node(i-1,j-1)__________________.Node(i,j-1)
     *                   |                |
     *                  \|/              \|/
     *                  a1+               a2-
     */
    

    
    
private:
    
    int nodeIDX;
    
    std::vector<int> indicesOfNeighborNodes;
    
    std::vector<int> neighborhoodIndices;// 8 for 2D
    
    std::vector<Subzone> subzones;
    
    std::map<std::string,ScalarVar> scalarVariables;
    
    std::map<std::string,VectorVar> vectorVariables;
    

public:
    
    Zone();
    Zone(std::map<std::string,ScalarVar>, std::vector<Subzone>, std::vector<int>,  std::vector<int>);
    ~Zone();
    std::vector<Subzone> getSubzones();
    ScalarVar getVariable(std::string);
    VectorVar getVectorVariable(std::string);
    std::vector<ScalarVar> getAllScalarVariables();
    std::vector<VectorVar> getAllVectorVariables();
    std::vector<int> getIndicesOfNeighborNodes();
    std::vector<int> getNeighborhoodIndices();
    void setIndicesOfNeighborNodes(std::vector<int>);
    void setNeighborhoodIndices(std::vector<int>);
    void setScalarVariable(ScalarVar);
    void setVectorVariable(VectorVar);
    void setSubzones(std::vector<Subzone> subzones);
    void setSubzoneScalarVariable(int, ScalarVar);
    void setSubzoneVectorVariable(int, VectorVar);
    void setSubzoneVertices(int, std::vector<std::vector<double>>);
    void setNodeIDX(int);
    int getNodeIDX();
};
#endif /* Zone_hpp */