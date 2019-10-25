#ifndef Subzone_hpp
#define Subzone_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "Node.hpp"
#include "../common/variables/ScalarVar.hpp"
#include "../common/variables/VectorVar.hpp"

class Subzone
{
    
    /*
     *            - - - - - - - -x Zone(i,j)
     *           |               |
     *           |               .
     *           |               |
     *           |   Subzone 1   .
     *           |               |
     *           |               .
     *           |               |
     *           .Node(i-1,j-1)___
     *                   
     *                  
     *                  
     */
    //for 2D quadrilateral zone

    /* corner = link to Node
     * Subzone 1 -> Node(i-1,j-1)
     * Subzone 2 -> Node(i  ,j-1)
     * Subzone 3 -> Node(i  ,j  )
     * Subzone 4 -> Node(i-1,j  )
     */

    
private:
   
    std::map<std::string,ScalarVar> scalarVariables;
    std::map<std::string,VectorVar> vectorVariables;
    std::vector<std::vector<double>> verticesCoordinates;
    
public:
    Subzone();
    int getNumber() const;
    
    std::vector<std::vector<double>> getVerticesCoordinates() const;
    
    void setVerticesCoordinates(std::vector<std::vector<double>>);
    
    ScalarVar getVariable(std::string);
    void setVariable(ScalarVar);
    std::vector<ScalarVar> getAllScalarVariables();
    
    VectorVar getVectorVariable(std::string);
    void setVectorVariable(VectorVar);
    std::vector<VectorVar> getAllVectorVariables();
};
#endif /* Subzone_hpp */