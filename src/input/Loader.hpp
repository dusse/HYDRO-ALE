#ifndef Loader_hpp
#define Loader_hpp

#include <Python.h>
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"

class Loader
{
private:
    Logger logger;
    double timeStep;
    int maxTimestepsNum;
    int timestepsNum2Write;
    std::string outputDir;
    std::string fileNameTemplate;
    
    PyObject *pInstance;
    
public:
    int resolution[3];
    int mpiDomains[3];
    int mpiCoords[3];
    const int SIMULATION_SIZE = 2;
    int dim;
    int subzoneNum;
    double boxCoordinates[3][2];
    double boxSizes[3];
    double spatialSteps[3];
    std::vector<int> neighbors2Send;//27
    std::vector<int> neighbors2Recv;//27
    Loader();
    ~Loader();
    void load();
    void init3Darr(double* ,char *);
    void initDensity(double*);
    void initVelocityX(double*);
    std::vector<double> getVelocity(double,double,double);
    double getDensity(double,double,double);
    double getTimeStep();
    int getMaxTimestepsNum();
    int getTimestepsNum2Write();
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    PyObject * getPythonClassInstance(std::string className);
};
#endif /* Loader_hpp */
