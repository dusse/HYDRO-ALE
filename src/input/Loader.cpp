#include "Loader.hpp"

using namespace std;

Loader::Loader(){
    logger.writeMsg("[Loader] Loader start!\n");
}

const char  *INIT_CLASS_NAME = "Initializer";
const string  BRACKETS ="()";
const string  BRACKETS_3DOUBLE = "(ddd)";
const string  GET_TIMESTEP = "getTimestep";
const string  GET_MAX_TIMESTEPS_NUM = "getMaxTimestepsNum";
const string  GET_TIMESTEP_WRITE = "getEachTimestepsNum2WriteFile";
const string  GET_OUTPUT_DIR = "getOutputDir";
const string  GET_FILENAME_TEMPLATE = "getOutputFilenameTemplate";
const string  GET_DENSITY = "getDensity";
const string  GET_VELOCITY = "getVelocity";
const string  GET_VELOCITY_X = "getVelocityX";
const string  GET_MPI_DOMAIN_NUM = "mpiDomainNum";
const string  dirs[] = {"X", "Y", "Z"};

void Loader::load()
{
        MPI_Comm com;
        int periods[3];
        periods[0] = 0;
        periods[1] = 0;
        periods[2] = 0;
    
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	int cs[3];
    	int ss[3];


    	string msg = "[Loader] initialization...for proc = " + to_string(rank);
    	logger.writeMsg(msg.c_str());
    	Py_Initialize();
    	this->pInstance = getPythonClassInstance(INIT_CLASS_NAME);


    	for (int n=0; n<3; n++){
    	     string domainNum ="get"+dirs[n]+GET_MPI_DOMAIN_NUM;
    	      this->mpiDomains[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance, strdup(domainNum.c_str()), strdup(BRACKETS.c_str())));
    	 }
    	MPI_Cart_create(MPI_COMM_WORLD, 3, mpiDomains, periods, 1, &com);
    	MPI_Cart_coords(com, rank, 3, cs);
    	for(int i = 0; i < 3; i++){
    		this->mpiCoords[i] = cs[i];
    	}
    	int wr;

    	for(int i = 0; i < 27; i++){
    		this->neighbors2Send.push_back(MPI_PROC_NULL);
    		this->neighbors2Recv.push_back(MPI_PROC_NULL);
    	}
    	for (int a = -1; a <= +1; a++){
    		for (int b = -1; b <= +1; b++){
    			for (int c = -1; c <= +1; c++){
    	                /* __ guess the coordinate of the neighbor __ */
    	                ss[0] = cs[0]+a;
    	                ss[1] = cs[1]+b;
    	                ss[2] = cs[2]+c;

    	                /* __ if coordinate out of range : no neighbor __ */
    	                if ((ss[0] >= 0 && ss[0] < mpiDomains[0]) && (ss[1] >= 0 && ss[1] < mpiDomains[1]) && (ss[2] >= 0 && ss[2] < mpiDomains[2]))
    	                {
    	                	MPI_Cart_rank(com, ss, &wr);

    	                    //-1 -> left, bottom, back
    	                    // 0 -> same position
    	                    //+1 -> right, top, front
    	                	this->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] = wr;
    	                	this->neighbors2Recv[(1-c)+3*((1-b)+3*(1-a))] = wr;
    	                }
    	            }
    	        }
    	    }


    	double boxSizePerDomain[3];
    	for (int n=0; n<3; n++){
        
    		string res ="get"+dirs[n]+"resolution";
    		string left_str ="get"+dirs[n]+"left";
    		string right_str="get"+dirs[n]+"right";
    		this->boxCoordinates[n][0] = PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(left_str.c_str()),strdup(BRACKETS.c_str())));
    		this->boxCoordinates[n][1] = PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(right_str.c_str()),strdup(BRACKETS.c_str())));
        	this->boxSizes[n]= boxCoordinates[n][1] - boxCoordinates[n][0];

        	boxSizePerDomain[n] = boxSizes[n]/mpiDomains[n];
        	this->boxCoordinates[n][0] = boxCoordinates[n][0] + cs[n]*boxSizePerDomain[n];
        	this->boxCoordinates[n][1] = boxCoordinates[n][0] + boxSizePerDomain[n];

        	this->resolution[n] = (int) PyInt_AsLong(PyObject_CallMethod(pInstance, strdup(res.c_str()),strdup(BRACKETS.c_str())));
        	if (resolution[n] == 1){
        		this->spatialSteps[n]=boxSizes[n];
        	}else{
        		this->spatialSteps[n]=boxSizePerDomain[n]/(double)(resolution[n]-1);
        	}
        
    }
    this->dim = SIMULATION_SIZE;
    this->subzoneNum = pow(2, dim);
    this->timeStep = PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(GET_TIMESTEP.c_str()), strdup(BRACKETS.c_str())));
    this->maxTimestepsNum = PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(GET_MAX_TIMESTEPS_NUM.c_str()), strdup(BRACKETS.c_str())));
    this->timestepsNum2Write = PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(GET_TIMESTEP_WRITE.c_str()), strdup(BRACKETS.c_str())));
    
    this->outputDir=PyString_AS_STRING(PyObject_CallMethod(pInstance, strdup(GET_OUTPUT_DIR.c_str()), strdup(BRACKETS.c_str())));
    
    this->fileNameTemplate=PyString_AS_STRING(PyObject_CallMethod(pInstance, strdup(GET_FILENAME_TEMPLATE.c_str()), strdup(BRACKETS.c_str())));
    
    if(rank == 0){


    	printf("[Loader]  Box size: \n");

    	for (int n=0; n<3; n++){
    		printf("[Loader]  [%s] [%1.5f, %1.5f] l = %1.5f res = %3i => step  = %1.6f\n", dirs[n].c_str(), boxCoordinates[n][0],boxCoordinates[n][1],boxSizes[n],resolution[n],spatialSteps[n]);
    	}
        msg = "[Loader] timeStep = "+to_string(timeStep);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] subzoneNum = "+to_string(subzoneNum);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] Max timeStep number = "+to_string(maxTimestepsNum);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  timeStep number to write to file = "+to_string(timestepsNum2Write);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  outputDir = "+outputDir;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  filename template = "+fileNameTemplate;
        logger.writeMsg(msg.c_str(), INFO);
        
    	logger.writeMsg("[Loader] initialization...OK!");
    	}

}


void Loader::initDensity(double* resultDensity){
     init3Darr(resultDensity, strdup(GET_DENSITY.c_str()));
}

void Loader::initVelocityX(double* resultVelocityX){
    init3Darr(resultVelocityX, strdup(GET_VELOCITY_X.c_str()));
}

double Loader::getTimeStep(){
    return timeStep;
}

int Loader::getMaxTimestepsNum(){
    return maxTimestepsNum;
}

int Loader::getTimestepsNum2Write(){
    return timestepsNum2Write;
}

string Loader::getOutputDir() const {
    return outputDir;
}

string Loader::getFilenameTemplate(){
    return fileNameTemplate;
}

void Loader::init3Darr(double* resultArr,char *varName){
    PyObject *pValue;
    int i,j,k;
    double x,y,z;
    for (i = 0 ; i < resolution[0] ; i++)
    {
        for (j = 0 ; j < resolution[1] ; j++)
        {
            for (k = 0 ; k < resolution[2] ; k++)
            {
                x = i*spatialSteps[0] + boxCoordinates[0][0];
                y = j*spatialSteps[1] + boxCoordinates[1][0];
                z = k*spatialSteps[1] + boxCoordinates[2][0];
                pValue = PyObject_CallMethod(pInstance, varName, strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
                resultArr[IDX(i, j, k, resolution[0], resolution[1], resolution[2])] = PyFloat_AsDouble(pValue);
            }
        }
    }
}

vector<double> Loader::getVelocity(double x,double y,double z){
    PyObject *pValue;
    vector<double> velocity;
    
    for (int n=0; n<3; n++){
        string varName =GET_VELOCITY+dirs[n];

        pValue = PyObject_CallMethod(pInstance, strdup(varName.c_str()), strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
        velocity.push_back(PyFloat_AsDouble(pValue));
    }
    return velocity;
}

double Loader::getDensity(double x,double y,double z){
    return PyFloat_AsDouble(PyObject_CallMethod(pInstance, strdup(GET_DENSITY.c_str()), strdup(BRACKETS_3DOUBLE.c_str()),x,y,z));
}

PyObject * Loader::getPythonClassInstance(string className){
    PyObject  *pName, *pModule, *pDict, *pClass, *pInstance;
    string msg = "[Loader] Start to instantiate the Python class " + className;
    logger.writeMsg(msg.c_str());
    pName = PyString_FromString(className.c_str());
    pModule = PyImport_Import(pName);
    pDict = PyModule_GetDict(pModule);
    pClass = PyDict_GetItemString(pDict, className.c_str());
    if (PyCallable_Check(pClass))
    {
        pInstance = PyObject_CallObject(pClass, NULL);
    }
    else
    {
        logger.writeMsg("[Loader] Cannot instantiate the Python class");
        pInstance = nullptr;
    }
    logger.writeMsg("[Loader] finish to instantiate the Python class ");
    
    return pInstance;
}

Loader::~Loader()
{
    Py_Finalize();
    logger.writeMsg("[Loader] FINALIZE...OK!");
}
