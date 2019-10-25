#include "Writer.hpp"

using namespace std;

const int   SIMULATION_SIZE = 2;

Writer::Writer(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr, shared_ptr<LaserManager> laserMnr, shared_ptr<HeatManager> heatMnr):
				loader(move(load)), gridMgr(move(gridMnr)),laserMgr(move(laserMnr)),heatMgr(move(heatMnr))
{
    logger.reset(new Logger());

    this->outputDir = loader->getOutputDir();
    this->fileNamePattern = loader->getFilenameTemplate();
    
    logger->writeMsg("[Writer] initialize...OK");
}

void Writer::write(int fileNum){
    logger->writeMsg("[Writer] writing...");
	int xRes = loader->resolution[0], yRes = loader->resolution[1];
    int totalNodeNum = gridMgr->getTotalNodeNumber();
    int totalZoneNum = gridMgr->getTotalZoneNumber();
    int ijZone, ijNode;
	vector<vector<VectorVar>> vectorVars = gridMgr->getVectorVariablesForAllNodes();
	vector<vector<ScalarVar>> scalarVars = gridMgr->getScalarVariablesForAllZones();
	vector<vector<VectorVar>> vectorVarsZone = gridMgr->getVectorVariablesForAllZones();
    vector<ScalarVar> temperature = heatMgr->getCurrentTemperature();
    vector<VectorVar> fluxQ = heatMgr->getCurrentFlux();
    vector<ScalarVar> divQ = heatMgr->getFluxDivergence();
    vector<ScalarVar> divI = laserMgr->getIntensityDivergence();
    vector<VectorVar> coordsExt = gridMgr->getExtendedGridCoordinates();
    
	MPI_Info info = MPI_INFO_NULL;
	hid_t access = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
	string fileName = outputDir + fileNamePattern + std::to_string(fileNum) + ".h5";
	hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
	const string groupname = "/vars";
	hid_t group = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

	for (int idx_var = 0; idx_var < vectorVars[0].size(); idx_var++) {
	    for (int dir = 0; dir < vectorVars[0][idx_var].getValue().size(); dir++) {
	        string varName = vectorVars[0][idx_var].getName()+"_"+to_string(dir);
	        double* var = new double[totalNodeNum];
	        for (ijNode = 0; ijNode < totalNodeNum; ijNode++) {
	            var[ijNode] = (vectorVars[ijNode][idx_var].getValue())[dir];
	        }
	        int size[3] = {xRes, yRes, 1};
	        writeParallel(fileID, group, dxpl_id, varName, var, size);
	        delete [] var;
	    }
	}

	for (int idx_var=0;idx_var<vectorVarsZone[0].size();idx_var++) {
	    for (int dir = 0; dir < vectorVarsZone[0][idx_var].getValue().size(); dir++) {
	        string varName = vectorVarsZone[0][idx_var].getName()+"_"+to_string(dir);
	        double* var = new double[totalZoneNum];
	        for (ijZone = 0; ijZone < totalZoneNum; ijZone++) {
	            var[ijZone] = (vectorVarsZone[ijZone][idx_var].getValue())[dir];
	        }
	        int size[3] = {xRes+1, yRes+1, 1};
	        writeParallel(fileID, group, dxpl_id, varName, var, size);
	        delete [] var;
	    }
	}

	for (int idx_var=0;idx_var<scalarVars[0].size();idx_var++) {
	    string varName = scalarVars[0][idx_var].getName();
	    double* var = new double[totalZoneNum];
	    for (ijZone = 0; ijZone < totalZoneNum; ijZone++) {
	        var[ijZone] = scalarVars[ijZone][idx_var].getValue();
	    }
	    int size[3] = {xRes+1, yRes+1, 1};
	    writeParallel(fileID, group, dxpl_id, varName, var, size);
	    delete [] var;
	}

	for (int dir = 0; dir < SIMULATION_SIZE; dir++) {
	    string varName = "crdext_"+to_string(dir);
	    double* var = new double[(xRes+2)*(yRes+2)];
	    for (int idx = 0; idx < (xRes+2)*(yRes+2); idx++) {
	        var[idx] = coordsExt[idx].getValue()[dir];
	    }
	    int size[3] = {xRes+2, yRes+2, 1};
        writeParallel(fileID, group, dxpl_id, varName, var, size);
	    delete [] var;
	}
    int size[3] = {xRes+1, yRes+1, 1};
    string varName = "temp";
    double* varT = new double[totalZoneNum];
    for (int idx = 0; idx < totalZoneNum; idx++) {
        varT[idx] = temperature[idx].getValue();
    }
    
    writeParallel(fileID, group, dxpl_id, varName, varT, size);
    delete [] varT;
    
    varName = "divq";
    double* varDivQ = new double[totalZoneNum];
    for (int idx = 0; idx < totalZoneNum; idx++) {
        varDivQ[idx] = divQ[idx].getValue();
    }
    writeParallel(fileID, group, dxpl_id, varName, varDivQ, size);
    delete [] varDivQ;
    
    varName = "divi";
    double* varDivI = new double[totalZoneNum];
    for (int idx = 0; idx < totalZoneNum; idx++) {
        varDivI[idx] = divI[idx].getValue();
    }
    writeParallel(fileID, group, dxpl_id, varName, varDivI, size);
    delete [] varDivI;
    
    
    for (int dir = 0; dir < SIMULATION_SIZE; dir++) {
        string varName = "flux_"+to_string(dir);
        double* var = new double[totalZoneNum];
        for (int idx = 0; idx < totalZoneNum; idx++) {
            var[idx] = fluxQ[idx].getValue()[dir];
        }
        writeParallel(fileID, group, dxpl_id, varName, var, size);
        delete [] var;
    }
    
    H5Gclose(group);
    H5Pclose(dxpl_id);
    H5Fclose(fileID);
    logger->writeMsg("[Writer] writing...OK");
}



inline bool doesFileExist(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

void Writer::writeParallel(hid_t fileID, hid_t group, hid_t dxpl_id, string dsetName , const double* data, int sizes[3]){

	int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	hsize_t offset[3];
	hsize_t nLoc[3];
	hsize_t nGlob[3];

	for (int dir = 0; dir < 3; dir++) {
		offset[dir] = sizes[dir]*loader->mpiCoords[dir]*(loader->mpiDomains[dir]-1);
		nLoc[dir] = sizes[dir];
		nGlob[dir] = sizes[dir]*loader->mpiDomains[dir];
	}

    hid_t memspace = H5Screate_simple(3, nLoc, NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);

    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);

    H5Sclose(memspace);
    H5Sclose(filespace);
}


