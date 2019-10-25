# need to set
# HDF5_PATH = ...
# MPI_PATH = ...
# PYTHON27_INC= ...
# PYTHON27_LIB= ...

LIBS=-lpython2.7 -lhdf5
INCLUDES=-I$(PYTHON27_INC) -I$(HDF5_PATH)/include/ -I$(MPI_PATH)/include/

DSRC = ./src
DEXE = ./

LD_LIBRARY_PATH=$(HDF5_PATH)/lib/:$(PYTHON27_LIB)


export LIBRARY_PATH=$LIBRARY_PATH:$(LD_LIBRARY_PATH)

CXX = $(MPI_PATH)/bin/mpicxx
CXXFLAGS  = -g3 -Wall -c -std=c++11 -Wno-sign-compare -Wno-unused-variable

_SRCS =  $(DSRC)/core/SimulationManager.cpp \
               $(DSRC)/grid/GridManager.cpp \
               $(DSRC)/grid/Node.cpp \
               $(DSRC)/grid/Subzone.cpp \
               $(DSRC)/grid/Zone.cpp \
               $(DSRC)/input/Loader.cpp \
               $(DSRC)/misc/Logger.cpp \
               $(DSRC)/misc/Misc.cpp \
               $(DSRC)/output/Writer.cpp \
               $(DSRC)/physics/heat/HeatManager.cpp \
               $(DSRC)/physics/laser/LaserManager.cpp \
               $(DSRC)/common/variables/ScalarVar.cpp \
               $(DSRC)/common/variables/VectorVar.cpp \
               $(DSRC)/solvers/LagrangianSolver.cpp \
               $(DSRC)/solvers/Remaper.cpp \
               $(DSRC)/HYDRO-ALE.cpp \

_OBJS            = $(_SRCS:.cpp=.o)

_EXEN            = $(DEXE)/hydro_ale.exe

all : $(_EXEN)


$(_EXEN) : $(_OBJS)
	@echo 'Building target: $@'
	$(CXX) -o $@ $^  $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o : %.cpp
	$(CXX) $(INCLUDES) -o $@ $< $(CXXFLAGS) 


clean :
	rm -f $(_OBJS)


