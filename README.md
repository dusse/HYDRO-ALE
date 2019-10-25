# HYDRO-ALE
Arbitrary Lagrangian-Eulerian code for laser plasma interaction. 
C++11, MPI, HDF5.

the work is based on the pretty understandable thesis work of Milan Kucharik and the associated papers:
Kucharık, Milan. Arbitrary Lagrangian-Eulerian (ALE) methods in plasma physics. Diss. Ph. D. Thesis, 2006


1. before 'make' need to set 
      HDF5_PATH= path to hdf5 lib (last well used 1.10.5)
      MPI_PATH= path to mpi lib (last well used openmpi 9.0.0)
      PYTHON27_INC= path to python2.7 include
      PYTHON27_LIB= path to python2.7 lib

2. for running default example from src/input/Initializer.py
      mpirun -n 2 hydro_ale.exe
      
3. normally need to set input python file
      mpirun -n 2 hydro_ale.exe PATH/TO/PYTHON/INPUT/FILE
      
4. also before running need to create output folder and set in the python file.
      
5. For visualization use python notebook in folder NOTEBOOK



TODO:
1. now it is a 2D version
2. fix remapping for parallel version (see Remaper.cpp)
3. switch on Heat Manager (see LagrangianSolver.cpp)

