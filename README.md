MolSim GroupB
===

## Build Instructions
### Run with CLion
- Open the local CMake configurations
- Edit the MolSim CMake Application:
``` 
  Target: MolSim
  Executable: MolSim
  Program arguments: input/eingabe-sonne.txt
  Working directory: $ProjectFileDir$ (the project root containing input/ and src/)
```
When using VTK output, follow these steps:
- Open CMake settings in File>Settings>Build, Execution, Deployment>CMake
- Add -DENABLE_VTK_OUTPUT=ON to the CMake options field
- Now you can run the MolSim Application using CMake

### Build and run (CLI)
From project root execute the following commands in that order:
``` 
mkdir build
cd build
cmake -DBUILD_DOC=ON ..
make -j 4 # compilation with (4) parallel jobs
./MolSim ../input/eingabe-sonne.txt # any appropriate input file 
```
When using VTK output, run this instead:
```
mkdir build
cd build
cmake -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 .. # insert your VTK directory
make -j 4 # compilation with (4) parallel jobs
./MolSim ../input/eingabe-sonne.txt # any appropriate input file 
```

### Simulation
```
./MolSim ../input/eingabe-sonne.txt t_end delta_t # choose appropriate t_end, delta_t like 1000, 0.014 
```

You will find the generated output files under build/output

### Utility
In scripts/ you can find a clang-format-project.sh, used run clang format on the entire project and rebuild.sh, which can be used to recompile and build the project code cleanly.