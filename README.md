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
cmake --build .
./MolSim ../input/eingabe-sonne.txt # any appropriate input file 
```
When using VTK output, run this instead:
```
mkdir build
cd build
cmake -DBUILD_DOC=ON -DENABLE_VTK_OUTPUT=ON -DVTK_DIR=/usr/local/vtk/lib/cmake/vtk-9.5 .. # insert your VTK directory
cmake --build .
./MolSim ../input/eingabe-collision.txt 5 0.0002 file info P:OFF
```

### Simulation
```
     "./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "P:ON]"
```

You will find the generated output files under build/output

### Utility
In scripts/ you can find a clang-format-project.sh, used run clang format on the entire project and rebuild.sh, which can be used to recompile and build the project code cleanly.
### Benchmarks
You can find optimized benchmark scripts in the scripts/ directory