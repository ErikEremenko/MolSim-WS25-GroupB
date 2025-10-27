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
- Now you can run the MolSim Application using CMake

### Build and run (CLI)
From project root execute the following commands in that order:
``` 
mkdir build
cd build
cmake ..
make -j 4 # compilation with (4) parallel jobs
./MolSim ../input/eingabe-sonne.txt # any appropriate input file 
```