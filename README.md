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
cmake CMakeLists.txt
make
./MolSim input/eingabe-sonne.txt # (you can use any appropriate input file here) 
```