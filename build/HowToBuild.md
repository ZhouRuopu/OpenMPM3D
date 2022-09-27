# Step 1: Install CMake
- Download CMake installer from https://cmake.org/download/ according to your own system.
- Install CMake as default.

# Step 2: Install VTK (Visualization Toolkit)
- Download the source codes of VTK from https://vtk.org/download/
- Configure and Generate the project of VTK by CMake
- Compile the project 

# Step 3: Build the project of MPM3D-CPP
- Configure and Generate the project of MPM3D-CPP by CMake. (You may need to add the path of VTK project manually)
- Compile the project of MPM3D-CPP

# Step 4: Run the simulation
- usage
```
$ [ExeFilePath]/MPM3D.exe [InputFilePath]/inputfile.xml
```
- Import the results into Paraview for visualization