# Fast_RNRR
This repository includes the source code the paper "Quasi-Newton Solver for Robust Non-Rigid Registration" (CVPR2020).

## Dependencies
1. OpenCV 
2. Eigen
3. OpenMesh
4. libigl

## Compilation
The code is compiled using [CMake](https://cmake.org/) and tested on Ubuntu 16.04 (gcc5.4.0) and on Windows with Visual Studio 2015. An executable `NonRigidreg` will be generated.

## Usage
The program is run with five input parameters:
1.　an input file storing the source mesh;
2.　an input file storing the target mesh or point cloud; 
3.　an output file storing the registered source mesh; 
4.　registration method (0 for N-ICP, 1 for RPTS, 2 for SVR-l0, 3 for Our robust non-rigid registration);
5.　an landmark file (two columns, first column includes the indexes in source file, second column includes the indexes in target file, each row is a pair correspondences separated by space).

The last two parameters can be ignored, our robust non-rigid registration method without landmarks will be used in this case.

Example:
```
$ NonRigidreg src.obj tar.obj source_reg_target.obj 3 landmark.txt
```
