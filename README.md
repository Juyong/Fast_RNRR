# Fast_FNRR
This repository includes the source code the paper "Quasi-Newton Solver for Robust Non-Rigid Registration" (CVPR2020).

Authors: Yuxin Yao, [Bailin Deng](http://www.bdeng.me/), [Weiwei Xu](http://www.cad.zju.edu.cn/home/weiweixu/) and [Juyong Zhang](http://staff.ustc.edu.cn/~juyong/) .

# Dependencies
1. [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)
2. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
3. [OpenCV](https://opencv.org/)
4. [libigl](https://github.com/libigl/libigl)

# Compilation
The code is compiled using [CMake](https://cmake.org/) and tested on Ubuntu 16.04 (gcc5.4.0) and on Windows with Visual Studio 2015. An executable `NonRigidreg` will be generated.

# Usage
The program is run with five input parameters:
```
$ NonRigidreg <srcFile> <tarFile> <outFile> <regMethod> <landmarkFile>
```
1.`<srcFile>`: an input file storing the source mesh;

2.`<tarFile>`: an input file storing the target mesh or point cloud; 

3.`<outFile>`: an output file storing the registered source mesh; 

4.`<regMethod>`: registration method (0 for N-ICP, 1 for RPTS, 2 for SVR-l0, 3 for Our robust non-rigid registration);

5.`<landmarkFile>`: an landmark file (2xn matrix, first column includes the indexes in source file, second column includes the indexes in target file, each row is a pair correspondences separated by space).

`<regMethod>` and `<landmarkFile>` can be ignored, our robust non-rigid registration method without landmarks will be used in this case.

## Notes
1. This code supports non-rigid registration from mesh to mesh or point cloud. If the calculation method of geodesics are changed to support point clouds, our method is also applicable when template is a point cloud.
2. This code contains the calculation of SHOT features and the diffusion pruning method, these methods can be used by setting `paras.corres_type = SHOT` and `paras.pruning_type = DIFFUSION` in `main.cpp`.

# Citation
Please cite the following papers if it helps your research:
```
@inproceedings{Yao2019Quasi,
      title={Quasi-Newton Solver for Robust Non-Rigid Registration},
      author={Yao, Yuxin and Deng, Bailin and Xu, Weiwei and Zhang, Juyong},
      booktitle={IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
      year={2020}
}
```

# Acknowledgement
The [SHOT implementation](https://github.com/fedassa/SHOT) provided by the author was used in the calculation part of the shot feature 
and the construction part of deformation graph in the source code of http://www.liuyebin.com/nonrigid.html was adopted.

# License
This code is protected under patent. It is for research purposes only at your university (research institution) only. If you are interested in business purposes/for-profit use, please contact Prof.Zhang(the corresponding author, email: juyong@ustc.edu.cn).