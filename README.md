# Fast_RNRR
This repository includes the source code of the paper "Quasi-Newton Solver for Robust Non-Rigid Registration", (CVPR2020), [https://arxiv.org/abs/2004.04322](https://arxiv.org/abs/2004.04322).

Authors: Yuxin Yao, [Bailin Deng](http://www.bdeng.me/), [Weiwei Xu](http://www.cad.zju.edu.cn/home/weiweixu/) and [Juyong Zhang](http://staff.ustc.edu.cn/~juyong/).

This code is protected under patent. It can be only used for research purposes. If you are interested in business purposes/for-profit use, please contact Juyong Zhang (the corresponding author, email: juyong@ustc.edu.cn).

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
1. This code supports non-rigid registration from mesh to mesh or point cloud. If the method of geodesics computation supports point clouds, our method is also applicable when template format is point cloud.
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
We adopt the [SHOT implementation](https://github.com/fedassa/SHOT) to compute the shot feature. To construct the deformation graph, we adopt the method in the source code of http://www.liuyebin.com/nonrigid.html.

This work was supported by the National Natural Science Foundation of China (No. 61672481), and Youth Innovation Promotion Association CAS (No. 2018495).
