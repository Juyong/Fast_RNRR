# Faust

##Datasets
Faust datasets [ http://faust.is.tue.mpg.de/ ]

##Registration
Register from tr_reg_0i0.obj to tr_reg_0ij.obj, where i from 0 to 9 and j from 1 to 9. 

##Initialization
These examples use SHOT and diffusion pruning to get the initial correspondences, we save the initial correspondences in `stepinfo_file*/shot.txt` and `stepinfo_file*/prune.txt`.

##Results
We save the deformed meshes in 0ij.ply and the iteration information in 0ij.txt, in which the first column is the time, 4-th is normalized RMSE, and 6-th is the normalized scaling.
