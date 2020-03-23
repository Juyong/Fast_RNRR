#supplementary_big_deformation 

##Datasets
"I-crane" and "swing" in Human-motion [http://people.csail.mit.edu/drdaniel/mesh_animation/#data/]

##Registration
"I-crane": 
Register from mesh_0000.obj to mesh_0019.obj; 
Register from mesh_0000.obj to mesh_0030.obj.

"swing": 
Register from mesh_0000.obj to  mesh_0022.obj;
Register from mesh_0000.obj to  mesh_0025.obj.

##Initialization
These examples use SHOT and diffusion pruning to get the initial correspondence, we save the initial correspondences in `stepinfo_file*/shot.txt` and `stepinfo_file*/prune.txt`.

##Results
We save the deformed meshes in `NICP/*.ply`, `RPTS/*.ply` and `ours/*.ply`