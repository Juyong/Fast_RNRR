#Noise

##Datasets
"Jumping" in Human-motion [http://people.csail.mit.edu/drdaniel/mesh_animation/#data/]

##Registration
0.3 dense noise (d03) and 0.7 dense noise (d07): 
Register from mesh_0000.obj to  mesh_0011.obj.

5% sparse noise (s005) and 50% sparse noise (s05): 
Register from mesh_0000.obj to  mesh_0018.obj.

##Initialization
These examples use SHOT and diffusion pruning(d03 & s005) and 60 manually labeled landmarks(d07 & s05) to get the initial correspondence, we save the initial correspondence in `init/landmark.txt`.


##Results
We save the deformed meshes in `NICP/*.ply`, `RPTS/*.ply` and `ours/*.ply`
