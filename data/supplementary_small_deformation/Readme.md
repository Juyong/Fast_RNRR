#supplementary_small_deformation 

##Datasets
Some models in Human-motion [http://people.csail.mit.edu/drdaniel/mesh_animation/#data/]

##Registration
"I-crane": 
Register from mesh_0022.obj to  mesh_0024.obj;

"march1": 
Register from mesh_0024.obj to  mesh_0026.obj;

"samba": 
Register from mesh_0015.obj to  mesh_0017.obj;

"squat1":
Register from mesh_0023.obj to  mesh_0025.obj;

"swing": 
Register from mesh_0020.obj to  mesh_0022.obj;

##Initialization
These examples use the nearest point as initial correspondence.

##Results
We save the deformed meshes in `NICP/*.ply`, `RPTS/*.ply` and `ours/*.ply`