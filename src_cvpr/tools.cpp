#include "tools.h"
#include "randutil.h"
#include <iostream>

#ifdef __linux__		
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

int GetBoundaryPoints(Mesh *polydata, bool* &boundaryPointsIds)
{
    boundaryPointsIds = new bool[polydata->n_vertices()];
    for(int po=0; po<polydata->n_vertices(); po++)
    {
        boundaryPointsIds[po] = false;
    }
	int num_boundaryPoints = 0;
	for (auto it = polydata->edges_begin(); it != polydata->edges_end(); it++)
	{
		if (polydata->is_boundary(*it))
		{
			int idx1 = polydata->halfedge_handle(*it, 0).idx();
			int idx2 = polydata->halfedge_handle(*it, 1).idx();
			if (!boundaryPointsIds[idx1])
			{
				boundaryPointsIds[idx1] = true;
				num_boundaryPoints++;
			}
			if (!boundaryPointsIds[idx2])
			{
				boundaryPointsIds[idx2] = true;
				num_boundaryPoints++;
			}
		}
	}
    return num_boundaryPoints;
}


double frand() { return rand() * 1.0 / RAND_MAX; }


double box_muller(double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
    static double y2;
    static int use_last = 0;

    if (use_last)		        /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * drand() - 1.0;
            x2 = 2.0 * drand() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );
}

double mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh)
{
    OpenMesh::Vec3d max(-1e12, -1e12, -1e12);
    OpenMesh::Vec3d min(1e12, 1e12, 1e12);
    for(auto it = src_mesh.vertices_begin(); it != src_mesh.vertices_end(); it++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(src_mesh.point(*it)[j] > max[j])
            {
                max[j] = src_mesh.point(*it)[j];
            }
            if(src_mesh.point(*it)[j] < min[j])
            {
                min[j] = src_mesh.point(*it)[j];
            }
        }
    }

    for(auto it = tar_mesh.vertices_begin(); it != tar_mesh.vertices_end(); it++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(tar_mesh.point(*it)[j] > max[j])
            {
                max[j] = tar_mesh.point(*it)[j];
            }
            if(tar_mesh.point(*it)[j] < min[j])
            {
                min[j] = tar_mesh.point(*it)[j];
            }
        }
    }
    double scale = (max-min).norm();

    for(auto it = src_mesh.vertices_begin(); it != src_mesh.vertices_end(); it++)
    {
        OpenMesh::Vec3d p = src_mesh.point(*it);
        p = p/scale;
        src_mesh.set_point(*it, p);
    }

    for(auto it = tar_mesh.vertices_begin(); it != tar_mesh.vertices_end(); it++)
    {
        OpenMesh::Vec3d p = tar_mesh.point(*it);
        p = p/scale;
        tar_mesh.set_point(*it, p);
    }

    return scale;
}

#ifdef __linux__
bool my_mkdir(std::string file_path)
{
    if(access(file_path.c_str(), 06))
   {
       std::cout << "file_path : (" << file_path << ") didn't exist or no write ability!!" << std::endl;
       if(mkdir(file_path.c_str(), S_IRWXU))
       {
           std::cout << "mkdir " << file_path << " is wrong! please check upper path " << std::endl;
           exit(0);
       }
       std::cout<< "mkdir " << file_path << " success!! " << std::endl;
   }
}
#endif