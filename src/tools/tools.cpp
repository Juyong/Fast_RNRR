#include "tools.h"
#include <iostream>

#ifdef __linux__		
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

const float SHRINKING_FACTOR = 7.5f;
const int NO_PROGRESS_STREAK_THRESHOLD = 100;

Vec3 Eigen2Vec(Vector3 s)
{
    return Vec3(s[0], s[1], s[2]);
}

Vector3 Vec2Eigen(Vec3 s)
{
    return Vector3(s[0], s[1], s[2]);
}

Scalar mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh)
{
    Vec3 max(-1e12, -1e12, -1e12);
    Vec3 min(1e12, 1e12, 1e12);
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
    Scalar scale = (max-min).norm();

    for(auto it = src_mesh.vertices_begin(); it != src_mesh.vertices_end(); it++)
    {
        Vec3 p = src_mesh.point(*it);
        p = p/scale;
        src_mesh.set_point(*it, p);
    }

    for(auto it = tar_mesh.vertices_begin(); it != tar_mesh.vertices_end(); it++)
    {
        Vec3 p = tar_mesh.point(*it);
        p = p/scale;
        tar_mesh.set_point(*it, p);
    }

    return scale;
}

void Mesh2VF(Mesh & mesh, MatrixXX& V, Eigen::MatrixXi& F)
{
    V.resize(mesh.n_vertices(),3);
    F.resize(mesh.n_faces(),3);
    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++)
    {
        V(it->idx(), 0) = mesh.point(*it)[0];
        V(it->idx(), 1) = mesh.point(*it)[1];
        V(it->idx(), 2) = mesh.point(*it)[2];
    }

    for (auto fit = mesh.faces_begin(); fit != mesh.faces_end(); fit++)
    {
        int i = 0;
        for (auto vit = mesh.fv_begin(*fit); vit != mesh.fv_end(*fit); vit++)
        {
            F(fit->idx(), i) = vit->idx();
            i++;
            if (i > 3)
            {
                std::cout << "Error!! one face has more than 3 points!" << std::endl;
                break;
            }
        }
    }
}

bool read_landmark(const char* filename, std::vector<int>& landmark_src, std::vector<int>& landmark_tar)
{
    std::ifstream in(filename);
    std::cout << "filename = " << filename << std::endl;
    if (!in)
    {
        std::cout << "Can't open the landmark file!!" << std::endl;
        return false;
    }
    int x, y;
    landmark_src.clear();
    landmark_tar.clear();
    while (!in.eof())
    {
        if (in >> x >> y) {
            landmark_src.push_back(x);
            landmark_tar.push_back(y);
        }
    }
    in.close();
    std::cout << "landmark_src = " << landmark_src.size() << " tar = " << landmark_tar.size() << std::endl;
    return true;
}

bool read_fixedvex(const char* filename, std::vector<int>& vertices_list)
{
    std::ifstream in(filename);
    std::cout << "filename = " << filename << std::endl;
    if (!in)
    {
        std::cout << "Can't open the landmark file!!" << std::endl;
        return false;
    }
    int x;
    vertices_list.clear();
    while (!in.eof())
    {
        if (in >> x) {
            vertices_list.push_back(x);
        }
    }
    in.close();
    std::cout << "the number of fixed vertices = " << vertices_list.size() << std::endl;
    return true;
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
