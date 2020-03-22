#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshDefinition.h"
#include <iostream>
#include <fstream>

enum { TRIANGLE = 0, QUAD, N_MESH_MODES };
int checkMeshMode(Mesh& mesh)
{
    int mesh_mode;
    Mesh::FaceIter fIt = mesh.faces_begin();
    Mesh::FaceIter fEnd = mesh.faces_end();
    Mesh::FaceEdgeIter fe_it;
    int count = 1;
    int meshType[3] = {0};
    for(fIt; fIt != fEnd; ++fIt)
    {
        fe_it = mesh.fe_iter(*fIt);
        while(--fe_it)
        {
            ++count;
        }
        if(count == 4)
        {
            meshType[1]++;
        }
        else if(count == 3)
        {
            meshType[0]++;
        }
        else
        {
            meshType[2]++;
        }
        count = 1;
    }
    size_t faceNum = mesh.n_faces();
    if(meshType[0] == faceNum)//triangle
    {
        mesh_mode = TRIANGLE;
    }
    else if(meshType[1] == faceNum)//no
    {
        mesh_mode = QUAD;
    }
    else
    {
        mesh_mode = N_MESH_MODES;
    }
    return mesh_mode;
}

void printBasicMeshInfo(Mesh& mesh, bool use_demean)
{
    if (mesh.n_vertices() == 0)
        printf("No Mesh\n");

    switch(checkMeshMode(mesh))
    {
    case TRIANGLE:
        printf("Triangle Mesh.\n");
        break;
    case QUAD:
        printf("Quadrilateral Mesh.\n");
        break;
    default:
        printf("General Mesh.\n");
        break;
    }

    printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n",
        mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());

    OpenMesh::Vec3d MeshScales;
    MeshScales[0] = 0.0; MeshScales[1] = 0.0; MeshScales[2] = 0.0;
    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        MeshScales += mesh.point(*v_it);
    }
    MeshScales /= mesh.n_vertices();
    mesh.data(mesh.vertex_handle(0)).mesh_trans = MeshScales;

    if(use_demean)
    {
        for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        {
            mesh.point(*v_it) = mesh.point(*v_it)- MeshScales;
        }
    }
}

bool read_data(const std::string filename, Mesh& mesh, bool use_demean)
{
    bool read_OK = OpenMesh::IO::read_mesh(mesh, filename);

	std::cout << "filename = " << filename << std::endl;
    if (read_OK)
    {
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();

        mesh.request_face_normals();
        mesh.request_vertex_normals();

        printBasicMeshInfo(mesh, use_demean);

        mesh.update_face_normals();
        mesh.update_vertex_normals();

        return true;
    }
    std::cout << "#vertices = " << mesh.n_vertices() << std::endl;
    return false;
}

bool write_data(const char* filename, Mesh& mesh, bool use_demean, double scale)
{
    OpenMesh::Vec3d MeshScales;
    MeshScales = mesh.data(mesh.vertex_handle(0)).mesh_trans;

    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        mesh.point(*v_it) = mesh.point(*v_it)*scale;
    }

    if(use_demean)
    {
        for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        {
            mesh.point(*v_it) = mesh.point(*v_it) + MeshScales;
        }
    }
    bool ok = OpenMesh::IO::write_mesh(mesh, filename);
	return ok;
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

std::string num2str(int num, const int size, bool is_add0)
{
    char s_id[10]={0};
    s_id[size] = '\0';
    int pos = 0;
    for (int i = size-1;i >= 0;--i)
    {
        s_id[i] = num % 10 + '0';
        num /= 10;
        if(!is_add0)
        {
            if(num==0)
            {
                pos = i;
                is_add0 = true;
            }
        }
    }
    std::string s_num(s_id);
    return s_num.substr(pos, size-pos);
}

