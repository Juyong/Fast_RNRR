#pragma once
#include <iostream>
#include <fstream>
#include "Types.h"

enum { TRIANGLE = 0, QUAD, N_MESH_MODES};
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

void printBasicMeshInfo(Mesh& mesh)
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

    printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n\n",
        mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());
}

bool read_data(const std::string filename, Mesh& mesh)
{
    OpenMesh::IO::Options opt_read = OpenMesh::IO::Options::VertexNormal;
    mesh.request_vertex_normals();
    bool read_OK = OpenMesh::IO::read_mesh(mesh, filename,opt_read);

	std::cout << "filename = " << filename << std::endl;
    if (read_OK)
    {
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();

        mesh.request_face_normals();
        printBasicMeshInfo(mesh);

        mesh.update_face_normals();
        if(mesh.n_faces()>0)
            mesh.update_vertex_normals();

        Vec3 MeshScales;
        MeshScales[0] = 0.0; MeshScales[1] = 0.0; MeshScales[2] = 0.0;
        for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        {
            MeshScales += mesh.point(*v_it);
        }
        MeshScales /= mesh.n_vertices();
        return true;
    }
    std::cout << "#vertices = " << mesh.n_vertices() << std::endl;
    return false;
}

bool write_data(const char* filename, Mesh& mesh, Scalar scale)
{
    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {

        mesh.point(*v_it) = mesh.point(*v_it)*scale;
    }
    bool ok = OpenMesh::IO::write_mesh(mesh, filename);
	return ok;
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

