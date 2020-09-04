#ifndef MESHDEFINITION_H
#define MESHDEFINITION_H

#include <OpenMesh/Core/Geometry/VectorT.hh>
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	FaceTraits
	{
	};

	EdgeTraits
	{
	};

	HalfedgeTraits
	{
		HalfedgeT() :face_he_var(-1)
		{
		};
	private:
		int face_he_var;
	public:
		int get_face_he_var() { return face_he_var; };
		void set_face_he_var(int fhe) { face_he_var = fhe; };
	};

	VertexTraits
	{
		VertexT() : new_pos_fixed(false)
		{
		};
	private:
		OpenMesh::Vec3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
	public:
		void set_New_Pos(const OpenMesh::Vec3d& n_p){ new_pos = n_p; };
		OpenMesh::Vec3d& get_New_Pos(){ return new_pos; };
		void set_new_pos_fixed(bool f){ new_pos_fixed = f; };
		bool get_new_pos_fixed(){ return new_pos_fixed; };
        OpenMesh::Vec3d mesh_trans;
        double dist;
        int prev;
		double dist_to_corres = 0.0;
		double angle_to_corres = 0.0;
	};

};

//typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> Mesh;
typedef OpenMesh::Vec3f								Vec3f;
typedef OpenMesh::Vec3d								Vec3d;

typedef Mesh::Point									Point;
typedef Mesh::Color									Color;
typedef Mesh::VertexHandle							VertexHandle;
typedef Mesh::FaceHandle							FaceHandle;
typedef Mesh::EdgeHandle							EdgeHandle;
typedef Mesh::HalfedgeHandle						HalfedgeHandle;

typedef Mesh::VertexIter							VertexIter;
typedef Mesh::EdgeIter								EdgeIter;
typedef Mesh::FaceIter								FaceIter;
typedef Mesh::HalfedgeIter							HalfedgeIter;

typedef Mesh::VertexVertexIter						VertexVertexIter;
typedef Mesh::VertexFaceIter						VertexFaceIter;
typedef Mesh::VertexIHalfedgeIter					VertexIHalfedgeIter;
typedef Mesh::VertexOHalfedgeIter					VertexOHalfedgeIter;

typedef Mesh::FaceVertexIter						FaceVertexIter;
typedef Mesh::ConstFaceVertexIter					ConstFaceVertexIter;
typedef Mesh::ConstHalfedgeIter						ConstHalfedgeIter;
typedef Mesh::ConstVertexOHalfedgeIter				ConstVertexOHalfedgeIter;
typedef Mesh::ConstVertexIHalfedgeIter				ConstVertexIHalfedgeIter;

bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);//just copy the code from openmesh
bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);

bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p);


#endif
