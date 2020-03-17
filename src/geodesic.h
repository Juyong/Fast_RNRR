#include "MeshDefinition.h"
#include <vector>
#include <Eigen/Sparse>
#include <igl/heat_geodesics.h>
#include <igl/read_triangle_mesh.h>
#include <igl/exact_geodesic.h>

#ifndef _GEODESIC_H_
#define _GEODESIC_H_

namespace svr
{
	class geodesic
	{
	public:
		geodesic() = delete;

		geodesic(Mesh& _Mesh);
		geodesic(const Mesh &_Mesh) : geodesic(const_cast<Mesh&>(_Mesh)) {}
		bool seed(int _seedIdx, Eigen::VectorXd &_GeoDist);

	private:

		double m_timeStep;
		int m_vertexNum;
		int m_faceNum;
		int m_edgeNum;
		
		const Mesh &m_Mesh;
		
		Eigen::SparseMatrix<double> A0;
		Eigen::SparseMatrix<double> A1;
		Eigen::SparseMatrix<double> A2;
		Eigen::SparseMatrix<double> A3;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> Solver0;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> Solver1;
	};

//////////////////////////////////////////////////////////////////////////
}	//	namespace DynamicReconstruction end
//////////////////////////////////////////////////////////////////////////

class hgeodesic
{
public:
    ~hgeodesic();
    hgeodesic(Mesh& src_mesh);
    hgeodesic(Mesh& src_mesh, Mesh& tar_mesh);
    igl::HeatGeodesicsData<double> data1, data2;
    int n_src_vertices;
    int n_tar_vertices;

private:
    Eigen::MatrixXd V1, V2;
    Eigen::MatrixXi F1, F2;
    void    Mesh2VF(Mesh & mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
};

#endif
