#include "geodesic.h"
#include <cmath>
#include <iostream>

#pragma warning(disable: 4267)

namespace svr
{
	//////////////////////////////////////////////////////////////////////////
	//	Define Geodesic constructor
	//////////////////////////////////////////////////////////////////////////
	geodesic::geodesic(Mesh &_Mesh) : m_Mesh(_Mesh)
	{
		m_vertexNum = m_Mesh.n_vertices();
		m_faceNum = m_Mesh.n_faces();
		m_edgeNum = m_Mesh.n_edges();

		//////////////////////////////////////////////////////////////////////////
		//	Resize A0, A1, A2, A3
		//////////////////////////////////////////////////////////////////////////
		A0.resize(m_vertexNum, m_vertexNum);
		A1.resize(3*m_faceNum, m_vertexNum);
		A2.resize(m_vertexNum, 3*m_faceNum);
		A3.resize(m_vertexNum, m_vertexNum);

		//////////////////////////////////////////////////////////////////////////
		//	Calculate area and normal for each face of the mesh
		//	Area and normals of faces can be reused, so it is better to calculate
		//	them separately.
		//////////////////////////////////////////////////////////////////////////
		Eigen::VectorXd FaceArea(m_faceNum);
		std::vector<Vec3d> FaceNormal(m_faceNum);
		ConstHalfedgeIter pHalfedge = m_Mesh.halfedges_begin();
		ConstHalfedgeIter pHalfedgeEnd = m_Mesh.halfedges_end();

		for (; pHalfedge != pHalfedgeEnd; ++pHalfedge)
		{
			HalfedgeHandle h0 = *pHalfedge;
			if (!m_Mesh.is_boundary(h0))
			{
				HalfedgeHandle h1 = m_Mesh.next_halfedge_handle(h0);
				HalfedgeHandle h2 = m_Mesh.next_halfedge_handle(h1);
				VertexHandle vh0 = m_Mesh.to_vertex_handle(h2);
				VertexHandle vh1 = m_Mesh.to_vertex_handle(h0);
				VertexHandle vh2 = m_Mesh.to_vertex_handle(h1);

				Point v0 = m_Mesh.point(vh0);
				Point v1 = m_Mesh.point(vh1);
				Point v2 = m_Mesh.point(vh2);

				/*Vec3d e1 = OpenMesh::vector_cast<Vec3d, Vec3f>(v1 - v0);
				Vec3d e2 = OpenMesh::vector_cast<Vec3d, Vec3f>(v2 - v0);*/
				Vec3d e1 = v1 - v0;
				Vec3d e2 = v2 - v0;

				Vec3d fN = OpenMesh::cross(e1, e2);
				double area = fN.norm() / 2.0;

				FaceHandle fh = m_Mesh.face_handle(h0);
				FaceArea(fh.idx()) = area;
				FaceNormal[fh.idx()] = fN.normalize();
			}
		}

		//////////////////////////////////////////////////////////////////////////
		//	Calculate average edge length and time step
		//	timeStep = averageEdgeLen * averageEdgeLen
		//
		//	Note that when dealing with nonuniform triangle mesh, you should choose
		//	maximum edge length as spacing length.
		//	timeStep = maximumEdgeLen * maximumEdgeLen
		//////////////////////////////////////////////////////////////////////////
#if 0
		double averageEdgeLen = 0.0;

		for (int i = 0; i < m_edgeNum; ++i)
		{
			EdgeHandle eh = m_Mesh.edge_handle(i);
			double edgeLen = m_Mesh.calc_edge_length(eh);
			averageEdgeLen += (double)edgeLen;
		}
		averageEdgeLen /= m_edgeNum;
		m_timeStep = 1.0 * averageEdgeLen * averageEdgeLen;
#else
		double maximumEdgeLen = 0.0;

		for (int i = 0; i < m_edgeNum; ++i)
		{
			EdgeHandle eh = m_Mesh.edge_handle(i);
			double edgeLen = m_Mesh.calc_edge_length(eh);
			if (edgeLen > maximumEdgeLen)
			{
				maximumEdgeLen = edgeLen;
			}
		}
		m_timeStep = 1.0 * maximumEdgeLen * maximumEdgeLen;
#endif

		//////////////////////////////////////////////////////////////////////////
		//	Construct VoronoiAreaMatrix A
		//	Construct CotangentWeightMatrix Lc
		//	
		//	A is a diagonal sparse matrix. Each element in diagonal represents
		//	Voronoi area for each vertex.
		//
		//	Lc is a symmetric sparse matrix. Each non-zero element represents cotangent
		//	weights for each vertex.
		//////////////////////////////////////////////////////////////////////////
		std::vector<Eigen::Triplet<double>> VoronoiAreaTripletList;
		std::vector<Eigen::Triplet<double>> CotangentWeightTripletList;
		Eigen::SparseMatrix<double> VoronoiAreaMatrix(m_vertexNum, m_vertexNum);
		Eigen::SparseMatrix<double> CotangentWeightMatrix(m_vertexNum, m_vertexNum);

		for (int i = 0; i < m_vertexNum; ++i)
		{
			VertexHandle vh0 = m_Mesh.vertex_handle(i);

			double vertexVoronoiArea = 0.0;
			double vertexCotangentWeight = 0.0;

			Mesh::ConstVertexOHalfedgeIter pVertexOHalfedge = m_Mesh.cvoh_iter(vh0);
			for (; pVertexOHalfedge.is_valid(); ++pVertexOHalfedge)
			{
				HalfedgeHandle h0 = *pVertexOHalfedge;

				if (!m_Mesh.is_boundary(h0))
				{
					HalfedgeHandle h1 = m_Mesh.next_halfedge_handle(h0);
					HalfedgeHandle t0 = m_Mesh.opposite_halfedge_handle(h0);
					HalfedgeHandle t1 = m_Mesh.next_halfedge_handle(t0);
					FaceHandle fh = m_Mesh.face_handle(h0);
					VertexHandle vh1 = m_Mesh.to_vertex_handle(h0);
					VertexHandle vh2 = m_Mesh.to_vertex_handle(h1);
					VertexHandle vh3 = m_Mesh.to_vertex_handle(t1);

					Point v0 = m_Mesh.point(vh0);
					Point v1 = m_Mesh.point(vh1);
					Point v2 = m_Mesh.point(vh2);
					Point v3 = m_Mesh.point(vh3);

					/*Vec3d e0 = OpenMesh::vector_cast<Vec3d, Vec3f>(v1 - v2);
					Vec3d e1 = OpenMesh::vector_cast<Vec3d, Vec3f>(v0 - v2);
					Vec3d e2 = OpenMesh::vector_cast<Vec3d, Vec3f>(v1 - v3);
					Vec3d e3 = OpenMesh::vector_cast<Vec3d, Vec3f>(v0 - v3);*/
					Vec3d e0 = (v1 - v2);
					Vec3d e1 = (v0 - v2);
					Vec3d e2 = (v1 - v3);
					Vec3d e3 = (v0 - v3);

					Vec3d fN = cross(e1, e0);

					double alpha = acos(OpenMesh::dot(e0, e1) / (e0.length() * e1.length()));
					double beta = acos(OpenMesh::dot(e2, e3) / (e2.length() * e3.length()));
					double fArea = fN.norm() / 2.0;
					double cotangentWeight = 0.5 * (1.0 / tan(alpha) + 1.0 / tan(beta));

					vertexCotangentWeight += -cotangentWeight;
					vertexVoronoiArea += fArea;

					CotangentWeightTripletList.push_back(Eigen::Triplet<double>(vh0.idx(), vh1.idx(), cotangentWeight));
				}
			}

			CotangentWeightTripletList.push_back(Eigen::Triplet<double>(vh0.idx(), vh0.idx(), vertexCotangentWeight));
			VoronoiAreaTripletList.push_back(Eigen::Triplet<double>(vh0.idx(), vh0.idx(), 1.0 / 3.0 * vertexVoronoiArea));
		}

		VoronoiAreaMatrix.setFromTriplets(VoronoiAreaTripletList.begin(), VoronoiAreaTripletList.end());
		CotangentWeightMatrix.setFromTriplets(CotangentWeightTripletList.begin(), CotangentWeightTripletList.end());

		//////////////////////////////////////////////////////////////////////////
		//	Setup A0 = A - t*Lc, t represents time step
		//////////////////////////////////////////////////////////////////////////
		A0 = VoronoiAreaMatrix - m_timeStep * CotangentWeightMatrix;
		
		//////////////////////////////////////////////////////////////////////////
		//	Setup A3 = Lc + diag(gamma)
		//	gamma is a tiny disturbance for diagonal element, and gamma = 1e-10.
		//////////////////////////////////////////////////////////////////////////
		A3 = CotangentWeightMatrix;

#if 0
		double gamma = 1e-10;
		for (int i = 0; i < m_vertexNum; ++i)
		{
			A3.coeffRef(i, i) += gamma;
		}
#endif	
		//////////////////////////////////////////////////////////////////////////
		//	Setup A1, Grad(U) = A1 * U
		//
		//	A1 is a sparse matrix, and not symmetric. Multiplying A1 to U is
		//	equivalent to calculating heat gradient of each face.
		//
		//	U is a vector of size vertexNum. Each element u is the heat kernel of
		//	each vertex.
		//
		//	Grad(U) is also a vector, but with size of 3*faceNum. The element j at
		//	location 3*j+0, 3*j+1 and 3*j+2 represents the gradient of heat kernel
		//	for face j.
		//	
		//	Note that in order to calculate geodesic distance at each vertex, I
		//	need to know gradient of distance at this vertex. Therefore, after I
		//	get Grad(U), I also need to normalize each vector element j, and negate
		//	it, since Eikonal equation ||phi|| = 1 requires unit vector and heat
		//	gradient is opposite to distance gradient. This is done separately in
		//	Seed(), for each seed.
		//////////////////////////////////////////////////////////////////////////
		std::vector<Eigen::Triplet<double>> A1TripletList;
		for (int i = 0; i < m_faceNum; ++i)
		{
			FaceHandle fh = m_Mesh.face_handle(i);
			ConstFaceVertexIter fvIter = m_Mesh.cfv_iter(fh);
			VertexHandle vh0 = *fvIter; ++fvIter;
			VertexHandle vh1 = *fvIter; ++fvIter;
			VertexHandle vh2 = *fvIter;

			Point v0 = m_Mesh.point(vh0);
			Point v1 = m_Mesh.point(vh1);
			Point v2 = m_Mesh.point(vh2);

			/*Vec3d e0 = OpenMesh::vector_cast<Vec3d, Vec3f>(v2 - v1);
			Vec3d e1 = OpenMesh::vector_cast<Vec3d, Vec3f>(v0 - v2);
			Vec3d e2 = OpenMesh::vector_cast<Vec3d, Vec3f>(v1 - v0);*/
			Vec3d e0 = (v2 - v1);
			Vec3d e1 = (v0 - v2);
			Vec3d e2 = (v1 - v0);

			int idx0 = vh0.idx();
			int idx1 = vh1.idx();
			int idx2 = vh2.idx();
			int fIdx = fh.idx();

			Vec3d fN = FaceNormal[fIdx];
			double area = 2*FaceArea[fIdx];

			Vec3d grad0 = OpenMesh::cross(fN, e0);
			Vec3d grad1 = OpenMesh::cross(fN, e1);
			Vec3d grad2 = OpenMesh::cross(fN, e2);

			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+0, idx0, grad0[0]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+0, idx1, grad1[0]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+0, idx2, grad2[0]/area));

			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+1, idx0, grad0[1]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+1, idx1, grad1[1]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+1, idx2, grad2[1]/area));

			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+2, idx0, grad0[2]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+2, idx1, grad1[2]/area));
			A1TripletList.push_back(Eigen::Triplet<double>(3*fIdx+2, idx2, grad2[2]/area));
		}

		A1.setFromTriplets(A1TripletList.begin(), A1TripletList.end());
		
		//////////////////////////////////////////////////////////////////////////
		//	Setup A2, Div(X) = A2 * X
		//
		//	Notation:
		//	G = Grad(U);
		//	gj0 = G(3*j+0), gj1 = G(3*j+1), gj2 = G(3*j+2);
		//	gj = [gj0, gj1, gj2].norm();
		//	[X(3*j+0), X(3*j+1), X(3*j+2)] = [-gj0, -gj1, -gj2] / gj;
		//
		//	A2 is a sparse matrix, and not symmetric. Multiply A2 to X is
		//	equivalent to calculating divergent of each vertex of the mesh.
		//////////////////////////////////////////////////////////////////////////
		std::vector<Eigen::Triplet<double>> A2TripletList;
		for (int i = 0; i < m_vertexNum; ++i)
		{
			VertexHandle vh0 = m_Mesh.vertex_handle(i);
			ConstVertexOHalfedgeIter pVertexOHalfedge = m_Mesh.cvoh_iter(vh0);
			for (; pVertexOHalfedge.is_valid(); ++pVertexOHalfedge)
			{
				HalfedgeHandle h0 = *pVertexOHalfedge;
				if (!m_Mesh.is_boundary(h0))
				{
					HalfedgeHandle h1 = m_Mesh.next_halfedge_handle(h0);
					VertexHandle vh1 = m_Mesh.to_vertex_handle(h0);
					VertexHandle vh2 = m_Mesh.to_vertex_handle(h1);
					Point v0 = m_Mesh.point(vh0);
					Point v1 = m_Mesh.point(vh1);
					Point v2 = m_Mesh.point(vh2);

					/*Vec3d e1 = OpenMesh::vector_cast<Vec3d, Vec3f>(v1 - v0);
					Vec3d e2 = OpenMesh::vector_cast<Vec3d, Vec3f>(v2 - v0);
					Vec3d e0 = OpenMesh::vector_cast<Vec3d, Vec3f>(v2 - v1);*/
					Vec3d e1 = (v1 - v0);
					Vec3d e2 = (v2 - v0);
					Vec3d e0 = (v2 - v1);

					double theta1 = acos(OpenMesh::dot(e0, e2) / (e0.length() * e2.length()));
					double theta2 = acos(OpenMesh::dot(-e0, e1) / (e0.length() * e1.length()));

					double cotTheta1 = 0.5 / tan(theta1);
					double cotTheta2 = 0.5 / tan(theta2);

					FaceHandle fh = m_Mesh.face_handle(h0);
					int fIdx = fh.idx();

					A2TripletList.push_back(Eigen::Triplet<double>(i, 3*fIdx+0, cotTheta1*e1[0]+cotTheta2*e2[0]));
					A2TripletList.push_back(Eigen::Triplet<double>(i, 3*fIdx+1, cotTheta1*e1[1]+cotTheta2*e2[1]));
					A2TripletList.push_back(Eigen::Triplet<double>(i, 3*fIdx+2, cotTheta1*e1[2]+cotTheta2*e2[2]));
				}
			}
		}

		A2.setFromTriplets(A2TripletList.begin(), A2TripletList.end());

		//////////////////////////////////////////////////////////////////////////
		//	Decompose A0 and A3 using SimplicialLDLT
		//////////////////////////////////////////////////////////////////////////
		Solver0.compute(A0);
		Solver1.compute(A3);
	}

	//////////////////////////////////////////////////////////////////////////
	//	Define Seed member function
	//////////////////////////////////////////////////////////////////////////
	bool geodesic::seed(int seedIdx, Eigen::VectorXd &_GeoDist)
	{	
		//	Check whether source vertex index is valid
		assert(seedIdx >= 0 && seedIdx < m_vertexNum);
		if (seedIdx < 0 || seedIdx >= m_vertexNum)
		{
			std::cerr << "ERROR: Invalid seed index !" << std::endl;
			return false;
		}

		//////////////////////////////////////////////////////////////////////////
		//	Setup working vectors
		//
		//	HeatKernel = U;
		//	HeatGrad = Grad(U);
		//	DistGrad = X;
		//	DistDiv = Div(X);
		//	HeatDelta = delta;
		//	GeoDist = phi
		//////////////////////////////////////////////////////////////////////////
		Eigen::VectorXd HeatKernel(m_vertexNum);
		Eigen::VectorXd HeatGrad(3*m_faceNum);
		Eigen::VectorXd DistGrad(3*m_faceNum);
		Eigen::VectorXd DistDiv(m_vertexNum);
		Eigen::VectorXd GeoDist(m_vertexNum);

		Eigen::VectorXd HeatDelta(m_vertexNum);
		memset((void*)HeatDelta.data(), 0, HeatDelta.size() * sizeof(double));
		HeatDelta(seedIdx) = 1.0;

		//////////////////////////////////////////////////////////////////////////
		//	Solve geodesic distance for each seed
		//////////////////////////////////////////////////////////////////////////
		HeatKernel = Solver0.solve(HeatDelta);
		HeatGrad = A1 * HeatKernel;

#pragma omp parallel for
		for (int j = 0; j < m_faceNum; ++j)
		{
			Eigen::Vector3d grad;
			Eigen::Vector3d gradAbs;

			grad(0) = HeatGrad(3 * j + 0);
			grad(1) = HeatGrad(3 * j + 1);
			grad(2) = HeatGrad(3 * j + 2);

			gradAbs(0) = fabs(grad(0));
			gradAbs(1) = fabs(grad(1));
			gradAbs(2) = fabs(grad(2));

			double maxAbs = gradAbs.maxCoeff();

			grad(0) /= maxAbs;
			grad(1) /= maxAbs;
			grad(2) /= maxAbs;

			grad /= grad.norm();

			DistGrad(3 * j + 0) = -grad(0);
			DistGrad(3 * j + 1) = -grad(1);
			DistGrad(3 * j + 2) = -grad(2);
		}

		DistDiv = A2 * DistGrad;
		GeoDist = Solver1.solve(DistDiv);

		double constShift = GeoDist(seedIdx);
		for (int i = 0; i < m_vertexNum; ++i)
		{
			GeoDist(i) -= constShift;
		}

		//////////////////////////////////////////////////////////////////////////
		//	Convert VectorXd to VectorXf
		//////////////////////////////////////////////////////////////////////////
		_GeoDist = GeoDist.cast<double>();

		return true;
	}
//////////////////////////////////////////////////////////////////////////
}	//	namespace DynamicReconstruction end
//////////////////////////////////////////////////////////////////////////
hgeodesic::hgeodesic(Mesh& src_mesh, Mesh& tar_mesh)
{

    if(src_mesh.n_faces())
    {
        Mesh2VF(src_mesh, V1, F1);
        data1.use_intrinsic_delaunay = true;
        double t1 = std::pow(igl::avg_edge_length(V1, F1), 2);
        igl::heat_geodesics_precompute(V1, F1, t1, data1);
    }
    else
        std::cout << "WARNING!! source mesh no face!" << std::endl;
    if(tar_mesh.n_faces())
    {
        Mesh2VF(tar_mesh, V2, F2);
        data2.use_intrinsic_delaunay = true;
        double t2 = std::pow(igl::avg_edge_length(V2, F2), 2);
        igl::heat_geodesics_precompute(V2, F2, t2, data2);
    }
    else
        std::cout << "WARNING!! target mesh no face!" << std::endl;

    n_src_vertices = src_mesh.n_vertices();
    n_tar_vertices = tar_mesh.n_vertices();
}

hgeodesic::hgeodesic(Mesh& src_mesh)
{
    Mesh2VF(src_mesh, V1, F1);

    data1.use_intrinsic_delaunay = true;

    double t1 = std::pow(igl::avg_edge_length(V1, F1), 2);

    igl::heat_geodesics_precompute(V1, F1, t1, data1);

    n_src_vertices = src_mesh.n_vertices();
    n_tar_vertices = 0;
}

hgeodesic::~hgeodesic()
{

}

void hgeodesic::Mesh2VF(Mesh & mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
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
