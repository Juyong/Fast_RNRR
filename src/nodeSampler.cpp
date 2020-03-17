//#pragma once
#include "nodeSampler.h"
#include "geodesic.h"
#include <memory>
#include <sstream>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <median.h>

//#include <igl/heat_geodesics.h>
//#include <igl/read_triangle_mesh.h>
//#include <igl/exact_geodesic.h>

namespace svr
{
    //	Define helper functions
    static auto square = [](const double argu) { return argu * argu; };
    static auto cube = [](const double argu) { return argu * argu * argu; };
    static auto max = [](const double lhs, const double rhs) { return lhs > rhs ? lhs : rhs; };

    //------------------------------------------------------------------------
    //	Node Sampling based on geodesic distance metric
    //
    //	Note that this member function samples nodes along some axis.
    //	Each node is not covered by any other node. And distance between each
    //	pair of nodes is at least sampling radius.
    //------------------------------------------------------------------------
    void nodeSampler::sample(const Mesh &mesh, igl::HeatGeodesicsData<double>& data,double sampleRadiusRatio, sampleAxis axis)
    {
        //	Save numbers of vertex and edge
        m_meshVertexNum = mesh.n_vertices();
        m_meshEdgeNum = mesh.n_edges();

        //	Calculate average edge length of bound mesh
        for (size_t i = 0; i < m_meshEdgeNum; ++i)
        {
            OpenMesh::EdgeHandle eh = mesh.edge_handle(i);
            double edgeLen = mesh.calc_edge_length(eh);
            m_averageEdgeLen += edgeLen;
        }
        m_averageEdgeLen /= m_meshEdgeNum;

        //	Sampling radius is calculated as averageEdgeLen multiplied by sampleRadiusRatio
        m_sampleRadius = sampleRadiusRatio * m_averageEdgeLen;


        //	Reorder mesh vertex along axis
        std::vector<size_t> vertexReorderedAlongAxis(m_meshVertexNum);
        size_t vertexIdx = 0;
        std::generate(vertexReorderedAlongAxis.begin(), vertexReorderedAlongAxis.end(), [&vertexIdx]() -> size_t { return vertexIdx++; });
        std::sort(vertexReorderedAlongAxis.begin(), vertexReorderedAlongAxis.end(), [&mesh, axis](const size_t &lhs, const size_t &rhs) -> bool {
            size_t lhsIdx = lhs;
            size_t rhsIdx = rhs;
            OpenMesh::VertexHandle vhl = mesh.vertex_handle(lhsIdx);
            OpenMesh::VertexHandle vhr = mesh.vertex_handle(rhsIdx);
            Mesh::Point vl = mesh.point(vhl);
            Mesh::Point vr = mesh.point(vhr);
            return vl[axis] > vr[axis];
        });

        //	Sample nodes using radius of m_sampleRadius
        size_t firstVertexIdx = vertexReorderedAlongAxis[0];
        std::unique_ptr<geodesic> pGeodesic = std::make_unique<geodesic>(mesh);
        Eigen::VectorXd geoDistVector(m_meshVertexNum);
        geoDistVector.setZero();
        pGeodesic->seed(firstVertexIdx, geoDistVector);
        m_geoDistContainer.push_back(geoDistVector);
        m_nodeContainer.emplace_back(0, firstVertexIdx);

        for (auto &vertexIdx : vertexReorderedAlongAxis)
        {
            bool IsNode = true;
            for (int k = 0; k < m_geoDistContainer.size(); ++k)
            {
                double dist = m_geoDistContainer.at(k)(vertexIdx);
                if (dist < m_sampleRadius)
                {
                    IsNode = false;
                    break;
                }
            }
            if (IsNode)
            {
                geoDistVector.setZero();
                pGeodesic->seed(vertexIdx, geoDistVector);
                m_geoDistContainer.push_back(geoDistVector);
                m_nodeContainer.emplace_back(m_geoDistContainer.size() - 1, vertexIdx);
            }
        }

        std::cout << "m_radius = " << m_sampleRadius << std::endl;
    }

    // type =1 fixed radius, type = 0, fixed number of nodes;
	void nodeSampler::farthest_point_sample(const Mesh &mesh, igl::HeatGeodesicsData<double>& data, double sampleRadiusRatio, int n_nodes, int type)
	{
		
        ////////////////////////////////////////////////////////////////////////////////////////////
		//	Save numbers of vertex and edge
		m_meshVertexNum = mesh.n_vertices();
		m_meshEdgeNum = mesh.n_edges();

		//	Calculate average edge length of bound mesh
		for (size_t i = 0; i < m_meshEdgeNum; ++i)
		{
			OpenMesh::EdgeHandle eh = mesh.edge_handle(i);
			double edgeLen = mesh.calc_edge_length(eh);
			m_averageEdgeLen += edgeLen;
		}
		m_averageEdgeLen /= m_meshEdgeNum;

		//	Sampling radius is calculated as averageEdgeLen multiplied by sampleRadiusRatio
		m_sampleRadius = sampleRadiusRatio * m_averageEdgeLen;

		std::cout << "sampleRadiusRatio = " << sampleRadiusRatio << " m_averageEdgeLen = "
			<< m_sampleRadius << std::endl;

		int firstVertexIdx = 0;
		Eigen::VectorXd geoDistVector;
		std::unique_ptr<geodesic> pGeodesic = std::make_unique<geodesic>(mesh);
		pGeodesic->seed(firstVertexIdx, geoDistVector);
		m_geoDistContainer.push_back(geoDistVector);

		//Eigen::MatrixXd all_geo_dists(m_meshVertexNum, m_meshVertexNum);
		//all_geo_dists.col(0) = geoDistVector;

		std::cout << "sample first vertex idx = " << firstVertexIdx << std::endl;

		m_nodeContainer.emplace_back(0, firstVertexIdx);

		bool stop = false;
		int num_cur_nodes = 1;
		int vertexIdx=0;
		double radius=0.0;
		std::vector<bool> is_nodes(m_meshVertexNum, false);
		is_nodes[firstVertexIdx] = true;
		Eigen::VectorXd min_geodist2set = geoDistVector;

		while (!stop)
		{
			// get next vertexidx
		    min_geodist2set.maxCoeff(&vertexIdx);		

			geoDistVector.setZero();
			pGeodesic->seed(vertexIdx, geoDistVector);
			m_geoDistContainer.push_back(geoDistVector);
			m_nodeContainer.emplace_back(m_geoDistContainer.size() - 1, vertexIdx);
			//all_geo_dists.col(num_cur_nodes) = geoDistVector;
			min_geodist2set = min_geodist2set.cwiseMin(geoDistVector);
			num_cur_nodes++;

			// get radius
			Eigen::VectorXd all_radius(num_cur_nodes);
			for (int i = 0; i < num_cur_nodes; i++)
			{
				Eigen::VectorXd temp(num_cur_nodes);
				temp.setZero();
				for (int j = 0; j < num_cur_nodes; j++)
				{
					if (j == i)
						temp[j] = 1e12;
					else
						temp[j] = m_geoDistContainer[m_nodeContainer[i].first][m_nodeContainer[j].second];
				}
				all_radius[i] = temp.minCoeff();
			}
			radius = all_radius.maxCoeff();
			// std::cout << "radius = " << radius << " num_cur_nodes = " << num_cur_nodes << std::endl;
			//

 			stop = type ? (radius <= m_sampleRadius) : (num_cur_nodes >= n_nodes);
 			if(num_cur_nodes>=m_meshVertexNum)
                break;
			
		}

		std::cout << "\nfarthest point sample: m_sample Radius = " << radius
			<< " num of nodes = " << num_cur_nodes << "\n" << std::endl;
	}

    //------------------------------------------------------------------------
    //	Node Sampling based on geodesic distance metric
    //
    //	Note that this member function samples nodes along some axis.
    //	Each node is not covered by any other node. And distance between each
    //	pair of nodes is at least sampling radius.
    //------------------------------------------------------------------------
    void nodeSampler::mysample(Mesh &mesh, igl::HeatGeodesicsData<double>& data, int n_nodes, int max_neighbors)
    {
        //	Save numbers of vertex and edge
        m_meshVertexNum = mesh.n_vertices();
        m_meshEdgeNum = mesh.n_edges();
        size_t m_sampleVertexNum = n_nodes;
        double lambda0 = 1;

        // Calculate Gaussian curvature
        Eigen::VectorXd all_curvature(m_meshVertexNum);
        std::vector<bool> is_node(m_meshVertexNum, false);
        Eigen::VectorXd all_min_dist2node(m_meshVertexNum);
        all_min_dist2node.setOnes(); all_min_dist2node *= 1e12;
        Eigen::MatrixXd all_geodist(m_meshVertexNum, m_sampleVertexNum);

        // Calculate Geodesic Distance
        std::vector<size_t> sample_nodes;

        // initial nodes
        size_t max_curva_idx;
        calcCurvature(mesh, all_curvature);
        all_curvature.maxCoeff(&max_curva_idx);
        sample_nodes.push_back(max_curva_idx);
        is_node[max_curva_idx] = true;

        // begin push nodes
        size_t cur_node_idx = max_curva_idx;
        size_t cur_num_nodes = sample_nodes.size() - 1;
        size_t prev_cur_node;
        std::unique_ptr<geodesic> pGeodesic = std::make_unique<geodesic>(mesh);
        Eigen::VectorXd geoDistVector(m_meshVertexNum);

        while (cur_num_nodes < m_sampleVertexNum)
        {
            geoDistVector.setZero();
            geoDistVector.setZero();
            pGeodesic->seed(cur_node_idx, geoDistVector);
            all_geodist.col(cur_num_nodes) = geoDistVector;
            // begin test
            double min_geodist = geoDistVector.minCoeff();
            if (min_geodist < 0)
            {
                std::cout << "min_geodist = " << min_geodist << "  cur_idx = " << cur_node_idx <</* " min_idx = " << min_idx <<*/ std::endl;
                for(int i = 0; i < m_meshVertexNum; i++)
                {
                    if(geoDistVector[i] < 0)
                        geoDistVector[i] = 0.0;
                }
            }
            all_geodist.col(cur_num_nodes) = geoDistVector;
            sample_nodes.push_back(cur_node_idx);
            is_node[cur_node_idx] = true;
            m_geoDistContainer.push_back(geoDistVector);
            m_nodeContainer.emplace_back(m_geoDistContainer.size() - 1, cur_node_idx);

            // update all_min_dist2node
            for (int i = 0; i < m_meshVertexNum; i++)
            {
                if (!is_node[i])
                {
                    all_min_dist2node[i] = all_min_dist2node[i] > geoDistVector[i] ? geoDistVector[i] : all_min_dist2node[i];
                }
                else
                {
                    all_min_dist2node[i] = 0.0;
                    all_curvature[i] = 0.0;
                }
            }

            double max_geo_dist = all_min_dist2node.maxCoeff();
            m_sampleRadius = max_geo_dist;
            double lambda = lambda0 * max_geo_dist / (2 * M_PI);
            Eigen::VectorXd temp = all_min_dist2node.cwiseProduct(lambda * all_curvature);

            prev_cur_node = cur_node_idx;
            temp.maxCoeff(&cur_node_idx);
            if (prev_cur_node == cur_node_idx || max_geo_dist <= 0)
                break;

            cur_num_nodes++;
        }
        non_unisamples_Radius.resize(m_meshVertexNum);
        int n_neighbor = max_neighbors;
        for (int i = 0; i < m_meshVertexNum; i++)
        {
            bool flag = true;
            // sort all_geodist(i,:) and selected n_neibor-th node, save the length
            for (int k = 0; k < n_neighbor; k++)
            {
                for (int j = k; j < m_sampleVertexNum; j++)
                {
                    if (all_geodist(i, j) < all_geodist(i, k))
                    {
                        double t = all_geodist(i, j);
                        all_geodist(i, j) = all_geodist(i, k);
                        all_geodist(i, k) = t;
                    }
                }
                if (all_geodist(i, k) > m_sampleRadius)
                {
                    non_unisamples_Radius[i] = (m_sampleRadius + 1e-3);
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                non_unisamples_Radius[i] = (all_geodist(i, n_neighbor - 1) + 1e-3);

            }
            if(non_unisamples_Radius[i]<=0)
            {
                non_unisamples_Radius[i] = 0.01;
                std::cout << "sample radius [ " << i << " ] = " << all_geodist.block(i, 0, 1, n_neighbor) << std::endl;
            }

        }
    }

    void nodeSampler::calcCurvature(Mesh & mesh, Eigen::VectorXd& curvas)
    {
        curvas.resize(mesh.n_vertices());
        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            std::vector<Point> pts;
            for (auto vv_it = mesh.vv_begin(v_it); vv_it != mesh.vv_end(v_it); vv_it++)
            {
                pts.push_back(mesh.point(*vv_it));
            }
            Point p0 = mesh.point(*v_it);
            double area = 0.0;
            double angle = 0.0;
            for (int i = 0; i < pts.size(); i++)
            {
                Point p1 = pts[i];
                Point p2 = pts[(i + 1) % pts.size()];
                double temp = (p1 - p0) | (p2 - p0);
                double ang = std::acos(temp / (p1 - p0).norm()/(p2 - p0).norm());
                angle += ang;

                Eigen::Matrix3d A;
                A(0, 0) = p0[1] * (p1[2] - p2[2]) - p0[2] * (p1[1] - p2[1]) + p1[1] * p2[2] - p1[2] * p2[1];
                A(0, 1) = -(p0[0] * (p1[2] - p2[2]) - p0[2] * (p1[0] - p2[0]) + p1[0] * p2[2] - p1[2] * p2[0]);
                A(0, 2) = p0[0] * (p1[1] - p2[1]) - p0[1] * (p1[0] - p2[0]) + p1[0] * p2[1] - p1[1] * p2[0];
                A(1, 0) = 2 * (p0[0] - p1[0]);
                A(1, 1) = 2 * (p0[1] - p1[1]);
                A(1, 2) = 2 * (p0[2] - p1[2]);
                A(2, 0) = 2 * (p0[0] - p2[0]);
                A(2, 1) = 2 * (p0[1] - p2[1]);
                A(2, 2) = 2 * (p0[2] - p2[2]);

                Eigen::Vector3d b;
                b[0] = p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0]);
                b[1] = pow(p0.norm(), 2) - pow(p1.norm(), 2);
                b[2] = pow(p0.norm(), 2) - pow(p2.norm(), 2);

                Eigen::Vector3d center = A.inverse() * b;
                Point p3(center[0], center[1], center[2]);

                area += ((p3 - p0) % (p1 - p0)).norm() / 4 + ((p3 - p0) % (p2 - p0)).norm() / 4;
            }
            // calculate area
            curvas[(*v_it).idx()] = (2 * M_PI - angle);
        }
    }

    //------------------------------------------------------------------------
    //	Define updateWeight member function
    //	This member function recalculates vertex-node and node-node weights.
    //	This member function does not change node topology.
    //	This member function may change node topology if mesh is no longer isometric.
    //------------------------------------------------------------------------
    void nodeSampler::updateWeight(Mesh &mesh)
    {
        //	Recalculate geodesic distance
        m_geoDistContainer.clear();

        std::unique_ptr<geodesic> pGeodesic = std::make_unique<geodesic>(mesh);
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
        {
            size_t vIdx = m_nodeContainer.at(nodeIdx).second;
            Eigen::VectorXd geoDistVector(mesh.n_vertices());
            pGeodesic->seed(vIdx, geoDistVector);
            m_geoDistContainer.push_back(geoDistVector);
        }

#if 1
        //	Update vertex-node weights
        for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
        {
            double totalWeight = 0.0f;
            for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
            {
                size_t nodeIdx = eachNeighbor.first;
                double dist = m_geoDistContainer.at(nodeIdx)(vertexIdx);
                double weight = max(0.0f, cube(1.0f - square(dist / (2.0f * m_sampleRadius))));
                m_vertexGraph.at(vertexIdx).at(nodeIdx) = weight;

                totalWeight += weight;
            }

            //	Normalize vertex-node weights
            for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
            {
                size_t nodeIdx = eachNeighbor.first;
                m_vertexGraph.at(vertexIdx).at(nodeIdx) /= totalWeight;
            }
        }

        //	Update node-node weights
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
        {
            double totalWeight = 0.0;
            for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
            {
                size_t neighborIdx = eachNeighbor.first;
                size_t neighborVertexIdx = m_nodeContainer.at(neighborIdx).second;
                double dist = m_geoDistContainer.at(nodeIdx)(neighborVertexIdx);
                double weight = max(0.0, cube(1.0 - square(dist / (2 * m_sampleRadius))));
                m_nodeGraph.at(nodeIdx).at(neighborIdx) = weight;

                totalWeight += weight;
            }

            //	Normalize node-node weights
            for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
            {
                size_t neighborIdx = eachNeighbor.first;
                m_nodeGraph.at(nodeIdx).at(neighborIdx) /= totalWeight;
            }
        }
#else
        //	Update node graph topology and weight based on reference mesh
        m_vertexGraph.clear();
        m_nodeGraph.clear();
        constructGraph();
#endif
    }

    //------------------------------------------------------------------------
    //	Construct vertex graph and node graph and calculate vertex-node weights
    //	and node-node weights
    //------------------------------------------------------------------------
    void nodeSampler::constructGraph(bool is_uniform)
    {
        //	Construct vertex graph and calculate vertex weights
        //	Vertex weights formula: w(vi, xj) = max(0, (1-(d(vi, xj)/r)^2)^3)
        //	r = 2*m_sampleRadius
        m_vertexGraph.resize(m_meshVertexNum);
        Eigen::VectorXd snum(m_meshVertexNum);
        for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
        {
            double totalWeight = 0.0f;
            int num = 0;
            for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
            {
                double dist = m_geoDistContainer.at(nodeIdx)(vertexIdx);
                double weight;
                if (is_uniform)
                {
                    weight = max(0.0f, cube(1.0f - square(dist / (2.0f * m_sampleRadius))));
                }
                else
                {
                    weight = max(0.0f, cube(1.0f - square(dist / (2.0f *non_unisamples_Radius[vertexIdx]))));
                }


                if (weight > 0.0f)
                {
                    m_vertexGraph.at(vertexIdx).emplace(nodeIdx, weight);
                    totalWeight += weight;
                }
            }

            // std::cout << "m_sampleRadius = " << m_sampleRadius << std::endl;
            //	Normalize vertex-node weights
            for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
            {
                num++;
                size_t neighborIdx = eachNeighbor.first;
                m_vertexGraph.at(vertexIdx).at(neighborIdx) /= totalWeight;
            }
            snum[vertexIdx] = num;
        }

        //	Construct node topology
        std::vector<std::set<size_t>> nodeTopo(m_nodeContainer.size());
        for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
        {
            for (auto &eachNeighbor0 : m_vertexGraph[vertexIdx])
            {
                size_t nodeIdx0 = eachNeighbor0.first;
                for (auto &eachNeighbor1 : m_vertexGraph[vertexIdx])
                {
                    size_t nodeIdx1 = eachNeighbor1.first;
                    if (nodeIdx0 != nodeIdx1)
                    {
                        nodeTopo[nodeIdx0].insert(nodeIdx1);
                    }
                }
            }
        }

        //	Calculate node weights
        //	Node weights formula: w(xi, xj) = max(0, (1-(d(xi, xj)/r)^2)^3)
        //	r = 2*m_sampleRadius
        m_nodeGraph.resize(m_nodeContainer.size());
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
        {
            double totalWeight = 0.0f;
            size_t vertexIdx = getNodeVertexIdx(nodeIdx);
            for (auto &eachNeighbor : nodeTopo[nodeIdx])
            {
                size_t neighborNodeVertexIdx = m_nodeContainer.at(eachNeighbor).second;
                double dist = m_geoDistContainer.at(nodeIdx)(neighborNodeVertexIdx);
                double weight;
                if (is_uniform)
                {
                    weight = max(0.0f, cube(1.0f - square(dist / (2.0f * m_sampleRadius))));
                }
                else
                {
                    weight = max(0.0f, cube(1.0f - square(dist / (2.0f * non_unisamples_Radius[vertexIdx]))));
                }
                if (weight > 0.0f)
                {
                    m_nodeGraph.at(nodeIdx).emplace(eachNeighbor, weight);
                    totalWeight += weight;
                }
            }

            //	Normalize node-node weights
            for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
            {
                size_t neighborIdx = eachNeighbor.first;
                m_nodeGraph.at(nodeIdx).at(neighborIdx) /= totalWeight;
            }
        }
    }


    void nodeSampler::initWeight(Eigen::SparseMatrix<double>& matPV, Eigen::MatrixXd & matP,
        Eigen::SparseMatrix<double>& matB, Eigen::MatrixXd& matD, Eigen::VectorXd& smoothw,
        const Mesh& mesh, bool all_vertices_smooth)
    {
        std::vector<Eigen::Triplet<double>> coeff;
        matP.setZero();
        // data coeff
        for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
        {
            //int test_vn = 0;
            for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
            {
                //test_vn++;
                size_t nodeIdx = eachNeighbor.first;
                double weight = m_vertexGraph.at(vertexIdx).at(nodeIdx);
                Point vi = mesh.point(mesh.vertex_handle(vertexIdx));
                Point pj = mesh.point(mesh.vertex_handle(nodeIdx));
                coeff.push_back(Eigen::Triplet<double>(vertexIdx, nodeIdx * 4, weight * (vi - pj)[0]));
                coeff.push_back(Eigen::Triplet<double>(vertexIdx, nodeIdx * 4 + 1, weight * (vi - pj)[1]));
                coeff.push_back(Eigen::Triplet<double>(vertexIdx, nodeIdx * 4 + 2, weight * (vi - pj)[2]));
                coeff.push_back(Eigen::Triplet<double>(vertexIdx, nodeIdx * 4 + 3, weight * 1.0));
                matP(vertexIdx, 0) += weight * pj[0];
                matP(vertexIdx, 1) += weight * pj[1];
                matP(vertexIdx, 2) += weight * pj[2];
            }
            //std::cout << "vertex_idx = " << vertexIdx << " test_vn = " << test_vn << std::endl;
        }
        matPV.setFromTriplets(coeff.begin(), coeff.end());

        // smooth coeff
        coeff.clear();
        matD.setZero();
        if (all_vertices_smooth)
        {
            for (size_t heIdx = 0; heIdx < mesh.n_halfedges(); ++heIdx)
            {
                size_t viIdx = mesh.from_vertex_handle(mesh.halfedge_handle(heIdx)).idx();
                size_t vjIdx = mesh.to_vertex_handle(mesh.halfedge_handle(heIdx)).idx();
                Point vi = mesh.point(mesh.vertex_handle(viIdx));
                double test = 0.0;
                for (auto &eachNeighbor : m_vertexGraph[vjIdx])
                {
                    size_t nodeIdx = eachNeighbor.first;
                    double weight = m_vertexGraph.at(vjIdx).at(nodeIdx);
                    Point pk = mesh.point(mesh.vertex_handle(nodeIdx));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4, weight * (vi - pk)[0]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 1, weight * (vi - pk)[1]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 2, weight * (vi - pk)[2]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 3, weight * 1.0));
                    matD(heIdx, 0) -= weight * pk[0];
                    matD(heIdx, 1) -= weight * pk[1];
                    matD(heIdx, 2) -= weight * pk[2];
                    test += weight;
                }
                for (auto &eachNeighbor : m_vertexGraph[viIdx])
                {
                    size_t nodeIdx = eachNeighbor.first;
                    double weight = m_vertexGraph.at(viIdx).at(nodeIdx);
                    Point pk = mesh.point(mesh.vertex_handle(nodeIdx));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4, -weight * (vi - pk)[0]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 1, -weight * (vi - pk)[1]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 2, -weight * (vi - pk)[2]));
                    coeff.push_back(Eigen::Triplet<double>(heIdx, nodeIdx * 4 + 3, -weight * 1.0));
                    matD(heIdx, 0) += weight * pk[0];
                    matD(heIdx, 1) += weight * pk[1];
                    matD(heIdx, 2) += weight * pk[2];
                }
            }
            smoothw.setOnes();
        }
        else
        {
            smoothw.setZero();
            //std::cout << "m_node container size = " << m_nodeContainer.size() << std::endl;
            for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
            {
                size_t vIdx0 = getNodeVertexIdx(nodeIdx);
                VertexHandle vh0 = mesh.vertex_handle(vIdx0);
                Point v0 = mesh.point(vh0);
                //std::cout << "node idx = " << nodeIdx << " vidx = " << vIdx0 << " neighbor = " << m_nodeGraph.size() << std::endl;
                //int test_n = 0;
                for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
                {
                    //test_n++;
                    size_t neighborIdx = eachNeighbor.first;
                    size_t vIdx1 = getNodeVertexIdx(neighborIdx);
                    Point v1 = mesh.point(mesh.vertex_handle(vIdx1));
                    Point dv = v0 - v1;
                    int k = nodeIdx * (m_nodeContainer.size() - 1) + neighborIdx;

                    coeff.push_back(Eigen::Triplet<double>(k, neighborIdx * 4, dv[0]));
                    coeff.push_back(Eigen::Triplet<double>(k, neighborIdx * 4 + 1, dv[1]));
                    coeff.push_back(Eigen::Triplet<double>(k, neighborIdx * 4 + 2, dv[2]));
                    coeff.push_back(Eigen::Triplet<double>(k, neighborIdx * 4 + 3, 1.0));
                    coeff.push_back(Eigen::Triplet<double>(k, nodeIdx * 4 + 3, -1.0));

                    //smoothw[k] = m_nodeGraph.at(nodeIdx).at(neighborIdx);
                    smoothw[k] = 1;
                    matD(k, 0) = dv[0];
                    matD(k, 1) = dv[1];
                    matD(k, 2) = dv[2];
                }
                //std::cout << "test_n = " << test_n << "\n\n" << std::endl;
            }
        }

        matB.setFromTriplets(coeff.begin(), coeff.end());
    }

    void nodeSampler::print_nodes(Mesh & mesh, std::string file_path)
    {
        std::string namev = file_path + "nodes.obj";
		std::ofstream out1(namev);
		//std::cout << "print nodes to " << name << std::endl;
		for (int i = 0; i < m_nodeContainer.size(); i++)
		{
			int vexid = m_nodeContainer[i].second;
			out1 << "v " << mesh.point(mesh.vertex_handle(vexid))[0] << " " << mesh.point(mesh.vertex_handle(vexid))[1]
				<< " " << mesh.point(mesh.vertex_handle(vexid))[2] << std::endl;
		}
		out1.close();
		std::string namee = file_path + "edges.txt";
		std::ofstream out2(namee);
		//std::cout << "print nodes to " << name << std::endl;
		for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
		{
			size_t vIdx0 = getNodeVertexIdx(nodeIdx);
			for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
			{
				size_t neighborIdx = eachNeighbor.first;
				size_t vIdx1 = getNodeVertexIdx(neighborIdx);
				out2 << nodeIdx << " " << neighborIdx << " " << vIdx0 << " " << vIdx1 << std::endl;
			}
		}
		out2.close();

    }

    double nodeSampler::get_graph_weights(Mesh & mesh)
    {
        double sumweights = 0;
        int num = 0;
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
		{
            size_t vIdx0 = getNodeVertexIdx(nodeIdx);
            VertexHandle vh0 = mesh.vertex_handle(vIdx0);
            Point v0 = mesh.point(vh0);
			for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
			{
				size_t neighborIdx = eachNeighbor.first;
				size_t vIdx1 = getNodeVertexIdx(neighborIdx);
                Point v1 = mesh.point(mesh.vertex_handle(vIdx1));
                Point dv = v0 - v1;

				sumweights += dv.norm();
                num++;
			}
		}
    
        return sumweights;
    }
}
//	//------------------------------------------------------------------------
//	//	Define writeToCacheFile member function
//	//	Write all member data to file
//	//------------------------------------------------------------------------
//	void nodeSampler::writeToCacheFile(std::ofstream &cacheStream) const
//	{
//		//	Write control parameters
//		cacheStream << "m_meshVertexNum = " << m_meshVertexNum << std::endl;
//		cacheStream << "m_meshEdgeNum = " << m_meshEdgeNum << std::endl;
//		cacheStream << "m_averageEdgeLen = " << m_averageEdgeLen << std::endl;
//		cacheStream << "m_sampleRadius = " << m_sampleRadius << std::endl;

//		//	Write all containers
//		cacheStream << "m_nodeContainer" << std::endl;
//		for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
//		{
//			cacheStream << m_nodeContainer[nodeIdx].first << " ";
//			cacheStream << m_nodeContainer[nodeIdx].second << std::endl;
//		}
//		cacheStream << "D" << std::endl;

//		cacheStream << "m_vertexGraph" << std::endl;
//		for (size_t vertexIdx = 0; vertexIdx < m_vertexGraph.size(); ++vertexIdx)
//		{
//			for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
//			{
//				cacheStream << eachNeighbor.first << " ";
//				cacheStream << eachNeighbor.second << std::endl;
//			}
//			cacheStream << "D" << std::endl;
//		}

//		cacheStream << "m_nodeGraph" << std::endl;
//		for (size_t nodeIdx = 0; nodeIdx < m_nodeGraph.size(); ++nodeIdx)
//		{
//			for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
//			{
//				cacheStream << eachNeighbor.first << " ";
//				cacheStream << eachNeighbor.second << std::endl;
//			}
//			cacheStream << "D" << std::endl;
//		}
//	}

//	//------------------------------------------------------------------------
//	//	Define loadFromCacheFile member function
//	//	Load all control parameters and container data from cache file
//	//	And calculate geodesic distance based on reference mesh
//	//------------------------------------------------------------------------
//    void nodeSampler::loadFromCacheFile(std::ifstream &cacheStream, Mesh &mesh)
//	{
//		//	Clear previous control parameters
//		m_meshVertexNum = 0;
//		m_meshEdgeNum = 0;
//		m_averageEdgeLen = 0.0f;
//		m_sampleRadius = 0.0f;

//		//	Load control parameters
//		std::string lineStr;
//		std::getline(cacheStream, lineStr);
//		size_t meshVertexNumPos = lineStr.find("m_meshVertexNum");
//		if (meshVertexNumPos != std::string::npos)
//		{
//			m_meshVertexNum = std::stoull(lineStr.substr(lineStr.find("=") + 2));
//		}

//		std::getline(cacheStream, lineStr);
//		size_t meshEdgeNumPos = lineStr.find("m_meshEdgeNum");
//		if (meshEdgeNumPos != std::string::npos)
//		{
//			m_meshEdgeNum = std::stoull(lineStr.substr(lineStr.find("=") + 2));
//		}

//		std::getline(cacheStream, lineStr);
//		size_t averageEdgeLenPos = lineStr.find("m_averageEdgeLen");
//		if (averageEdgeLenPos != std::string::npos)
//		{
//			m_averageEdgeLen = std::stof(lineStr.substr(lineStr.find("=") + 2));
//		}

//		std::getline(cacheStream, lineStr);
//		size_t sampleRadiusPos = lineStr.find("m_sampleRadius");
//		if (sampleRadiusPos != std::string::npos)
//		{
//			m_sampleRadius = std::stof(lineStr.substr(lineStr.find("=") + 2));
//		}

//		//	Clear all containers
//		m_nodeContainer.clear();
//		m_vertexGraph.clear();
//		m_nodeGraph.clear();

//		//	Load container data
//		std::getline(cacheStream, lineStr);
//		if (lineStr == "m_nodeContainer")
//		{
//			while (std::getline(cacheStream, lineStr) && lineStr != "D")
//			{
//				std::istringstream parseLine(lineStr);
//				std::pair<size_t, size_t> nodePair;
//				parseLine >> nodePair.first >> nodePair.second;
//				m_nodeContainer.push_back(nodePair);
//			}
//		}

//		std::getline(cacheStream, lineStr);
//		if (lineStr == "m_vertexGraph")
//		{
//			for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
//			{
//				std::map<size_t, double> vertexNeighbors;
//				while (std::getline(cacheStream, lineStr) && lineStr != "D")
//				{
//					std::istringstream parseLine(lineStr);
//					std::pair<size_t, double> nodePair;
//					parseLine >> nodePair.first >> nodePair.second;
//					vertexNeighbors.emplace(nodePair);
//				}
//				m_vertexGraph.push_back(vertexNeighbors);
//			}
//		}

//		std::getline(cacheStream, lineStr);
//		if (lineStr == "m_nodeGraph")
//		{
//			for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
//			{
//				std::map<size_t, double> nodeNeighbors;
//				while (std::getline(cacheStream, lineStr) && lineStr != "D")
//				{
//					std::istringstream parseLine(lineStr);
//					std::pair<size_t, double> nodePair;
//					parseLine >> nodePair.first >> nodePair.second;
//					nodeNeighbors.emplace(nodePair);
//				}
//				m_nodeGraph.push_back(nodeNeighbors);
//			}
//		}

//		//	Calculate geodesic distance based on reference mesh
//		updateWeight(mesh);
//	}
//}
