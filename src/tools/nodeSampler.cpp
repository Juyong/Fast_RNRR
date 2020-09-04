//#pragma once
#include "nodeSampler.h"
#include "tools.h"

#include <igl/heat_geodesics.h>
#include <igl/read_triangle_mesh.h>
#include "geodesic/geodesic_algorithm_exact.h"

namespace svr
{
    //	Define helper functions
    static auto square = [](const Scalar argu) { return argu * argu; };
    static auto cube = [](const Scalar argu) { return argu * argu * argu; };
    static auto max = [](const Scalar lhs, const Scalar rhs) { return lhs > rhs ? lhs : rhs; };

    //------------------------------------------------------------------------
    //	Node Sampling based on geodesic distance metric
    //
    //	Note that this member function samples nodes along some axis.
    //	Each node is not covered by any other node. And distance between each
    //	pair of nodes is at least sampling radius.
    //------------------------------------------------------------------------
    // heat method
   Scalar nodeSampler::sample(Mesh &mesh, Scalar sampleRadiusRatio, sampleAxis axis)
    {
        //	Save numbers of vertex and edge
        m_meshVertexNum = mesh.n_vertices();
        m_meshEdgeNum = mesh.n_edges();
        m_mesh = & mesh;

        //	Calculate average edge length of bound mesh
        for (size_t i = 0; i < m_meshEdgeNum; ++i)
        {
            OpenMesh::EdgeHandle eh = mesh.edge_handle(i);
            Scalar edgeLen = mesh.calc_edge_length(eh);
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
        VectorX geoDistVector(m_meshVertexNum);
        geoDistVector.setZero();

        igl::HeatGeodesicsData<Scalar> data;
        MatrixXX V1;
        Eigen::MatrixXi F1;
        data.use_intrinsic_delaunay = true;
        Mesh2VF(mesh, V1, F1);
        igl::heat_geodesics_precompute(V1, F1, data);

        igl::heat_geodesics_solve(data, (Eigen::MatrixXi(1,1)<<firstVertexIdx).finished(), geoDistVector);
        m_geoDistContainer.push_back(geoDistVector/2.0);
        m_nodeContainer.emplace_back(0, firstVertexIdx);
        VertexNodeIdx.resize(m_meshVertexNum, -1);

        for (auto &vertexIdx : vertexReorderedAlongAxis)
        {
            bool IsNode = true;
            for (size_t k = 0; k < m_geoDistContainer.size(); ++k)
            {
                Scalar dist = m_geoDistContainer.at(k)(vertexIdx);
                if (dist < m_sampleRadius)
                {
                    IsNode = false;
                    break;
                }
            }
            if (IsNode)
            {
                geoDistVector.setZero();
                igl::heat_geodesics_solve(data, (Eigen::MatrixXi(1,1)<<vertexIdx).finished(), geoDistVector);
                m_geoDistContainer.push_back(geoDistVector/2.0);
                m_nodeContainer.emplace_back(m_geoDistContainer.size() - 1, vertexIdx);
                VertexNodeIdx[vertexIdx] = m_geoDistContainer.size()-1;
            }
        }

        std::cout << "m_radius = " << m_sampleRadius << std::endl;
        return m_sampleRadius;
    }


    // Local geodesic calculation
    Scalar nodeSampler::sampleAndconstuct(Mesh &mesh, Scalar sampleRadiusRatio, sampleAxis axis)
    {
        //	Save numbers of vertex and edge
        m_meshVertexNum = mesh.n_vertices();
        m_meshEdgeNum = mesh.n_edges();
        m_mesh = & mesh;

        //	Calculate average edge length of bound mesh
        for (size_t i = 0; i < m_meshEdgeNum; ++i)
        {
            OpenMesh::EdgeHandle eh = mesh.edge_handle(i);
            Scalar edgeLen = mesh.calc_edge_length(eh);
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
        VertexNodeIdx.resize(m_meshVertexNum);
        VertexNodeIdx.setConstant(-1);
        VertexNodeIdx[firstVertexIdx] = 0;
        size_t cur_node_idx = 0;

        m_vertexGraph.resize(m_meshVertexNum);
        VectorX weight_sum = VectorX::Zero(m_meshVertexNum);

        for (auto &vertexIdx : vertexReorderedAlongAxis)
        {
            if(VertexNodeIdx[vertexIdx] < 0 && m_vertexGraph.at(vertexIdx).empty())
            {
                m_nodeContainer.emplace_back(cur_node_idx, vertexIdx);
                VertexNodeIdx[vertexIdx] = cur_node_idx;

                std::vector<size_t> neighbor_verts;
                geodesic::GeodesicAlgorithmExact geoalg(&mesh, vertexIdx, m_sampleRadius);
                geoalg.propagate(vertexIdx, neighbor_verts);
                for(size_t i = 0; i < neighbor_verts.size(); i++)
                {
                    int neighIdx = neighbor_verts[i];
                    Scalar geodist = mesh.data(mesh.vertex_handle(neighIdx)).geodesic_distance;
                    if(geodist < m_sampleRadius)
                    {
                        Scalar weight = std::pow(1-std::pow(geodist/m_sampleRadius, 2), 3);
                        m_vertexGraph.at(neighIdx).emplace(std::pair<int, Scalar>(cur_node_idx, weight));
                        weight_sum[neighIdx] += weight;
                    }
                }
                cur_node_idx++;
            }
        }

        m_nodeGraph.resize(cur_node_idx);
        for (auto &vertexIdx : vertexReorderedAlongAxis)
        {
            for(auto &node: m_vertexGraph[vertexIdx])
            {
                size_t nodeIdx = node.first;
                for(auto &neighNode: m_vertexGraph[vertexIdx])
                {
                    size_t neighNodeIdx = neighNode.first;
                    if(nodeIdx != neighNodeIdx)
                    {
                        m_nodeGraph.at(nodeIdx).emplace(std::pair<int, Scalar>(neighNodeIdx, 1.0));
                    }
                }
                m_vertexGraph.at(vertexIdx).at(nodeIdx) /= weight_sum[vertexIdx];
            }
        }
        return m_sampleRadius;
    }



    void nodeSampler::initWeight(Eigen::SparseMatrix<Scalar>& matPV, MatrixXX & matP,
        Eigen::SparseMatrix<Scalar>& matB, MatrixXX& matD, VectorX& smoothw)
    {
        std::vector<Eigen::Triplet<Scalar>> coeff;
        matP.setZero();
        Eigen::VectorXi nonzero_num = Eigen::VectorXi::Zero(m_mesh->n_vertices());
        // data coeff
        for (size_t vertexIdx = 0; vertexIdx < m_meshVertexNum; ++vertexIdx)
        {
            Mesh::Point vi = m_mesh->point(m_mesh->vertex_handle(vertexIdx));
            for (auto &eachNeighbor : m_vertexGraph[vertexIdx])
            {
                size_t nodeIdx = eachNeighbor.first;
                Scalar weight = m_vertexGraph.at(vertexIdx).at(nodeIdx);
                Mesh::Point pj = m_mesh->point(m_mesh->vertex_handle(getNodeVertexIdx(nodeIdx)));
                coeff.push_back(Eigen::Triplet<Scalar>(vertexIdx, nodeIdx * 4, weight * (vi[0] - pj[0])));
                coeff.push_back(Eigen::Triplet<Scalar>(vertexIdx, nodeIdx * 4 + 1, weight * (vi[1] - pj[1])));
                coeff.push_back(Eigen::Triplet<Scalar>(vertexIdx, nodeIdx * 4 + 2, weight * (vi[2] - pj[2])));
                coeff.push_back(Eigen::Triplet<Scalar>(vertexIdx, nodeIdx * 4 + 3, weight * 1.0));
                matP(vertexIdx, 0) += weight * pj[0];
                matP(vertexIdx, 1) += weight * pj[1];
                matP(vertexIdx, 2) += weight * pj[2];
            }
            nonzero_num[vertexIdx] = m_vertexGraph[vertexIdx].size();
        }
        matPV.setFromTriplets(coeff.begin(), coeff.end());

        // smooth coeff
        coeff.clear();
        int max_edge_num = nodeSize() * (nodeSize()-1);
        matB.resize(max_edge_num, 4 * nodeSize());
        matD.resize(max_edge_num, 3);
        smoothw.resize(max_edge_num);
        matD.setZero();
        int edge_id = 0;
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
        {
            size_t vIdx0 = getNodeVertexIdx(nodeIdx);
            Mesh::VertexHandle vh0 = m_mesh->vertex_handle(vIdx0);
            Mesh::Point v0 = m_mesh->point(vh0);
            for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
            {
                size_t neighborIdx = eachNeighbor.first;
                size_t vIdx1 = getNodeVertexIdx(neighborIdx);
                Mesh::Point v1 = m_mesh->point(m_mesh->vertex_handle(vIdx1));
                Mesh::Point dv = v0 - v1;
                int k = edge_id;

                coeff.push_back(Eigen::Triplet<Scalar>(k, neighborIdx * 4, dv[0]));
                coeff.push_back(Eigen::Triplet<Scalar>(k, neighborIdx * 4 + 1, dv[1]));
                coeff.push_back(Eigen::Triplet<Scalar>(k, neighborIdx * 4 + 2, dv[2]));
                coeff.push_back(Eigen::Triplet<Scalar>(k, neighborIdx * 4 + 3, 1.0));
                coeff.push_back(Eigen::Triplet<Scalar>(k, nodeIdx * 4 + 3, -1.0));

                Scalar dist = dv.norm();
                if(dist > 0)
                {
                    smoothw[k] = 1.0/ dist;
                }
                else
                {
                    smoothw[k] = 0.0;
                    std::cout << "node repeat";
                    exit(1);
                }
                matD(k, 0) = dv[0];
                matD(k, 1) = dv[1];
                matD(k, 2) = dv[2];
                edge_id++;
            }
        }
        matB.setFromTriplets(coeff.begin(), coeff.end());
        matD.conservativeResize(edge_id, 3);
        matB.conservativeResize(edge_id, matPV.cols());
        smoothw.conservativeResize(edge_id);
        smoothw *= edge_id/(smoothw.sum());
    }


    void nodeSampler::print_nodes(Mesh & mesh, std::string file_path)
    {
        std::string namev = file_path + "nodes.obj";
        std::ofstream out1(namev);
        //std::cout << "print nodes to " << name << std::endl;
        for (size_t i = 0; i < m_nodeContainer.size(); i++)
        {
            int vexid = m_nodeContainer[i].second;
            out1 << "v " << mesh.point(mesh.vertex_handle(vexid))[0] << " " << mesh.point(mesh.vertex_handle(vexid))[1]
                << " " << mesh.point(mesh.vertex_handle(vexid))[2] << std::endl;
        }
        Eigen::VectorXi nonzero_num = Eigen::VectorXi::Zero(m_nodeContainer.size());
        for (size_t nodeIdx = 0; nodeIdx < m_nodeContainer.size(); ++nodeIdx)
        {
            for (auto &eachNeighbor : m_nodeGraph[nodeIdx])
            {
                size_t neighborIdx = eachNeighbor.first;
                out1 << "l " << nodeIdx+1 << " " << neighborIdx+1 << std::endl;
            }
            nonzero_num[nodeIdx] = m_nodeGraph[nodeIdx].size();
        }
        std::cout << "node neighbor min = " << nonzero_num.minCoeff() << " max = "
                  << nonzero_num.maxCoeff() << " average = " << nonzero_num.mean() << std::endl;
        out1.close();
        std::string namee = file_path + "edges.txt";
        std::ofstream out2(namee);
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
}
