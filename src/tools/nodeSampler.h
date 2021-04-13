#ifndef NODESAMPLER_H
#define NODESAMPLER_H
#include "Types.h"

namespace svr
{
    //------------------------------------------------------------------------
    //	Define neighborIter class
    //------------------------------------------------------------------------
    class neighborIter
    {
    public:
        neighborIter(const std::map<size_t, Scalar> &nodeNeighbors)
        {
            m_neighborIter = nodeNeighbors.begin();
            m_neighborEnd = nodeNeighbors.end();
        }

        neighborIter& operator++()
        {
            if (m_neighborIter != m_neighborEnd)
                ++m_neighborIter;
            return *this;
        }

        neighborIter operator++(int)
        {
            neighborIter tempIter(*this);
            ++*this;
            return tempIter;
        }

        const std::pair<const size_t, Scalar>& operator*() { return *m_neighborIter; }
        std::map<size_t, Scalar>::const_iterator operator->() { return m_neighborIter; }
        bool is_valid() { return m_neighborIter != m_neighborEnd; }
        size_t getIndex() { return m_neighborIter->first; }
        Scalar getWeight() { return m_neighborIter->second; }

    private:
        std::map<size_t, Scalar>::const_iterator m_neighborIter;
        std::map<size_t, Scalar>::const_iterator m_neighborEnd;
    };

    //------------------------------------------------------------------------
    //	Define node sampling class
    //------------------------------------------------------------------------
    class nodeSampler
    {
    public:
        enum sampleAxis { X_AXIS, Y_AXIS, Z_AXIS };
        nodeSampler() {};

        // return sample radius
//        Scalar sample(Mesh &mesh, Scalar sampleRadiusRatio, sampleAxis axis);
        Scalar sampleAndconstuct(Mesh &mesh, Scalar sampleRadiusRatio, sampleAxis axis);

        void updateWeight(Mesh &mesh);
        void constructGraph(bool is_uniform);

        size_t nodeSize() const { return m_nodeContainer.size(); }

        neighborIter getNodeNodeIter(size_t nodeIdx) const { return neighborIter(m_nodeGraph[nodeIdx]); }
        neighborIter getVertexNodeIter(size_t vertexIdx) const { return neighborIter(m_vertexGraph[vertexIdx]); }
        size_t getNodeVertexIdx(size_t nodeIdx) const { return m_nodeContainer.at(nodeIdx).second; }

        size_t getVertexNeighborSize(size_t vertexIdx) const { return m_vertexGraph.at(vertexIdx).size(); }
        size_t getNodeNeighborSize(size_t nodeIdx) const { return m_nodeGraph.at(nodeIdx).size(); }
        void initWeight(Eigen::SparseMatrix<Scalar>& matPV, MatrixXX & matP,
                        Eigen::SparseMatrix<Scalar>& matB, MatrixXX& matD, VectorX& smoothw);
        void print_nodes(Mesh & mesh, std::string file_path);

    private:
        size_t m_meshVertexNum = 0;
        size_t m_meshEdgeNum = 0;
        Scalar m_averageEdgeLen = 0.0f;
        Scalar m_sampleRadius = 0.0f;
        std::vector<Scalar> non_unisamples_Radius;
        Mesh * m_mesh;

        std::vector<std::pair<size_t, size_t>> m_nodeContainer;
        std::vector<std::map<size_t, Scalar>> m_vertexGraph; // vertex node graph
        std::vector<std::map<size_t, Scalar>> m_nodeGraph;
        std::vector<std::map<size_t, Scalar>> m_nodeVertexGraph;
        std::vector<VectorX> m_geoDistContainer;
        Eigen::VectorXi      VertexNodeIdx;


    };
}
#endif
