#ifndef NODESAMPLER_H
#define NODESAMPLER_H
#include "MeshDefinition.h"
#include <Eigen/Eigen>
#include "geodesic.h"

namespace svr
{
    //------------------------------------------------------------------------
    //	Define neighborIter class
    //------------------------------------------------------------------------
    class neighborIter
    {
    public:
        neighborIter(const std::map<size_t, double> &nodeNeighbors)
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

        const std::pair<const size_t, double>& operator*() { return *m_neighborIter; }
        std::map<size_t, double>::const_iterator operator->() { return m_neighborIter; }
        bool is_valid() { return m_neighborIter != m_neighborEnd; }
        size_t getIndex() { return m_neighborIter->first; }
        double getWeight() { return m_neighborIter->second; }

    private:
        std::map<size_t, double>::const_iterator m_neighborIter;
        std::map<size_t, double>::const_iterator m_neighborEnd;
    };

    //------------------------------------------------------------------------
    //	Define node sampling class
    //------------------------------------------------------------------------
    class nodeSampler
    {
    public:
        enum sampleAxis { X_AXIS, Y_AXIS, Z_AXIS };
        nodeSampler() {};

        void sample(const Mesh &mesh,igl::HeatGeodesicsData<double>& data, double sampleRadiusRatio, sampleAxis axis);
        void mysample(Mesh &mesh, igl::HeatGeodesicsData<double>& data, int n_nodes, int max_neighbors);
        void farthest_point_sample(const Mesh &mesh, igl::HeatGeodesicsData<double>& data, double sampleRadiusRatio, int n_nodes, int type);

        void updateWeight(Mesh &mesh);
        void constructGraph(bool is_uniform);
        void calcCurvature(Mesh & mesh, Eigen::VectorXd& curvas);

        size_t nodeSize() const { return m_nodeContainer.size(); }

        neighborIter getNodeNodeIter(size_t nodeIdx) const { return neighborIter(m_nodeGraph[nodeIdx]); }
        neighborIter getVertexNodeIter(size_t vertexIdx) const { return neighborIter(m_vertexGraph[vertexIdx]); }
        size_t getNodeVertexIdx(size_t nodeIdx) const { return m_nodeContainer.at(nodeIdx).second; }

        size_t getVertexNeighborSize(size_t vertexIdx) const { return m_vertexGraph.at(vertexIdx).size(); }
        size_t getNodeNeighborSize(size_t nodeIdx) const { return m_nodeGraph.at(nodeIdx).size(); }

//		void drawNodes(const Mesh &mesh, const std::string &nodeFile) const;
//		void drawNodeGraph(const Mesh &mesh, const std::string &nodeGraphFile) const;

//		void drawGeoDist(const Mesh &mesh, const std::string &geoDistFolder) const;
//		void outputGeoDistInTxt(const Mesh &mesh, const std::string &geoDistTxtFolder) const;

//		void outputNodeInfo(const std::string &nodeInfoFile) const;
//		void outputGraphInfo(const std::string &graphInfoFile) const;

//		void writeToCacheFile(std::ofstream &cacheStream) const;
//        void loadFromCacheFile(std::ifstream &cacheStream, Mesh &mesh);

        void initWeight(Eigen::SparseMatrix<double>& matPV, Eigen::MatrixXd & matP, Eigen::SparseMatrix<double>& matB, Eigen::MatrixXd& matD, Eigen::VectorXd& smoothw, const Mesh& mesh, bool all_vertices_smooth = false);
        void print_nodes(Mesh & mesh, std::string file_path);

        double get_graph_weights(Mesh & mesh);

    private:
        size_t m_meshVertexNum = 0;
        size_t m_meshEdgeNum = 0;
        double m_averageEdgeLen = 0.0f;
        double m_sampleRadius = 0.0f;
        std::vector<double> non_unisamples_Radius;

        std::vector<std::pair<size_t, size_t>> m_nodeContainer;
        std::vector<std::map<size_t, double>> m_vertexGraph;
        std::vector<std::map<size_t, double>> m_nodeGraph;
        std::vector<Eigen::VectorXd> m_geoDistContainer;
    };
}
#endif
