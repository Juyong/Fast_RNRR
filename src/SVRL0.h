#pragma once
#include "nodeSampler.h"
#include "Registration.h"


struct gnParams
{
    gnParams(double alphaRigid = 1000.0f, double alphaSmooth = 1000.0f,
        double alphaPoint = 0.1f, double alphaPlane = 1.0f) :
        m_alphaRigid(alphaRigid),
        m_alphaSmooth(alphaSmooth),
        m_alphaPoint(alphaPoint),
        m_alphaPlane(alphaPlane) {}

    double m_alphaRigid = 0.0f;									///<  Rigid control parameter
    double m_alphaSmooth = 0.0f;									///<  Smooth control parameter
    double m_alphaPoint = 0.0f;									///<  Point fitting control parameter
    double m_alphaPlane = 0.0f;									///<  Plane fitting control parameter

    void relaxParams(void)
    {
        m_alphaRigid /= 2.0f;
        m_alphaSmooth /= 2.0f;
    }
};

class SVRL0 : public Registration
{
public:
    SVRL0() {};
    ~SVRL0() {};
    virtual void Initialize();
    virtual double DoNonRigid();
    double gaussNewton(const Mesh &currentNonRigidMesh,
    const Mesh &currentScanMesh,
    const std::vector<std::pair<size_t, size_t>> &vertexPair,
    const svr::nodeSampler& sampler,
    Eigen::VectorXd &affineVector,
    const gnParams &controlParams);
    Eigen::VectorXd initAffineVector(size_t nodeSize) const;
    /*Eigen::VectorXd innerICP(const Mesh &currentNonRigidMesh,
        const Mesh &currentScanMesh,
        const svr::nodeSampler &sampler,
        const gnParams &controlParams);*/
    void accumulateAffine(const Eigen::VectorXd &affineVector,
        Eigen::VectorXd &accumulatedAffineVector,
        size_t nodeSize) const;
    double updateNonRigidMesh(Mesh &currentNonRigidMesh,
        const Eigen::VectorXd &affineVector, const svr::nodeSampler &sampler, double& gt_mean_err, double& gt_max_err);


    void initTempK(std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
        const svr::nodeSampler &sampler) const;

    ///  @brief		Solve sub-problem 1 of L0 optimization
    ///				SubProblem1 returns pruned tempK based on threshold lambda/beta
    ///  @return	Pruned tempK
    ///  @note		Input affineVector should be the return from non-rigid icp
    void subProblem1(const MatrixXX &corresMesh,
        std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
        const Eigen::VectorXd &affineVector, double beta, double lambda,
        const svr::nodeSampler &anchorSampler);

    ///  @brief		Solve sub-problem 2 of L0 optimization
    ///  @return	Sparse affineVector_L0
    double subProblem2(const Mesh &anchorMesh,
        const std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
        Eigen::VectorXd &affineVector_L0,
        const Eigen::VectorXd &accumAffineVector, double beta,
        const svr::nodeSampler &anchorSamplerr);

private:
    size_t m_templateVertexNum = 0;

    //-------------------------------------
    //	Define control parameters
    //-------------------------------------
    double m_cosLimit = 0.0f;									///<  Maximum tolerated theta between vertex pair
    double m_distLimit = 0.0f;									///<  Maximum tolerated distance between vertex pair
    double m_sampleRadiusCtrl = 0.0f;							///<  Node sampling control parameter
    double m_lambdaInit = 0.0f;									///<  Lambda initial parameter
    double m_nonRigidThreshold = 0.0f;							///<  A node has significant non-rigidity once its non-rigidity is above this threshold
    double m_nonRigidRatioThreshold = 0.0f;						///<  Execute L0 when ratio of significant non-rigid nodes is above this threshold
    static const double c_betaMax;								///<  Maximum beta controls L0 optimization
    gnParams m_gnParams;

    std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> tempK;
    svr::nodeSampler src_sample_nodes;
    int num_sample_nodes;

};
