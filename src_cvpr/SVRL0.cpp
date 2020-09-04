#include "SVRL0.h"
#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>

//	Define control parameters for rigid term and data term of subproblem 2
static double g_alphaRigid = 10.0;
static double g_alphaData = 1.0;

void SVRL0::Initialize()
{
    g_alphaData = pars_.alpha;
    g_alphaRigid = pars_.beta;
    nonrigid_init();
    //	Initialize temporary variable K
    tempK.resize(n_src_vertex_);

    // Put source mesh into tempK
    if(pars_.is_unifom_sample)
            src_sample_nodes.sample(*src_mesh_, hgeodesic_->data1, pars_.uni_sample_radio, svr::nodeSampler::X_AXIS);
        else
            src_sample_nodes.mysample(*src_mesh_, hgeodesic_->data1, n_src_vertex_ * pars_.n_uni_sample_radio, pars_.n_uni_sample_neighbors);

    src_sample_nodes.constructGraph(pars_.is_unifom_sample);
    if(pars_.print_each_step_info)
    {
        std::string out_obj = pars_.out_each_step_info + "init_node.obj";
        src_sample_nodes.print_nodes(*src_mesh_, out_obj);//init sample nodes
    }

    initTempK(tempK, src_sample_nodes);
    std::cout << "tempk.size = " << tempK.size() << std::endl;

    num_sample_nodes = src_sample_nodes.nodeSize();
    std::cout << "sample node size = " << num_sample_nodes << std::endl;

    mat_U0_.resize(3, n_src_vertex_);

    if(correspondence_pairs_.size() < n_src_vertex_)
    {
        VPairs temp_corres;
        FindClosestPoints(temp_corres);
        for(int i = 0; i < correspondence_pairs_.size(); i++)
        {
            temp_corres[correspondence_pairs_[i].first] =
                    std::pair<int,int>(correspondence_pairs_[i].first, correspondence_pairs_[i].second);
        }
        correspondence_pairs_ = temp_corres;
    }
}

//double SVRL0::DoNonRigid()
//{
//    // paras
//    double beta = 1;
//    double lambda = 0.02;

//    //	Initialize temporary variable K
//    std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> tempK(n_src_vertex_);

//    // Put source mesh into tempK
//    svr::nodeSampler src_sample_nodes;

//    if(pars_.is_unifom_sample)
//            src_sample_nodes.sample(*src_mesh_, pgeodesic_->data1, pars_.uni_sample_radio, svr::nodeSampler::X_AXIS);
//        else
//            src_sample_nodes.mysample(*src_mesh_, pgeodesic_->data1, n_src_vertex_ * pars_.n_uni_sample_radio, pars_.n_uni_sample_neighbors);
//    src_sample_nodes.constructGraph(pars_.is_unifom_sample);

//    initTempK(tempK, src_sample_nodes);

//    int num_sample_nodes = src_sample_nodes.nodeSize();
//    std::cout << "sample node size = " << num_sample_nodes << std::endl;

//    Eigen::VectorXd affineVector_L0 = initAffineVector(num_sample_nodes);
//    Eigen::VectorXd accumAffineVector = affineVector_L0;

//    mat_U0_.resize(3, n_src_vertex_);
//    FindClosestPoints(correspondence_pairs_);

//    //	Solve L0 minimization
//    for (int iter = 0; iter < pars_.max_outer_iters; iter++)
//    {

//        std::cout << "out_iter = " << iter << std::endl;
//        //src_sample_nodes.print_nodes(*src_mesh_, iter);
//        // according correspondence_pairs to update mat_U0_;
//        mat_U0_.setZero();
//        for (int i = 0; i < correspondence_pairs_.size(); i++)
//        {
//            mat_U0_.col(correspondence_pairs_[i].first) = tar_points_.col(correspondence_pairs_[i].second);
//        }

//        // get auxiliary variable k_ij
//        subProblem1(mat_U0_, tempK, affineVector_L0, beta, lambda, src_sample_nodes);
//        subProblem2(*src_mesh_, tempK, affineVector_L0, accumAffineVector, beta, src_sample_nodes);
//        //std::cout << "affine vector L0 = " << affineVector_L0.head(10).transpose() << std::endl;

//        double e1, e2;
//        updateNonRigidMesh(*src_mesh_, affineVector_L0, src_sample_nodes, e1, e2);
//        FindClosestPoints(correspondence_pairs_);

//        if(pars_.print_each_step_info)
//        {
//            std::string out_node = pars_.out_each_step_info + "/iter_nodes_" + std::to_string(iter) + ".obj";
//            src_sample_nodes.print_nodes(*src_mesh_, out_node);
//            std::string out_obj = pars_.out_each_step_info + "/iter_obj_" + std::to_string(iter) + ".ply";
//            OpenMesh::IO::write_mesh( *src_mesh_, out_obj.c_str());
//        }
//        src_sample_nodes.updateWeight(*src_mesh_);
//    }
//    return 0.0;
//}

double SVRL0::DoNonRigid()
{
    // paras
//    double beta = 1;
//    double lambda = 0.02;
    Eigen::VectorXd affineVector_L0 = initAffineVector(num_sample_nodes);
    Eigen::VectorXd accumAffineVector = affineVector_L0;
    Eigen::VectorXd prev_affine = affineVector_L0;

    Timer time;
    Timer::EventID begin_time, run_time;
    pars_.each_energys.clear();
    pars_.each_gt_mean_errs.clear();
    pars_.each_gt_max_errs.clear();
    pars_.each_times.clear();
    pars_.each_iters.clear();

    // print initial run results
    pars_.each_energys.push_back(0.0);
    pars_.each_gt_max_errs.push_back(pars_.init_gt_max_errs);
    pars_.each_gt_mean_errs.push_back(pars_.init_gt_mean_errs);
    pars_.each_iters.push_back(0);
    pars_.each_times.push_back(pars_.non_rigid_init_time);

    begin_time = time.get_time();

    //	Solve L0 minimization
    for (int iter = 0; iter < pars_.max_outer_iters; iter++)
    {

        //std::cout << "out_iter = " << iter << std::endl;
        // according correspondence_pairs to update mat_U0_;
        mat_U0_.setZero();
        corres_pair_ids_.setZero();
        for (int i = 0; i < correspondence_pairs_.size(); i++)
        {
            mat_U0_.col(correspondence_pairs_[i].first) = tar_points_.col(correspondence_pairs_[i].second);
            corres_pair_ids_[i] = 1;
        }

//        for (int i = 0; i < n_src_vertex_; i++)
//        {
//            mat_U0_.col(i) = tar_points_.col(i);
//            corres_pair_ids_[i] = 1;
//        }
        //std::cout << "corres_pair_ids_ size = " << corres_pair_ids_.sum() << std::endl;

        // begin test

        double lambda = 0.02;
        double beta = 50 * lambda;
        unsigned int optCnt = 0;
//        while (beta < 1e6 && optCnt < pars_.max_inner_iters)
//        {
            //mainLogger << std::endl << "beta = " << beta << ", lambda = " << lambda << std::endl;

            subProblem1(mat_U0_, tempK, affineVector_L0, beta, lambda, src_sample_nodes);
            subProblem2(*src_mesh_, tempK, affineVector_L0, accumAffineVector, beta, src_sample_nodes);
//            beta *= 2.0;
//            ++optCnt;
//            std::cout << "beta = " << beta << " optcnt = " << optCnt << std::endl;
//        }

        // end test

       /* // get auxiliary variable k_ij
        subProblem1(mat_U0_, tempK, affineVector_L0, beta, lambda, src_sample_nodes);
        subProblem2(*src_mesh_, tempK, affineVector_L0, accumAffineVector, beta, src_sample_nodes);
        //std::cout << "affine vector L0 = " << affineVector_L0.head(10).transpose() << std::endl*/;


        double gt_mean_err, gt_max_err;
        updateNonRigidMesh(*src_mesh_, affineVector_L0, src_sample_nodes, gt_mean_err, gt_max_err);

        run_time = time.get_time();
        pars_.each_gt_mean_errs.push_back(gt_mean_err);
        pars_.each_gt_max_errs.push_back(gt_max_err);
        double energy = 0.0;
        pars_.each_energys.push_back(energy);
        double eps_time = time.elapsed_time(begin_time, run_time);
        pars_.each_times.push_back(eps_time);
        pars_.each_iters.push_back(1);

        if(pars_.print_each_step_info)
        {
            std::cout << "iter = " << iter << " time = " << eps_time << " energy = " << energy << " gt_err = " << gt_mean_err << std::endl;
            std::string out_node = pars_.out_each_step_info + "/iter_nodes_" + std::to_string(iter) + ".obj";
            src_sample_nodes.print_nodes(*src_mesh_, out_node);
            std::string out_obj = pars_.out_each_step_info + "/iter_obj_" + std::to_string(iter) + ".ply";
            OpenMesh::IO::write_mesh( *src_mesh_, out_obj.c_str());
        }

        FindClosestPoints(correspondence_pairs_);

        if((affineVector_L0 - prev_affine).norm() < pars_.stop)
        {
            break;
        }
        prev_affine = affineVector_L0;
        //src_sample_nodes.updateWeight(*src_mesh_);
    }
    return 0.0;
}


//------------------------------------------------------------------------
//	Define initTempK member function
//------------------------------------------------------------------------
void SVRL0::initTempK(std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
    const svr::nodeSampler &sampler) const
{
    for (size_t nodeIdx0 = 0; nodeIdx0 < sampler.nodeSize(); ++nodeIdx0)
    {
        svr::neighborIter nIter = sampler.getNodeNodeIter(nodeIdx0);
        for (; nIter.is_valid(); ++nIter)
        {
            size_t nodeIdx1 = nIter.getIndex();
            tempK[nodeIdx0][nodeIdx1].first = false;
            tempK[nodeIdx0][nodeIdx1].second = Eigen::Vector3d::Zero();
        }
    }
}


Eigen::VectorXd SVRL0::initAffineVector(size_t nodeSize) const
{
    Eigen::VectorXd affineVector(12 * nodeSize);

#pragma omp parallel for
    for (int nodeIdx = 0; nodeIdx < nodeSize; ++nodeIdx)
    {
        affineVector(12 * nodeIdx + 0) = 1.0f;
        affineVector(12 * nodeIdx + 1) = 0.0f;
        affineVector(12 * nodeIdx + 2) = 0.0f;

        affineVector(12 * nodeIdx + 3) = 0.0f;
        affineVector(12 * nodeIdx + 4) = 1.0f;
        affineVector(12 * nodeIdx + 5) = 0.0f;

        affineVector(12 * nodeIdx + 6) = 0.0f;
        affineVector(12 * nodeIdx + 7) = 0.0f;
        affineVector(12 * nodeIdx + 8) = 1.0f;

        affineVector(12 * nodeIdx + 9) = 0.0f;
        affineVector(12 * nodeIdx + 10) = 0.0f;
        affineVector(12 * nodeIdx + 11) = 0.0f;
    }

    return affineVector;
}


double SVRL0::gaussNewton(const Mesh &currentNonRigidMesh,
    const Mesh &currentScanMesh,
    const std::vector<std::pair<size_t, size_t>> &vertexPair,
    const svr::nodeSampler& sampler,
    Eigen::VectorXd &affineVector,
    const gnParams &controlParams)
{
    double gaussNewtonEnergy = 0.0;
    double gaussNewtonEnergyPrev = 0.0;
    size_t gaussNewtonIterCnt = 0;

    while (true)
    {
        //------------------------------------------------------------------------
        //	Initialize numerical containers
        //------------------------------------------------------------------------
        std::vector<Eigen::Triplet<double>> JrList;
        std::vector<Eigen::Triplet<double>> AsList;
        std::vector<Eigen::Triplet<double>> ApList;
        std::vector<Eigen::Triplet<double>> AqList;

        Eigen::SparseMatrix<double> Jr;
        Eigen::SparseMatrix<double> As;
        Eigen::SparseMatrix<double> Ap;
        Eigen::SparseMatrix<double> Aq;

        Eigen::VectorXd fr;
        Eigen::VectorXd bs;
        Eigen::VectorXd bp;
        Eigen::VectorXd bq;

        Eigen::SparseMatrix<double> A0(12 * sampler.nodeSize(), 12 * sampler.nodeSize());
        Eigen::SparseMatrix<double> A1(12 * sampler.nodeSize(), 12 * sampler.nodeSize());
        Eigen::SparseMatrix<double> A2(12 * sampler.nodeSize(), 12 * sampler.nodeSize());
        Eigen::SparseMatrix<double> A3(12 * sampler.nodeSize(), 12 * sampler.nodeSize());

        Eigen::VectorXd b0(12 * sampler.nodeSize());
        Eigen::VectorXd b1(12 * sampler.nodeSize());
        Eigen::VectorXd b2(12 * sampler.nodeSize());
        Eigen::VectorXd b3(12 * sampler.nodeSize());

#pragma omp sections
        {
            //	First independent section
            //	Construct Jr, fr, A0, and b0
#pragma omp section
        {
            //------------------------------------------------------------------------
            //	Construct fr
            //------------------------------------------------------------------------
            std::vector<double> frVector;
            for (int i = 0; i < sampler.nodeSize(); ++i)
            {
                double rigid[6] = { 0.0 };
                rigid[0] += affineVector(12 * i + 0) * affineVector(12 * i + 3);
                rigid[0] += affineVector(12 * i + 1) * affineVector(12 * i + 4);
                rigid[0] += affineVector(12 * i + 2) * affineVector(12 * i + 5);

                rigid[1] += affineVector(12 * i + 3) * affineVector(12 * i + 6);
                rigid[1] += affineVector(12 * i + 4) * affineVector(12 * i + 7);
                rigid[1] += affineVector(12 * i + 5) * affineVector(12 * i + 8);

                rigid[2] += affineVector(12 * i + 0) * affineVector(12 * i + 6);
                rigid[2] += affineVector(12 * i + 1) * affineVector(12 * i + 7);
                rigid[2] += affineVector(12 * i + 2) * affineVector(12 * i + 8);

                rigid[3] += affineVector(12 * i + 0) * affineVector(12 * i + 0);
                rigid[3] += affineVector(12 * i + 1) * affineVector(12 * i + 1);
                rigid[3] += affineVector(12 * i + 2) * affineVector(12 * i + 2);
                rigid[3] -= 1.0;

                rigid[4] += affineVector(12 * i + 3) * affineVector(12 * i + 3);
                rigid[4] += affineVector(12 * i + 4) * affineVector(12 * i + 4);
                rigid[4] += affineVector(12 * i + 5) * affineVector(12 * i + 5);
                rigid[4] -= 1.0;

                rigid[5] += affineVector(12 * i + 6) * affineVector(12 * i + 6);
                rigid[5] += affineVector(12 * i + 7) * affineVector(12 * i + 7);
                rigid[5] += affineVector(12 * i + 8) * affineVector(12 * i + 8);
                rigid[5] -= 1.0;

                frVector.push_back(rigid[0]);
                frVector.push_back(rigid[1]);
                frVector.push_back(rigid[2]);
                frVector.push_back(rigid[3]);
                frVector.push_back(rigid[4]);
                frVector.push_back(rigid[5]);
            }

            fr.resize(frVector.size());
            for (int i = 0; i < frVector.size(); ++i)
            {
                fr(i) = frVector[i];
            }

            //------------------------------------------------------------------------
            //	Construct Jr
            //------------------------------------------------------------------------
            size_t jOffset = 0;
            for (int i = 0; i < sampler.nodeSize(); ++i)
            {
                double x[9] = { 0.0 };
                x[0] = affineVector(12 * i + 0);
                x[1] = affineVector(12 * i + 1);
                x[2] = affineVector(12 * i + 2);
                x[3] = affineVector(12 * i + 3);
                x[4] = affineVector(12 * i + 4);
                x[5] = affineVector(12 * i + 5);
                x[6] = affineVector(12 * i + 6);
                x[7] = affineVector(12 * i + 7);
                x[8] = affineVector(12 * i + 8);

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[3]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[4]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[5]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[0]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[1]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[2]));
                ++jOffset;

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[6]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[7]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[8]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[3]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[4]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[5]));
                ++jOffset;

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[6]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[7]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[8]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[0]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[1]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[2]));
                ++jOffset;

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, 2 * x[0]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, 2 * x[1]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, 2 * x[2]));
                ++jOffset;

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, 2 * x[3]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, 2 * x[4]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, 2 * x[5]));
                ++jOffset;

                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, 2 * x[6]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, 2 * x[7]));
                JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, 2 * x[8]));
                ++jOffset;
            }

            Jr.resize(jOffset, 12 * sampler.nodeSize());
            Jr.setFromTriplets(JrList.begin(), JrList.end());

            Eigen::SparseMatrix<double> JrTJr = Eigen::SparseMatrix<double>(Jr.transpose()) * Jr;
            A0 = controlParams.m_alphaRigid * JrTJr;
            b0 = controlParams.m_alphaRigid * JrTJr * affineVector.cast<double>();
            b0 -= controlParams.m_alphaRigid * Jr.transpose() * fr;
        }

        //	Second independent section
        //	Construct As, bs, A1 and b1
#pragma omp section
        {
            //------------------------------------------------------------------------
            //	Construct As and bs
            //------------------------------------------------------------------------
            std::vector<double> bsVector;
            size_t sOffset = 0;

            for (int idx0 = 0; idx0 < sampler.nodeSize(); ++idx0)
            {
                size_t nodeIdx0 = sampler.getNodeVertexIdx(idx0);
                VertexHandle vh0 = currentNonRigidMesh.vertex_handle(nodeIdx0);
                Point v0 = currentNonRigidMesh.point(vh0);

                svr::neighborIter nIter = sampler.getNodeNodeIter(idx0);
                for (; nIter.is_valid(); ++nIter)
                {
                    size_t idx1 = nIter.getIndex();
                    double weight = nIter.getWeight();
                    double weightRoot = sqrt(weight);
                    size_t nodeIdx1 = sampler.getNodeVertexIdx(idx1);
                    VertexHandle vh1 = currentNonRigidMesh.vertex_handle(nodeIdx1);
                    Point v1 = currentNonRigidMesh.point(vh1);

                    double vec[3] = { 0.0 };
                    vec[0] = weightRoot * (v1[0] - v0[0]);
                    vec[1] = weightRoot * (v1[1] - v0[1]);
                    vec[2] = weightRoot * (v1[2] - v0[2]);

                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 0, vec[0]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 3, vec[1]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 6, vec[2]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 9, weightRoot));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx1 + 9, -weightRoot));
                    ++sOffset;

                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 1, vec[0]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 4, vec[1]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 7, vec[2]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 10, weightRoot));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx1 + 10, -weightRoot));
                    ++sOffset;

                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 2, vec[0]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 5, vec[1]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 8, vec[2]));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx0 + 11, weightRoot));
                    AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * idx1 + 11, -weightRoot));
                    ++sOffset;

                    bsVector.push_back(vec[0]);
                    bsVector.push_back(vec[1]);
                    bsVector.push_back(vec[2]);
                }
            }

            As.resize(sOffset, 12 * sampler.nodeSize());
            As.setFromTriplets(AsList.begin(), AsList.end());

            bs.resize(bsVector.size());
            for (int i = 0; i < bsVector.size(); ++i)
            {
                bs(i) = bsVector[i];
            }

            A1 = controlParams.m_alphaSmooth * Eigen::SparseMatrix<double>(As.transpose()) * As;
            b1 = controlParams.m_alphaSmooth * As.transpose() * bs;
        }

        //	Third independent section
        //	Construct Ap, bp, A2 and b2
#pragma omp section
        {
            //------------------------------------------------------------------------
            //	Construct Ap and bp
            //------------------------------------------------------------------------
            std::vector<double> bpVector;
            size_t pOffset = 0;

            for (int i = 0; i < vertexPair.size(); ++i)
            {
                size_t vIdx = vertexPair[i].first;
                size_t cIdx = vertexPair[i].second;
                VertexHandle vh = currentNonRigidMesh.vertex_handle(vIdx);
                VertexHandle ch = currentScanMesh.vertex_handle(cIdx);
                Point v = currentNonRigidMesh.point(vh);
                Point c = currentScanMesh.point(ch);

                double rhs[3] = { 0.0 };

                svr::neighborIter nIter = sampler.getVertexNodeIter(vIdx);
                for (; nIter.is_valid(); ++nIter)
                {
                    size_t nIdx = nIter.getIndex();
                    double weight = nIter.getWeight();
                    size_t nodeIdx = sampler.getNodeVertexIdx(nIdx);
                    VertexHandle nh = currentNonRigidMesh.vertex_handle(nodeIdx);
                    Point n = currentNonRigidMesh.point(nh);

                    double vec[3] = { 0.0 };
                    vec[0] = weight * (v[0] - n[0]);
                    vec[1] = weight * (v[1] - n[1]);
                    vec[2] = weight * (v[2] - n[2]);

                    ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nIdx + 0, vec[0]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nIdx + 3, vec[1]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nIdx + 6, vec[2]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nIdx + 9, weight));

                    ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nIdx + 1, vec[0]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nIdx + 4, vec[1]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nIdx + 7, vec[2]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nIdx + 10, weight));

                    ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nIdx + 2, vec[0]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nIdx + 5, vec[1]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nIdx + 8, vec[2]));
                    ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nIdx + 11, weight));

                    rhs[0] -= weight * n[0];
                    rhs[1] -= weight * n[1];
                    rhs[2] -= weight * n[2];
                }

                rhs[0] += c[0];
                rhs[1] += c[1];
                rhs[2] += c[2];

                bpVector.push_back(rhs[0]);
                bpVector.push_back(rhs[1]);
                bpVector.push_back(rhs[2]);

                pOffset += 3;
            }

            Ap.resize(pOffset, 12 * sampler.nodeSize());
            Ap.setFromTriplets(ApList.begin(), ApList.end());

            bp.resize(bpVector.size());
            for (int i = 0; i < bpVector.size(); ++i)
            {
                bp(i) = bpVector[i];
            }

            A2 = controlParams.m_alphaPoint * Eigen::SparseMatrix<double>(Ap.transpose()) * Ap;
            b2 = controlParams.m_alphaPoint * Ap.transpose() * bp;
        }

        //	Fourth independent section
        //	Construct Aq, bq, A3 and b3
#pragma omp section
        {
            //------------------------------------------------------------------------
            //	Construct Aq and bq
            //------------------------------------------------------------------------
            std::vector<double> bqVector;
            size_t qOffset = 0;

            for (int i = 0; i < vertexPair.size(); ++i)
            {
                size_t vIdx = vertexPair[i].first;
                size_t cIdx = vertexPair[i].second;
                VertexHandle vh = currentNonRigidMesh.vertex_handle(vIdx);
                VertexHandle ch = currentScanMesh.vertex_handle(cIdx);
                Point v = currentNonRigidMesh.point(vh);
                Point c = currentScanMesh.point(ch);
                //Vec3f cN = currentScanMesh.normal(ch);
                Vec3d cN = currentNonRigidMesh.normal(vh);

                double rhs = 0.0;

                svr::neighborIter nIter = sampler.getVertexNodeIter(vIdx);
                for (; nIter.is_valid(); ++nIter)
                {
                    size_t nIdx = nIter.getIndex();
                    double weight = nIter.getWeight();
                    size_t nodeIdx = sampler.getNodeVertexIdx(nIdx);
                    VertexHandle nh = currentNonRigidMesh.vertex_handle(nodeIdx);
                    Point n = currentNonRigidMesh.point(nh);

                    double vec[3] = { 0.0 };
                    vec[0] = weight * (v[0] - n[0]);
                    vec[1] = weight * (v[1] - n[1]);
                    vec[2] = weight * (v[2] - n[2]);

                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 0, cN[0] * vec[0]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 1, cN[1] * vec[0]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 2, cN[2] * vec[0]));

                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 3, cN[0] * vec[1]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 4, cN[1] * vec[1]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 5, cN[2] * vec[1]));

                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 6, cN[0] * vec[2]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 7, cN[1] * vec[2]));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 8, cN[2] * vec[2]));

                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 9, cN[0] * weight));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 10, cN[1] * weight));
                    AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nIdx + 11, cN[2] * weight));

                    rhs -= cN[0] * weight * n[0];
                    rhs -= cN[1] * weight * n[1];
                    rhs -= cN[2] * weight * n[2];
                }

                rhs += cN[0] * c[0];
                rhs += cN[1] * c[1];
                rhs += cN[2] * c[2];

                bqVector.push_back(rhs);

                ++qOffset;
            }

            Aq.resize(qOffset, 12 * sampler.nodeSize());
            Aq.setFromTriplets(AqList.begin(), AqList.end());

            bq.resize(bqVector.size());
            for (int i = 0; i < bqVector.size(); ++i)
            {
                bq(i) = bqVector[i];
            }

            A3 = controlParams.m_alphaPlane * Eigen::SparseMatrix<double>(Aq.transpose()) * Aq;
            b3 = controlParams.m_alphaPlane * Aq.transpose() * bq;
        }
        }

        //------------------------------------------------------------------------
        //	Construct linear equation A * x = b
        //	A = m_alphaRigid * JrT * Jr +
        //		m_alphaSmooth * AsT * As +
        //		m_alphaPoint * ApT * Ap +
        //		m_alphaPlane * AqT * Aq +
        //
        //	b = - m_alphaRigid * JrT * fr +
        //		  m_alphaSmooth * AsT * bs +
        //		  m_alphaPoint * ApT * bp +
        //		  m_alphaPlane * AqT * bq +
        //		  m_alphaRigid * JrT * Jr * affineVector
        //
        //	x = nextAffineVector
        //------------------------------------------------------------------------
        Eigen::SparseMatrix<double> A = A0 + A1 + A2 + A3;
        Eigen::VectorXd b = b0 + b1 + b2 + b3;

        /*Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> CholmodSolver;
        CholmodSolver.compute(A);
        Eigen::VectorXd nextAffineVector = CholmodSolver.solve(b);*/
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> SimplicialLDLTSolver;
        SimplicialLDLTSolver.compute(A);
        Eigen::VectorXd nextAffineVector = SimplicialLDLTSolver.solve(b);

        //------------------------------------------------------------------------
        //	Calculate total energy energyGaussNewton
        //	Break the loop if relative energy change < 1e-6f,
        //	or total iteration is greater than 50
        //------------------------------------------------------------------------
        Eigen::VectorXd RigidVector(12 * sampler.nodeSize());
        Eigen::VectorXd SmoothVector(12 * sampler.nodeSize());
        Eigen::VectorXd PointVector(12 * sampler.nodeSize());
        Eigen::VectorXd PlaneVector(12 * sampler.nodeSize());

        double RigidEnergy = 0.0;
        double SmoothEnergy = 0.0;
        double PointEnergy = 0.0;
        double PlaneEnergy = 0.0;

        //	Parallelize energy calculation
#pragma omp sections
        {
            //	First section calculates rigid energy
#pragma omp section
        {
            RigidVector = fr + Jr * (nextAffineVector - affineVector.cast<double>());
            RigidEnergy = controlParams.m_alphaRigid * RigidVector.transpose() * RigidVector;
        }
        //	Second section calculates smooth energy
#pragma omp section
        {
            SmoothVector = As * nextAffineVector - bs;
            SmoothEnergy = controlParams.m_alphaSmooth * SmoothVector.transpose() * SmoothVector;
        }
        //	Third section calculates point-to-point energy
#pragma omp section
        {
            PointVector = Ap * nextAffineVector - bp;
            PointEnergy = controlParams.m_alphaPoint * PointVector.transpose() * PointVector;
        }
        //	Fourth section calculates point-to-plane energy
#pragma omp section
        {
            PlaneVector = Aq * nextAffineVector - bq;
            PlaneEnergy = controlParams.m_alphaPlane * PlaneVector.transpose() * PlaneVector;
        }
        }

        //	Sum up all four energy
        gaussNewtonEnergy = RigidEnergy + SmoothEnergy + PointEnergy + PlaneEnergy;

        //	Break Gauss-Newton loop if relative energy change is less than 1e-6,
        //	or total iteration is greater than 50
        double gaussNewtonEnergyChange = fabs(gaussNewtonEnergy - gaussNewtonEnergyPrev) / (gaussNewtonEnergyPrev + 1.0);
        //if (verbose)
        //{
        //	mainLogger << "GN iteration " << gaussNewtonIterCnt << ", "
        //		<< "gnEnergy = " << gaussNewtonEnergy << ", gnEnergyChange = " << gaussNewtonEnergyChange << std::endl;
        //}

        affineVector = nextAffineVector.cast<double>();

        if (gaussNewtonEnergyChange < 1e-6 || gaussNewtonIterCnt >= 50) break;
        gaussNewtonEnergyPrev = gaussNewtonEnergy;

        ++gaussNewtonIterCnt;
    }

    return gaussNewtonEnergy;
}

//------------------------------------------------------------------------
//	Define accumulatedAffine member function
//------------------------------------------------------------------------
void SVRL0::accumulateAffine(const Eigen::VectorXd &affineVector,
    Eigen::VectorXd &accumulatedAffineVector,
    size_t nodeSize) const
{
#pragma omp parallel for
    for (int nodeIdx = 0; nodeIdx < nodeSize; ++nodeIdx)
    {
        Eigen::Matrix3d rotMat;
        Eigen::Vector3d shiftVec;
        Eigen::Matrix3d accumRotMat;
        Eigen::Vector3d accumShiftVec;

        //	Initialize rotMat
        rotMat(0, 0) = affineVector(12 * nodeIdx + 0);
        rotMat(1, 0) = affineVector(12 * nodeIdx + 1);
        rotMat(2, 0) = affineVector(12 * nodeIdx + 2);
        rotMat(0, 1) = affineVector(12 * nodeIdx + 3);
        rotMat(1, 1) = affineVector(12 * nodeIdx + 4);
        rotMat(2, 1) = affineVector(12 * nodeIdx + 5);
        rotMat(0, 2) = affineVector(12 * nodeIdx + 6);
        rotMat(1, 2) = affineVector(12 * nodeIdx + 7);
        rotMat(2, 2) = affineVector(12 * nodeIdx + 8);

        //	Initialize shiftVec
        shiftVec(0) = affineVector(12 * nodeIdx + 9);
        shiftVec(1) = affineVector(12 * nodeIdx + 10);
        shiftVec(2) = affineVector(12 * nodeIdx + 11);

        //	Initialize accumRotMat
        accumRotMat(0, 0) = accumulatedAffineVector(12 * nodeIdx + 0);
        accumRotMat(1, 0) = accumulatedAffineVector(12 * nodeIdx + 1);
        accumRotMat(2, 0) = accumulatedAffineVector(12 * nodeIdx + 2);
        accumRotMat(0, 1) = accumulatedAffineVector(12 * nodeIdx + 3);
        accumRotMat(1, 1) = accumulatedAffineVector(12 * nodeIdx + 4);
        accumRotMat(2, 1) = accumulatedAffineVector(12 * nodeIdx + 5);
        accumRotMat(0, 2) = accumulatedAffineVector(12 * nodeIdx + 6);
        accumRotMat(1, 2) = accumulatedAffineVector(12 * nodeIdx + 7);
        accumRotMat(2, 2) = accumulatedAffineVector(12 * nodeIdx + 8);

        //	Initialize accumShiftVec
        accumShiftVec(0) = accumulatedAffineVector(12 * nodeIdx + 9);
        accumShiftVec(1) = accumulatedAffineVector(12 * nodeIdx + 10);
        accumShiftVec(2) = accumulatedAffineVector(12 * nodeIdx + 11);

        //	Accumulate rotation and translation
        accumRotMat = rotMat * accumRotMat;
        accumShiftVec = shiftVec + accumShiftVec;

        //	Assign accumulated results back to accumulatedAffineVector
        accumulatedAffineVector(12 * nodeIdx + 0) = accumRotMat(0, 0);
        accumulatedAffineVector(12 * nodeIdx + 1) = accumRotMat(1, 0);
        accumulatedAffineVector(12 * nodeIdx + 2) = accumRotMat(2, 0);
        accumulatedAffineVector(12 * nodeIdx + 3) = accumRotMat(0, 1);
        accumulatedAffineVector(12 * nodeIdx + 4) = accumRotMat(1, 1);
        accumulatedAffineVector(12 * nodeIdx + 5) = accumRotMat(2, 1);
        accumulatedAffineVector(12 * nodeIdx + 6) = accumRotMat(0, 2);
        accumulatedAffineVector(12 * nodeIdx + 7) = accumRotMat(1, 2);
        accumulatedAffineVector(12 * nodeIdx + 8) = accumRotMat(2, 2);

        accumulatedAffineVector(12 * nodeIdx + 9) = accumShiftVec(0);
        accumulatedAffineVector(12 * nodeIdx + 10) = accumShiftVec(1);
        accumulatedAffineVector(12 * nodeIdx + 11) = accumShiftVec(2);
    }
}

//------------------------------------------------------------------------
//	Define updateNonRigidMesh member function
//------------------------------------------------------------------------
double SVRL0::updateNonRigidMesh(Mesh &currentNonRigidMesh,
    const Eigen::VectorXd &affineVector,
    const svr::nodeSampler &sampler, double& gt_mean_err, double& gt_max_err)
{
    //	Replicate a copy of node coordinates
    std::vector<Point> NodeCoords(sampler.nodeSize());
#pragma omp parallel for
    for (int i = 0; i < sampler.nodeSize(); ++i)
    {
        int nodeIdx = sampler.getNodeVertexIdx(i);
        VertexHandle nh = currentNonRigidMesh.vertex_handle(nodeIdx);
        Point vNode = currentNonRigidMesh.point(nh);
        NodeCoords[i] = vNode;
    }

    //	Update vertex coordinates of CurrentNonRigidMesh
    Eigen::VectorXd gt_errs = Eigen::VectorXd::Zero(n_src_vertex_);
#pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; ++i)
    {
        VertexHandle vh = currentNonRigidMesh.vertex_handle(i);
        Point v = currentNonRigidMesh.point(vh);
        Eigen::Vector3d vUpdated = Eigen::Vector3d::Zero();

        svr::neighborIter nIter = sampler.getVertexNodeIter(i);
        for (; nIter.is_valid(); ++nIter)
        {
            int idx = nIter.getIndex();
            double weight = nIter.getWeight();
            Point &vNode = NodeCoords[idx];

            Eigen::Matrix3d Aj = Eigen::Matrix3d::Identity();
            Eigen::Vector3d bj = Eigen::Vector3d::Zero();
            Eigen::Vector3d vi = Eigen::Vector3d::Zero();
            Eigen::Vector3d uj = Eigen::Vector3d::Zero();

            Aj(0, 0) = affineVector(12 * idx + 0);
            Aj(1, 0) = affineVector(12 * idx + 1);
            Aj(2, 0) = affineVector(12 * idx + 2);
            Aj(0, 1) = affineVector(12 * idx + 3);
            Aj(1, 1) = affineVector(12 * idx + 4);
            Aj(2, 1) = affineVector(12 * idx + 5);
            Aj(0, 2) = affineVector(12 * idx + 6);
            Aj(1, 2) = affineVector(12 * idx + 7);
            Aj(2, 2) = affineVector(12 * idx + 8);

            bj(0) = affineVector(12 * idx + 9);
            bj(1) = affineVector(12 * idx + 10);
            bj(2) = affineVector(12 * idx + 11);

            vi(0) = v[0];
            vi(1) = v[1];
            vi(2) = v[2];

            uj(0) = vNode[0];
            uj(1) = vNode[1];
            uj(2) = vNode[2];

            vUpdated += weight * (Aj * (vi - uj) + uj + bj);
        }
        currentNonRigidMesh.set_point(vh, Point(vUpdated(0), vUpdated(1), vUpdated(2)));
        if(pars_.calc_gt_err)
            gt_errs[i] = (vUpdated - tar_points_.col(i)).squaredNorm();
    }
    gt_mean_err = sqrt(gt_errs.sum()/n_src_vertex_);
    gt_max_err = sqrt(gt_errs.maxCoeff());
    return gt_mean_err;
}

//------------------------------------------------------------------------
    //	Define subProblem1 member function
    //	Prune tempK
    //------------------------------------------------------------------------
    void SVRL0::subProblem1(const MatrixXX &corresMesh,
        std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
        const Eigen::VectorXd &affineVector, double beta, double lambda,
        const svr::nodeSampler &sampler)
{
#pragma omp parallel for
    for (int nodeIdx0 = 0; nodeIdx0 < sampler.nodeSize(); ++nodeIdx0)
    {
        size_t vIdx0 = sampler.getNodeVertexIdx(nodeIdx0);
        Point v0(corresMesh(0,vIdx0), corresMesh(1, vIdx0), corresMesh(2, vIdx0));

        svr::neighborIter nIter = sampler.getNodeNodeIter(nodeIdx0);
        for (; nIter.is_valid(); ++nIter)
        {
            size_t nodeIdx1 = nIter.getIndex();

            // v1 is the neighbor of v0 (node)
            size_t vIdx1 = sampler.getNodeVertexIdx(nodeIdx1);
            Point v1(corresMesh(0, vIdx1), corresMesh(1, vIdx1), corresMesh(2, vIdx1));

            double delta[3] = { 0 };
            double sqSum = 0.0f;

            // reg-perm Dx_ij
            delta[0] += affineVector(12 * nodeIdx0 + 0) * (v1[0] - v0[0]);
            delta[0] += affineVector(12 * nodeIdx0 + 3) * (v1[1] - v0[1]);
            delta[0] += affineVector(12 * nodeIdx0 + 6) * (v1[2] - v0[2]);
            delta[0] += v0[0] + affineVector(12 * nodeIdx0 + 9);
            delta[0] -= v1[0] + affineVector(12 * nodeIdx1 + 9);

            delta[1] += affineVector(12 * nodeIdx0 + 1) * (v1[0] - v0[0]);
            delta[1] += affineVector(12 * nodeIdx0 + 4) * (v1[1] - v0[1]);
            delta[1] += affineVector(12 * nodeIdx0 + 7) * (v1[2] - v0[2]);
            delta[1] += v0[1] + affineVector(12 * nodeIdx0 + 10);
            delta[1] -= v1[1] + affineVector(12 * nodeIdx1 + 10);

            delta[2] += affineVector(12 * nodeIdx0 + 2) * (v1[0] - v0[0]);
            delta[2] += affineVector(12 * nodeIdx0 + 5) * (v1[1] - v0[1]);
            delta[2] += affineVector(12 * nodeIdx0 + 8) * (v1[2] - v0[2]);
            delta[2] += v0[2] + affineVector(12 * nodeIdx0 + 11);
            delta[2] -= v1[2] + affineVector(12 * nodeIdx1 + 11);

            sqSum += delta[0] * delta[0];
            sqSum += delta[1] * delta[1];
            sqSum += delta[2] * delta[2];

#pragma omp critical
            {
                if (sqSum < lambda / beta)
                {
                    tempK[nodeIdx0][nodeIdx1].first = false;
                    tempK[nodeIdx0][nodeIdx1].second = Eigen::Vector3d::Zero();
                }
                else
                {
                    tempK[nodeIdx0][nodeIdx1].first = true;
                    tempK[nodeIdx0][nodeIdx1].second = Eigen::Vector3d(delta[0], delta[1], delta[2]);
                }
            }
        }
    }
}



    //------------------------------------------------------------------------
    //	Define subProblem1 member function
    //	Solve energy function of subProblem2
    //------------------------------------------------------------------------
    double SVRL0::subProblem2(const Mesh &anchorMesh,
        const std::vector<std::map<size_t, std::pair<bool, Eigen::Vector3d>>> &tempK,
        Eigen::VectorXd &affineVector_L0,
        const Eigen::VectorXd &accumAffineVector, double beta,
        const svr::nodeSampler &anchorSampler)
    {
        double gaussNewtonEnergy = 0.0;
        double gaussNewtonEnergyPrev = 0.0;
        size_t gaussNewtonIterCnt = 0;

        //bool run_once = true;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> CholmodSolver;

        while (true)
        {
            //------------------------------------------------------------------------
            //	Initialize numerical containers
            //------------------------------------------------------------------------
            std::vector<Eigen::Triplet<double>> JrList;
            std::vector<Eigen::Triplet<double>> AdList;
            std::vector<Eigen::Triplet<double>> AkList;

            Eigen::SparseMatrix<double> Jr;
            Eigen::SparseMatrix<double> Ad;
            Eigen::SparseMatrix<double> Ak;

            Eigen::VectorXd fr;
            Eigen::VectorXd bd;
            Eigen::VectorXd bk;

            Eigen::SparseMatrix<double> A0(12 * anchorSampler.nodeSize(), 12 * anchorSampler.nodeSize());
            Eigen::SparseMatrix<double> A1(12 * anchorSampler.nodeSize(), 12 * anchorSampler.nodeSize());
            Eigen::SparseMatrix<double> A2(12 * anchorSampler.nodeSize(), 12 * anchorSampler.nodeSize());

            Eigen::VectorXd b0(12 * anchorSampler.nodeSize());
            Eigen::VectorXd b1(12 * anchorSampler.nodeSize());
            Eigen::VectorXd b2(12 * anchorSampler.nodeSize());

#pragma omp sections
            {
                //	First independent section
                //	Construct Jr, fr, A0 and b0
#pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct fr
                //------------------------------------------------------------------------
                std::vector<double> frVector;
                for (int i = 0; i < anchorSampler.nodeSize(); ++i)
                {
                    double rigid[6] = { 0.0 };
                    rigid[0] += affineVector_L0(12 * i + 0) * affineVector_L0(12 * i + 3);
                    rigid[0] += affineVector_L0(12 * i + 1) * affineVector_L0(12 * i + 4);
                    rigid[0] += affineVector_L0(12 * i + 2) * affineVector_L0(12 * i + 5);

                    rigid[1] += affineVector_L0(12 * i + 3) * affineVector_L0(12 * i + 6);
                    rigid[1] += affineVector_L0(12 * i + 4) * affineVector_L0(12 * i + 7);
                    rigid[1] += affineVector_L0(12 * i + 5) * affineVector_L0(12 * i + 8);

                    rigid[2] += affineVector_L0(12 * i + 0) * affineVector_L0(12 * i + 6);
                    rigid[2] += affineVector_L0(12 * i + 1) * affineVector_L0(12 * i + 7);
                    rigid[2] += affineVector_L0(12 * i + 2) * affineVector_L0(12 * i + 8);

                    rigid[3] += affineVector_L0(12 * i + 0) * affineVector_L0(12 * i + 0);
                    rigid[3] += affineVector_L0(12 * i + 1) * affineVector_L0(12 * i + 1);
                    rigid[3] += affineVector_L0(12 * i + 2) * affineVector_L0(12 * i + 2);
                    rigid[3] -= 1.0;

                    rigid[4] += affineVector_L0(12 * i + 3) * affineVector_L0(12 * i + 3);
                    rigid[4] += affineVector_L0(12 * i + 4) * affineVector_L0(12 * i + 4);
                    rigid[4] += affineVector_L0(12 * i + 5) * affineVector_L0(12 * i + 5);
                    rigid[4] -= 1.0;

                    rigid[5] += affineVector_L0(12 * i + 6) * affineVector_L0(12 * i + 6);
                    rigid[5] += affineVector_L0(12 * i + 7) * affineVector_L0(12 * i + 7);
                    rigid[5] += affineVector_L0(12 * i + 8) * affineVector_L0(12 * i + 8);
                    rigid[5] -= 1.0;

                    frVector.push_back(rigid[0]);
                    frVector.push_back(rigid[1]);
                    frVector.push_back(rigid[2]);
                    frVector.push_back(rigid[3]);
                    frVector.push_back(rigid[4]);
                    frVector.push_back(rigid[5]);
                }

                fr.resize(frVector.size());
                for (int i = 0; i < frVector.size(); ++i)
                {
                    fr(i) = frVector[i];
                }

                //------------------------------------------------------------------------
                //	Construct Jr
                //------------------------------------------------------------------------
                size_t jOffset = 0;
                for (int i = 0; i < anchorSampler.nodeSize(); ++i)
                {
                    double x[9] = { 0.0 };
                    x[0] = affineVector_L0(12 * i + 0);
                    x[1] = affineVector_L0(12 * i + 1);
                    x[2] = affineVector_L0(12 * i + 2);
                    x[3] = affineVector_L0(12 * i + 3);
                    x[4] = affineVector_L0(12 * i + 4);
                    x[5] = affineVector_L0(12 * i + 5);
                    x[6] = affineVector_L0(12 * i + 6);
                    x[7] = affineVector_L0(12 * i + 7);
                    x[8] = affineVector_L0(12 * i + 8);

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[5]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[8]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[5]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[8]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, 2 * x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, 2 * x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, 2 * x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, 2 * x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, 2 * x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, 2 * x[5]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, 2 * x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, 2 * x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, 2 * x[8]));
                    ++jOffset;
                }

                Jr.resize(jOffset, 12 * anchorSampler.nodeSize());
                Jr.setFromTriplets(JrList.begin(), JrList.end());

                Eigen::SparseMatrix<double> JrTJr = Eigen::SparseMatrix<double>(Jr.transpose()) * Jr;
                A0 = g_alphaRigid * JrTJr;
                b0 = g_alphaRigid * JrTJr * affineVector_L0.cast<double>();
                b0 -= g_alphaRigid * Jr.transpose() * fr;
            }

            //	Second independent section
            //	Construct Ad, bd, A1 and b1
#pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct Ad and bd
                //------------------------------------------------------------------------
                std::vector<double> bdVector;
                size_t dOffset = 0;

                for (size_t vertexIdx = 0; vertexIdx < anchorMesh.n_vertices(); ++vertexIdx)
                {
                    VertexHandle vh = anchorMesh.vertex_handle(vertexIdx);
                    Point v = anchorMesh.point(vh);

                    double rhs[3] = { 0.0 };

                    svr::neighborIter nIter = anchorSampler.getVertexNodeIter(vertexIdx);
                    for (; nIter.is_valid(); ++nIter)
                    {
                        size_t nodeIdx = nIter.getIndex();
                        double weight = nIter.getWeight();
                        VertexHandle nh = anchorMesh.vertex_handle(anchorSampler.getNodeVertexIdx(nodeIdx));
                        Point n = anchorMesh.point(nh);

                        double vec[3] = { 0 };
                        vec[0] = weight * (v[0] - n[0]);
                        vec[1] = weight * (v[1] - n[1]);
                        vec[2] = weight * (v[2] - n[2]);

                        AdList.push_back(Eigen::Triplet<double>(dOffset + 0, 12 * nodeIdx + 0, vec[0]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 0, 12 * nodeIdx + 3, vec[1]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 0, 12 * nodeIdx + 6, vec[2]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 0, 12 * nodeIdx + 9, weight));

                        AdList.push_back(Eigen::Triplet<double>(dOffset + 1, 12 * nodeIdx + 1, vec[0]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 1, 12 * nodeIdx + 4, vec[1]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 1, 12 * nodeIdx + 7, vec[2]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 1, 12 * nodeIdx + 10, weight));

                        AdList.push_back(Eigen::Triplet<double>(dOffset + 2, 12 * nodeIdx + 2, vec[0]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 2, 12 * nodeIdx + 5, vec[1]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 2, 12 * nodeIdx + 8, vec[2]));
                        AdList.push_back(Eigen::Triplet<double>(dOffset + 2, 12 * nodeIdx + 11, weight));

                        rhs[0] -= weight * n[0];
                        rhs[1] -= weight * n[1];
                        rhs[2] -= weight * n[2];

                        ///////////////////////////////////////////////////
                        ////\B4?\A6\D0\E8?\D0??\AC\D6\F7?\CA\C7\D5\E2\B8\F6\B4\FA\C2\EB\CA\C7\D4\DA?\B5\C0??\B5?»Œ\B5?\F9\B4\A1\C9\CFß’\B5?\AC\D0\E8?\D0?\C4?\B6\D4?\B5\C4?\D7\EE\BD\FC\B5\E3\B5\C4\D7\F8\B1\EA\D0\CE?
                        //rhs[0] += accumAffineVector(12 * nodeIdx + 0) * vec[0];
                        //rhs[0] += accumAffineVector(12 * nodeIdx + 3) * vec[1];
                        //rhs[0] += accumAffineVector(12 * nodeIdx + 6) * vec[2];
                        //rhs[0] += weight * (n[0] + accumAffineVector(12 * nodeIdx + 9));

                        //rhs[1] += accumAffineVector(12 * nodeIdx + 1) * vec[0];
                        //rhs[1] += accumAffineVector(12 * nodeIdx + 4) * vec[1];
                        //rhs[1] += accumAffineVector(12 * nodeIdx + 7) * vec[2];
                        //rhs[1] += weight * (n[1] + accumAffineVector(12 * nodeIdx + 10));

                        //rhs[2] += accumAffineVector(12 * nodeIdx + 2) * vec[0];
                        //rhs[2] += accumAffineVector(12 * nodeIdx + 5) * vec[1];
                        //rhs[2] += accumAffineVector(12 * nodeIdx + 8) * vec[2];
                        //rhs[2] += weight * (n[2] + accumAffineVector(12 * nodeIdx + 11));
                    }

                    rhs[0] += mat_U0_(0, vertexIdx);
                    rhs[1] += mat_U0_(1, vertexIdx);
                    rhs[2] += mat_U0_(2, vertexIdx);

                    bdVector.push_back(rhs[0]);
                    bdVector.push_back(rhs[1]);
                    bdVector.push_back(rhs[2]);

                    dOffset += 3;
                }

                Ad.resize(dOffset, 12 * anchorSampler.nodeSize());
                Ad.setFromTriplets(AdList.begin(), AdList.end());

                bd.resize(dOffset);
                for (int i = 0; i < dOffset; ++i)
                {
                    bd(i) = bdVector[i];
                }

                A1 = g_alphaData * Eigen::SparseMatrix<double>(Ad.transpose()) * Ad;
                b1 = g_alphaData * Ad.transpose() * bd;
            }

            //	Third independent section
            //	Construct Ak, bk, A2 and b2
#pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct Ak and bk
                //------------------------------------------------------------------------
                std::vector<double> bkVector;
                size_t kOffset = 0;

                for (int idx0 = 0; idx0 < anchorSampler.nodeSize(); ++idx0)
                {
                    size_t nodeIdx0 = anchorSampler.getNodeVertexIdx(idx0);
                    VertexHandle vh0 = anchorMesh.vertex_handle(nodeIdx0);
                    Point v0 = anchorMesh.point(vh0);

                    svr::neighborIter nIter = anchorSampler.getNodeNodeIter(idx0);
                    for (; nIter.is_valid(); ++nIter)
                    {
                        size_t idx1 = nIter.getIndex();
                        size_t nodeIdx1 = anchorSampler.getNodeVertexIdx(idx1);
                        VertexHandle vh1 = anchorMesh.vertex_handle(nodeIdx1);
                        Point v1 = anchorMesh.point(vh1);

                        const Eigen::Vector3d &k = tempK[idx0].at(idx1).second;

                        double vec[3] = { 0.0 };
                        vec[0] = v1[0] - v0[0];
                        vec[1] = v1[1] - v0[1];
                        vec[2] = v1[2] - v0[2];

                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 0, vec[0]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 3, vec[1]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 6, vec[2]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 9, 1.0));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx1 + 9, -1.0));
                        ++kOffset;

                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 1, vec[0]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 4, vec[1]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 7, vec[2]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 10, 1.0));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx1 + 10, -1.0));
                        ++kOffset;

                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 2, vec[0]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 5, vec[1]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 8, vec[2]));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx0 + 11, 1.0));
                        AkList.push_back(Eigen::Triplet<double>(kOffset, 12 * idx1 + 11, -1.0));
                        ++kOffset;

                        bkVector.push_back(vec[0] + k(0));
                        bkVector.push_back(vec[1] + k(1));
                        bkVector.push_back(vec[2] + k(2));
                    }
                }

                Ak.resize(kOffset, 12 * anchorSampler.nodeSize());
                Ak.setFromTriplets(AkList.begin(), AkList.end());

                bk.resize(bkVector.size());
                for (int i = 0; i < bkVector.size(); ++i)
                {
                    bk(i) = bkVector[i];
                }

                A2 = beta * Eigen::SparseMatrix<double>(Ak.transpose()) * Ak;
                b2 = beta * Ak.transpose() * bk;
            }
            }

            //------------------------------------------------------------------------
            //	Construct linear equation A * x = b
            //	A = g_alphaRigid * JrT * Jr +
            //		g_alphaData * AdT * Ad +
            //		beta * AkT * Ak
            //
            //	b = - g_alphaRigid * JrT * fr +
            //		  g_alphaData * AdT * bd +
            //		  beta * AkT * bk +
            //		  g_alphaRigid * JrT * Jr * affineVector_L0
            //
            //	x = nextAffineVector_L0
            //------------------------------------------------------------------------
            Eigen::SparseMatrix<double> A = A0 + A1 + A2;
            Eigen::VectorXd b = b0 + b1 + b2;

            //Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> CholmodSolver;

            CholmodSolver.compute(A);
            Eigen::VectorXd nextAffineVector_L0 = CholmodSolver.solve(b);

            //------------------------------------------------------------------------
            //	Calculate total energy energyGaussNewton
            //	Break the loop if relative energy change < 1e-6f,
            //	or total iteration is greater than 50
            //------------------------------------------------------------------------
            Eigen::VectorXd RigidVector(12 * anchorSampler.nodeSize());
            Eigen::VectorXd DataVector(12 * anchorSampler.nodeSize());
            Eigen::VectorXd L0Vector(12 * anchorSampler.nodeSize());

            double RigidEnergy = 0.0;
            double DataEnergy = 0.0;
            double L0Energy = 0.0;

            //	Parallelize energy calculation
#pragma omp sections
            {
                //	First section calculates rigid energy
#pragma omp section
            {
                RigidVector = fr + Jr * (nextAffineVector_L0 - affineVector_L0.cast<double>());
                RigidEnergy = g_alphaRigid * RigidVector.transpose() * RigidVector;
            }
            //	Second section calculates data energy
#pragma omp section
            {
                DataVector = Ad * nextAffineVector_L0 - bd;
                DataEnergy = g_alphaData* DataVector.transpose() * DataVector;
            }
            //	Third section calculates L0 energy
#pragma omp section
            {
                L0Vector = Ak * nextAffineVector_L0 - bk;
                L0Energy = beta * L0Vector.transpose() * L0Vector;
            }
            }

            //	Sum up all three energy
            gaussNewtonEnergy = RigidEnergy + DataEnergy + L0Energy;

            //	Break Gauss-Newton loop if relative energy change is less than 1e-6f,
            //	or total iteration is greater than 50
            double gaussNewtonEnergyChange = fabs(gaussNewtonEnergy - gaussNewtonEnergyPrev) / (gaussNewtonEnergyPrev + 1.0);
            //std::cout << "GN iteration " << gaussNewtonIterCnt << ", "
            //	<< "gnEnergy = " << gaussNewtonEnergy << ", gnEnergyChange = " << gaussNewtonEnergyChange << std::endl;

            affineVector_L0 = nextAffineVector_L0.cast<double>();

            if (gaussNewtonEnergyChange < 1e-6 || gaussNewtonIterCnt >= 50) break;
            gaussNewtonEnergyPrev = gaussNewtonEnergy;

            ++gaussNewtonIterCnt;
        }

        return gaussNewtonEnergy;
    }
