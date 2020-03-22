#include "Registration.h"
#include <fstream>
#include <iostream>
#include "shot.h"
#include "median.h"
#include <opencv2/opencv.hpp>
#include "Prune.h"
#include <memory>
#if (__cplusplus >= 201402L) || (defined(_MSC_VER) && _MSC_VER >= 1800)
#define MAKE_UNIQUE std::make_unique
#else
#define MAKE_UNIQUE company::make_unique
#endif


Registration::Registration() {
    target_tree = NULL;
    src_mesh_ = NULL;
    tar_mesh_ = NULL;
    use_cholesky_solver_ = true;
    use_pardiso_ = true;
    update_tarGeotree = true;
    //geo_kdtree = NULL;
    //pgeodesic_ = NULL;
};

Registration::~Registration()
{
    if (target_tree != NULL)
    {
        delete target_tree;
        target_tree = NULL;
    }
}

// initialize before rigid transformation
void Registration::rigid_init(Mesh & src_mesh, Mesh & tar_mesh, RegParas& paras)
{
    src_mesh_ = new Mesh;
    tar_mesh_ = new Mesh;
    src_mesh_ = &src_mesh;
    tar_mesh_ = &tar_mesh;
    pars_ = paras;
    n_src_vertex_ = src_mesh_->n_vertices();
    n_tar_vertex_ = tar_mesh_->n_vertices();
    if(pars_.pruning_type == DIFFUSION)
        hgeodesic_ = new hgeodesic(src_mesh, tar_mesh);
    else
         hgeodesic_ = new hgeodesic(src_mesh);

    corres_pair_ids_.resize(n_src_vertex_);

    // construct kd Tree
    tar_points_.resize(3, n_tar_vertex_);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n_tar_vertex_; i++)
    {
        tar_points_(0, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[0];
        tar_points_(1, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[1];
        tar_points_(2, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[2];
    }
    target_tree = new KDtree(tar_points_);
    //FindClosestPoints(correspondence_pairs_);
    InitCorrespondence(correspondence_pairs_);
}


void Registration::nonrigid_init()
{
    mat_U0_.resize(3, n_src_vertex_);

    // welsch parameters
    weight_d_.resize(n_src_vertex_);
    weight_s_.resize(src_mesh_->n_halfedges());

    // Initialize correspondences
    InitCorrespondence(correspondence_pairs_);

   //double init_nu = 0.0;
    Eigen::VectorXd init_nus(correspondence_pairs_.size());
    for(int i = 0; i < correspondence_pairs_.size(); i++)
    {
        init_nus[i] = (src_mesh_->point(src_mesh_->vertex_handle(correspondence_pairs_[i].first))
                    - tar_mesh_->point(tar_mesh_->vertex_handle(correspondence_pairs_[i].second))).norm();
    }
    igl::median(init_nus, pars_.Data_nu);
     //pars_.Data_nu = init_nus.sum() / correspondence_pairs_.size();

    if(pars_.calc_gt_err&&n_src_vertex_ == n_tar_vertex_)
    {
        Eigen::VectorXd gt_err(n_src_vertex_);
        for(int i = 0; i < n_src_vertex_; i++)
        {
            gt_err[i] = (src_mesh_->point(src_mesh_->vertex_handle(i)) - tar_mesh_->point(tar_mesh_->vertex_handle(i))).norm();
        }
        pars_.init_gt_mean_errs = std::sqrt(gt_err.squaredNorm()/n_src_vertex_);
        pars_.init_gt_max_errs = gt_err.maxCoeff();
    }
}


// Rigid Registration
double Registration::DoRigid()
{
    Eigen::Matrix3Xd rig_tar_v = Eigen::Matrix3Xd::Zero(3, n_src_vertex_);
    Eigen::Matrix3Xd rig_src_v = Eigen::Matrix3Xd::Zero(3, n_src_vertex_);
    Eigen::Affine3d old_T;

    //std::cout << "before for loop done !" << std::endl;
	corres_pair_ids_.setZero();
    for(int iter = 0; iter < pars_.rigid_iters; iter++)
    {
        for (int i = 0; i < correspondence_pairs_.size(); i++)
        {
            rig_src_v.col(i) = Eigen::Map<Eigen::Vector3d>(src_mesh_->point(src_mesh_->vertex_handle(correspondence_pairs_[i].first)).data(), 3, 1);
            rig_tar_v.col(i) = Eigen::Map<Eigen::Vector3d>(tar_mesh_->point(tar_mesh_->vertex_handle(correspondence_pairs_[i].second)).data(), 3, 1);
            corres_pair_ids_[correspondence_pairs_[i].first] = 1;
        }
//        if(pars_.calc_gt_err && tar_points_.cols() == n_src_vertex_)
//        {
//            std::cout << "iters = " << iter << " gt_err = " << (rig_src_v - tar_points_).norm()/sqrt(n_src_vertex_);
//        }
//        std::cout << "Error = " << (rig_src_v - rig_tar_v).norm()/ sqrt(n_src_vertex_) << std::endl;

        //std::cout << "rigid_src_v = " << rig_src_v << " \nrigid_tar_v = " << rig_tar_v << std::endl;
        old_T = rigid_T_;
        rigid_T_ = point_to_point(rig_src_v, rig_tar_v, corres_pair_ids_);
       // std::cout << "rigid_T = " << rigid_T_.matrix() << std::endl;

        if((old_T.matrix() - rigid_T_.matrix()).norm() < 1e-3)
        {
            break;
        }

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
        for (int i = 0; i < n_src_vertex_; i++)
        {
            OpenMesh::Vec3d p = src_mesh_->point(src_mesh_->vertex_handle(i));
            //Eigen::Map<Eigen::Vector3d>(p.data(), 3)
            Eigen::Vector3d temp = rigid_T_ * Eigen::Map<Eigen::Vector3d>(p.data(), 3);
            p[0] = temp[0];
            p[1] = temp[1];
            p[2] = temp[2];
            src_mesh_->set_point(src_mesh_->vertex_handle(i), p);
        }

        // Find correspondence
        FindClosestPoints(correspondence_pairs_);
        //std::cout << "Find closest points done !" << std::endl;
        SimplePruning(correspondence_pairs_, pars_.use_distance_reject, pars_.use_normal_reject);
        //std::cout << "simple pruning done !" << std::endl;
    }

    return 0;
}


void Registration::FindClosestPoints(VPairs & corres)
{
    corres.resize(n_src_vertex_);
    //src_mesh_->update_vertex_normals();
    Eigen::VectorXd gt_errs = Eigen::VectorXd::Zero(n_src_vertex_);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n_src_vertex_; i++)
    {
        double mini_dist;
        int idx = target_tree->closest(src_mesh_->point(src_mesh_->vertex_handle(i)).data(), mini_dist);
//        OpenMesh::Vec3d src_normal = src_mesh_->normal(src_mesh_->vertex_handle(i));
//        OpenMesh::Vec3d tar_normal = tar_mesh_->normal(tar_mesh_->vertex_handle(idx));
//        double angle = acos(src_normal | tar_normal / (src_normal.norm()*tar_normal.norm()));
//        src_mesh_->data(src_mesh_->vertex_handle(i)).dist_to_corres = mini_dist;
//        src_mesh_->data(src_mesh_->vertex_handle(i)).angle_to_corres = abs(angle);
        std::pair<int, int> pair(i, idx);
        corres[i] = pair;
        gt_errs[i] = mini_dist*mini_dist;
    }
    //std::cout << "out test closest distance = " << sqrt(gt_errs.sum()/n_src_vertex_) << std::endl;
}

void Registration::SimplePruning(VPairs & corres, bool use_distance = true, bool use_normal = true) const
{
    // Distance and normal
    Eigen::VectorXd tar_min_dists(n_tar_vertex_);
    tar_min_dists.setConstant(1e10);
    Eigen::VectorXi min_idxs(n_tar_vertex_);
    min_idxs.setConstant(-1);
    src_mesh_->update_vertex_normals();

    Eigen::VectorXd corres_idx = Eigen::VectorXd::Zero(n_src_vertex_);
    for(int i = 0; i < corres.size(); i++)
    {
        double dist = (src_mesh_->point(src_mesh_->vertex_handle(corres[i].first))
                       - tar_mesh_->point(tar_mesh_->vertex_handle(corres[i].second))).norm();
        OpenMesh::Vec3d src_normal = src_mesh_->normal(src_mesh_->vertex_handle(corres[i].first));
        OpenMesh::Vec3d tar_normal = tar_mesh_->normal(tar_mesh_->vertex_handle(corres[i].second));
        double angle = acos(src_normal | tar_normal / (src_normal.norm()*tar_normal.norm()));
        if((!use_distance || dist < pars_.distance_threshold)
            && (!use_normal || src_mesh_->n_faces() == 0 || angle < pars_.normal_threshold))
        {
            corres_idx[i] = 1;
        }

//        if (src_mesh_->data(src_mesh_->vertex_handle(corres[i].first)).dist_to_corres < tar_min_dists[corres[i].second])
//        {
//            min_idxs[corres[i].second] = corres[i].first;
//            tar_min_dists[corres[i].second] = src_mesh_->data(src_mesh_->vertex_handle(corres[i].first)).dist_to_corres;
//        }
    }

    // picky
    /*for (int i = 0; i < corres.size(); i++)
    {
        if (min_idxs[corres[i].second] == corres[i].first)
        {
            corres_pair_ids_[corres[i].first] = 1;
        }
        else
        {
            corres_pair_ids_[corres[i].first] = 0;
        }
    }*/
    VPairs corres2;
    for (auto it = corres.begin(); it != corres.end(); it++)
    {
        if (corres_idx[(*it).first] == 1)
        {
            corres2.push_back(*it);
        }
    }
    corres.clear();
    corres = corres2;
}

void Registration::LandMarkCorres(VPairs & corres)
{
    corres.clear();
    if (pars_.landmark_src.size() != pars_.landmark_tar.size())
    {
        std::cout << "Error: landmark data wrong!!" << std::endl;
    }
    n_landmark_nodes_ = pars_.landmark_tar.size();
    for (int i = 0; i < n_landmark_nodes_; i++)
    {
        //corres_pair_ids_[pars_.landmark_src[i]] = 1;
        std::pair<int, int> pair(pars_.landmark_src[i], pars_.landmark_tar[i]);
        corres.push_back(pair);
		if (pair.first > n_src_vertex_ || pair.first < 0)
			std::cout << "Error: source index in Landmark is out of range!" << std::endl;
		if (pair.second > n_src_vertex_ || pair.second < 0)
			std::cout << "Error: target index in Landmark is out of range!" << std::endl;
	}
    std::cout << " use landmark and landmark is ... " << pars_.landmark_src.size() << std::endl;
}

void Registration::find_good_init()
{
    InitCorrespondence(correspondence_pairs_);
}


void Registration::InitCorrespondence(VPairs & corres)
{
    switch(pars_.corres_type)
    {
    case SHOT:
    {
        FindSHOTCorrespondence(corres);
        if(!corres.size())
        {
            FindClosestPoints(corres);
        }
        break;
    }
    case CLOSEST:
    {
        FindClosestPoints(corres);
        break;
    }
    case LANDMARK:
    {
        LandMarkCorres(corres);
        break;
    }
    }


    // test correct points
    int nn = 0;
    for (int i = 0; i < corres.size(); i++)
    {
        if (corres[i].first == corres[i].second)
        {
            nn++;
        }
    }
    pars_.n_correct_shot_pair = nn;
    pars_.n_total_shot_pair = corres.size();

    if(pars_.print_each_step_info)
    {
        std::string shot_file = pars_.out_each_step_info + "shot_corres.txt";
        std::ofstream out_corres1(shot_file);
        for (int i = 0; i < corres.size(); i++)
        {
            out_corres1 << corres[i].first << " " << corres[i].second << std::endl;
        }
        out_corres1.close();
    }


    switch(pars_.pruning_type)
    {
    case DIFFUSION:
    {
        Prune prune;
        double prune_t = omp_get_wtime();
        prune.Initialize(hgeodesic_, corres,  pars_.diffusion_para);
        prune.RunPruning(corres);
        double prune_t2 = omp_get_wtime();
        std::cout << "pruning time = " << prune_t2 - prune_t << std::endl;
        if(!corres.size())
        {
           FindClosestPoints(corres);
           SimplePruning(corres, pars_.use_distance_reject, pars_.use_normal_reject);
        }
        break;
    }
    case SIMPLE:
    {
        SimplePruning(corres);
        break;
    }
    default:
    {
        break;
    }
    }


    //test corres
    nn = 0;
    for (int i = 0; i < corres.size(); i++)
    {
        if (corres[i].first == corres[i].second)
        {
            nn++;
        }
    }
    //std::cout << "after pruning !!, the correct pair is " << nn << " / " << corres.size() << std::endl;
    pars_.n_correct_prune_pair = nn;
    pars_.n_total_prune_pair = corres.size();

    if(pars_.print_each_step_info)
    {
        std::string file_diffprune = pars_.out_each_step_info + "prune_corres.txt";
        std::ofstream out_corres2(file_diffprune);
        for (int i = 0; i < corres.size(); i++)
        {
            out_corres2 << corres[i].first << " " << corres[i].second << std::endl;
        }
        out_corres2.close();
    }

    // end test output
}



void Registration::FindSHOTCorrespondence(VPairs & correspondence_pairs)
{
    CmdLineParams params;
    cv::Mat src_feat, tar_feat;
    int src_nActualFeat, tar_nActualFeat;
    Eigen::Matrix3Xd src_points(3,n_src_vertex_);
    for (auto it = src_mesh_->vertices_begin(); it != src_mesh_->vertices_end(); it++)
    {
        src_points(0, (*it).idx()) = src_mesh_->point(*it)[0];
        src_points(1, (*it).idx()) = src_mesh_->point(*it)[1];
        src_points(2, (*it).idx()) = src_mesh_->point(*it)[2];
    }
    KDtree *src_tree = new KDtree(src_points);
    src_nActualFeat = CalcSHOTfeature(src_mesh_, src_tree, n_src_vertex_, src_feat);
    tar_nActualFeat = CalcSHOTfeature(tar_mesh_, target_tree, n_tar_vertex_, tar_feat);

    cv::flann::Index kdtree(src_feat, cv::flann::KDTreeIndexParams());
    std::vector<float> dists;
    dists.resize(2);
    std::vector<int> knn;
    knn.resize(2);
    int correctMatch = 0, totalMatch = 0;

    correspondence_pairs.clear();
    for (int i = 0; i < tar_nActualFeat; i++)
    {
 //       std::vector<float> query;
//        for (int j = 0; j < tar_feat.cols(); j++)
//            query.push_back(tar_feat.at(i, j));

        kdtree.knnSearch(tar_feat.row(i), knn, dists, 2, cv::flann::SearchParams());
        //std::cout << "dists[0] = " << dists[0] << std::endl;

        assert(dists[0] <= dists[1]);

        if (dists[0] <= params.matchTh * params.matchTh * dists[1])
        {
            //printf("Match %d: %d\n", i, knn[0]);
            std::pair<int, int> pair(knn[0], i);
            correspondence_pairs.push_back(pair);
            if (i == knn[0])
            {
                correctMatch++;
            }

            totalMatch++;
        }

    }
    std::cout << "correct match is " << correctMatch << " / " << totalMatch << std::endl;
    delete src_tree;
}

int Registration::CalcSHOTfeature(Mesh* mesh, KDtree* kdtree, int nFeats, cv::Mat& Features)
{
    CmdLineParams params;
    double meshRes = CalcEdgelength(mesh, 1);
    std::cout << "meshRes = " << meshRes << std::endl;

    std::cout << "radius = " << meshRes * params.radiusMR << std::endl;
    Random3DDetector detector(nFeats, false, meshRes * params.radiusMR, params.minNeighbors);

    Feature3D* feat;
    int nActualFeat = detector.extract(mesh, kdtree, feat);

    // set shot params;
    SHOTParams shotParams;
    shotParams.radius = meshRes * params.radiusMR;
    shotParams.localRFradius = meshRes * params.radiusMR;
    shotParams.minNeighbors = params.minNeighbors;
    shotParams.shapeBins = params.shotShapeBins;
    shotParams.colorBins = params.shotColorBins;
    shotParams.describeColor = params.describeColor;
    shotParams.describeShape = params.describeShape;
    shotParams.nThreads = params.nThreads;
    shotParams.describeColor = false;
    shotParams.shapeBins = 10;

    SHOTDescriptor descriptor(shotParams);

    double** Desc;
    descriptor.describe(mesh, kdtree, feat, Desc, nActualFeat);
    std::cout << "src describe done!" << std::endl;
    cv::Mat features(nFeats, descriptor.getDescriptorLength(), CV_32FC1);
    for (int i = 0; i < nActualFeat; i++)
    {
        for (int j = 0; j < descriptor.getDescriptorLength(); j++)
        {
            features.at<float>(i, j) = Desc[i][j];
        }
    }
    Features = features;
    return nActualFeat;
}


// *type: 0 :median, 1: average
double Registration::CalcEdgelength(Mesh* mesh, int type)
{
    double med;
    if(mesh->n_faces()!=0)
    {
        Eigen::VectorXd edges_length(mesh->n_edges());
        for(int i = 0; i < mesh->n_edges();i++)
        {
            OpenMesh::VertexHandle vi = mesh->from_vertex_handle(mesh->halfedge_handle(mesh->edge_handle(i),0));
            OpenMesh::VertexHandle vj = mesh->to_vertex_handle(mesh->halfedge_handle(mesh->edge_handle(i),0));
            edges_length[i] = (mesh->point(vi) - mesh->point(vj)).norm();
        }
        if (type == 0)
            igl::median(edges_length, med);
        else
            med = edges_length.mean();
    }
    else
    {
        // source is mesh, target may be point cloud.
        Eigen::VectorXd edges_length(mesh->n_vertices());
        int nk = 7;
        for(int i = 0; i<mesh->n_vertices(); i++)
        {
            int* id = new int[nk];
            double *dist = new double[nk];
            target_tree->query(mesh->point(mesh->vertex_handle(i)).data(), nk, id, dist);
            Eigen::VectorXd k_dist = Eigen::Map<Eigen::VectorXd>(dist, nk);
            if (type == 0)
                igl::median(k_dist.tail(nk - 1), edges_length[i]);
            else
                edges_length[i] = k_dist.tail(nk - 1).mean();
            delete[]id;
            delete[]dist;
        }
        if (type == 0)
            igl::median(edges_length, med);
        else
            med = edges_length.mean();
    }
    return med;
}

/// @param Source (one 3D point per column)
/// @param Target (one 3D point per column)
/// @param Confidence weights
template <typename Derived1, typename Derived2, typename Derived3>
Eigen::Affine3d Registration::point_to_point(Eigen::MatrixBase<Derived1>& X,
    Eigen::MatrixBase<Derived2>& Y,
    const Eigen::MatrixBase<Derived3>& w) {
    /// Normalize weight vector
    Eigen::VectorXd w_normalized = w / w.sum();
    /// De-mean
    Eigen::Vector3d X_mean, Y_mean;
    for (int i = 0; i<3; ++i) {
        X_mean(i) = (X.row(i).array()*w_normalized.transpose().array()).sum();
        Y_mean(i) = (Y.row(i).array()*w_normalized.transpose().array()).sum();
    }
    X.colwise() -= X_mean;
    Y.colwise() -= Y_mean;

    /// Compute transformation
    Eigen::Affine3d transformation;
    Eigen::Matrix3d sigma = X * w_normalized.asDiagonal() * Y.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
        Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
        transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
    }
    else {
        transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();
    }
    transformation.translation().noalias() = Y_mean - transformation.linear()*X_mean;
    /// Re-apply mean
    X.colwise() += X_mean;
    Y.colwise() += Y_mean;
    /// Return transformation
    return transformation;
}
