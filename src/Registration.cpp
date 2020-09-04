#include "Registration.h"
#include <igl/median.h>


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
    corres_pair_ids_.resize(n_src_vertex_);

    tar_points_.resize(3, n_tar_vertex_);
    #pragma omp parallel for
    for (int i = 0; i < n_tar_vertex_; i++)
    {
        tar_points_(0, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[0];
        tar_points_(1, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[1];
        tar_points_(2, i) = tar_mesh_->point(tar_mesh_->vertex_handle(i))[2];
    }

        // construct kd Tree
    target_tree = new KDtree(tar_points_);
    InitCorrespondence(correspondence_pairs_);
}

/// Find self edge median of point cloud
template<typename Derived1>
Scalar Registration::FindKnearestMed(Eigen::MatrixBase<Derived1>& X, int nk)
{
    nanoflann::KDTreeAdaptor<Eigen::MatrixBase<Derived1>, 3, nanoflann::metric_L2_Simple> kdtree(X);
    VectorX X_nearest(X.cols());
#pragma omp parallel for
    for (int i = 0; i<X.cols(); i++)
    {
        int* id = new int[nk];
        Scalar *dist = new Scalar[nk];
        kdtree.query(X.col(i).data(), nk, id, dist);
        VectorX k_dist = Eigen::Map<VectorX>(dist, nk);
        igl::median(k_dist.tail(nk - 1), X_nearest[i]);
        delete[]id;
        delete[]dist;
    }
    Scalar med;
    igl::median(X_nearest, med);
    return med;
}


void Registration::nonrigid_init()
{
    mat_U0_.resize(3, n_src_vertex_);

    // welsch parameters
    weight_d_.resize(n_src_vertex_);
    weight_s_.resize(src_mesh_->n_halfedges());

    // Initialize correspondences
    InitCorrespondence(correspondence_pairs_);

    VectorX init_nus(correspondence_pairs_.size());
    for(size_t i = 0; i < correspondence_pairs_.size(); i++)
    {
        Vector3 closet = correspondence_pairs_[i].position;
        init_nus[i] = (src_mesh_->point(src_mesh_->vertex_handle(correspondence_pairs_[i].src_idx))
                    - Vec3(closet[0], closet[1], closet[2])).norm();
    }
    igl::median(init_nus, pars_.Data_nu);

    if(pars_.calc_gt_err&&n_src_vertex_ == n_tar_vertex_)
    {
        VectorX gt_err(n_src_vertex_);
        for(int i = 0; i < n_src_vertex_; i++)
        {
            gt_err[i] = (src_mesh_->point(src_mesh_->vertex_handle(i)) - tar_mesh_->point(tar_mesh_->vertex_handle(i))).norm();
        }
        pars_.init_gt_mean_errs = std::sqrt(gt_err.squaredNorm()/n_src_vertex_);
        pars_.init_gt_max_errs = gt_err.maxCoeff();
    }
}

// Rigid Registration
Scalar Registration::DoRigid()
{
    Matrix3X rig_tar_v = Matrix3X::Zero(3, n_src_vertex_);
    Matrix3X rig_src_v = Matrix3X::Zero(3, n_src_vertex_);
    Affine3 old_T;

    //std::cout << "before for loop done !" << std::endl;
	corres_pair_ids_.setZero();
    for(int iter = 0; iter < pars_.rigid_iters; iter++)
    {
        for (size_t i = 0; i < correspondence_pairs_.size(); i++)
        {
            rig_src_v.col(i) = Eigen::Map<Vector3>(src_mesh_->point(src_mesh_->vertex_handle(correspondence_pairs_[i].src_idx)).data(), 3, 1);
            rig_tar_v.col(i) = correspondence_pairs_[i].position;
            corres_pair_ids_[correspondence_pairs_[i].src_idx] = 1;
        }
        old_T = rigid_T_;
        rigid_T_ = point_to_point(rig_src_v, rig_tar_v, corres_pair_ids_);

        if((old_T.matrix() - rigid_T_.matrix()).norm() < 1e-3)
        {
            break;
        }

        #pragma omp parallel for
        for (int i = 0; i < n_src_vertex_; i++)
        {
            Vec3 p = src_mesh_->point(src_mesh_->vertex_handle(i));
            Vector3 temp = rigid_T_ * Eigen::Map<Vector3>(p.data(), 3);
            p[0] = temp[0];
            p[1] = temp[1];
            p[2] = temp[2];
            src_mesh_->set_point(src_mesh_->vertex_handle(i), p);
        }

        // Find correspondence
        FindClosestPoints(correspondence_pairs_);
        SimplePruning(correspondence_pairs_, pars_.use_distance_reject, pars_.use_normal_reject);
    }

    return 0;
}

void Registration::FindClosestPoints(VPairs & corres)
{
    corres.resize(n_src_vertex_);

    #pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; i++)
    {
        Scalar mini_dist;
        int idx = target_tree->closest(src_mesh_->point(src_mesh_->vertex_handle(i)).data(), mini_dist);
        Closest c;
        c.src_idx = i;
        c.position = tar_points_.col(idx);
        c.normal = Vec2Eigen(tar_mesh_->normal(tar_mesh_->vertex_handle(idx)));
        corres[i] = c;
    }
}

void Registration::SimplePruning(VPairs & corres, bool use_distance = true, bool use_normal = true)
{
    // Distance and normal
    VectorX tar_min_dists(n_tar_vertex_);
    tar_min_dists.setConstant(1e10);
    Eigen::VectorXi min_idxs(n_tar_vertex_);
    min_idxs.setConstant(-1);
    src_mesh_->update_vertex_normals();

    VectorX corres_idx = VectorX::Zero(n_src_vertex_);
    for(size_t i = 0; i < corres.size(); i++)
    {
        Vector3 closet = corres[i].position;
        Scalar dist = (src_mesh_->point(src_mesh_->vertex_handle(corres[i].src_idx))
                       - Eigen2Vec(closet)).norm();
        Vec3 src_normal = src_mesh_->normal(src_mesh_->vertex_handle(corres[i].src_idx));
        Vec3 tar_normal = Eigen2Vec(corres[i].normal);

        Scalar angle = acos(src_normal | tar_normal / (src_normal.norm()*tar_normal.norm()));
        if((!use_distance || dist < pars_.distance_threshold)
            && (!use_normal || src_mesh_->n_faces() == 0 || angle < pars_.normal_threshold))
        {

            corres_idx[i] = 1;
        }
    }
    if(pars_.use_fixedvex)
    {
        for(size_t i = 0; i < pars_.fixed_vertices.size(); i++)
        {
            int idx = pars_.fixed_vertices[i];
            corres_idx[idx] = 0;
        }
    }
    VPairs corres2;
    for (auto it = corres.begin(); it != corres.end(); it++)
    {
        if (corres_idx[(*it).src_idx] == 1)
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
        Closest c;
        c.src_idx = pars_.landmark_src[i];
        OpenMesh::VertexHandle vh = tar_mesh_->vertex_handle(pars_.landmark_tar[i]);

        if (c.src_idx > n_src_vertex_ || c.src_idx < 0)
            std::cout << "Error: source index in Landmark is out of range!" << std::endl;
        if (vh.idx() < 0)
            std::cout << "Error: target index in Landmark is out of range!" << std::endl;

        c.position = Vec2Eigen(tar_mesh_->point(vh));
        c.normal = Vec2Eigen(tar_mesh_->normal(vh));
        corres.push_back(c);
	}
    std::cout << " use landmark and landmark is ... " << pars_.landmark_src.size() << std::endl;
}

void Registration::InitCorrespondence(VPairs & corres)
{
    FindClosestPoints(corres);
    SimplePruning(corres);

    if(pars_.use_landmark)
    {
        corres.clear();
        for(int i = 0; i < pars_.landmark_src.size(); i++)
        {
            Closest c;
            c.src_idx = pars_.landmark_src[i];
            c.tar_idx = pars_.landmark_tar[i];
            c.position = tar_points_.col(c.tar_idx);
            c.normal = Vec2Eigen(tar_mesh_->normal(tar_mesh_->vertex_handle(c.tar_idx)));
            corres.push_back(c);
        }
    }
}

// *type: 0 :median, 1: average
Scalar Registration::CalcEdgelength(Mesh* mesh, int type)
{
    Scalar med;
    if(mesh->n_faces() > 0)
    {
        VectorX edges_length(mesh->n_edges());
        for(size_t i = 0; i < mesh->n_edges();i++)
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
        VectorX edges_length(mesh->n_vertices());
        int nk = 7;
        for(size_t i = 0; i<mesh->n_vertices(); i++)
        {
            int* id = new int[nk];
            Scalar *dist = new Scalar[nk];
            target_tree->query(mesh->point(mesh->vertex_handle(i)).data(), nk, id, dist);
            VectorX k_dist = Eigen::Map<VectorX>(dist, nk);
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
Affine3 Registration::point_to_point(Eigen::MatrixBase<Derived1>& X,
    Eigen::MatrixBase<Derived2>& Y,
    const Eigen::MatrixBase<Derived3>& w) {
    /// Normalize weight vector
    VectorX w_normalized = w / w.sum();
    /// De-mean
    Vector3 X_mean, Y_mean;
    for (int i = 0; i<3; ++i) {
        X_mean(i) = (X.row(i).array()*w_normalized.transpose().array()).sum();
        Y_mean(i) = (Y.row(i).array()*w_normalized.transpose().array()).sum();
    }
    X.colwise() -= X_mean;
    Y.colwise() -= Y_mean;

    /// Compute transformation
    Affine3 transformation;
    Matrix33 sigma = X * w_normalized.asDiagonal() * Y.transpose();
    Eigen::JacobiSVD<Matrix33> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
        Vector3 S = Vector3::Ones(); S(2) = -1.0;
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
