#include "NonrigidICP.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

NonrigidICP::NonrigidICP()
{

}

NonrigidICP::~NonrigidICP()
{

}

void NonrigidICP::Initialize()
{
    nonrigid_init();
    alpha = pars_.alpha;
    gamma = 1.0;

    alpha_M_G.clear();
    for (int i = 0; i < src_mesh_->n_edges(); ++i)
    {
        int a = src_mesh_->from_vertex_handle(src_mesh_->halfedge_handle(src_mesh_->edge_handle(i),0)).idx();
        int b = src_mesh_->to_vertex_handle(src_mesh_->halfedge_handle(src_mesh_->edge_handle(i),0)).idx();

        for (int j = 0; j < 3; j++) alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + j, a*4 + j, alpha));
        alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + 3, a*4 + 3, alpha * gamma));

        for (int j = 0; j < 3; j++) alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + j, b*4 + j, -alpha));
        alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + 3, b*4 + 3, -alpha * gamma));
    }
    pars_.alpha = pars_.alpha / src_mesh_->n_edges() * n_src_vertex_;
}

double NonrigidICP::DoNonRigid()
{
    int n = n_src_vertex_;
    int m = src_mesh_->n_edges();

    Timer time;
    Timer::EventID begin_time, run_time;
    pars_.each_energys.clear();
    pars_.each_gt_mean_errs.clear();
    pars_.each_gt_max_errs.clear();
    pars_.each_times.clear();
    pars_.each_iters.clear();
    Eigen::MatrixX3d prev_X = Eigen::MatrixX3d::Zero(4*n, 3);
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;

    // print initial run results
    pars_.each_energys.push_back(0.0);
    pars_.each_gt_max_errs.push_back(pars_.init_gt_max_errs);
    pars_.each_gt_mean_errs.push_back(pars_.init_gt_mean_errs);
    pars_.each_iters.push_back(0);
    pars_.each_times.push_back(pars_.non_rigid_init_time);

    bool run_once = true;
    begin_time = time.get_time();
    for(int out_iter = 0; out_iter < pars_.max_outer_iters; out_iter++)
    {
        // according correspondence_pairs to update mat_U0_;
        mat_U0_.setZero();
        corres_pair_ids_.setZero();
        Eigen::SparseMatrix<double> A(4*m + n, 4*n);
        std::vector< Eigen::Triplet<double> > W_D(4*n);
        Eigen::MatrixX3d B = Eigen::MatrixX3d::Zero(4*m + n, 3);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < correspondence_pairs_.size(); i++)
        {
            mat_U0_.col(correspondence_pairs_[i].first) = tar_points_.col(correspondence_pairs_[i].second);
            corres_pair_ids_[correspondence_pairs_[i].first] = 1;
        }


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                W_D[i*4+j] = (Eigen::Triplet<double>(4*m + i, i*4 + j, corres_pair_ids_[i]
                                                     * src_mesh_->point(src_mesh_->vertex_handle(i))[j]));
                B(4*m + i, j) = corres_pair_ids_[i] * mat_U0_(j, i);
            }
            W_D[i*4+3] = (Eigen::Triplet<double>(4*m + i, i*4 + 3, corres_pair_ids_[i]));
        }

        std::vector< Eigen::Triplet<double> > _A = alpha_M_G;
        _A.insert(_A.end(), W_D.begin(), W_D.end());
        A.setFromTriplets(_A.begin(), _A.end());


        Eigen::SparseMatrix<double> ATA = Eigen::SparseMatrix<double>(A.transpose()) * A;
        Eigen::MatrixX3d ATB = Eigen::SparseMatrix<double>(A.transpose()) * B;

        if(run_once)
        {
            solver.analyzePattern(ATA);
            run_once = false;
        }
        solver.factorize(ATA);
        if (solver.info()!=Eigen::Success)
        {
            std::cerr << "Decomposition failed" << std::endl;
            return 1;
        }

        Eigen::MatrixX3d X = solver.solve(ATB);
        Eigen::Matrix3Xd XT = X.transpose();
        Eigen::VectorXd  gt_errs = Eigen::VectorXd(n);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; ++i)
        {
            OpenMesh::Vec3d p = src_mesh_->point(src_mesh_->vertex_handle(i));
            Eigen::Vector4d point(p[0], p[1], p[2], 1.0);
            Eigen::Vector3d point_transformed = XT.block<3, 4>(0, 4*i) * point;
            p[0] = point_transformed[0];
            p[1] = point_transformed[1];
            p[2] = point_transformed[2];
            src_mesh_->set_point(src_mesh_->vertex_handle(i), p);
            if(pars_.calc_gt_err)
                gt_errs[i] = (point_transformed - tar_points_.col(i)).squaredNorm();
        }

        run_time = time.get_time();
        double gt_err = std::sqrt(gt_errs.sum()/n);
        pars_.each_gt_mean_errs.push_back(gt_err);
        pars_.each_gt_max_errs.push_back(sqrt(gt_errs.maxCoeff()));
        double energy = ( A * X - B).squaredNorm();
        pars_.each_energys.push_back(energy);
        double eps_time = time.elapsed_time(begin_time, run_time);
        pars_.each_times.push_back(eps_time);
        pars_.each_iters.push_back(1);

        if(pars_.print_each_step_info)
        {
            std::string out_obj = pars_.out_each_step_info + "/iter_obj_" + std::to_string(out_iter) + ".ply";
            OpenMesh::IO::write_mesh( *src_mesh_, out_obj.c_str());
            std::cout << "iter = " << out_iter << " time = " << eps_time << " energy = " << energy << " gt_err = " << gt_err << std::endl;
        }

        FindClosestPoints(correspondence_pairs_);

        if((X - prev_X).norm() < pars_.stop)
        {
            break;
        }
        prev_X = X;
    }
    return 0.0;
}
