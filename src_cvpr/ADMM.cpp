#include "ADMM.h"
#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>

ADMM::ADMM() {
};

ADMM::~ADMM()
{
}

void ADMM::Initialize()
{
    nonrigid_init();
	// update matrix V, L, J
	std::vector<Triplet> coeffv, coeffL, coeffJ;
	MatrixXX block_matI(4, 3);
	block_matI.setIdentity();

	//initialize V, X, R, K
	mat_V0_.resize(n_src_vertex_, 4 * n_src_vertex_);
	mat_V_.resize(n_src_vertex_, 4 * n_src_vertex_);
	mat_U_.resize(3, n_src_vertex_);

	mat_X_.resize(4 * n_src_vertex_, 3); mat_X_.setZero();
    mat_R_.resize(3 * n_src_vertex_, 3); mat_R_.setZero();
	mat_L_.resize(4 * n_src_vertex_, 4 * n_src_vertex_);
	mat_J_.resize(4 * n_src_vertex_, 3 * n_src_vertex_);
	mat_A0_.resize(4 * n_src_vertex_, 4 * n_src_vertex_);
	mat_VU_.resize(4 * n_src_vertex_, 3);
	direction_.resize(4 * n_src_vertex_, 3);
	grad_X_.resize(4 * n_src_vertex_, 3); grad_X_.setZero();
	coeffv.reserve(4 * n_src_vertex_);
	coeffL.reserve(3 * n_src_vertex_);
	coeffJ.reserve(3 * n_src_vertex_);

	for (int i = 0; i < n_src_vertex_; i++)
	{
		// fill u and L
		for (int j = 0; j < 3; j++)
		{
			coeffv.push_back(Triplet(i, 4 * i + j, src_mesh_->point(src_mesh_->vertex_handle(i))[j]));
			coeffL.push_back(Triplet(4 * i + j, 4 * i + j, 1.0));
			coeffJ.push_back(Triplet(4 * i + j, 3 * i + j, 1.0));
		}
		coeffv.push_back(Triplet(i, 4 * i + 3, 1.0));
		// fill X
		mat_X_.block(4 * i, 0, 4, 3) = block_matI;
		// fill R
		mat_R_.block(3 * i, 0, 3, 3) = block_matI.block(0, 0, 3, 3);
	}
	mat_V0_.setFromTriplets(coeffv.begin(), coeffv.end());
	mat_L_.setFromTriplets(coeffL.begin(), coeffL.end());
	mat_J_.setFromTriplets(coeffJ.begin(), coeffJ.end());
	coeffv.clear();
	coeffL.clear();
	coeffJ.clear();

	// update B and A0 by initial V
	Sparse_Kronecker_Product(mat_B_);
	mat_B0_ = mat_B_;

    mat_ADMM_A_.resize(n_src_vertex_, 3);
    mat_ADMM_C_.resize(n_src_vertex_, 3);
    weight_d_.resize(n_src_vertex_);
    weight_s_.resize(src_mesh_->n_halfedges());
}
double ADMM::DoNonRigid()
{
    //MatrixXX prev_X, prev_grad;
    double new_err = 0.0;
    MatrixXX curV =  MatrixXX::Zero(n_src_vertex_, 3);
    MatrixXX prevV = MatrixXX::Zero(n_src_vertex_, 3);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    bool run_once = true;

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
    for(int out_iters = 0; out_iters < pars_.max_outer_iters; out_iters++)
    {
        //ADMM coeff
        mat_ADMM_Y1_.setZero(n_src_vertex_, 3);
        mat_ADMM_Y2_.setZero(src_mesh_->n_halfedges(), 3);
        mu1 = 0.5;
        mu2 = 0.5;
        raw1 = 2.0;
        raw2 = 2.0;

		mat_U0_.setZero();
        corres_pair_ids_.setZero();
        #ifdef USE_OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < correspondence_pairs_.size(); i++)
        {
            mat_U0_.col(correspondence_pairs_[i].first) = tar_points_.col(correspondence_pairs_[i].second);
            corres_pair_ids_[correspondence_pairs_[i].first] = 1;
        }

        //Establishing WS
        #ifdef USE_OPENMP
        #pragma omp parallel for
        #endif
        for (int e_it = 0; e_it < src_mesh_->n_halfedges(); e_it++)
        {
            Eigen::Vector4d vi_4D;
            int i = src_mesh_->from_vertex_handle(src_mesh_->halfedge_handle(e_it)).idx();
            int j = src_mesh_->to_vertex_handle(src_mesh_->halfedge_handle(e_it)).idx();
            vi_4D = mat_V0_.block(i, 4 * i, 1, 4).transpose();
            weight_s_[e_it] = 1.0 / (((vi_4D.transpose()*mat_X_.block(4 * i, 0, 4, 3) - vi_4D.transpose()*mat_X_.block(4 * j, 0, 4, 3))).lpNorm<1>() + ipsilonS);
        }

        // Establishing WD
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
        for (int v_idx = 0; v_idx < n_src_vertex_; v_idx++)
        {
            Eigen::Vector4d v_4D;
            v_4D = mat_V0_.block(v_idx, 4 * v_idx, 1, 4).transpose();
            if(abs(corres_pair_ids_[v_idx]) < 1e-3)
            {
                weight_d_[v_idx] = 0;
            }
            else
                weight_d_[v_idx] = 1.0 / ((v_4D.transpose()*(mat_X_.block(4 * v_idx, 0, 4, 3))
                                       - mat_U0_.transpose().block(v_idx, 0, 1, 3)).lpNorm<1>() + ipsilonD);
        }

        mat_V_ = weight_d_.asDiagonal() * mat_V0_;
        mat_B_ = weight_s_.asDiagonal() * mat_B0_;
        mat_U_ = mat_U0_ * weight_d_.asDiagonal();

        int iter;
        for (iter = 0; iter < pars_.max_inner_iters; iter++)
        {
            MatrixXX old_X = mat_X_;

            //ADMM_C
            mat_ADMM_C_ =  (mat_V_ * mat_X_ - mat_U_.transpose()) - (1 / mu1)*mat_ADMM_Y1_;
            #ifdef USE_OPENMP
            #pragma omp parallel for
            #endif
            for (int v_idx1 = 0; v_idx1 < n_src_vertex_; v_idx1++)
            {
                mat_ADMM_C_(v_idx1,0) = Sign(mat_ADMM_C_(v_idx1, 0))*std::max(fabs(mat_ADMM_C_(v_idx1, 0)) - 1 / mu1, 0.0);
                mat_ADMM_C_(v_idx1, 1) = Sign(mat_ADMM_C_(v_idx1, 1))*std::max(fabs(mat_ADMM_C_(v_idx1, 1)) - 1 / mu1, 0.0);
                mat_ADMM_C_(v_idx1, 2) = Sign(mat_ADMM_C_(v_idx1, 2))*std::max(fabs(mat_ADMM_C_(v_idx1, 2)) - 1 / mu1, 0.0);
            }

            // ADMM_A
            mat_ADMM_A_ = mat_B_ * mat_X_ - (1 / mu2)*mat_ADMM_Y2_;
            #ifdef USE_OPENMP
            #pragma omp parallel for
            #endif
            for (int v_idx1 = 0; v_idx1 < src_mesh_->n_halfedges(); v_idx1++)
            {
                for (int v_idx2 = 0; v_idx2 < 3; v_idx2++)
                {
                    mat_ADMM_A_(v_idx1, v_idx2) = Sign(mat_ADMM_A_(v_idx1, v_idx2))
                            *std::max(fabs(mat_ADMM_A_(v_idx1, v_idx2)) - pars_.alpha / mu2, 0.0);
                }
            }


            UpdateR();
            Eigen::SparseMatrix<double> LDL = mu1*mat_V_.transpose()*mat_V_ + mu2*mat_B_.transpose()*mat_B_ + 2 * pars_.beta *mat_L_;
            if(run_once)
            {
                solver.analyzePattern(LDL);
            }
            solver.factorize(LDL);
            mat_X_ = solver.solve(mat_V_.transpose()*(mat_ADMM_Y1_ + mu1*(mat_ADMM_C_ + mat_U_.transpose()))
                 + mat_B_.transpose()*(mat_ADMM_Y2_ + mu2*mat_ADMM_A_) + 2 * pars_.beta*mat_J_*mat_R_);

            mat_ADMM_Y1_ = mat_ADMM_Y1_ + mu1*(mat_ADMM_C_ - (mat_V_*mat_X_ - mat_U_.transpose()));
            mat_ADMM_Y2_ = mat_ADMM_Y2_ + mu2*(mat_ADMM_A_ - mat_B_*mat_X_);
            mu1 = raw1*mu1;
            mu2 = raw2*mu2;

            if ((old_X - mat_X_).norm() < pars_.stop)
            {
                break;
            }
        }

        Eigen::VectorXd gt_errs= Eigen::VectorXd::Zero(n_src_vertex_);
        // according new transformation X to update V and src_mesh_
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n_src_vertex_; i++)
        {
            MatrixXX tar_p = mat_V0_.block(i, 4 * i, 1, 4) * mat_X_.block(4 * i, 0, 4, 3);
            //std::cout << "mat V0.block = " << mat_V0_.block(i, 4 * i, 1, 4) << std::endl;
            OpenMesh::Vec3d p(tar_p(0, 0), tar_p(0, 1), tar_p(0, 2));
            src_mesh_->set_point(src_mesh_->vertex_handle(i), p);
            curV.row(i) = tar_p;
            if(pars_.calc_gt_err)
                gt_errs[i] = (tar_p.transpose() - tar_points_.col(i)).squaredNorm();
        }

        run_time = time.get_time();
        double gt_err = std::sqrt(gt_errs.sum()/n_src_vertex_);
        pars_.each_gt_mean_errs.push_back(gt_err);
        pars_.each_gt_max_errs.push_back(sqrt(gt_errs.maxCoeff()));
        double energy = 0.0;
        pars_.each_energys.push_back(energy);
        double eps_time = time.elapsed_time(begin_time, run_time);
        pars_.each_times.push_back(eps_time);
        pars_.each_iters.push_back(iter);

        if(pars_.print_each_step_info)
        {
            std::cout << "iter = " << out_iters << " time = " << eps_time << " inner iter = " << iter << " energy = " << energy << " gt_err = " << gt_err << std::endl;
            std::string out_obj = pars_.out_each_step_info + "/iter_obj_" + std::to_string(out_iters) + ".ply";
            OpenMesh::IO::write_mesh( *src_mesh_, out_obj.c_str());
        }

        // Find clost points
        FindClosestPoints(correspondence_pairs_);
        SimplePruning(correspondence_pairs_, pars_.use_distance_reject, pars_.use_normal_reject);

        // stop condition should be revised
        if((curV - prevV).rowwise().norm().maxCoeff() < pars_.stop)
        {
            break;
        }
        prevV = curV;
    }
	return new_err;
}


int ADMM::Sign(double n)
{
	if (n > 0)
		return 1;
	else if (n < 0)
		return -1;
	else
		return 0;
}

double ADMM::EvalError(double& data_err, double& smooth_err, double& orth_err)
{
	// data error
    data_err = (mat_V_ * mat_X_ - mat_U_.transpose()).lpNorm<1>();
    smooth_err = (mat_B_ * mat_X_).lpNorm<1>();
	Eigen::VectorXd orth_errs(n_src_vertex_);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < n_src_vertex_; i++)
	{
		orth_errs[i] = (mat_X_.block(4 * i, 0, 3, 3) - mat_R_.block(3 * i, 0, 3, 3)).squaredNorm();
	}
	orth_err = orth_errs.sum();
	return data_err + pars_.alpha * smooth_err + pars_.beta * orth_err;
}


void ADMM::EvalGradient()
{
	//Compute the gradient.
	grad_X_ = 2 * (mat_A0_ * mat_X_ - mat_VU_ - pars_.beta * mat_J_ * mat_R_);
}

void ADMM::UpdateR()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < n_src_vertex_; i++)
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat_X_.block(4 * i, 0, 3, 3), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d V1 = svd.matrixV(), U1 = svd.matrixU();
		mat_R_.block(3 * i, 0, 3, 3) = U1 * (V1.transpose());
        double deg_R = mat_R_.block(3*i, 0, 3,3).determinant();
        if(fabs(deg_R -1.0) > 1e-3)
        {
            Eigen::Vector3d s(1,1,-1);
            mat_R_.block(3 * i, 0, 3, 3) = U1 * s.asDiagonal() * (V1.transpose());
        }
        deg_R = mat_R_.block(3*i, 0, 3,3).determinant();
	}
}

void ADMM::Sparse_Kronecker_Product(Eigen::SparseMatrix<double> & B)
{
	B.resize(src_mesh_->n_halfedges(), 4 * n_src_vertex_);
    std::vector<Triplet> coefficients(src_mesh_->n_halfedges() * 8);
	double temp = 1.0;
	int k = 0;
	for (auto e_it = src_mesh_->halfedges_begin(); e_it != src_mesh_->halfedges_end(); ++e_it)
	{
		int i = src_mesh_->from_vertex_handle(*e_it).idx();
		int j = src_mesh_->to_vertex_handle(*e_it).idx();
		coefficients[8 * k + 0] = (Triplet((*e_it).idx(), i * 4 + 0, temp*(src_mesh_->point(src_mesh_->vertex_handle(i))[0])));
		coefficients[8 * k + 1] = (Triplet((*e_it).idx(), i * 4 + 1, temp*(src_mesh_->point(src_mesh_->vertex_handle(i))[1])));
		coefficients[8 * k + 2] = (Triplet((*e_it).idx(), i * 4 + 2, temp*(src_mesh_->point(src_mesh_->vertex_handle(i))[2])));
		coefficients[8 * k + 3] = (Triplet((*e_it).idx(), i * 4 + 3, temp));
		coefficients[8 * k + 4] = (Triplet((*e_it).idx(), j * 4 + 0, -temp * (src_mesh_->point(src_mesh_->vertex_handle(i))[0])));
		coefficients[8 * k + 5] = (Triplet((*e_it).idx(), j * 4 + 1, -temp * (src_mesh_->point(src_mesh_->vertex_handle(i))[1])));
		coefficients[8 * k + 6] = (Triplet((*e_it).idx(), j * 4 + 2, -temp * (src_mesh_->point(src_mesh_->vertex_handle(i))[2])));
		coefficients[8 * k + 7] = (Triplet((*e_it).idx(), j * 4 + 3, -temp));
		k++;
	}
	B.setFromTriplets(coefficients.begin(), coefficients.end());
}
