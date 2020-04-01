#pragma once
#include "QNwelsch.h"
#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
typedef Eigen::Triplet<Scalar> Triplet;



QNwelsch::QNwelsch() {
};

QNwelsch::~QNwelsch()
{
}

void QNwelsch::Initialize()
{
	ldlt_ = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>;
	nonrigid_init();
	// Initialize BFGS paras
	iter_ = 0;
	col_idx_ = 0;

	weight_d_.setOnes();
	weight_s_.setOnes();

	// our sample method
	// construct node graph
	double time1 = clock();
	if (pars_.is_unifom_sample)
	{
		src_sample_nodes.sample(*src_mesh_, hgeodesic_->data1, pars_.uni_sample_radio, svr::nodeSampler::X_AXIS);
		// src_sample_nodes.farthest_point_sample(*src_mesh_, hgeodesic_->data1, pars_.uni_sample_radio, 0, 1);
	}
	else
		src_sample_nodes.mysample(*src_mesh_, hgeodesic_->data1, n_src_vertex_ * pars_.n_uni_sample_radio, pars_.n_uni_sample_neighbors);

	src_sample_nodes.constructGraph(pars_.is_unifom_sample);
	if (pars_.print_each_step_info)
	{
		std::string out_node = pars_.out_each_step_info + "/init_";
		src_sample_nodes.print_nodes(*src_mesh_, out_node);//init sample nodes
	}

	num_sample_nodes = src_sample_nodes.nodeSize();
	pars_.num_sample_nodes = num_sample_nodes;
	int num_sample_edges;
	if (pars_.is_allvertices_smooth)
		num_sample_edges = src_mesh_->n_halfedges();
	else
		num_sample_edges = num_sample_nodes*(num_sample_nodes - 1);
	//std::cout << "num_sample_nodes = " << num_sample_nodes << std::endl;

	all_s_.resize(4 * num_sample_nodes * 3, pars_.lbfgs_m); all_s_.setZero();
	all_t_.resize(4 * num_sample_nodes * 3, pars_.lbfgs_m); all_t_.setZero();

	Smat_X_.resize(4 * num_sample_nodes, 3); Smat_X_.setZero();
	Weight_PV_.resize(n_src_vertex_, 4 * num_sample_nodes);
	Smat_P_.resize(n_src_vertex_, 3);

	Smat_B_.resize(num_sample_edges, 4 * num_sample_nodes);
	Smat_D_.resize(num_sample_edges, 3);
	Sweight_s_.resize(num_sample_edges);

	Smat_R_.resize(3 * num_sample_nodes, 3); Smat_R_.setZero();
	Smat_L_.resize(4 * num_sample_nodes, 4 * num_sample_nodes);
	Smat_J_.resize(4 * num_sample_nodes, 3 * num_sample_nodes);

	std::vector<Triplet> coeffv(4 * num_sample_nodes);
	std::vector<Triplet> coeffL(3 * num_sample_nodes);
	std::vector<Triplet> coeffJ(3 * num_sample_nodes);
	for (int i = 0; i < num_sample_nodes; i++)
	{
		size_t vIdx0 = src_sample_nodes.getNodeVertexIdx(i);
		VertexHandle vh0 = src_mesh_->vertex_handle(vIdx0);
		Point v0 = src_mesh_->point(vh0);

		// Smat_X_
		Smat_X_.block(4 * i, 0, 3, 3) = MatrixXX::Identity(3, 3);

		// Smat_R_
		Smat_R_.block(3 * i, 0, 3, 3) = MatrixXX::Identity(3, 3);

		// Smat_L_
		coeffL[3 * i] = Triplet(4 * i, 4 * i, 1.0);
		coeffL[3 * i + 1] = Triplet(4 * i + 1, 4 * i + 1, 1.0);
		coeffL[3 * i + 2] = Triplet(4 * i + 2, 4 * i + 2, 1.0);

		// Smat_J_
		coeffJ[3 * i] = Triplet(4 * i, 3 * i, 1.0);
		coeffJ[3 * i + 1] = Triplet(4 * i + 1, 3 * i + 1, 1.0);
		coeffJ[3 * i + 2] = Triplet(4 * i + 2, 3 * i + 2, 1.0);

	}
	//Smat_P_.setFromTriplets(coeffv.begin(), coeffv.end());
	Smat_L_.setFromTriplets(coeffL.begin(), coeffL.end());
	Smat_J_.setFromTriplets(coeffJ.begin(), coeffJ.end());
	direction_.resize(4 * num_sample_nodes, 3);

	// Weight_PV_
	src_sample_nodes.initWeight(Weight_PV_, Smat_P_, Smat_B_, Smat_D_, Sweight_s_, *src_mesh_, pars_.is_allvertices_smooth);
	//std::cout << "sweight_s_ = " << Sweight_s_.transpose() << std::endl;

	int real_edge = 0;
	for (int i = 0; i < Smat_B_.rows(); i++)
	{
		if (Smat_B_.row(i).norm()>0)
		{
			real_edge++;
		}
	}
	pars_.alpha = pars_.alpha * n_src_vertex_ / real_edge;
	pars_.beta = pars_.beta * n_src_vertex_ / num_sample_nodes;

}

double QNwelsch::DoNonRigid()
{
	double data_err, smooth_err, orth_err;
	MatrixXX curV = MatrixXX::Zero(n_src_vertex_, 3);
	MatrixXX prevV = MatrixXX::Zero(n_src_vertex_, 3);

	SparseMatrix Weight_PV0 = Weight_PV_;
	SparseMatrix Smat_B0 = Smat_B_;
	MatrixXX	 Smat_D0 = Smat_D_;
	// welsch_sweight
	Eigen::VectorXd welsch_weight_s = Eigen::VectorXd::Ones(Sweight_s_.rows(), 1);
	bool run_once = true;

	// Data perm parameters
	double nu1 = pars_.Data_initk * pars_.Data_nu;
	double average_len = std::min(CalcEdgelength(src_mesh_, 1), CalcEdgelength(tar_mesh_, 1));
	double end_nu1 = pars_.Data_endk * average_len;
	double nu2 = pars_.Smooth_nu * average_len;

	// Smooth perm parameters
	ori_alpha = pars_.alpha;
	ori_beta = pars_.beta;
	pars_.alpha = ori_alpha * nu1 * nu1 / (nu2 * nu2);
	pars_.beta = ori_beta * 2.0 * nu1 * nu1;

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

	bool dynamic_stop = false;
	begin_time = time.get_time();
	while (!dynamic_stop)
	{
		for (int out_iter = 0; out_iter < pars_.max_outer_iters; out_iter++)
		{
			// according correspondence_pairs to update mat_U0_;
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

			// int welsch_iter;
			int total_inner_iters = 0;
			MatrixXX old_X = Smat_X_;


			// update V,U and D
			weight_d_ = weight_d_.cwiseSqrt().cwiseProduct(corres_pair_ids_);
			Weight_PV_ = weight_d_.asDiagonal() * Weight_PV0;
			welsch_weight_s = welsch_weight_s.cwiseSqrt().cwiseProduct(Sweight_s_);
			Smat_B_ = welsch_weight_s.asDiagonal() * Smat_B0;

			Smat_D_ = welsch_weight_s.asDiagonal() * Smat_D0;
			Smat_UP_ = weight_d_.asDiagonal() * (mat_U0_.transpose() - Smat_P_);

			// construct matrix A0 and pre-decompose
			mat_A0_ = (Weight_PV_).transpose() * Weight_PV_
				+ pars_.alpha * (Smat_B_).transpose() * Smat_B_
				+ pars_.beta * Smat_L_;

			if (run_once)
			{
				ldlt_->analyzePattern(mat_A0_);
				run_once = false;
			}
			ldlt_->factorize(mat_A0_);

			// auxiliary variable 4m * 3
			mat_VU_ = (Weight_PV_).transpose() * Smat_UP_
				+ pars_.alpha * (Smat_B_).transpose() * Smat_D_;

			total_inner_iters += QNSolver(data_err, smooth_err, orth_err);

			Eigen::VectorXd gt_errs = Eigen::VectorXd::Zero(n_src_vertex_);
			MatrixXX target = Weight_PV0 * Smat_X_ + Smat_P_;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
			for (int i = 0; i < n_src_vertex_; i++)
			{
				MatrixXX tar_p = target.block(i, 0, 1, 3);
				OpenMesh::Vec3d p(tar_p(0, 0), tar_p(0, 1), tar_p(0, 2));
				src_mesh_->set_point(src_mesh_->vertex_handle(i), p);
				curV.row(i) = tar_p;
				if (pars_.calc_gt_err)
					gt_errs[i] = (tar_p.transpose() - tar_points_.col(i)).squaredNorm();
			}


			double energy = 0.0;
			// update weight
			if (pars_.use_welsch_function)
			{
				weight_d_ = (Weight_PV0 * Smat_X_ + Smat_P_ - mat_U0_.transpose()).rowwise().norm();
				welsch_weight(weight_d_, nu1);
				welsch_weight_s = ((Smat_B0 * Smat_X_ - Smat_D0)).rowwise().norm();
				welsch_weight(welsch_weight_s, nu2);
			}

			run_time = time.get_time();
			double gt_err = std::sqrt(gt_errs.sum() / n_src_vertex_);
			pars_.each_gt_mean_errs.push_back(gt_err);
			pars_.each_gt_max_errs.push_back(sqrt(gt_errs.maxCoeff()));

			pars_.each_energys.push_back(energy);
			double eps_time = time.elapsed_time(begin_time, run_time);
			pars_.each_times.push_back(eps_time);
			pars_.each_iters.push_back(total_inner_iters);

			if (pars_.print_each_step_info)
			{
				std::cout << "iter = " << out_iter << " time = " << eps_time << " energy = " << energy << " gt_err = " << gt_err
					<< "  inner iter = " << total_inner_iters << std::endl;
			}

			// Find clost points
			FindClosestPoints(correspondence_pairs_);
			//        SimplePruning(correspondence_pairs_, pars_.use_distance_reject, pars_.use_normal_reject);

			// stop condition should be revised
			if ((curV - prevV).rowwise().norm().maxCoeff() < pars_.stop)
			{
				break;
			}
			prevV = curV;
		}
		if (fabs(nu1 - end_nu1)<1e-8 || !pars_.use_Dynamic_nu || !pars_.use_welsch_function)
			dynamic_stop = true;
		nu1 = std::max(0.5*nu1, end_nu1);
		nu2 *= 0.5;
		pars_.alpha = ori_alpha * nu1 * nu1 / (nu2 * nu2);
		pars_.beta = ori_beta * 2 * nu1 * nu1;
	}
	return 0;
}

double QNwelsch::sample_energy(double& data_err, double& smooth_err, double& orth_err)
{
	data_err = (Weight_PV_ * Smat_X_ - Smat_UP_).squaredNorm();
	smooth_err = ((Smat_B_ * Smat_X_ - Smat_D_)).squaredNorm();
	Eigen::VectorXd orth_errs(num_sample_nodes);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_sample_nodes; i++)
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(Smat_X_.block(4 * i, 0, 3, 3), Eigen::ComputeThinU | Eigen::ComputeThinV);
		// Eigen::Matrix3d V1 = svd.matrixV(), U1 = svd.matrixU();

		if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
			Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
			// Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
			Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();
		}
		else {
			// Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixV()*svd.matrixU().transpose();
			Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixU()*svd.matrixV().transpose();
		}
		orth_errs[i] = (Smat_X_.block(4 * i, 0, 3, 3) - Smat_R_.block(3 * i, 0, 3, 3)).squaredNorm();
	}
	orth_err = orth_errs.sum();
	return data_err + pars_.alpha * smooth_err + pars_.beta * orth_err;
}

void QNwelsch::update_R()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_sample_nodes; i++)
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(Smat_X_.block(4 * i, 0, 3, 3), Eigen::ComputeThinU | Eigen::ComputeThinV);
		if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
			Eigen::Vector3d S = Eigen::Vector3d::Ones(); S(2) = -1.0;
			// Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
			Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();
		}
		else {
			// Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixV()*svd.matrixU().transpose();
			Smat_R_.block(3 * i, 0, 3, 3) = svd.matrixU()*svd.matrixV().transpose();
		}
	}
}

void QNwelsch::sample_gradient()
{
	grad_X_ = 2 * (mat_A0_ * Smat_X_ - mat_VU_ - pars_.beta * Smat_J_ * Smat_R_);
}

// col vector
void QNwelsch::LBFGS(int iter, MatrixXX & dir) const
{
	VectorX rho(pars_.lbfgs_m);
	VectorX kersi(pars_.lbfgs_m);
	MatrixXX q(4 * num_sample_nodes, 3);
	MatrixXX temp(4 * num_sample_nodes, 3);
	int k = iter;
	q.setZero();
	dir = q;
	q = -grad_X_;
	int m_k = std::min(k, pars_.lbfgs_m);
	for (int i = k - 1; i > k - m_k - 1; i--)
	{
		int col = (pars_.lbfgs_m + col_idx_ - (k - 1 - i)) % pars_.lbfgs_m;
		rho(k - 1 - i) = all_t_.col(col).transpose().dot(all_s_.col(col));
		Scalar lbfgs_err_scalar = Eigen::Map<VectorX>(q.data(), q.size()).dot(all_s_.col(col));
		kersi(k - 1 - i) = lbfgs_err_scalar / rho(k - 1 - i);
		Eigen::Map<VectorX>(q.data(), q.size()) -= kersi(k - 1 - i) * all_t_.col(col);
	}
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int cid = 0; cid < 3; cid++)
	{
		dir.col(cid) = ldlt_->solve(q.col(cid));
	}
	//    dir = ldlt_->solve(q);

	for (int i = k - m_k; i < k; i++)
	{
		int col = (pars_.lbfgs_m + col_idx_ - (k - 1 - i)) % pars_.lbfgs_m;
		Scalar lbfgs_err_scalar = all_t_.col(col).dot(Eigen::Map<VectorX>(dir.data(), dir.size()));
		Scalar eta = kersi(k - 1 - i) - lbfgs_err_scalar / rho(k - 1 - i);
		Eigen::Map<VectorX>(dir.data(), dir.size()) += all_s_.col(col) * eta;
	}

	rho.resize(0);
	kersi.resize(0);
	q.resize(0, 0);
	temp.resize(0, 0);
	return;
}

int QNwelsch::QNSolver(double& data_err, double& smooth_err, double& orth_err)
{
	MatrixXX prev_X;
	int count_linesearch = 0;
	//double new_err = EvalError(data_err, smooth_err, orth_err);
	double new_err = sample_energy(data_err, smooth_err, orth_err);
	int iter;
	for (iter = 0; iter <= pars_.max_inner_iters; iter++)
	{
		//EvalGradient();
		sample_gradient();
		// update decent direction
		if (iter == 0)
		{
			MatrixXX temp = -grad_X_;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
			for (int cid = 0; cid < 3; cid++)
			{
				direction_.col(cid) = ldlt_->solve(temp.col(cid));
			}

			col_idx_ = 0;
			all_s_.col(col_idx_) = -Eigen::Map<Eigen::VectorXd>(Smat_X_.data(), 4 * num_sample_nodes * 3);
			all_t_.col(col_idx_) = -Eigen::Map<Eigen::VectorXd>(grad_X_.data(), 4 * num_sample_nodes * 3);
		}
		else
		{
			all_s_.col(col_idx_) += Eigen::Map<Eigen::VectorXd>(Smat_X_.data(), 4 * num_sample_nodes * 3);
			all_t_.col(col_idx_) += Eigen::Map<Eigen::VectorXd>(grad_X_.data(), 4 * num_sample_nodes * 3);

			// Get descent direction
			LBFGS(iter, direction_);
			col_idx_ = (col_idx_ + 1) % pars_.lbfgs_m;
			all_s_.col(col_idx_) = -Eigen::Map<Eigen::VectorXd>(Smat_X_.data(), 4 * num_sample_nodes * 3);
			all_t_.col(col_idx_) = -Eigen::Map<Eigen::VectorXd>(grad_X_.data(), 4 * num_sample_nodes * 3);
		}
		double alpha = 2.0;
		prev_X = Smat_X_;
		double prev_err = new_err;
		double gamma = 0.3;
		double x = (grad_X_.transpose() * direction_).trace();
		do
		{
			alpha /= 2;
			Smat_X_ = prev_X + alpha * direction_;
			//new_err = EvalError(data_err, smooth_err, orth_err);
			new_err = sample_energy(data_err, smooth_err, orth_err);
			count_linesearch++;
		} while (new_err > prev_err + gamma * alpha * x);

		// std::cout << "iter = " << iter << " err = " << new_err << std::endl;
		if (fabs(new_err - prev_err)< pars_.stop)
		{
			break;
		}
		iter_++;
	}
	return iter;
}

double QNwelsch::welsch_error(double nu1, double nu2)
{
	Eigen::VectorXd w_data = (Weight_PV_ * Smat_X_ - Smat_UP_).rowwise().norm();
	double data_err = welsch_energy(w_data, nu1);
	Eigen::VectorXd s_data = ((Smat_B_ * Smat_X_ - Smat_D_)).rowwise().norm();
	double smooth_err = welsch_energy(s_data, nu2);
	Eigen::VectorXd orth_errs(n_src_vertex_);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_sample_nodes; i++)
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(Smat_X_.block(4 * i, 0, 3, 3), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d V1 = svd.matrixV(), U1 = svd.matrixU();
		double deg_R = Smat_R_.block(3 * i, 0, 3, 3).determinant();
		if (abs(deg_R + 1) < 1e-3)
		{
			Eigen::Vector3d s(1, 1, -1);
			Smat_R_.block(3 * i, 0, 3, 3) = U1 * s.asDiagonal() * (V1.transpose());
		}
		orth_errs[i] = (Smat_X_.block(4 * i, 0, 3, 3) - Smat_R_.block(3 * i, 0, 3, 3)).squaredNorm();
	}
	return data_err + pars_.alpha * smooth_err + pars_.beta * orth_errs.sum();
}

double QNwelsch::welsch_energy(Eigen::VectorXd& r, double p) {
	VectorX energy(r.rows());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i<r.rows(); ++i) {
		energy[i] = 1.0 - std::exp(-r(i)*r(i) / (2 * p*p));
	}
	return energy.sum();
}

void QNwelsch::welsch_weight(Eigen::VectorXd& r, double p) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i<r.rows(); ++i) {
		r(i) = std::exp(-r(i)*r(i) / (2 * p*p));
	}
}
