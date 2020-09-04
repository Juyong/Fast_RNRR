#pragma once
#include "Registration.h"
typedef Eigen::Triplet<Scalar> Triplet;
#define ipsilonD 1
#define ipsilonS 1

class ADMM : public Registration
{
public:
	ADMM();
    ~ADMM();

private:
    Eigen::SparseMatrix<double> mat_VX_; // the set of deformed point clouds n * 4n = diag(v1,v2,..,vn)
	int count = 0;
    Eigen::SparseMatrix<double> mat_V_; // the set of source point clouds n * 4n = diag(v1,v2,..,vn)
    Eigen::SparseMatrix<double> mat_V0_; // the set of initial source point clouds n * 4n
	MatrixXX mat_X_;	 // the set of transformation matrices 4n * 3 = [X1,X2,...,Xn]'
	MatrixXX mat_U_;	 // the set of correspondence point in target point clouds 3 * n = [u1, u2, ..., un]
    Eigen::SparseMatrix<double> mat_B_; // the set of the cross of coefficient and vertice B = {b_ij} = {k_ij * v_i'}  |E| * 4n;
    Eigen::SparseMatrix<double> mat_L_; // the set of selected matrix 4n * 4n = \sum_i Gi'Si'SiGi G = 4 * 4n S = 3 * 4
    Eigen::SparseMatrix<double> mat_J_; // the set of selected matrix 4n * 3n = \sum_i Gi'Si'Ti T = 3 * 3n
	MatrixXX mat_R_;	 // the set of rotation matrix 3n * 3 = [R1, R2,...,Rn]'
    Eigen::SparseMatrix<double> mat_B0_; // the initial mat_B_


	// ADMM paras
	MatrixXX mat_ADMM_A_; // Matrix A in ADMM
	MatrixXX mat_ADMM_C_; // Matrix C in ADMM
	MatrixXX mat_ADMM_Y1_; // Multiplier Y1 in ADMM
	MatrixXX mat_ADMM_Y2_; // Multiplier Y2 in ADMM
	double mu1;
	double mu2;
	double raw1;
	double raw2;

public:
    virtual void Initialize();
    virtual double DoNonRigid();

private:
	void EvalGradient();
	double EvalError(double& data_err, double& smooth_err, double& orth_err);
	void UpdateR();
    void Sparse_Kronecker_Product(Eigen::SparseMatrix<double> & B);

	//ADMM
	int Sign(double n);

};
