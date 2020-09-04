#ifndef QN_WELSCH_H_
#define QN_WELSCH_H_
#include "Registration.h"
#include "nodeSampler.h"

typedef Eigen::SparseMatrix<Scalar> SparseMatrix;

class QNwelsch : public Registration
{
public:
    QNwelsch();
    ~QNwelsch();
    virtual double DoNonRigid();
    virtual void Initialize();

private:
    double welsch_error(double nu1, double nu2);
    double welsch_energy(Eigen::VectorXd& r, double p);
    void welsch_weight(Eigen::VectorXd& r, double p);

    void LBFGS(int iter, MatrixXX & dir) const;
    int QNSolver(double& data_err, double& smooth_err, double& orth_err);

	double sample_energy(double& data_err, double& smooth_err, double& orth_err);
	void   sample_gradient();
    void   update_R();

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>* ldlt_;

private:
    // BFGS paras
    MatrixXX all_s_; // si = X_i+1 - X_i; all_s_ = 4n * 3lbfgs_m_
    MatrixXX all_t_; // ti = Grad(X_i+1) - Grad(X_i); all_t_ = 4n * 3lbfgs_m_
    int iter_;
    int col_idx_;

	// Sample paras
	int				num_sample_nodes; // (m) the number of sample nodes
	// simplily nodes storage structure
    svr::nodeSampler src_sample_nodes;
	// variable
	MatrixXX		Smat_X_;	// (4m * 3) transformations of sample nodes
	// data perm
    Eigen::SparseMatrix<double>	Weight_PV_; // (n * 4m) the weighting matrix of the sampling points affecting all vertices
	MatrixXX		Smat_P_;    // (n * 3) all sample nodes' coordinates
	// smooth perm
    Eigen::SparseMatrix<double>	Smat_B_;	// (|e| * 4m) the smooth between nodes, |e| is the edges' number of sample node graph;
	MatrixXX		Smat_D_;	// (|e| * 3) the different coordinate between xi and xj
	VectorX			Sweight_s_; // (|e|) the smooth weight;
	// orth perm
	MatrixXX		Smat_R_;	// (3m * 3)
    Eigen::SparseMatrix<double>	Smat_L_;	// (4m * 4m)
    Eigen::SparseMatrix<double>	Smat_J_;	// (4m * 3m)
    MatrixXX		Smat_UP_;   // aux matrix
    MatrixXX        mat_N0_;

    double          ori_alpha;
    double          ori_beta;
};
#endif
