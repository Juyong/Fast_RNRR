#ifndef REGISTRATION_H_
#define REGISTRATION_H_
#include <time.h>
#include "nanoflann.h"
#include "tools.h"
#include <opencv2/opencv.hpp>
#include "geodesic.h"
//#include <Eigen/PardisoSupport>
#include "Types.h"
#include "OmpHelper.h"
//#include "PardisoSolver.h"

//#define EIGEN_USE_MKL_ALL


class Registration
{
public:
    Registration();
    virtual ~Registration();

    Mesh* src_mesh_;
    Mesh* tar_mesh_;
    int n_src_vertex_;
    int n_tar_vertex_;
    int n_landmark_nodes_;

protected:
    // non-rigid Energy function paras
    VectorX weight_d_;	 // E_data weight n * n diag(w1, w2, ... wn)
    VectorX weight_s_;	 // E_smooth weight |E| * |E| diag;
    VectorX weight_4o_;	 // E_orth weight

    MatrixXX grad_X_;	 // gradient of E about X; 4n * 3
    Eigen::SparseMatrix<double> mat_A0_;	 // intial Hessian approximation 4n * 4n
    MatrixXX direction_;  // descent direction 4n * 3
    MatrixXX tar_points_; // target_mesh  3*m;
    MatrixXX mat_VU_;         // U'V  4n*3
    MatrixXX mat_U0_;      // initial mat_U_; 3 * n

    KDtree* target_tree; // correspondence paras

    // Rigid paras
    Eigen::Affine3d rigid_T_;  // rigid registration transform matrix

    // Check correspondence points
    Eigen::VectorXd corres_pair_ids_;
    VPairs correspondence_pairs_;
    Eigen::SparseMatrix<double> sub_V_;
    MatrixXX sub_U_;
    int current_n_;

    // dynamic welsch parasmeters
    bool init_nu;
    double end_nu;
    double nu;

    bool update_tarGeotree;


public:
    // adjusted paras
    bool use_cholesky_solver_;
    bool use_pardiso_;
    RegParas pars_;

public:
    void nonrigid_init();
    virtual double DoNonRigid() { return 0.0; }
    double DoRigid();
    void rigid_init(Mesh& src_mesh, Mesh& tar_mesh, RegParas& paras);
    virtual void Initialize(){}

    // just to output corres file
    void find_good_init();

private:
    cv::Mat tar_geodist_;
    Eigen::VectorXi  init_geo_pairs;


protected:
    hgeodesic*  hgeodesic_;

    //point to point rigid registration
    template <typename Derived1, typename Derived2, typename Derived3>
    Eigen::Affine3d point_to_point(Eigen::MatrixBase<Derived1>& X,
        Eigen::MatrixBase<Derived2>& Y, const Eigen::MatrixBase<Derived3>& w);

    // Find correct correspondences
    void InitCorrespondence(VPairs & corres);
    void FindClosestPoints(VPairs & corres);
    void FindSHOTCorrespondence(VPairs & correspondence_pairs);
    int  CalcSHOTfeature(Mesh* mesh, KDtree* kdtree, int nFeats, cv::Mat& Features);

    // Pruning method
    void SimplePruning(VPairs & corres, bool use_distance, bool use_normal) const;
    void DiffusionPruning(VPairs & src_corres);

    // Use landmark;
    void LandMarkCorres(VPairs & correspondence_pairs);

    // Aux_tool function
    double CalcEdgelength(Mesh* mesh, int type);

};
#endif
