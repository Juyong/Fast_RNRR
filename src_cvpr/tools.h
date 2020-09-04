#ifndef TOOL_H_
#define TOOL_H_
#include <vector>
#include <string>
#include "MeshDefinition.h"
#include <Eigen/Dense>

enum CorresType {CLOSEST, SHOT, LANDMARK};
enum PruningType {SIMPLE, DIFFUSION, NONE};

struct RegParas
{
    int		max_outer_iters;    // nonrigid max iters
    int		max_inner_iters;    // nonrigid max quasi-newton iters
    double	alpha;              // smooth
    double	beta;               // orth
    bool	use_lbfgs;          // use lbfgs speed up or not
    int		lbfgs_m;            // lbfgs parameters
    bool	use_normal_reject;  // use normal reject or not
    bool	use_distance_reject;// use distance reject or not
    double	normal_threshold;
    double	distance_threshold;
    bool	use_welsch_function; // use robust welsch function as energy function or just use L2-norm
    int         rigid_iters;         // rigid registration max iterations
    bool	use_landmark;
    bool        calc_gt_err;         // calculate ground truth error (DEBUG)

    CorresType		 corres_type; // select initial correspondences type [SHOT, CLOSEST]
    PruningType		 pruning_type; // select initial reject incorrect correspondences type [DIFFUSION, SIMPLE, or not]
    std::vector<int> landmark_src;
    std::vector<int> landmark_tar;

    // dynamic welsch method parameters
    bool    use_Dynamic_nu;
    double  Data_nu;
    double  Smooth_nu;
    double  Data_initk;
    double  Data_endk;
    double  stop;


    // diffusion parameters
    double  diffusion_para;

	// Sample para
    double  uni_sample_radio;       // uniform sample radio
    double  n_uni_sample_radio;     // non-uniform sample radio
    int     n_uni_sample_neighbors; // non-uniform sample maximum neighbors
    bool    is_unifom_sample;       // use uniform or non-uniform sample
    bool    is_allvertices_smooth;  // smooth constraint on all vertices or just sample nodes
    bool    print_each_step_info;   // debug output : each step nodes, correspondences

    // output path
    std::string out_gt_file;
    std::string out_each_step_info;
    int     n_correct_shot_pair;
    int     n_correct_prune_pair;
    int     n_total_shot_pair;
    int     n_total_prune_pair;
    int     num_sample_nodes;

    std::vector<double> each_times;
    std::vector<double> each_gt_mean_errs;
    std::vector<double> each_gt_max_errs;
    std::vector<double> each_energys;
    std::vector<double> each_iters;
    double  non_rigid_init_time;
    double  init_gt_mean_errs;
    double  init_gt_max_errs;


    RegParas() // default
    {
        max_outer_iters = 100;
        max_inner_iters = 20;
        alpha = 1.0;  // smooth
        beta = 1.0;   // orth
        use_lbfgs = true;
        lbfgs_m = 5;
        use_normal_reject = false;
        use_distance_reject = false;
        normal_threshold = M_PI / 3;
        distance_threshold = 0.3;
        use_welsch_function = true;
        rigid_iters = 15;
        use_landmark = false;
        calc_gt_err = false;

        // dynamic welsch method parameters
        use_Dynamic_nu = true;
        Data_nu = 0.0;
        Smooth_nu = 40;
        Data_initk =10;
        Data_endk = 0.5;
        stop = 1e-3;

        // Sample para
        uni_sample_radio = 5;
        n_uni_sample_radio = 0.05;
        n_uni_sample_neighbors = 7;
        is_unifom_sample = true;
        is_allvertices_smooth = false;
        print_each_step_info = false;

        corres_type = CLOSEST;
        pruning_type = NONE;

        diffusion_para = 0.5;
        n_correct_shot_pair = 0;
        n_correct_prune_pair = 0;
        n_total_shot_pair = 0;
        n_total_prune_pair = 0;

        non_rigid_init_time = .0;
        init_gt_mean_errs = .0;
        init_gt_max_errs = .0;
    }

};

// SHOT feature parameters
struct CmdLineParams
{
	double matchTh;
    double radiusMR;
    double sigmaNoiseMR;
    int nFeat;
    int minNeighbors;
	double rotationAngle;
	double rotationAxis[3];

    int nThreads;

    bool describeColor;
    bool describeShape;

    std::string datapath;
    std::string outputFile;

    int shotShapeBins;
    int shotColorBins;

    CmdLineParams() //default
    {
        shotShapeBins = 10;
        shotColorBins = 30;
        matchTh = 0.9f;
        radiusMR = 5;
        sigmaNoiseMR = 0.3;
        nFeat = 3000;
        minNeighbors = 5;
        describeShape = true;
        describeColor = false;
        rotationAngle = 60.0f;
        rotationAxis[0] = 0.75f;
        rotationAxis[1] = 0.1f;
        rotationAxis[2] = 1-0.75f*0.75f-0.1f*0.1f;
    }
};


inline float EuclideanDistance2(float* point1, float* point2) { return (point1[0] - point2[0])*(point1[0] - point2[0]) + (point1[1] - point2[1])*(point1[1] - point2[1]) + (point1[2] - point2[2])*(point1[2] - point2[2]); };

inline double EuclideanDistance2(double* point1, double* point2) { return (point1[0] - point2[0])*(point1[0] - point2[0]) + (point1[1] - point2[1])*(point1[1] - point2[1]) + (point1[2] - point2[2])*(point1[2] - point2[2]); };

inline float EuclideanDistance(float* point1, float* point2) { return sqrt((point1[0] - point2[0])*(point1[0] - point2[0]) + (point1[1] - point2[1])*(point1[1] - point2[1]) + (point1[2] - point2[2])*(point1[2] - point2[2])); };

inline double EuclideanDistance(double* point1, double* point2) { return sqrt((point1[0] - point2[0])*(point1[0] - point2[0]) + (point1[1] - point2[1])*(point1[1] - point2[1]) + (point1[2] - point2[2])*(point1[2] - point2[2])); };

//extern void Mesh2VF(const Mesh & mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

int GetBoundaryPoints(Mesh *polydata, bool* &boundaryPointsIds);


double frand();


double box_muller(double m, double s);	/* normal random variate generator */

double mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh);

#ifdef __linux__
bool my_mkdir(std::string file_path);
#endif

#endif
