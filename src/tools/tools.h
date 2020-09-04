#ifndef TOOL_H_
#define TOOL_H_
#include "Types.h"

enum CorresType {CLOSEST, LANDMARK};
enum PruningType {SIMPLE, NONE};

struct RegParas
{
    int		max_outer_iters;    // nonrigid max iters
    int		max_inner_iters;    // nonrigid max quasi-newton iters
    Scalar	alpha;              // smooth
    Scalar	beta;               // orth
    Scalar  gamma;              // landmarks
    bool	use_lbfgs;          // use lbfgs speed up or not
    int		lbfgs_m;            // lbfgs parameters
    bool	use_normal_reject;  // use normal reject or not
    bool	use_distance_reject;// use distance reject or not
    Scalar	normal_threshold;
    Scalar	distance_threshold;
    int     rigid_iters;         // rigid registration max iterations
    bool	use_landmark;
    bool    use_fixedvex;        // Set the point which haven't correspondences
    bool    calc_gt_err;         // calculate ground truth error (DEBUG)
    bool    data_use_welsch;  // use robust welsch function as energy function or just use L2-norm
    bool    smooth_use_welsch; // use robust welsch function as energy function or just use L2-norm

    std::vector<int> landmark_src;
    std::vector<int> landmark_tar;
    std::vector<int> fixed_vertices;

    // dynamic welsch method parameters
    bool    use_Dynamic_nu;
    Scalar  Data_nu;
    Scalar  Smooth_nu;
    Scalar  Data_initk;
    Scalar  Data_endk;
    Scalar  stop;

	// Sample para
    Scalar  uni_sample_radio;       // uniform sample radio
    bool    print_each_step_info;   // debug output : each step nodes, correspondences

    // output path
    std::string out_gt_file;
    std::string out_each_step_info;
    int         num_sample_nodes;

    std::vector<Scalar> each_times;
    std::vector<Scalar> each_gt_mean_errs;
    std::vector<Scalar> each_gt_max_errs;
    std::vector<Scalar> each_energys;
    std::vector<Scalar> each_iters;
    std::vector<Vector3> each_term_energy;
    Scalar  non_rigid_init_time;
    Scalar  init_gt_mean_errs;
    Scalar  init_gt_max_errs;


    RegParas() // default
    {
        max_outer_iters = 100;
        max_inner_iters = 20;
        alpha = 100.0;  // smooth
        beta = 100.0;   // orth
        gamma = 1e6;
        use_lbfgs = true;
        lbfgs_m = 5;
        use_normal_reject = false;
        use_distance_reject = false;
        normal_threshold = M_PI / 3;
        distance_threshold = 0.3;
        rigid_iters = 0;
        use_landmark = false;
        use_fixedvex = false;
        calc_gt_err = false;
        data_use_welsch = true;
        smooth_use_welsch = true;

        // dynamic welsch method parameters
        use_Dynamic_nu = true;
        Data_nu = 0.0;
        Smooth_nu = 40;
        Data_initk =10;
        Data_endk = 0.5;
        stop = 1e-3;

        // Sample para
        uni_sample_radio = 5;
        print_each_step_info = false;

        non_rigid_init_time = .0;
        init_gt_mean_errs = .0;
    }

};

// normalize mesh
Scalar mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh);
// Convert Mesh to libigl format to calculate geodesic distance
void Mesh2VF(Mesh & mesh, MatrixXX& V, Eigen::MatrixXi& F);
Vec3 Eigen2Vec(Vector3 s);
Vector3 Vec2Eigen(Vec3 s);
// read landmark points into landmark_src and landmark_tar if they exist
bool read_landmark(const char* filename, std::vector<int>& landmark_src, std::vector<int>& landmark_tar);
// read fixed points into vertices_list if they exist
bool read_fixedvex(const char* filename, std::vector<int>& vertices_list);

#ifdef __linux__
bool my_mkdir(std::string file_path);
#endif

#endif
