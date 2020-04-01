#include "io_mesh.h"
#include "ADMM.h"
#include "QNwelsch.h"
#include "SVRL0.h"
#include "NonrigidICP.h"
#include "OmpHelper.h"


int main(int argc, char **argv)
{
    Mesh src_mesh;
    Mesh tar_mesh;
    std::string src_file;
    std::string tar_file;
    std::string out_file;
    std::string landmark_file;
    RegParas paras;

    paras.rigid_iters = 15;

    // our method paras
    paras.use_welsch_function = true;
    paras.corres_type = CLOSEST;
    paras.pruning_type = NONE;

    // sample
    paras.print_each_step_info = false;
    paras.use_distance_reject = false;
    paras.use_normal_reject = false;

    bool use_demean = true;

    enum METHOD {NICP, RPTS, SVR, QNWELSCH} method=QNWELSCH;
    std::string res_folder;

    if(argc==4)
    {
        src_file = argv[1];
        tar_file = argv[2];
        out_file = argv[3];
    }
    else if(argc==5)
    {
        src_file = argv[1];
        tar_file = argv[2];
        out_file = argv[3];
        method = METHOD(std::stoi(argv[4]));
    }
    else if(argc==6)
    {
        src_file = argv[1];
        tar_file = argv[2];
        out_file = argv[3];
        method = METHOD(std::stoi(argv[4]));
        landmark_file = argv[5];
        paras.corres_type = LANDMARK;
    }
    else
    {
        std::cout << "Usage: <srcFile> <tarFile> <outFile>\n    or <srcFile> <tarFile> <outFile> <regMethod> <landmarkFile>" << std::endl;
        exit(0);
    }

    // Setting paras
    paras.alpha = 1;
    paras.beta = 1000.0;

    paras.max_inner_iters = 20;
    paras.max_outer_iters = 100;


    read_data(src_file, src_mesh, use_demean);
    read_data(tar_file, tar_mesh, use_demean);
    if(src_mesh.n_vertices()==0 || tar_mesh.n_vertices()==0)
        exit(0);

    if(paras.corres_type==LANDMARK)
        read_landmark(landmark_file.c_str(), paras.landmark_src, paras.landmark_tar);
    double scale = mesh_scaling(src_mesh, tar_mesh);

    Registration* reg;
    switch (method)
    {
	case NICP:
	{
		reg = new NonrigidICP;

		break;
	}
    case RPTS:
    {
        reg = new ADMM;
        break;
    }
    case SVR:
    {
        reg = new SVRL0;
        break;
    }
	case QNWELSCH:
	{
		reg = new QNwelsch;
		break;
	}
    }

    Timer time;
    std::cout << "\nrigid registration to initial..." << std::endl;
    Timer::EventID time1 = time.get_time();
    reg->rigid_init(src_mesh, tar_mesh, paras);

    Timer::EventID time2 = time.get_time();
    std::cout << "rgid registration... " << std::endl;
    reg->DoRigid();
    // non-rigid initialize
    std::cout << "use method " << method << " to non-rigid registration..." << std::endl;
    Timer::EventID time3 = time.get_time();
    reg->Initialize();
    Timer::EventID time4 = time.get_time();
    reg->pars_.non_rigid_init_time = time.elapsed_time(time3, time4);
    std::cout << "non-rigid registration... " << std::endl;
    reg->DoNonRigid();
    Timer::EventID time5 = time.get_time();

    std::cout << "Registration done!\nrigid_init time : "
              << time.elapsed_time(time1, time2) << " s \trigid-reg run time = " << time.elapsed_time(time2, time3)
              << " s \nnon-rigid init time = "
              << time.elapsed_time(time3, time4) << " s \tnon-rigid run time = "
              << time.elapsed_time(time4, time5) << " s\n" << std::endl;

    if(use_demean)
        src_mesh.data(src_mesh.vertex_handle(0)).mesh_trans =  tar_mesh.data(tar_mesh.vertex_handle(0)).mesh_trans;

    write_data(out_file.c_str(), src_mesh, use_demean, scale);
    std::cout<< "write result to " << out_file << std::endl;
    delete reg;

    return 0;
}
