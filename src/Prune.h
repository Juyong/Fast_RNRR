#ifndef PRUNE_H_
#define PRUNE_H_
#include <MeshDefinition.h>
#include <vector>
#include "Types.h"
#include "geodesic.h"


class Prune
{
public:
    Prune();
    ~Prune();

    double delta;
    double c0;
    void RunPruning(VPairs& corres);
    void Initialize(hgeodesic* geo, VPairs& corres, double c_0);

private:
    Eigen::VectorXd  confidence_scores;
    Eigen::SparseMatrix<double>     K_;
    VPairs new_corres;
    Eigen::MatrixXd  k_ab; // upper save source_corres, lower save tar_corres_geodesic distance
    std::vector<int> Newseeds;
    std::vector<int> Delay;
    int nvertex1, nvertex2;
    double D1, D2;

    double  calc_diameter(const igl::HeatGeodesicsData<double>& data, int num);

};

#endif
