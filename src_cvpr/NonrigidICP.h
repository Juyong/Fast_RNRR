#ifndef NICP_H_
#define NICP_H_
#include "Registration.h"

//typedef std::pair<int,int> Edge;
//typedef boost::shared_ptr< std::vector<Edge> > Edges;
//typedef vtkSmartPointer<vtkPoints> Vertices;
//typedef boost::shared_ptr< std::vector<float> > Weights;

class NonrigidICP: public Registration
{
public:
    NonrigidICP();
    ~NonrigidICP();

    virtual double DoNonRigid();
    virtual void Initialize();

private:
    double gamma = 1.0; // smooth parameters
    double alpha = 10.0; // smooth coefficient
    std::vector< Eigen::Triplet<double> > alpha_M_G;
};

#endif
