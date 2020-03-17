#include <Prune.h>
#include "tools.h"

Prune::Prune()
{
}

Prune::~Prune()
{
}



void Prune::Initialize(hgeodesic* geo, VPairs& corres, double c_0)
{
    delta = 0.05;
    c0 = c_0;

    // calculate D (the diameter of the mesh) using the furthest sampling scheme and dijkstra's algorithm
    D1 = calc_diameter(geo->data1, geo->n_src_vertices);
    D2 = calc_diameter(geo->data2, geo->n_tar_vertices);

    // calculate all K
    K_.resize(corres.size(), corres.size());
    std::vector<Eigen::Triplet<double>> coeff;
    k_ab.resize(corres.size(), corres.size());
    k_ab.setZero();
    for(int i = 0; i < corres.size(); i++)
    {
        double t1 = omp_get_wtime();

        Eigen::VectorXd dist1(geo->n_src_vertices);
        Eigen::VectorXd dist2(geo->n_tar_vertices);
        Eigen::MatrixXi gamma(1,1);
        gamma(0,0) = corres[i].first;
        igl::heat_geodesics_solve(geo->data1, gamma, dist1);
        gamma(0,0) = corres[i].second;
        igl::heat_geodesics_solve(geo->data2, gamma, dist2);

        for(int j = i+1; j < corres.size(); j++)
        {
            double d1 = dist1[corres[j].first] /2.0;
            double d2 = dist2[corres[j].second]/2.0;

            k_ab(i, j) = d1;
            k_ab(j, i) = d2;

            if(d1 < delta * D1 && d2 < delta * D2)
            {
                double k_ab = std::min(d2/d1, d1/d2);
                if(k_ab > c0)
                {
                    double t = (k_ab - c0) * (k_ab - c0)/( (1-c0) * (1-c0));
                    coeff.push_back(Eigen::Triplet<double>(i, j, t));
                    coeff.push_back(Eigen::Triplet<double>(j, i, t));
                }
            }
        }
    }
    K_.setFromTriplets(coeff.begin(), coeff.end());
    // calculate all scores pi
    double total_dist = K_.sum();
    confidence_scores.resize(K_.cols());
    for(int i = 0; i < K_.cols(); i++)
    {
        confidence_scores[i] = K_.col(i).sum()/total_dist;
    }
    //std::cout << "confidence score = " << confidence_scores.sum() << std::endl;
}

void Prune::RunPruning(std::vector<std::pair<int, int>>& corres)
{
    // calculate
    int idx;
    confidence_scores.maxCoeff(&idx);
    Eigen::VectorXd is_in_corres = Eigen::VectorXd::Ones(corres.size(),1);
    std::vector<int> new_corres_idx;
    new_corres.push_back(corres[idx]);
    new_corres_idx.push_back(idx);
    Newseeds.push_back(idx);
    is_in_corres[idx] = 0;


    // calc global geodesic distance from s_i and t_u
    Eigen::VectorXd temp = Eigen::VectorXd::Ones(corres.size());

    int iter = 0;
    while(temp.sum())
    {
        temp = confidence_scores.cwiseProduct(is_in_corres);
        temp.maxCoeff(&idx);
        is_in_corres[idx] = 0;

        // update B1, B2 and B
        Eigen::VectorXi B1(new_corres.size()), B2(new_corres.size()), B(new_corres.size());
        B1.setZero();
        B2.setZero();
        B.setZero();

        // save calculated geodesic in new_corres
        std::vector<double> geodist_in_new_corres(new_corres.size());
        for(int i = 0; i < new_corres.size(); i++)
        {
            double d1, d2;
            if(idx > new_corres_idx[i])
            {
                d1 = k_ab(idx, new_corres_idx[i]);
                d2= k_ab(new_corres_idx[i], idx);
            }
            else
            {
                d2 = k_ab(idx, new_corres_idx[i]);
                d1 = k_ab(new_corres_idx[i], idx);
            }


            geodist_in_new_corres[i] = d1/d2;

            if(d1 < delta * D1)
            {
                B1[i] = 1;
            }
            if(d2 < delta * D2)
            {
                B2[i] = 1;
            }
            if(d1 < delta * D1 && d2< delta * D2)
            {
                B[i] = 1;
            }
        }

        // judge C1 and C2 and C3
        // C1
        if(B.sum() > 0 && 1.0 * B.sum()/std::max(B1.sum(), B2.sum()) > 0.5)
        {
            int flag = 0;
            for(int i = 0; i < new_corres.size(); i++)
            {
                if(B[i])
                {
                    if(std::min(geodist_in_new_corres[i], 1.0/ geodist_in_new_corres[i])>c0)
                    {
                        flag++;
                    }
                }
            }
            if(flag == B.sum())
            {
                new_corres.push_back(corres[idx]);
                new_corres_idx.push_back(idx);

                if(Delay.size()>0)
                {
                    for(int j = 0; j < Delay.size(); j++)
                        is_in_corres[Delay[j]] = 1;
                    Delay.clear();
                }
            }
        }
        // C2
        else if(B1.sum() == 0 && B2.sum() == 0)
        {
            int flag = 0;
            for(int i = 0; i < Newseeds.size(); i++)
            {
                double d1, d2;
                if(idx > Newseeds[i])
                {
                    d1 = k_ab(idx, Newseeds[i]);
                    d2 = k_ab(Newseeds[i], idx);
                }
                else
                {
                    d2 = k_ab(idx, Newseeds[i]);
                    d1 = k_ab(Newseeds[i], idx);
                }

                if(std::min(d1/d2, d2/d1) >= 0.7)
                {
                    flag++;
                }

            }
            if(flag == Newseeds.size())
            {
                new_corres.push_back(corres[idx]);
                new_corres_idx.push_back(idx);

                Newseeds.push_back(idx);
                if(Delay.size()>0)
                {
                    for(int j = 0; j < Delay.size(); j++)
                        is_in_corres[Delay[j]] = 1;
                    Delay.clear();
                }
            }
        }
        // C3
        else
        {
            Delay.push_back(idx);
        }
    }
    corres.clear();
    corres = new_corres;
    return;
}

double Prune::calc_diameter(const igl::HeatGeodesicsData<double>& data, int num)
{
    int cur_id = int(1.0 * rand()/RAND_MAX * num);
    int v0 = cur_id, v1 = cur_id;
    double d0 = -1, d1 = -1, d2 = 0;
    int n = 0;
    while(d2 > d0 && d2 > d0)
    {
        Eigen::VectorXd all_dijkstra_dist(num);
        Eigen::MatrixXi gamma(1,1);
        gamma(0,0) = cur_id;
        igl::heat_geodesics_solve(data, gamma, all_dijkstra_dist);
        int max_idx;
        d0 = d1;
        d1 = d2;
        d2 = all_dijkstra_dist.maxCoeff(&max_idx);
//        std::cout << "n = " << n++ << " cur_idx = " << cur_id << " max_idx = " << max_idx
//                  <<" d0 = " << d1 <<  " max_dist = " << d2 << std::endl;
        cur_id = max_idx;
    }
    return d1/2.0;
}
