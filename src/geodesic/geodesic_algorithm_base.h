#ifndef GEODESIC_ALGORITHM_BASE
#define GEODESIC_ALGORITHM_BASE

#include "geodesic_algorithm_exact_elements.h"
#include "geodesic_constants_and_simple_functions.h"
#define hfid0 0
#define hfid1 1

namespace geodesic{
class GeodesicAlgorithmBase
{
public:
    vertex_pointer opposite_vertex(edge_pointer e, vertex_pointer v)
    {
        halfedge_handle hf0 = this->mesh()->halfedge_handle(e, hfid0);
        if(this->mesh()->from_vertex_handle(hf0).idx()==v.idx())
            return this->mesh()->to_vertex_handle(hf0);
        else
            return this->mesh()->from_vertex_handle(hf0);
    };

    vertex_pointer opposite_vertex(face_pointer f, edge_pointer e)
    {
        halfedge_handle hf = this->mesh()->halfedge_handle(e, hfid0);
        hf = this->mesh()->face_handle(hf)==f? hf : this->mesh()->opposite_halfedge_handle(hf);
        return this->mesh()->to_vertex_handle(this->mesh()->next_halfedge_handle(hf));
    };

    bool belongs_v(edge_pointer e, vertex_pointer v)
    {
        halfedge_handle hf = this->mesh()->halfedge_handle(e, hfid0);
        return this->mesh()->from_vertex_handle(hf) == v ||
                this->mesh()->to_vertex_handle(hf) == v;
    }

    edge_pointer next_edge(face_pointer f, edge_pointer e, vertex_pointer v)
    {

        halfedge_handle hf = this->mesh()->halfedge_handle(e, hfid0);
        hf = this->mesh()->face_handle(hf)==f? hf : this->mesh()->opposite_halfedge_handle(hf);
        halfedge_handle next_hf;
        if(this->mesh()->from_vertex_handle(hf)==v)
        {
            next_hf = this->mesh()->next_halfedge_handle(this->mesh()->next_halfedge_handle(hf));
        }
        else if(this->mesh()->to_vertex_handle(hf)==v)
            next_hf = this->mesh()->next_halfedge_handle(hf);
        else
        {
            std::cout << "next_edge e and v has no connection" << std::endl;
            exit(1);
        }
        return this->mesh()->edge_handle(next_hf);
    };

    Scalar vertex_angle(face_pointer f, vertex_pointer v)
    {
        halfedge_handle hf0 = this->mesh()->halfedge_handle(f);
        vertex_pointer v0 = this->mesh()->from_vertex_handle(hf0);
        vertex_pointer v1 = this->mesh()->to_vertex_handle(hf0);
        vertex_pointer v2 = this->mesh()->to_vertex_handle(this->mesh()->next_halfedge_handle(hf0));
        Scalar l1 = (this->mesh()->point(v0) - this->mesh()->point(v1)).norm();
        Scalar l2 = (this->mesh()->point(v1) - this->mesh()->point(v2)).norm();
        Scalar l3 = (this->mesh()->point(v2) - this->mesh()->point(v0)).norm();
        if(v0.idx()==v.idx())
        {
            return angle_from_edges(l2, l3, l1);
        }
        if(v1.idx()==v.idx())
        {
            return angle_from_edges(l3, l1, l2);
        }
        if(v2.idx()==v.idx())
        {
            return angle_from_edges(l1, l2, l3);
        }
    };

    void update_edgelen()
    {
        for(auto eit = this->mesh()->edges_begin(); eit != this->mesh()->edges_end(); eit++)
        {
            halfedge_handle hf = this->mesh()->halfedge_handle(*eit, hfid0);
            Vec3 v1 = this->mesh()->point(this->mesh()->from_vertex_handle(hf));
            Vec3 v2 = this->mesh()->point(this->mesh()->to_vertex_handle(hf));
            this->mesh()->data(*eit).length = (v1-v2).norm();
        }
    };

    void build_adjacencies();

    edge_pointer opposite_edge(face_pointer f, vertex_pointer v)
    {
        for(auto e_it = this->mesh()->fe_begin(f); e_it != this->mesh()->fe_end(f); ++e_it)
        {
            edge_pointer e = *e_it;
            if(!belongs_v(e, v))
            {
                return e;
            }
        }
    };

    face_pointer opposite_face(edge_pointer e, face_pointer f)
    {
        halfedge_handle hf = this->mesh()->halfedge_handle(e, hfid0);
        return this->mesh()->face_handle(hf) == f?
                    this->mesh()->face_handle(this->mesh()->opposite_halfedge_handle(hf))
                  :this->mesh()->face_handle(hf);
    };

    void initialize(Mesh* mesh);

    GeodesicAlgorithmBase(Mesh* mesh)
    {
       initialize(mesh);
    };

    GeodesicAlgorithmBase(){m_mesh = NULL;};

    virtual ~GeodesicAlgorithmBase(){};

    virtual void print_statistics()		//print info about timing and memory usage in the propagation step of the algorithm
    {
        std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;
    };

    Mesh* mesh(){return m_mesh;};

    // propagate a window
    bool compute_propagated_parameters(Scalar pseudo_x,
                                       Scalar pseudo_y,
                                       Scalar start,
                                       Scalar end,		//start/end of the interval
                                       Scalar alpha,	//corner angle
                                       Scalar L,		//length of the new edge
                                       interval_pointer candidates,
                                       Scalar d);		//if it is the last interval on the edge

    // intersection point on an edge
    Scalar compute_positive_intersection(Scalar start,
                                         Scalar pseudo_x,
                                         Scalar pseudo_y,
                                         Scalar sin_alpha,
                                         Scalar cos_alpha);

    inline bool calculate_triangle_parameters(list_pointer &list, Triangle &Tri); // calculate the parameters of the triangle to be propagated

    list_pointer interval_list_0(edge_pointer e)
    {
        return &m_edge_interval_lists_0[e.idx()];
    };

    list_pointer interval_list_1(edge_pointer e)
    {
        return &m_edge_interval_lists_1[e.idx()];
    };


protected:

    Mesh* m_mesh;

    Triangle Tri; // the triangle to be propagated

    IntervalList wl_left, wl_right;

    std::vector<IntervalList> m_edge_interval_lists_0;		// windows propagated from adjacent_face[0] of the edge
    std::vector<IntervalList> m_edge_interval_lists_1;		// windows propagated from adjacent_face[1] of the edge

};

inline void GeodesicAlgorithmBase::initialize(Mesh* mesh)
{
    m_mesh = mesh;
    m_edge_interval_lists_0.resize(mesh->n_edges());
    m_edge_interval_lists_1.resize(mesh->n_edges());

    // initialize statistics
    m_queue_max_size      = 0;
    m_windows_propagation = 0;
    m_windows_wavefront   = 0;
    m_windows_peak        = 0;

    // initialize window lists, similar to half-edge structure
    for (unsigned i = 0; i < m_edge_interval_lists_0.size(); ++i)
    {
        edge_pointer edge = mesh->edge_handle(i);
        m_edge_interval_lists_0[i].initialize(edge);
        m_edge_interval_lists_1[i].initialize(edge);
        interval_list_0(edge)->start_vertex() = mesh->from_vertex_handle(mesh->halfedge_handle(edge, hfid0));
        interval_list_1(edge)->start_vertex() = mesh->from_vertex_handle(mesh->halfedge_handle(edge, hfid1));
    }
    // verify list links
    for (unsigned i = 0; i < mesh->n_faces(); ++i)
    {
        face_pointer f = (mesh->face_handle(i));
        vertex_pointer v[3];
        size_t j = 0;
        for (auto e_it = mesh->fe_begin(f); e_it != mesh->fe_end(f); ++e_it)
        {
            edge_pointer e = *e_it;
            if (mesh->face_handle(mesh->halfedge_handle(e, hfid0)) == f)
                v[j] = interval_list_0(e)->start_vertex();
            else
                v[j] = interval_list_1(e)->start_vertex();

            if ((interval_list_0(e)->start_vertex().idx() < 0) || (interval_list_1(e)->start_vertex().idx() < 0))
            {
                std::cout << "list link error" << std::endl;
                exit(1);
            }

            if (interval_list_0(e)->start_vertex() == interval_list_1(e)->start_vertex())
            {
                std::cout << "list link error" << std::endl;
                exit(1);
            }

            if (!((belongs_v(e, interval_list_0(e)->start_vertex())) &&
                  (belongs_v(e, interval_list_1(e)->start_vertex()))))
            {
                std::cout << "list link error" << std::endl;
                exit(1);
            }
            j++;
        }
        if ((v[0].idx() >-1 && v[0] == v[1]) || (v[0].idx() >-1 && v[0] == v[2]) || (v[1].idx() >-1 && v[1] == v[2]))
        {
            std::cout << "list link error" << std::endl;
            exit(1);
        }
    }
    update_edgelen();
    build_adjacencies();
};

inline Scalar GeodesicAlgorithmBase::compute_positive_intersection(Scalar start,
                                                                   Scalar pseudo_x,
                                                                   Scalar pseudo_y,
                                                                   Scalar sin_alpha,
                                                                   Scalar cos_alpha)
{
    //assert(pseudo_y < 0);
    assert(pseudo_y <= 0);

    Scalar denominator = sin_alpha*(pseudo_x - start) - cos_alpha*pseudo_y;
    if (denominator < 0.0)
    {
        return -1.0;
    }

    Scalar numerator = -pseudo_y*start;

    if (numerator < 1e-30)
    {
        return 0.0;
    }

    if (denominator < 1e-30)
    {
        return -1.0;
    }

    return numerator / denominator;
}

inline bool GeodesicAlgorithmBase::compute_propagated_parameters(Scalar pseudo_x,
                                                                 Scalar pseudo_y,
                                                                 Scalar begin,
                                                                 Scalar end,		//start/end of the interval
                                                                 Scalar alpha,	//corner angle
                                                                 Scalar L,		//length of the new edge
                                                                 interval_pointer candidates,
                                                                 Scalar d)
{
    assert(pseudo_y <= 0.0);
    assert(begin <= end);
    assert(begin >= 0);

    ++m_windows_propagation; // Statistics

    interval_pointer p = candidates;

    Scalar sin_alpha = sin(alpha);
    Scalar cos_alpha = cos(alpha);

    //important: for the first_interval, this function returns zero only if the new edge is "visible" from the source
    //if the new edge can be covered only after turn_over, the value is negative (-1.0)
    Scalar L1 = compute_positive_intersection(begin,
                                              pseudo_x,
                                              pseudo_y,
                                              sin_alpha,
                                              cos_alpha);

    if (L1 < 0 || L1 >= L) // Does not produce a window on the edge
        return false;

    Scalar L2 = compute_positive_intersection(end,
                                              pseudo_x,
                                              pseudo_y,
                                              sin_alpha,
                                              cos_alpha);

    if (L2 < 0 || L2 >= L) // Covers vertex
    {
        p->start() = L1;
        p->stop() = L;
        p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
        p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;
        assert(p->pseudo_y() <= 0.0);

        return true;
    }
    else
    {
        // Does not cover vertex
        p->start() = L1;
        p->stop() = L2;
        p->pseudo_x() = cos_alpha*pseudo_x + sin_alpha*pseudo_y;
        p->pseudo_y() = -sin_alpha*pseudo_x + cos_alpha*pseudo_y;
        assert(p->pseudo_y() <= 0.0);

        return true;
    }
}

inline bool GeodesicAlgorithmBase::calculate_triangle_parameters(list_pointer &list, Triangle &Tri) // Calculate the parameters of the triangle to be propagated
{
    OpenMesh::HalfedgeHandle hf0 = this->mesh()->halfedge_handle(list->edge(), hfid0);
    OpenMesh::HalfedgeHandle hf1 = this->mesh()->halfedge_handle(list->edge(), hfid1);
    size_t adjface_size=0;
    if(this->mesh()->face_handle(hf0).idx()>-1)
        adjface_size++;
    if(this->mesh()->face_handle(hf1).idx()>-1)
        adjface_size++;

    if (adjface_size > 1)
    {
        Tri.bottom_edge = list->edge();

        if (list == interval_list_0(Tri.bottom_edge))
            Tri.face = this->mesh()->face_handle(hf1);
        else
            Tri.face = this->mesh()->face_handle(hf0);

        Tri.top_vertex = opposite_vertex(Tri.face, Tri.bottom_edge);
        Tri.left_vertex = list->start_vertex();
        Tri.right_vertex = opposite_vertex(Tri.bottom_edge, Tri.left_vertex);

        Tri.left_edge = next_edge(Tri.face, Tri.bottom_edge, Tri.left_vertex);
        Tri.right_edge = next_edge(Tri.face, Tri.bottom_edge, Tri.right_vertex);

        Tri.top_alpha = vertex_angle(Tri.face, Tri.top_vertex);
        Tri.left_alpha = vertex_angle(Tri.face, Tri.left_vertex);
        Tri.right_alpha = vertex_angle(Tri.face, Tri.right_vertex);

        if (this->mesh()->face_handle(this->mesh()->halfedge_handle(Tri.left_edge, hfid0)) == Tri.face)
            Tri.left_list = interval_list_0(Tri.left_edge);
        else
            Tri.left_list = interval_list_1(Tri.left_edge);

        if (this->mesh()->face_handle(this->mesh()->halfedge_handle(Tri.right_edge, hfid0)) == Tri.face)
            Tri.right_list = interval_list_0(Tri.right_edge);
        else
            Tri.right_list = interval_list_1(Tri.right_edge);

        return false;
    }
    else
    {
        return true;
    }
}


inline void GeodesicAlgorithmBase::build_adjacencies()
{
    // define m_turn_around_flag for vertices
    std::vector<Scalar> total_vertex_angle(this->mesh()->n_vertices(), 0);
    for(auto f_it = this->mesh()->faces_begin(); f_it != this->mesh()->faces_end(); ++f_it)
    {
        halfedge_handle hf0 = this->mesh()->halfedge_handle(*f_it);
        vertex_pointer v0 = this->mesh()->from_vertex_handle(hf0);
        vertex_pointer v1 = this->mesh()->to_vertex_handle(hf0);
        vertex_pointer v2 = this->mesh()->to_vertex_handle(this->mesh()->next_halfedge_handle(hf0));
        Scalar l1 = (this->mesh()->point(v0) - this->mesh()->point(v1)).norm();
        Scalar l2 = (this->mesh()->point(v1) - this->mesh()->point(v2)).norm();
        Scalar l3 = (this->mesh()->point(v2) - this->mesh()->point(v0)).norm();

        total_vertex_angle[v0.idx()] += angle_from_edges(l2, l3, l1);
        total_vertex_angle[v1.idx()] += angle_from_edges(l3, l1, l2);
        total_vertex_angle[v2.idx()] += angle_from_edges(l1, l2, l3);
    }

    for(auto v_it = this->mesh()->vertices_begin(); v_it != this->mesh()->vertices_end(); ++v_it)
    {
        vertex_pointer v = *v_it;
        this->mesh()->data(v).saddle_or_boundary = (total_vertex_angle[v.idx()] > 2.0*M_PI - 1e-5);
    }

    for(auto e_it = this->mesh()->edges_begin(); e_it != this->mesh()->edges_end(); ++e_it)
    {
        edge_pointer e = *e_it;
        if(this->mesh()->is_boundary(e))
        {
            halfedge_handle hf = this->mesh()->halfedge_handle(e, hfid0);
            this->mesh()->data(this->mesh()->from_vertex_handle(hf)).saddle_or_boundary = true;
            this->mesh()->data(this->mesh()->to_vertex_handle(hf)).saddle_or_boundary = true;
        }
    }
}


}//geodesic

#endif
