#ifndef GEODESIC_ALGORITHM_EXACT
#define GEODESIC_ALGORITHM_EXACT

#include "geodesic_memory.h"
#include "geodesic_algorithm_base.h"
#include "geodesic_algorithm_exact_elements.h"
#include <string.h>

namespace geodesic {

	class GeodesicAlgorithmExact : public GeodesicAlgorithmBase
	{
	public:

		// basic functions related to class
        GeodesicAlgorithmExact(Mesh* mesh) :
			GeodesicAlgorithmBase(mesh),
            m_memory_allocator(mesh->n_edges(), mesh->n_edges())
        {ori_mesh = mesh;};

        // construct neighboring mesh around src within m_radius, if no mesh is constructed,
        // increase m_radius with increase_radio
        GeodesicAlgorithmExact(Mesh* mesh, size_t src, Scalar m_radius)
        {
            ori_mesh = mesh;
            Mesh* sub_mesh = new Mesh;
            bool ok = construct_submesh(sub_mesh, src, m_radius);
            if(ok)
            {
                GeodesicAlgorithmBase::initialize(sub_mesh);
                m_memory_allocator.reset(sub_mesh->n_edges(), sub_mesh->n_edges());
            }
            else
            {
                std::cerr << "Error:Some points cannot be covered under the specified radius, please increase the radius" << std::endl;
                exit(1);
            }
        };
		~GeodesicAlgorithmExact() {};
        void clear() {
            m_memory_allocator.clear();
            for(auto v_it = this->mesh()->vertices_begin();v_it != this->mesh()->vertices_end(); ++v_it)
            {
                this->mesh()->data(*v_it).geodesic_distance = GEODESIC_INF;
            }
       };

        // main entry
        void propagate(unsigned source, std::vector<size_t>& idxs);

        // print the resulting statistics
        void print_statistics();
	private:

        // simple functions
        void initialize_propagation_data();
        void create_pseudo_source_windows(vertex_pointer &v, bool UpdateFIFOQueue);
        void erase_from_queue(vertex_pointer& v);

        // propagate a windows list (Rule 1)
        void find_separating_point(list_pointer &list); // find the separating point of the windows and the list
        void propagate_windows_to_two_edges(list_pointer &list); // propagates windows to two edges accross a triangle face

        // pairwise windows checking (Rule 2)
        void check_with_vertices(list_pointer &list);
        windows_state check_between_two_windows(interval_pointer &w1, interval_pointer &w2); // Check two neighbouring crossing windows on same edge
        void pairwise_windows_checking(list_pointer &list); // Check crossing windows on same edge

        // main operation
        void propagate_one_windows_list(list_pointer &list);

        // construct neighboring mesh
        bool construct_submesh(Mesh* sub_mesh, size_t source_idx, Scalar radius);

		// member variables
        std::set<vertex_pointer> m_vertex_queue;
		std::queue<list_pointer> m_list_queue;                // FIFO queue for lists
		MemoryAllocator<Interval> m_memory_allocator;		  // quickly allocate and deallocate intervals 
        Scalar neighbor_radius;

        Eigen::VectorXi SubVidxfromMesh;
        std::vector<int> MeshVidxfromSub;

        Mesh* ori_mesh;
		unsigned m_source;
	};


	//----------------- simple functions ---------------------
	inline void GeodesicAlgorithmExact::initialize_propagation_data()
	{
		clear();

		// initialize source's parameters
        vertex_pointer source = (this->mesh()->vertex_handle(m_source));
        this->mesh()->data(source).geodesic_distance = 0;
        this->mesh()->data(source).state = VertexState::INSIDE;

		// initialize windows around source
		create_pseudo_source_windows(source, false);
	}

    inline void GeodesicAlgorithmExact::erase_from_queue(vertex_pointer& v)
	{
        assert(m_vertex_queue.count(v) <= 1);

        std::multiset<vertex_pointer>::iterator it = m_vertex_queue.find(v);
        if (it != m_vertex_queue.end())
            m_vertex_queue.erase(it);
	}

	inline void GeodesicAlgorithmExact::create_pseudo_source_windows(vertex_pointer &pseudo_source, bool inside_traversed_area)
	{
		// update vertices around pseudo_source
        for (auto e_it = this->mesh()->ve_begin(pseudo_source); e_it != this->mesh()->ve_end(pseudo_source); ++e_it)
		{
            edge_pointer   edge_it = *e_it;
            vertex_pointer vert_it = opposite_vertex(edge_it, pseudo_source);

            Scalar distance = this->mesh()->data(pseudo_source).geodesic_distance
                    + this->mesh()->data(edge_it).length;

            if (distance < this->mesh()->data(vert_it).geodesic_distance)
			{
				m_vertex_queue.erase(vert_it);

                this->mesh()->data(vert_it).geodesic_distance = distance;
                if (this->mesh()->data(vert_it).state == VertexState::OUTSIDE)
                    this->mesh()->data(vert_it).state = VertexState::FRONT;

                this->mesh()->data(vert_it).incident_face = this->mesh()->face_handle(this->mesh()->halfedge_handle(edge_it, hfid0));
                edge_pointer next_edge = geodesic::GeodesicAlgorithmBase::next_edge(
                            this->mesh()->data(vert_it).incident_face,edge_it, pseudo_source);
                this->mesh()->data(vert_it).incident_point =
                        (this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(next_edge, hfid0)) == pseudo_source) ?
                            0 : this->mesh()->data(next_edge).length;

                m_vertex_queue.insert(vert_it);
			}
		}

        // update pseudo_source windows around pseudo_source
        for(auto f_it = this->mesh()->vf_begin(pseudo_source); f_it != this->mesh()->vf_end(pseudo_source); ++f_it)
		{
            face_pointer face_it = *f_it;
            edge_pointer edge_it = geodesic::GeodesicAlgorithmBase::opposite_edge(face_it, pseudo_source);
            list_pointer list = (this->mesh()->face_handle(this->mesh()->halfedge_handle(edge_it, hfid0))==face_it)?
                    interval_list_0(edge_it) : interval_list_1(edge_it);

			// create a window
			interval_pointer candidate = new Interval;

			candidate->start() = 0;
            candidate->stop() = this->mesh()->data(edge_it).length;
            candidate->d() = this->mesh()->data(pseudo_source).geodesic_distance;
            Scalar angle = geodesic::GeodesicAlgorithmBase::vertex_angle(face_it, list->start_vertex());
            Scalar length = this->mesh()->data(geodesic::GeodesicAlgorithmBase::next_edge
                                               (face_it, edge_it,list->start_vertex())).length;
			candidate->pseudo_x() = cos(angle) * length;
			candidate->pseudo_y() = -sin(angle) * length;

			// insert into list
			list->push_back(candidate);

			// push into M_LIST_QUEUE if inside traversed area
            vertex_pointer v0 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(edge_it, hfid0));
            vertex_pointer v1 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(edge_it, hfid1));
            if ((inside_traversed_area) &&
                    ((this->mesh()->data(v0).state != VertexState::FRONT)
                     || (this->mesh()->data(v1).state != VertexState::FRONT)))
				m_list_queue.push(list);

			// Statistics
			++m_windows_wavefront;
			if (m_windows_peak < m_windows_wavefront)
				m_windows_peak = m_windows_wavefront;

		}
    }

	//----------------- propagate a windows list (Rule 1) ---------------------
	inline void GeodesicAlgorithmExact::find_separating_point(list_pointer &list)
    {
        const Scalar LOCAL_EPSILON = 1e-20 * this->mesh()->data(list->edge()).length; // numerical issue

        Scalar L = this->mesh()->data(Tri.left_edge).length;
        Scalar top_x = L * cos(Tri.left_alpha);
        Scalar top_y = L * sin(Tri.left_alpha);

        Scalar temp_geodesic = GEODESIC_INF;
        face_pointer temp_face_handle = this->mesh()->data(Tri.top_vertex).incident_face;
        Scalar temp_incident_point = this->mesh()->data(Tri.top_vertex).incident_point;

		interval_pointer iter = list->begin();

        Scalar wlist_sp = 0;
        Scalar wlist_pseudo_x = 0;
        Scalar wlist_pseudo_y = 0;

		while (iter != NULL)
		{
            interval_pointer &w = iter;

            Scalar w_sp = w->pseudo_x() - w->pseudo_y() * ((top_x - w->pseudo_x()) / (top_y - w->pseudo_y()));
            Scalar distance = GEODESIC_INF;

			// shortest path from the window
			if ((w_sp - w->start() > LOCAL_EPSILON) && (w_sp - w->stop() < -LOCAL_EPSILON))
			{
				distance = w->d() + sqrt((top_x - w->pseudo_x()) * (top_x - w->pseudo_x()) + (top_y - w->pseudo_y()) * (top_y - w->pseudo_y()));
				w->shortest_distance() = distance;
			}
			else if (w_sp - w->start() <= LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->start()) * (top_x - w->start()) + top_y * top_y) + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->start();
			}
			else if (w_sp - w->stop() >= -LOCAL_EPSILON)
			{
				distance = w->d() + sqrt((top_x - w->stop()) * (top_x - w->stop()) + top_y * top_y) + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
				w->shortest_distance() = distance;
				w_sp = w->stop();
			}

			// update information at top_t
            if (distance < temp_geodesic)
			{
                temp_geodesic = distance;
                temp_face_handle = Tri.face;
                vertex_pointer v0 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(list->edge(), hfid0));
                temp_incident_point = (list->start_vertex() == v0) ?
                            w_sp : this->mesh()->data(list->edge()).length - w_sp;
				wlist_sp = w_sp;
				wlist_pseudo_x = w->pseudo_x();
				wlist_pseudo_y = w->pseudo_y();
			}
			w->sp() = w_sp;

			iter = iter->next();
		}

        // update top_vertex and M_VERTEX_QUEUE
        if (temp_geodesic < this->mesh()->data(Tri.top_vertex).geodesic_distance)
		{
            if (this->mesh()->data(Tri.top_vertex).state == VertexState::FRONT) erase_from_queue(Tri.top_vertex);
            this->mesh()->data(Tri.top_vertex).geodesic_distance = temp_geodesic;
            this->mesh()->data(Tri.top_vertex).incident_face = temp_face_handle;
            this->mesh()->data(Tri.top_vertex).incident_point = temp_incident_point;
            if (this->mesh()->data(Tri.top_vertex).state == VertexState::FRONT)
            {
                m_vertex_queue.insert(Tri.top_vertex);
            }

            if ((this->mesh()->data(Tri.top_vertex).state == VertexState::INSIDE)
                    && (this->mesh()->data(Tri.top_vertex).saddle_or_boundary))
                create_pseudo_source_windows(Tri.top_vertex, true); // handle saddle vertex
		}

		list->sp() = wlist_sp;
		list->pseudo_x() = wlist_pseudo_x;
		list->pseudo_y() = wlist_pseudo_y;
	}

	inline void GeodesicAlgorithmExact::propagate_windows_to_two_edges(list_pointer &list)
	{
        const Scalar LOCAL_EPSILON = 1e-8 * this->mesh()->data(list->edge()).length; // numerical issue

		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		enum PropagationDirection
		{
			LEFT,
			RIGHT,
			BOTH
		};

		PropagationDirection direction;

		while (!list->empty() && (iter != NULL))
		{
			interval_pointer &w = iter;
			assert(w->start() <= w->stop());

			if (w->sp() < list->sp() - LOCAL_EPSILON)
			{
				// only propagate to left edge
                Scalar Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->stop(), 0, Intersect_X, Intersect_Y);
				if ((w->stop() < list->sp()) || ((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::LEFT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}				
			}
			else if (w->sp() > list->sp() + LOCAL_EPSILON)
			{
				// only propagate to right edge
                Scalar Intersect_X, Intersect_Y;

				// judge the positions of the two windows
				CalculateIntersectionPoint(list->pseudo_x(), list->pseudo_y(), list->sp(), 0, w->pseudo_x(), w->pseudo_y(), w->start(), 0, Intersect_X, Intersect_Y);
				if ((w->start() > list->sp())||((Intersect_Y <= 0) && (Intersect_Y >= list->pseudo_y()) && (Intersect_Y >= w->pseudo_y())))
				{
					direction = PropagationDirection::RIGHT;
				}
				else
				{
					direction = PropagationDirection::BOTH;
				}	
			}
			else
			{
				// propagate to both edges
				direction = PropagationDirection::BOTH;
			}

			bool ValidPropagation;
			interval_pointer right_w;

			switch (direction) {
			case PropagationDirection::LEFT:
				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
					Tri.left_alpha,
                    this->mesh()->data(Tri.left_edge).length,
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case PropagationDirection::RIGHT:
                ValidPropagation = compute_propagated_parameters(this->mesh()->data(Tri.bottom_edge).length - w->pseudo_x(),
					w->pseudo_y(),
                    this->mesh()->data(Tri.bottom_edge).length - w->stop(),
                    this->mesh()->data(Tri.bottom_edge).length - w->start(),
					Tri.right_alpha,
                    this->mesh()->data(Tri.right_edge).length,
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
                    Scalar length = this->mesh()->data(Tri.right_edge).length; // invert window
                    Scalar start = length - w->stop();
					w->stop() = length - w->start();
					w->start() = start;
					w->pseudo_x() = length - w->pseudo_x();

					list->erase(w);
					wl_right.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;
				break;

			case PropagationDirection:: BOTH:
				right_w = new Interval;
				memcpy(right_w, w, sizeof(Interval));

				ValidPropagation = compute_propagated_parameters(w->pseudo_x(),
					w->pseudo_y(),
					w->start(),
					w->stop(),
                    geodesic::GeodesicAlgorithmBase::vertex_angle(Tri.face, Tri.left_vertex),
                    this->mesh()->data(Tri.left_edge).length,
					w,
					w->d());

				iter_t = iter->next();
				if (ValidPropagation)
				{
					list->erase(w);
					wl_left.push_back(w);
				}
				else
				{
					list->erase(w);
					delete w;
					--m_windows_wavefront;
				}
				iter = iter_t;

                ValidPropagation = compute_propagated_parameters(this->mesh()->data(Tri.bottom_edge).length - right_w->pseudo_x(),
					right_w->pseudo_y(),
                    this->mesh()->data(Tri.bottom_edge).length - right_w->stop(),
                    this->mesh()->data(Tri.bottom_edge).length - right_w->start(),
                    geodesic::GeodesicAlgorithmBase::vertex_angle(Tri.face, Tri.right_vertex),
                    this->mesh()->data(Tri.right_edge).length,
					right_w,
					right_w->d());

				if (ValidPropagation)
				{
					// invert window
                    Scalar length = this->mesh()->data(Tri.right_edge).length;
                    Scalar start = length - right_w->stop();
					right_w->stop() = length - right_w->start();
					right_w->start() = start;
					right_w->pseudo_x() = length - right_w->pseudo_x();

					wl_right.push_back(right_w);

					++m_windows_wavefront;
					if (m_windows_peak < m_windows_wavefront)
						m_windows_peak = m_windows_wavefront;
				}
				else
				{
					delete right_w;
				}
				break;

			default:
				break;
			}
		}
	}

	//----------------- pairwise windows checking (Rule 2) ----------------------
	inline void GeodesicAlgorithmExact::check_with_vertices(list_pointer &list)
    {
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer iter_t;

		while ((!list->empty()) && (iter != NULL))
		{
			interval_pointer &w = iter;
			bool w_survive = true;

            edge_pointer   e = list->edge();
			vertex_pointer v1 = list->start_vertex();
            vertex_pointer v2 = opposite_vertex(e, v1);
            Scalar d1 = GEODESIC_INF;

			d1 = w->d() + sqrt((w->stop() - w->pseudo_x()) * (w->stop() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
            if (this->mesh()->data(v1).geodesic_distance + w->stop() < d1)
				w_survive = false;

			d1 = w->d() + sqrt((w->start() - w->pseudo_x()) * (w->start() - w->pseudo_x()) + w->pseudo_y() * w->pseudo_y());
            if (this->mesh()->data(v2).geodesic_distance + this->mesh()->data(e).length - w->start() < d1)
				w_survive = false;


			iter_t = iter;
			iter = iter->next();

			if (!w_survive)
			{
				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
			}
		}
	}

	inline windows_state GeodesicAlgorithmExact::check_between_two_windows(interval_pointer &w1, interval_pointer &w2)
	{
        Scalar NUMERCIAL_EPSILON = 1 - 1e-12;
		// we implement the discussed 6 cases as follows for simplicity

		if ((w1->start() >= w2->start()) && (w1->start() <= w2->stop())) // w1->start
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((w1->start() - w1->pseudo_x()) * (w1->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->start() - w2->pseudo_x()) * (w1->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->start() = w1->start();
			}
		}

		if ((w1->stop() >= w2->start()) && (w1->stop() <= w2->stop())) // w1->stop
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w1->stop(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((w1->stop() - w1->pseudo_x()) * (w1->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w1->stop() - w2->pseudo_x()) * (w1->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
				if (d1 < d2 * NUMERCIAL_EPSILON)
					w2->stop() = w1->stop();
			}
		}

		if ((w2->start() >= w1->start()) && (w2->start() <= w1->stop())) // w2->start
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((w2->start() - w1->pseudo_x()) * (w2->start() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->start() - w2->pseudo_x()) * (w2->start() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->start() = w2->start();
			}
		}

		if ((w2->stop() >= w1->start()) && (w2->stop() <= w1->stop())) // w2->stop
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w2->stop(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, Intersect_X, Intersect_Y);

			if ((Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((w2->stop() - w1->pseudo_x()) * (w2->stop() - w1->pseudo_x()) + (w1->pseudo_y()) * (w1->pseudo_y()));
				d2 = w2->d() + sqrt((w2->stop() - w2->pseudo_x()) * (w2->stop() - w2->pseudo_x()) + (w2->pseudo_y()) * (w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					w1->stop() = w2->stop();
			}
		}

		if (w1->start() >= w2->stop())
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w1->pseudo_x(), w1->pseudo_y(), w1->start(), 0, w2->pseudo_x(), w2->pseudo_y(), w2->stop(), 0, Intersect_X, Intersect_Y);

            face_pointer f = opposite_face(Tri.bottom_edge, Tri.face);
            edge_pointer e = next_edge(f, Tri.bottom_edge, Tri.left_vertex);
            Scalar angle = vertex_angle(f, Tri.left_vertex);
            Scalar Cx = this->mesh()->data(e).length * cos(angle);
            Scalar Cy = this->mesh()->data(e).length * -sin(angle);

            if ((PointInTriangle(Intersect_X, Intersect_Y, this->mesh()->data(Tri.bottom_edge).length, Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 * NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 * NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		if (w2->start() >= w1->stop())
		{
            Scalar Intersect_X, Intersect_Y;

			// judge the previous order of the two windows
			CalculateIntersectionPoint(w2->pseudo_x(), w2->pseudo_y(), w2->start(), 0, w1->pseudo_x(), w1->pseudo_y(), w1->stop(), 0, Intersect_X, Intersect_Y);

            face_pointer f = opposite_face(Tri.bottom_edge, Tri.face);
            edge_pointer e = next_edge(f, Tri.bottom_edge, Tri.left_vertex);
            Scalar angle = vertex_angle(f, Tri.left_vertex);
            Scalar Cx = this->mesh()->data(e).length * cos(angle);
            Scalar Cy = this->mesh()->data(e).length * -sin(angle);

            if ((PointInTriangle(Intersect_X, Intersect_Y, this->mesh()->data(Tri.bottom_edge).length, Cx, Cy))
				&& (Intersect_Y <= 0) && (Intersect_Y >= w1->pseudo_y()) && (Intersect_Y >= w2->pseudo_y()))
			{
                Scalar d1, d2;
				d1 = w1->d() + sqrt((Intersect_X - w1->pseudo_x()) * (Intersect_X - w1->pseudo_x()) + (Intersect_Y - w1->pseudo_y()) * (Intersect_Y - w1->pseudo_y()));
				d2 = w2->d() + sqrt((Intersect_X - w2->pseudo_x()) * (Intersect_X - w2->pseudo_x()) + (Intersect_Y - w2->pseudo_y()) * (Intersect_Y - w2->pseudo_y()));

				if (d1 < d2 - NUMERCIAL_EPSILON)
					return w2_invalid;
				if (d2 < d1 - NUMERCIAL_EPSILON)
					return w1_invalid;
			}
		}

		return both_valid;
	}

	inline void GeodesicAlgorithmExact::pairwise_windows_checking(list_pointer &list)
	{
		if (list->empty()) return;

		interval_pointer iter = list->begin();
		interval_pointer next, iter_t;

		next = iter->next();

		// traverse successive pairs of windows
		while ((!list->empty()) && (next != NULL))
		{
			windows_state ws = check_between_two_windows(iter, next);

			switch (ws)
			{
			case geodesic::w1_invalid:
				iter_t = iter;
				if (iter == list->begin())
				{
					iter = iter->next();
				}
				else
				{
					iter = iter->previous();
				}

				list->erase(iter_t);
				delete iter_t;
				--m_windows_wavefront;
				break;

			case geodesic::w2_invalid:
				list->erase(next);
				delete next;
				--m_windows_wavefront;
				break;

			case geodesic::both_valid:
				iter = iter->next();
				break;

			default:
				break;
			}

			next = iter->next();
		}
	}

	//------------------------- main operation ----------------------------
	inline void GeodesicAlgorithmExact::propagate_one_windows_list(list_pointer &list)
	{
		if (list->empty()) return;
        OpenMesh::HalfedgeHandle hf0 = this->mesh()->halfedge_handle(list->edge(), hfid0);
        OpenMesh::HalfedgeHandle hf1 = this->mesh()->halfedge_handle(list->edge(), hfid1);
        if (this->mesh()->face_handle(hf0).idx()>-1 && this->mesh()->face_handle(hf1).idx()>-1)
        {
			// Rule 2: pairwise windows checking
            check_with_vertices(list);
			pairwise_windows_checking(list);

            // Rule 1: "One Angle Two Sides"
            find_separating_point(list);
            propagate_windows_to_two_edges(list);
		}
	}

	//-------------------------- main entry --------------------------
    inline void GeodesicAlgorithmExact::propagate(unsigned source, std::vector<size_t>& idxs)
	{
		// initialization
        m_source = SubVidxfromMesh[source];
		initialize_propagation_data();
		while (!m_vertex_queue.empty())
		{
			// (1) pop a vertex from M_VERTEX_QUEUE
            vertex_pointer vert = *m_vertex_queue.begin();
			m_vertex_queue.erase(m_vertex_queue.begin());

			// (2) update wavefront
            this->mesh()->data(vert).state = VertexState::INSIDE;
            for(auto e_it = this->mesh()->ve_begin(vert); e_it != this->mesh()->ve_end(vert); ++e_it)
			{
                vertex_pointer vert_it = opposite_vertex(*e_it, vert);
                if (this->mesh()->data(vert_it).state == VertexState::OUTSIDE)
                    this->mesh()->data(vert_it).state = VertexState::FRONT;
			}

			// (3) handle saddle vertex
            if (this->mesh()->data(vert).saddle_or_boundary) create_pseudo_source_windows(vert, false);

            // (4) push window lists on the wavefront incident to v into M_LIST_QUEUE
            for(auto e_it = this->mesh()->ve_begin(vert); e_it != this->mesh()->ve_end(vert); ++e_it)
			{
                edge_pointer edge_it = *e_it;
                if (!interval_list_0(edge_it)->empty())
                {
                    m_list_queue.push(interval_list_0(edge_it));
                }
                if (!interval_list_1(edge_it)->empty())
                {
                    m_list_queue.push(interval_list_1(edge_it));
                }
			}


            for(auto f_it = this->mesh()->vf_begin(vert); f_it != this->mesh()->vf_end(vert); ++f_it)
			{
                edge_pointer   edge_it = opposite_edge(*f_it, vert);
                bool two_adjface = (this->mesh()->face_handle(this->mesh()->halfedge_handle(edge_it, hfid0)).idx()>-1)
                        && (this->mesh()->face_handle(this->mesh()->halfedge_handle(edge_it, hfid1)).idx()>-1);
                vertex_pointer vert_it;
                if(two_adjface)
                {
                    face_pointer faceid = opposite_face(edge_it, *f_it);
                    vert_it = opposite_vertex(faceid, edge_it);
                }
                if (!two_adjface || (this->mesh()->data(vert_it).state != VertexState::OUTSIDE))
                {
                    if (!interval_list_0(edge_it)->empty())
                    {
                        m_list_queue.push(interval_list_0(edge_it));
                    }
                    if (!interval_list_1(edge_it)->empty())
                    {
                        m_list_queue.push(interval_list_1(edge_it));
                    }
                }
			}


			// (5) propagate window lists in a FIFO order
			while (!m_list_queue.empty())
			{
				// pop an list from M_LIST_QUEUE
                list_pointer list = m_list_queue.front();

				m_list_queue.pop();

                bool is_boundary = calculate_triangle_parameters(list, Tri);
				if (!is_boundary)
				{
					// propagate the window list using Rule 1 and 2
					wl_left.clear(); wl_right.clear();
					propagate_one_windows_list(list);

					// merge windows lists
					if (!wl_left.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.left_list->empty())
						{
							Tri.left_list->begin()->previous() = wl_left.end();
							wl_left.end()->next() = Tri.left_list->begin();
                            Tri.left_list->begin() = wl_left.begin();
						}
						else
						{
							Tri.left_list->begin() = wl_left.begin();
                            Tri.left_list->end() = wl_left.end();
						}

						// push updated list into M_LIST_QUEUE
                        vertex_pointer v0 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(Tri.left_edge, hfid0));
                        vertex_pointer v1 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(Tri.left_edge, hfid1));
                        if (((this->mesh()->data(v0).state == VertexState::INSIDE)
                             || (this->mesh()->data(v1).state == VertexState::INSIDE))
                                && (!Tri.left_list->empty()))
                        {
							m_list_queue.push(Tri.left_list);
                        }
					}

					if (!wl_right.empty())
					{
						// in VTP, both "PrimeMerge" and "SecondMerge" connect window lists in an order-free way
						if (!Tri.right_list->empty())
						{
							Tri.right_list->end()->next() = wl_right.begin();
							wl_right.begin()->previous() = Tri.right_list->end();
							Tri.right_list->end() = wl_right.end();
						}
						else
						{
							Tri.right_list->begin() = wl_right.begin();
							Tri.right_list->end() = wl_right.end();
						}

						// push updated list into M_LIST_QUEUE
                        vertex_pointer v0 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(Tri.right_edge, hfid0));
                        vertex_pointer v1 = this->mesh()->from_vertex_handle(this->mesh()->halfedge_handle(Tri.right_edge, hfid1));
                        if (((this->mesh()->data(v0).state == VertexState::INSIDE)
                             || (this->mesh()->data(v1).state == VertexState::INSIDE)) && (!Tri.right_list->empty()))
							m_list_queue.push(Tri.right_list);
					}
				}

				list->clear();
			}
			// statistics
			if (m_vertex_queue.size() > m_queue_max_size)
				m_queue_max_size = m_vertex_queue.size();
		}
        idxs.clear();
        for(auto v_it = this->mesh()->vertices_begin(); v_it != this->mesh()->vertices_end(); ++v_it)
        {
            idxs.push_back(MeshVidxfromSub[v_it->idx()]);
            this->ori_mesh->data(this->ori_mesh->vertex_handle(MeshVidxfromSub[v_it->idx()])).geodesic_distance
                    = this->mesh()->data(*v_it).geodesic_distance;
        }
	}

    // construct sub mesh
    inline bool GeodesicAlgorithmExact::construct_submesh(Mesh* sub_mesh, size_t source_idx, Scalar radius)
    {
        std::queue<size_t> vertexlist;
        vertexlist.push(source_idx);
        Vec3 srcp = ori_mesh->point(ori_mesh->vertex_handle(source_idx));
        std::vector<bool> visited(ori_mesh->n_vertices(), false);
        std::vector<bool> added_face(ori_mesh->n_faces(), false);
        SubVidxfromMesh.resize(ori_mesh->n_vertices());
        SubVidxfromMesh.setConstant(-1);
        MeshVidxfromSub.clear();
        visited[source_idx] = true;
        while(!vertexlist.empty())
        {
            size_t vidx = vertexlist.front();
            vertexlist.pop();
            OpenMesh::VertexHandle vh = ori_mesh->vertex_handle(vidx);
            Vec3 vp = ori_mesh->point(vh);
            if((srcp - vp).norm() < radius)
            {
                vertex_pointer new_v = sub_mesh->add_vertex(vp);
                SubVidxfromMesh[vh.idx()] = new_v.idx();
                MeshVidxfromSub.push_back(vh.idx());
                for(auto vv_it = ori_mesh->vv_begin(vh); vv_it != ori_mesh->vv_end(vh); vv_it++)
                {
                    if(!visited[vv_it->idx()])
                    {
                        vertexlist.push(vv_it->idx());
                        visited[vv_it->idx()] = true;
                    }
                }
                for(auto vf_it = ori_mesh->vf_begin(vh); vf_it != ori_mesh->vf_end(vh); vf_it++)
                {
                    halfedge_handle hf = ori_mesh->halfedge_handle(*vf_it);
                    if(!added_face[vf_it->idx()])
                    {
                        vertex_pointer vh = ori_mesh->from_vertex_handle(hf);
                        vertex_pointer nextv = ori_mesh->to_vertex_handle(hf);
                        vertex_pointer thirdv = ori_mesh->to_vertex_handle(ori_mesh->next_halfedge_handle(hf));
                        if(SubVidxfromMesh[vh.idx()] >= 0
                            && SubVidxfromMesh[nextv.idx()] >= 0
                            && SubVidxfromMesh[thirdv.idx()] >= 0)
                        {
                            std::vector<vertex_pointer> vertices;
                            vertices.push_back(sub_mesh->vertex_handle(SubVidxfromMesh[vh.idx()]));
                            vertices.push_back(sub_mesh->vertex_handle(SubVidxfromMesh[nextv.idx()]));
                            vertices.push_back(sub_mesh->vertex_handle(SubVidxfromMesh[thirdv.idx()]));
                            sub_mesh->add_face(vertices);
                            added_face[vf_it->idx()] = true;
                        }
                    }
                }
            }
        }

        sub_mesh->delete_isolated_vertices();
        sub_mesh->garbage_collection();

        if(sub_mesh->n_vertices() > 0)
            return true;
        else
            return false;
    }

	//---------------------- print statistics --------------------------
	inline void GeodesicAlgorithmExact::print_statistics()
	{
		GeodesicAlgorithmBase::print_statistics();

        Scalar memory = sizeof(Interval);

		//std::cout << std::endl;
		std::cout << "Peak number of intervals on wave-front " << m_windows_peak << std::endl;
		std::cout << "uses about " << memory * m_windows_peak / 1e6 << "MB of memory" << std::endl;
		std::cout << "total interval propagation number " << m_windows_propagation << std::endl;
		std::cout << "maximum interval queue size is " << m_queue_max_size << std::endl;
	}
}		//geodesic
#endif
