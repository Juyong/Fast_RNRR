#ifndef GEODESIC_ALGORITHM_EXACT_ELEMENTS
#define GEODESIC_ALGORITHM_EXACT_ELEMENTS

#include "../tools/Types.h"

namespace geodesic {

	class Interval;
	class IntervalList;
	typedef Interval* interval_pointer;
    typedef IntervalList* list_pointer;
    typedef OpenMesh::FaceHandle face_pointer;
    typedef OpenMesh::EdgeHandle edge_pointer;
    typedef OpenMesh::VertexHandle vertex_pointer;
    typedef OpenMesh::HalfedgeHandle halfedge_handle;

	struct Triangle // Components of a face to be propagated
	{
        face_pointer face; // Current Face

        edge_pointer bottom_edge, // Edges
			left_edge,
			right_edge;

        vertex_pointer top_vertex, // Vertices
			left_vertex,
			right_vertex;

        Scalar top_alpha,
			left_alpha,
			right_alpha; // Angles

		list_pointer left_list,
			right_list; // Lists
	};

	class Interval						//interval of the edge
	{
	public:

		Interval() {};
		~Interval() {};

        Scalar& start() { return m_start; };
        Scalar& stop() { return m_stop; };
        Scalar& d() { return m_d; };
        Scalar& pseudo_x() { return m_pseudo_x; };
        Scalar& pseudo_y() { return m_pseudo_y; };

        Scalar& sp() { return m_sp; };

        Scalar& shortest_distance() { return m_shortest_distance; }

		interval_pointer& next() { return m_next; };
		interval_pointer& previous() { return m_previous; };

	private:
        Scalar m_start;						//initial point of the interval on the edge
        Scalar m_stop;
        Scalar m_d;							//distance from the source to the pseudo-source
        Scalar m_pseudo_x;					//coordinates of the pseudo-source in the local coordinate system
        Scalar m_pseudo_y;					//y-coordinate should be always negative

        Scalar m_sp;                        //separating point

        Scalar m_shortest_distance;         //shortest distance from the interval to top_vertex, for numerical precision issue

		interval_pointer m_next;			//pointer to the next interval in the list	
		interval_pointer m_previous;        //pointer to the previous interval in the list
	};

	class IntervalList						//list of the of intervals of the given edge
	{
	public:
        IntervalList() {  /*m_start = NULL; m_edge = NULL;*/ m_sp = -1; m_begin = m_end = NULL; }
		~IntervalList() {};

		void clear() { m_begin = m_end = NULL; }
        void initialize(edge_pointer e) { m_edge = e; }

        vertex_pointer& start_vertex() { return m_start; }
        edge_pointer& edge() { return m_edge; }

        Scalar& sp() { return m_sp; };

        Scalar& pseudo_x() { return m_pseudo_x; };
        Scalar& pseudo_y() { return m_pseudo_y; };

		// List operation
		interval_pointer& begin() { return m_begin; }

		interval_pointer& end() { return m_end; }

		bool empty() { return m_begin == NULL; }

		void push_back(interval_pointer & w)
		{
			if (!m_end)
			{
				w->previous() = NULL;
				w->next() = NULL;
				m_begin = m_end = w;
			}
			else
			{
				w->next() = NULL;
				w->previous() = m_end;
				m_end->next() = w;
				m_end = w;
			}
		}

		void erase(interval_pointer & w)
		{
			if ((w == m_begin) && (w == m_end))
			{
				m_begin = m_end = NULL;
			}
			else if (w == m_begin)
			{
				m_begin = m_begin->next();
				m_begin->previous() = NULL;
			}
			else if (w == m_end)
			{
				m_end = m_end->previous();
				m_end->next() = NULL;
			}
			else
			{
				w->previous()->next() = w->next();
				w->next()->previous() = w->previous();
			}
		}

	private:

        edge_pointer     m_edge;		    //edge that owns this list
        vertex_pointer   m_start;           //vertex from which the interval list starts

		interval_pointer m_begin;
		interval_pointer m_end;

        Scalar m_pseudo_x;
        Scalar m_pseudo_y;

        Scalar m_sp;                        //separating point
	};

}		//geodesic

#endif
