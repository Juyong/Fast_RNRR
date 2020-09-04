#pragma once
#include <nanoflann.hpp>
#include "Types.h"

namespace nanoflann {
	// KD-tree adaptor for working with data directly stored in an Eigen Matrix, without duplicating the data storage.
	// This code is adapted from the KDTreeEigenMatrixAdaptor class of nanoflann.hpp
	template <class MatrixType, int DIM = 3, class Distance = nanoflann::metric_L2, typename IndexType = int>
	struct KDTreeAdaptor {
		typedef KDTreeAdaptor<MatrixType, DIM, Distance> self_t;
		typedef typename MatrixType::Scalar              num_t;
		typedef typename Distance::template traits<num_t, self_t>::distance_t metric_t;
		typedef KDTreeSingleIndexAdaptor< metric_t, self_t, DIM, IndexType>  index_t;
		index_t* index;
		KDTreeAdaptor(const MatrixType &mat, const int leaf_max_size = 10) : m_data_matrix(mat) {
			const size_t dims = mat.rows();
            index = new index_t(dims, *this, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size, dims));
			index->buildIndex();
		}
		~KDTreeAdaptor() { delete index; }
		const MatrixType &m_data_matrix;
		/// Query for the num_closest closest points to a given point (entered as query_point[0:dim-1]).
		inline void query(const num_t *query_point, const size_t num_closest, IndexType *out_indices, num_t *out_distances_sq) const {
			nanoflann::KNNResultSet<typename MatrixType::Scalar, IndexType> resultSet(num_closest);
			resultSet.init(out_indices, out_distances_sq);
			index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
		}
		/// Query for the closest points to a given point (entered as query_point[0:dim-1]).
		inline IndexType closest(const num_t *query_point, num_t& out_distance_sq) const {
			IndexType out_indices;
			num_t out_distances_sq;
			query(query_point, 1, &out_indices, &out_distances_sq);
			out_distance_sq = sqrt(out_distances_sq);
			return out_indices;
		}
        // this function has problem
		inline IndexType findInRadius(const num_t *query_point, const num_t radius, std::vector<std::pair<IndexType, num_t>>& out_idxdist) const
		{
			out_idxdist.clear();
            index->radiusSearch(query_point, radius, out_idxdist, nanoflann::SearchParams());
			std::cout << "nanoflann out_idxdist.size = " << out_idxdist.size() << std::endl;
			return out_idxdist.size();
		};
		const self_t & derived() const { return *this; }
		self_t & derived() { return *this; }
		inline size_t kdtree_get_point_count() const { return m_data_matrix.cols(); }
		/// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
		inline num_t kdtree_distance(const num_t *p1, const size_t idx_p2, size_t size) const {
			num_t s = 0;
			for (size_t i = 0; i<size; i++) {
				const num_t d = p1[i] - m_data_matrix.coeff(i, idx_p2);
				s += d*d;
			}
			return s;
		}
		/// Returns the dim'th component of the idx'th point in the class:
		inline num_t kdtree_get_pt(const size_t idx, int dim) const {
			return m_data_matrix.coeff(dim, idx);
		}
		/// Optional bounding-box computation: return false to default to a standard bbox computation loop.
		template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
	};


    struct StdVecOfEigenVector3dRefAdapter
    {
    public:
        StdVecOfEigenVector3dRefAdapter(const std::vector<Vector3>& pps) : pts(pps) {};

        const std::vector<Vector3>& pts;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return pts.size(); }

        // Returns the dim'th component of the idx'th point in the class:
        // Since this is inlined and the "dim" argument is typically an immediate value, the
        //  "if/else's" are actually solved at compile time.
        inline Scalar kdtree_get_pt(const size_t idx, int dim) const
        {
            if (dim == 0) return pts[idx].x();
            else if (dim == 1) return pts[idx].y();
            else return pts[idx].z();
        }

        // Optional bounding-box computation: return false to default to a standard bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
        //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

    };
}


typedef nanoflann::KDTreeAdaptor<MatrixXX, -1, nanoflann::metric_L2_Simple> KDtree;
