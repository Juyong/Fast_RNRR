#pragma once
// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <Eigen/Dense>
#include <vector>
#include <nanoflann.h>
#include <igl/median.h>

/// Find self edge median of point cloud
template<typename Derived1>
double FindKnearestMed(Eigen::MatrixBase<Derived1>& X, int nk)
{
	nanoflann::KDTreeAdaptor<Eigen::MatrixBase<Derived1>, 3, nanoflann::metric_L2_Simple> kdtree(X);
	Eigen::VectorXd X_nearest(X.cols());
#pragma omp parallel for
	for (int i = 0; i<X.cols(); i++)
	{
		int* id = new int[nk];
		double *dist = new double[nk];
		kdtree.query(X.col(i).data(), nk, id, dist);
		Eigen::VectorXd k_dist = Eigen::Map<Eigen::VectorXd>(dist, nk);
		igl::median(k_dist.tail(nk - 1), X_nearest[i]);
		delete[]id;
		delete[]dist;
	}
	double med;
	igl::median(X_nearest, med);
	return med;
}
