/*
Copyright (C) 2010 Samuele Salti, Federico Tombari, all rights reserved.

This file is part of SHOT. SHOT has been developed by the
Computer Vision Laboratory of the University of Bologna
(http://www.vision.deis.unibo.it)

SHOT is an implementation of the work described in
F. Tombari, S. Salti and L. Di Stefano
"Unique Signatures of Histograms for Local Surface Description"
The 11th IEEE European Conference on Computer Vision (ECCV) 2010

Contacts:
Samuele Salti mailto:samuele.salti@unibo.it
Federico Tombari mailto:federico.tombari@unibo.it


SHOT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SHOT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SHOT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "shot.h"
#include "tools.h"
#include "assert.h"
#include "time.h"
#include <set>
#include <algorithm>
#include "randutil.h"


int findpointsInradius(double* query_point, Mesh* cloud, double radius, std::vector<std::pair<int, double>>& indices)
{
	indices.clear();
	for (auto it = cloud->vertices_begin(); it != cloud->vertices_end(); it++)
	{
		double dist = (query_point[0] - cloud->point(*it)[0]) *(query_point[0] - cloud->point(*it)[0])
			+ (query_point[1] - cloud->point(*it)[1]) *(query_point[1] - cloud->point(*it)[1])
			+ (query_point[2] - cloud->point(*it)[2]) *(query_point[2] - cloud->point(*it)[2]);
		if (dist < radius*radius)
		{
			indices.push_back(std::pair<int, double>((*it).idx(), dist));
		}
	}
	return indices.size();
}

Random3DDetector::Random3DDetector(int nPoints, bool random, double radius, int minNeigh, double borderDistance)
{
	m_featCapacity = nPoints;

	m_numFeat = nPoints;
	m_requestedNumFeat = nPoints;
	m_feat = new Feature3D[nPoints];
	m_minNeigh = minNeigh;
	m_radius = radius;
	m_borderDistance = borderDistance;


	if (random)
        srand(randreseed());
	else
		srand(0xFFFFFFFF);

}


Random3DDetector::~Random3DDetector()
{
	if (m_feat != NULL)
	{
		delete[] m_feat;
		m_feat = NULL;
	}
}



int Random3DDetector::extract(Mesh* cloud, KDtree* index, Feature3D* & feat)
{

	double point[3];
	int randomNum;

	int n_cloud_points = cloud->n_vertices();
	int detectablePointsLeft = n_cloud_points;

	bool *edgePointList;
	if (m_borderDistance > 0.0)
		GetBoundaryPoints(cloud, edgePointList);

	int numFeatures = Min(m_requestedNumFeat, n_cloud_points);
	// 更新feat大小为numFeatures，初始化操作
	updateFeatCapacity(feat, m_featCapacity, numFeatures);

	/* MeshPointsIndex index(cloud);
	vtkIdList* NNpoints;*/

	int extractedPoints = 0;

	std::set<int> ids;

	while (extractedPoints < numFeatures) {


		numFeatures = Min(m_requestedNumFeat, detectablePointsLeft);

		randomNum = randrange(n_cloud_points);
         //std::cout << "extracted points = " << randomNum << " random Num = " << randomNum << std::endl;

		// 取第randomNum个点
		//cloud->GetPoint(randomNum, point);
		point[0] = cloud->point(cloud->vertex_handle(randomNum))[0];
		point[1] = cloud->point(cloud->vertex_handle(randomNum))[1];
		point[2] = cloud->point(cloud->vertex_handle(randomNum))[2];

		if (ids.find(randomNum) != ids.end()) // found
			continue;
		ids.insert(randomNum);


		//Check the minimum number of neighbors (if radius != borderDistance)
		if (m_radius > 0.0 && m_radius != m_borderDistance) {

			// 找到与point距离为m_radius以内的邻域点
			std::vector<std::pair<int, double>> NNpoints;
			//int nNeighbors = index->findInRadius(point, m_radius, NNpoints);
			int nNeighbors = findpointsInradius(point, cloud, m_radius, NNpoints);

			if (nNeighbors - 1 < m_minNeigh) {
				detectablePointsLeft--;
				continue;
			}
		}

		//Check if the point is not too close to the border
		if (m_borderDistance > 0.0) {

			std::vector<std::pair<int, double>> NNpoints;
			//int nNeighbors = index->findInRadius(point, m_radius, NNpoints);
			int nNeighbors = findpointsInradius(point, cloud, m_radius, NNpoints);

			for (int j = 0; j<nNeighbors; j++) {
				if (edgePointList[NNpoints[j].first]) {
					detectablePointsLeft--;
					continue;
				}
			}

			//Check the minimum number of neighbors (if radius = borderDistance)
			if (nNeighbors < m_minNeigh && m_radius == m_borderDistance) {
				detectablePointsLeft--;
				continue;
			}
		}


		m_feat[extractedPoints].x = point[0];
		m_feat[extractedPoints].y = point[1];
		m_feat[extractedPoints].z = point[2];


		m_feat[extractedPoints].scale = -1;
		m_feat[extractedPoints].index = randomNum;

		extractedPoints++;

	}

	m_numFeat = extractedPoints;
	// for test
	std::cout << "extracted points = " << extractedPoints << std::endl;
//    std::ofstream out_file("E:/Non_Rigid_Registration/Toolcode/merge_mesh/test_data/feat_points.txt");
//	for (auto it = ids.begin(); it != ids.end(); it++)
//	{
//		out_file << *it << " " << *it << std::endl;
//	}
//    out_file.close();

	if (m_borderDistance > 0.0)
		delete[] edgePointList;

	feat = m_feat;
	return m_numFeat;
}


void Random3DDetector::updateFeatCapacity(Feature3D * & feat, int & capacity, const int nFeat)
{
	if (capacity < nFeat)
	{
		if (feat != NULL)
		{
			delete[] feat;
			feat = NULL;
		}
	}

	if (feat == NULL && nFeat > 0) {
		feat = new Feature3D[nFeat];
		capacity = nFeat;
	}

}





void getSHOTLocalRF(Mesh *cloud, std::vector<std::pair<int, double>>& NNpoints, double radius, int index, double *rfc)
{
	double originDouble[3];
	originDouble[0] = cloud->point(cloud->vertex_handle(index))[0];
	originDouble[1] = cloud->point(cloud->vertex_handle(index))[1];
	originDouble[2] = cloud->point(cloud->vertex_handle(index))[2];

	int nNeighbours = NNpoints.size();
	//std::cout << "localrf 1 " << std::endl;
	//double V_Vt[9];
	double currPoint[3];
	double* vij = new double[nNeighbours * 3];
	Eigen::Matrix3d covM = Eigen::Matrix3d::Zero();

	// Initialize covariance matrix
	double distance = 0.0;
	double sum = 0.0;

	int validNNpoints = 0;
	for (int ne = 0; ne < nNeighbours; ne++)
	{


		if (NNpoints[ne].first != index) { // perch� il KdTree restituisce anche il punto origine

			currPoint[0] = cloud->point(cloud->vertex_handle(NNpoints[ne].first))[0];
			currPoint[1] = cloud->point(cloud->vertex_handle(NNpoints[ne].first))[1];
			currPoint[2] = cloud->point(cloud->vertex_handle(NNpoints[ne].first))[2];

			// Difference between current point and origin
			vij[validNNpoints * 3 + 0] = currPoint[0] - originDouble[0];
			vij[validNNpoints * 3 + 1] = currPoint[1] - originDouble[1];
			vij[validNNpoints * 3 + 2] = currPoint[2] - originDouble[2];

			distance = radius - sqrt(vij[validNNpoints * 3 + 0] * vij[validNNpoints * 3 + 0] + vij[validNNpoints * 3 + 1] * vij[validNNpoints * 3 + 1] + vij[validNNpoints * 3 + 2] * vij[validNNpoints * 3 + 2]);

			// Multiply vij * vij'
			covM(0, 0) += distance * vij[validNNpoints * 3] * vij[validNNpoints * 3];
			covM(1, 1) += distance * vij[validNNpoints * 3 + 1] * vij[validNNpoints * 3 + 1];
			covM(2, 2) += distance * vij[validNNpoints * 3 + 2] * vij[validNNpoints * 3 + 2];

			double temp = distance * vij[validNNpoints * 3] * vij[validNNpoints * 3 + 1];
			covM(0, 1) += temp;
			covM(1, 0) += temp;
			temp = distance * vij[validNNpoints * 3] * vij[validNNpoints * 3 + 2];
			covM(0, 2) += temp;
			covM(2, 0) += temp;
			temp = distance * vij[validNNpoints * 3 + 1] * vij[validNNpoints * 3 + 2];
			covM(1, 2) += temp;
			covM(2, 1) += temp;

			sum += distance;
			validNNpoints++;

		}
	}
	if (validNNpoints < 5)
	{
		printf("Warning! Neighborhood has less than 5 vertexes. Aborting Local RF computation of feature point with index %d\n", index);
		rfc[0] = 1;
		rfc[1] = 0;
		rfc[2] = 0;

		rfc[3] = 0;
		rfc[4] = 1;
		rfc[5] = 0;

		rfc[6] = 0;
		rfc[7] = 0;
		rfc[8] = 1;

		return;
	}
	covM(0, 0) /= sum;
	covM(0, 1) /= sum;
	covM(0, 2) /= sum;

	covM(1, 0) /= sum;
	covM(1, 1) /= sum;
	covM(1, 2) /= sum;

	covM(2, 0) /= sum;
	covM(2, 1) /= sum;
	covM(2, 2) /= sum;



	Eigen::EigenSolver<Eigen::Matrix3d> es(covM);
	Eigen::Matrix3d eval = es.pseudoEigenvalueMatrix();
	Eigen::Matrix3d evect = es.pseudoEigenvectors();
	for (int i = 0; i < 3; i++)
	{
		for (int j = i + 1; j < 3; j++)
		{
			if (eval(j, j) > eval(i, i))
			{
				double t = eval(j, j);
				eval(j, j) = eval(i, i);
				eval(i, i) = t;
				Eigen::Vector3d temp = evect.col(j);
				evect.col(j) = evect.col(i);
				evect.col(i) = temp;
			}
		}
	}

	int plusNormal = 0, plusTangentDirection1 = 0;
	for (int ne = 0; ne < validNNpoints; ne++)
	{
		double dotProduct = vij[ne * 3] * evect(0, 0) + vij[ne * 3 + 1] * evect(1, 0) + vij[ne * 3 + 2] * evect(2, 0);
		if (dotProduct >= 0)
		{
			plusTangentDirection1++;
		}
		dotProduct = vij[ne * 3] * evect(0, 2) + vij[ne * 3 + 1] * evect(1, 2) + vij[ne * 3 + 2] * evect(2, 2);
		if (dotProduct >= 0)
		{
			plusNormal++;
		}
	}
	// ######## PATCH: directions might still be ambiguous if point density is the same on both emispheres ##########

	bool isNormalStillAmbiguous = false, isTgDirStillAmbiguous = false;

	if (std::abs(plusTangentDirection1 - validNNpoints + plusTangentDirection1) == 0)
		isTgDirStillAmbiguous = true;
	if (std::abs(plusNormal - validNNpoints + plusNormal) == 0)
		isNormalStillAmbiguous = true;

	std::vector<std::pair<double, int>> vij_sorted;
	if (isNormalStillAmbiguous || isTgDirStillAmbiguous) {

		std::pair<double, int> tempPair;

		for (int i = 0; i<validNNpoints; i++) {

			tempPair.first = sqrt(vij[i * 3] * vij[i * 3] + vij[i * 3 + 1] * vij[i * 3 + 1] + vij[i * 3 + 2] * vij[i * 3 + 2]);
			tempPair.second = i;
			vij_sorted.push_back(tempPair);
		}

		std::sort(vij_sorted.begin(), vij_sorted.end());
	}

	if (!isTgDirStillAmbiguous) {

		if (plusTangentDirection1 < validNNpoints - plusTangentDirection1)
		{
			evect(0, 0) *= -1;
			evect(1, 0) *= -1;
			evect(2, 0) *= -1;
		}
	}
	else
	{
		plusTangentDirection1 = 0;
		int pointsToDis = 5; ///std::min(valid_nn_points*2/2+1, 11);
		int nPointsToDis_half = 3;
		int indexToDis = validNNpoints / 2;
		double dotProduct;

		for (int i = -pointsToDis / 2; i <= pointsToDis / 2; i++) {

			dotProduct = vij[vij_sorted[indexToDis - i].second * 3] * evect(0, 0) + vij[vij_sorted[indexToDis - i].second * 3 + 1] * evect(1, 0) + vij[vij_sorted[indexToDis - i].second * 3 + 2] * evect(2, 0);
			if (dotProduct > 0)
				plusTangentDirection1++;
		}

		if (plusTangentDirection1 < nPointsToDis_half) {
			evect(0, 0) *= -1;
			evect(1, 0) *= -1;
			evect(2, 0) *= -1;
		}
	}



	if (!isNormalStillAmbiguous) {
		if (plusNormal < validNNpoints - plusNormal)
		{
			evect(0, 2) *= -1;
			evect(1, 2) *= -1;
			evect(2, 2) *= -1;
		}
	}
	else {
		plusNormal = 0;
		int pointsToDis = 5; ///std::min(valid_nn_points*2/2+1, 11);
		int nPointsToDis_half = 3;
		int indexToDis = validNNpoints / 2;
		double dotProduct;

		for (int i = -pointsToDis / 2; i <= pointsToDis / 2; i++) {
			dotProduct = vij[vij_sorted[indexToDis - i].second * 3] * evect(0, 2) + vij[vij_sorted[indexToDis - i].second * 3 + 1] * evect(1, 2) + vij[vij_sorted[indexToDis - i].second * 3 + 2] * evect(2, 2);
			if (dotProduct > 0)
				plusNormal++;
		}

		if (plusNormal < nPointsToDis_half) {
			evect(0, 2) *= -1;
			evect(1, 2) *= -1;
			evect(2, 2) *= -1;
		}
	}
	rfc[6] = evect(0, 2);
	rfc[7] = evect(1, 2);
	rfc[8] = evect(2, 2);

	rfc[0] = evect(0, 0);
	rfc[1] = evect(1, 0);
	rfc[2] = evect(2, 0);

	// The last reference direction "n3" is the cross product of the other two
    // vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
    Eigen::Vector3d temp = Eigen::Vector3d(rfc[6], rfc[7], rfc[8]).cwiseProduct(Eigen::Vector3d(rfc[0], rfc[1],rfc[2]));
    rfc[6] = temp[0];
    rfc[7] = temp[1];
    rfc[8] = temp[2];


	delete[] vij;
	return;
}


void getSHOTLocalRF(Mesh *cloud, KDtree * pointsIndex, double radius, int index, double *rfc)
{
	double originDouble[3];
	originDouble[0] = cloud->point(cloud->vertex_handle(index))[0];
	originDouble[1] = cloud->point(cloud->vertex_handle(index))[1];
	originDouble[2] = cloud->point(cloud->vertex_handle(index))[2];

	std::vector<std::pair<int, double>> NNpoints;
	// pointsIndex->findInRadius(originDouble, radius, NNpoints);
	findpointsInradius(originDouble, cloud, radius, NNpoints);
	getSHOTLocalRF(cloud, NNpoints, radius, index, rfc);
}




double SHOTDescriptor::sRGB_LUT[256] = { -1 };

double SHOTDescriptor::sXYZ_LUT[4000] = { -1 };


void SHOTDescriptor::RGB2CIELAB(unsigned char R, unsigned char G, unsigned char B, double &L, double &A, double &B2) {

	if (sRGB_LUT[0] < 0)
	{
		for (int i = 0; i < 256; i++)
		{
			double f = i / 255.0f;
			if (f > 0.04045)
				sRGB_LUT[i] = (double)pow((f + 0.055) / 1.055, 2.4);
			else
				sRGB_LUT[i] = f / 12.92;
		}

		for (int i = 0; i < 4000; i++)
		{
			double f = i / 4000.0;
			if (f > 0.008856)
				sXYZ_LUT[i] = pow(f, (double)0.3333);
			else
				sXYZ_LUT[i] = (7.787 * f) + (16.0 / 116.0);
		}
	}

	double fr = sRGB_LUT[R];
	double fg = sRGB_LUT[G];
	double fb = sRGB_LUT[B];

	/*if (fr > 0.04045)
	fr = (float) pow((fr + 0.055) / 1.055, 2.4);
	else
	fr = fr / 12.92;

	if (fg > 0.04045)
	fg = (float) pow((fg + 0.055) / 1.055, 2.4);
	else
	fg = fg / 12.92;

	if (fb > 0.04045)
	fb = (float) pow((fb + 0.055) / 1.055, 2.4);
	else
	fb = fb / 12.92;*/

	// Use white = D65
	const double x = fr * 0.412453 + fg * 0.357580 + fb * 0.180423;
	const double y = fr * 0.212671 + fg * 0.715160 + fb * 0.072169;
	const double z = fr * 0.019334 + fg * 0.119193 + fb * 0.950227;

	double vx = x / 0.95047;
	double vy = y;
	double vz = z / 1.08883;

	//printf("vx:%f vy:%f vz:%f\n", vx, vy, vz);

	vx = sXYZ_LUT[int(vx * 4000)];
	vy = sXYZ_LUT[int(vy * 4000)];
	vz = sXYZ_LUT[int(vz * 4000)];
	/*if (vx > 0.008856)
	vx = pow(vx, (float)0.3333);
	else
	vx = (7.787 * vx) + (16.0 / 116.0);

	if (vy > 0.008856)
	vy = (float) pow(vy, (float)0.3333);
	else
	vy = (7.787 * vy) + (16.0 / 116.0);

	if (vz > 0.008856)
	vz = (float) pow(vz, (float)0.3333);
	else
	vz = (7.787 * vz) + (16.0 / 116.0);*/

	L = 116.0 * vy - 16.0;
	if (L>100)
		L = 100;

	A = 500.0 * (vx - vy);
	if (A>120)
		A = 120;
	else if (A<-120)
		A = -120;

	B2 = 200.0 * (vy - vz);
	if (B2>120)
		B2 = 120;
	else if (B2<-120)
		B2 = -120;




}

SHOTDescriptor::SHOTDescriptor(SHOTParams params)
	: m_maxAngularSectors(32) //19.9.17 set to 32
	, m_nSect(32)
	, m_descCapacity(0)
	, m_desc(NULL)
{
	setParams(params, 1000);
}

SHOTDescriptor::~SHOTDescriptor()
{
	if (m_desc != NULL)
	{
		delete m_desc[0];
		delete m_desc;
		m_desc = NULL;
	}
}

void SHOTDescriptor::setParams(SHOTParams params, int nPoints)
{
	m_params = params;

	m_k = 32;	//always set to 32 "husks"

	int descLength = (params.describeShape) ? m_k*(params.shapeBins + 1) : 0;
	descLength += (params.describeColor) ? m_k*(params.colorBins + 1) : 0;

	// m_desc = ndesc * descLength
	updateDescCapacity(m_desc, m_descCapacity, nPoints, descLength);

	m_radius3_4 = (params.radius * 3) / 4;
	m_radius1_4 = params.radius / 4;
	m_radius1_2 = params.radius / 2;

}

void SHOTDescriptor::updateDescCapacity(double ** & desc, int & capacity, int nDesc, int descLength)
{
	//std::cout << "desclength = " << descLength << " m_desclen = " << m_descLength << std::endl;
	if ((capacity < nDesc) || (descLength > m_descLength))
	{
		if (desc != NULL)
		{
			delete[] desc[0];
			delete[] desc;
			desc = NULL;
		}
	}


	if (desc == NULL && nDesc > 0)
	{
		desc = new double*[nDesc];
		desc[0] = new double[nDesc*descLength];

		for (int i = 1; i < nDesc; ++i)
			desc[i] = desc[i - 1] + descLength;
		capacity = nDesc;

	}

	for (int i = 0; i < nDesc; i++)
		for (int j = 0; j < descLength; j++)
			desc[i][j] = 0.0;

	m_descLength = descLength;
}


void SHOTDescriptor::interpolateSingleChannel(Mesh* cloud, std::vector<std::pair<int, double>>& NNpoints, const std::vector<double> & distances, const double centralPoint[3], double rf[9], std::vector<double> & binDistance, int nr_bins, double* shot)
{

	for (int i_idx = 0; i_idx < NNpoints.size(); ++i_idx)
	{
		double point[3];
		point[0] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[0];
		point[1] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[1];
		point[2] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[2];

		double distance = sqrt(distances[i_idx]); //sqrt(distance_sqr);

		double delta[3] = { point[0] - centralPoint[0], point[1] - centralPoint[1], point[2] - centralPoint[2] };

		// Compute the Euclidean norm
		if (areEquals(distance, 0.0))
			continue;

		double xInFeatRef = (delta[0] * rf[0] + delta[1] * rf[1] + delta[2] * rf[2]);
		double yInFeatRef = (delta[0] * rf[3] + delta[1] * rf[4] + delta[2] * rf[5]);
		double zInFeatRef = (delta[0] * rf[6] + delta[1] * rf[7] + delta[2] * rf[8]);


		// To avoid numerical problems afterwards
		if (std::abs(yInFeatRef) < 1E-30)
			yInFeatRef = 0;
		if (std::abs(xInFeatRef) < 1E-30)
			xInFeatRef = 0;
		if (std::abs(zInFeatRef) < 1E-30)
			zInFeatRef = 0;



		unsigned char bit4 = ((yInFeatRef > 0) || ((yInFeatRef == 0.0) && (xInFeatRef < 0))) ? 1 : 0;
		unsigned char bit3 = ((xInFeatRef > 0) || ((xInFeatRef == 0.0) && (yInFeatRef > 0))) ? !bit4 : bit4;

		assert(bit3 == 0 || bit3 == 1);

		int desc_index = (bit4 << 3) + (bit3 << 2);

		desc_index = desc_index << 1;

		if ((xInFeatRef * yInFeatRef > 0) || (xInFeatRef == 0.0))
			desc_index += (std::abs(xInFeatRef) >= std::abs(yInFeatRef)) ? 0 : 4;
		else
			desc_index += (std::abs(xInFeatRef) > std::abs(yInFeatRef)) ? 4 : 0;

		desc_index += zInFeatRef > 0 ? 1 : 0;

		// 2 RADII
		desc_index += (distance > m_radius1_2) ? 2 : 0;

		int step_index = static_cast<int>(floor(binDistance[i_idx] + 0.5));
		int volume_index = desc_index * (nr_bins + 1);

		//Interpolation on the cosine (adjacent bins in the histogram)
		binDistance[i_idx] -= step_index;
		double intWeight = (1 - std::abs(binDistance[i_idx]));

		//19.9.17: (nr_bins + 1) used instead of nr_bins
		if (binDistance[i_idx] > 0)
			shot[volume_index + ((step_index + 1) % (nr_bins + 1))] += binDistance[i_idx];
		else
			shot[volume_index + ((step_index + nr_bins) % (nr_bins + 1))] += -binDistance[i_idx];

		//Interpolation on the distance (adjacent husks)
		if (distance > m_radius1_2) {	//external sphere

			double radiusDistance = (distance - m_radius3_4) / m_radius1_2;

			if (distance > m_radius3_4)	//most external sector, votes only for itself
				intWeight += 1 - radiusDistance;	//peso=1-d
			else {	//3/4 of radius, votes also for the internal sphere
				intWeight += 1 + radiusDistance;
				shot[(desc_index - 2) * (nr_bins + 1) + step_index] -= radiusDistance;
			}
		}
		else {	//internal sphere

			double radiusDistance = (distance - m_radius1_4) / m_radius1_2;

			if (distance < m_radius1_4)	//most internal sector, votes only for itself
				intWeight += 1 + radiusDistance;	//weight=1-d
			else {	//3/4 of radius, votes also for the external sphere
				intWeight += 1 - radiusDistance;
				shot[(desc_index + 2) * (nr_bins + 1) + step_index] += radiusDistance;
			}

		}

		//Interpolation on the inclination (adjacent vertical volumes)
		double inclinationCos = zInFeatRef / distance;
		if (inclinationCos < -1.0) inclinationCos = -1.0;
		if (inclinationCos > 1.0) inclinationCos = 1.0;

		double inclination = acos(inclinationCos);

		assert(inclination >= 0.0 && inclination <= DEG_180_TO_RAD);

		if (inclination > DEG_90_TO_RAD || (std::abs(inclination - DEG_90_TO_RAD) < 1e-30 && zInFeatRef <= 0)) {

			double inclinationDistance = (inclination - DEG_135_TO_RAD) / DEG_90_TO_RAD;
			if (inclination > DEG_135_TO_RAD)
				intWeight += 1 - inclinationDistance;
			else {
				intWeight += 1 + inclinationDistance;
				assert((desc_index + 1) * (nr_bins + 1) + step_index >= 0 && (desc_index + 1) * (nr_bins + 1) + step_index < m_descLength);
				shot[(desc_index + 1) * (nr_bins + 1) + step_index] -= inclinationDistance;
			}

		}
		else {

			double inclinationDistance = (inclination - DEG_45_TO_RAD) / DEG_90_TO_RAD;
			if (inclination < DEG_45_TO_RAD)
				intWeight += 1 + inclinationDistance;
			else {
				intWeight += 1 - inclinationDistance;
				assert((desc_index - 1) * (nr_bins + 1) + step_index >= 0 && (desc_index - 1) * (nr_bins + 1) + step_index < m_descLength);
				shot[(desc_index - 1) * (nr_bins + 1) + step_index] += inclinationDistance;
			}
		}

		if (yInFeatRef != 0.0 || xInFeatRef != 0.0)
		{
			//Interpolation on the azimuth (adjacent horizontal volumes)
			double azimuth = atan2(yInFeatRef, xInFeatRef);

			int sel = desc_index >> 2;
			double angularSectorSpan = DEG_45_TO_RAD;
			double angularSectorStart = -DEG_168_TO_RAD;

			double azimuthDistance = (azimuth - (angularSectorStart + angularSectorSpan*sel)) / angularSectorSpan;

			azimuthDistance = std::max(-0.5, std::min(azimuthDistance, 0.5));

			assert((azimuthDistance < 0.5 || areEquals(azimuthDistance, 0.5)) && (azimuthDistance > -0.5 || areEquals(azimuthDistance, -0.5)));

			if (azimuthDistance > 0) {
				intWeight += 1 - azimuthDistance;
				int interp_index = (desc_index + 4) % m_maxAngularSectors;
				assert(interp_index * (nr_bins + 1) + step_index >= 0 && interp_index * (nr_bins + 1) + step_index < m_descLength);
				shot[interp_index * (nr_bins + 1) + step_index] += azimuthDistance;
			}
			else {
				int interp_index = (desc_index - 4 + m_maxAngularSectors) % m_maxAngularSectors;
				assert(interp_index * (nr_bins + 1) + step_index >= 0 && interp_index * (nr_bins + 1) + step_index < m_descLength);
				intWeight += 1 + azimuthDistance;
				shot[interp_index * (nr_bins + 1) + step_index] -= azimuthDistance;
			}

		}

		assert(volume_index + step_index >= 0 && volume_index + step_index < m_descLength);
		shot[volume_index + step_index] += intWeight;
	}
}

//Quadrilinear interpolation; used when color and shape descriptions are both activated
void SHOTDescriptor::interpolateDoubleChannel(Mesh* cloud, std::vector<std::pair<int, double>>& NNpoints, const std::vector<double> & distances, const double centralPoint[3], double rf[9], std::vector<double> & binDistanceShape, std::vector<double> & binDistanceColor, double* shot)
{
	int shapeToColorStride = 32 * (m_params.shapeBins + 1);

	for (int i_idx = 0; i_idx < NNpoints.size(); ++i_idx)
	{
		double point[3];
		point[0] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[0];
		point[1] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[1];
		point[2] = cloud->point(cloud->vertex_handle(NNpoints[i_idx].first))[2];

		double distance = sqrt(distances[i_idx]); //sqrt(distance_sqr);

		double delta[3] = { point[0] - centralPoint[0], point[1] - centralPoint[1], point[2] - centralPoint[2] };

		// Compute the Euclidean norm

		if (areEquals(distance, 0.0))
			continue;

		double xInFeatRef = (delta[0] * rf[0] + delta[1] * rf[1] + delta[2] * rf[2]);
		double yInFeatRef = (delta[0] * rf[3] + delta[1] * rf[4] + delta[2] * rf[5]);
		double zInFeatRef = (delta[0] * rf[6] + delta[1] * rf[7] + delta[2] * rf[8]);

		// To avoid numerical problems afterwards
		if (std::abs(yInFeatRef) < 1E-30)
			yInFeatRef = 0;
		if (std::abs(xInFeatRef) < 1E-30)
			xInFeatRef = 0;
		if (std::abs(zInFeatRef) < 1E-30)
			zInFeatRef = 0;



		unsigned char bit4 = ((yInFeatRef > 0) || ((yInFeatRef == 0.0) && (xInFeatRef < 0))) ? 1 : 0;
		unsigned char bit3 = ((xInFeatRef > 0) || ((xInFeatRef == 0.0) && (yInFeatRef > 0))) ? !bit4 : bit4;

		assert(bit3 == 0 || bit3 == 1);

		int desc_index = (bit4 << 3) + (bit3 << 2);

		desc_index = desc_index << 1;

		if ((xInFeatRef * yInFeatRef > 0) || (xInFeatRef == 0.0))
			desc_index += (std::abs(xInFeatRef) >= std::abs(yInFeatRef)) ? 0 : 4;
		else
			desc_index += (std::abs(xInFeatRef) > std::abs(yInFeatRef)) ? 4 : 0;

		desc_index += zInFeatRef > 0 ? 1 : 0;

		// 2 RADII
		desc_index += (distance > m_radius1_2) ? 2 : 0;

		int step_index_shape = static_cast<int>(floor(binDistanceShape[i_idx] + 0.5));
		int step_index_color = static_cast<int>(floor(binDistanceColor[i_idx] + 0.5));

		int volume_index_shape = desc_index * (m_params.shapeBins + 1);
		int volume_index_color = shapeToColorStride + desc_index * (m_params.colorBins + 1);

		//Interpolation on the cosine (adjacent bins in the histrogram)
		binDistanceShape[i_idx] -= step_index_shape;
		binDistanceColor[i_idx] -= step_index_color;

		double intWeightShape = (1 - std::abs(binDistanceShape[i_idx]));
		double intWeightColor = (1 - std::abs(binDistanceColor[i_idx]));

		//19.9.17: (m_params.shapeBins + 1) used instead of m_params.shapeBins (same for colorBins)
		if (binDistanceShape[i_idx] > 0)
			shot[volume_index_shape + ((step_index_shape + 1) % (m_params.shapeBins + 1))] += binDistanceShape[i_idx];
		else
			shot[volume_index_shape + ((step_index_shape + m_params.shapeBins) % (m_params.shapeBins + 1))] -= binDistanceShape[i_idx];

		if (binDistanceColor[i_idx] > 0)
			shot[volume_index_color + ((step_index_color + 1) % (m_params.colorBins + 1))] += binDistanceColor[i_idx];
		else
			shot[volume_index_color + ((step_index_color + m_params.colorBins) % (m_params.colorBins + 1))] -= binDistanceColor[i_idx];

		//Interpolation on the distance (adjacent husks)
		if (distance > m_radius1_2) {	//external sphere

			double radiusDistance = (distance - m_radius3_4) / m_radius1_2;

			if (distance > m_radius3_4)	//most external sector, votes only for itself
			{
				intWeightShape += 1 - radiusDistance;	//weight=1-d
				intWeightColor += 1 - radiusDistance;	//weight=1-d
			}
			else {	//3/4 of radius, votes also for the internal sphere
				intWeightShape += 1 + radiusDistance;
				intWeightColor += 1 + radiusDistance;
				shot[(desc_index - 2) * (m_params.shapeBins + 1) + step_index_shape] -= radiusDistance;
				shot[shapeToColorStride + (desc_index - 2) * (m_params.colorBins + 1) + step_index_color] -= radiusDistance;
			}
		}
		else {	//internal sphere

			double radiusDistance = (distance - m_radius1_4) / m_radius1_2;

			if (distance < m_radius1_4)	//most internal sector, votes only for itself
			{
				intWeightShape += 1 + radiusDistance;
				intWeightColor += 1 + radiusDistance;	//weight=1-d
			}
			else {	//3/4 of radius, votes also for the external sphere
				intWeightShape += 1 - radiusDistance;	//weight=1-d
				intWeightColor += 1 - radiusDistance;	//weight=1-d
				shot[(desc_index + 2) * (m_params.shapeBins + 1) + step_index_shape] += radiusDistance;
				shot[shapeToColorStride + (desc_index + 2) * (m_params.colorBins + 1) + step_index_color] += radiusDistance;
			}

		}

		//Interpolation on the inclination (adjacent vertical volumes)
		double inclinationCos = zInFeatRef / distance;
		if (inclinationCos < -1.0) inclinationCos = -1.0;
		if (inclinationCos > 1.0) inclinationCos = 1.0;

		double inclination = acos(inclinationCos);

		assert(inclination >= 0.0 && inclination <= DEG_180_TO_RAD);

		if (inclination > DEG_90_TO_RAD || (std::abs(inclination - DEG_90_TO_RAD) < 1e-30 && zInFeatRef <= 0)) {

			double inclinationDistance = (inclination - DEG_135_TO_RAD) / DEG_90_TO_RAD;
			if (inclination > DEG_135_TO_RAD)
			{
				intWeightShape += 1 - inclinationDistance;
				intWeightColor += 1 - inclinationDistance;
			}
			else {
				intWeightShape += 1 + inclinationDistance;
				intWeightColor += 1 + inclinationDistance;
				assert((desc_index + 1) * (m_params.shapeBins + 1) + step_index_shape >= 0 && (desc_index + 1) * (m_params.shapeBins + 1) + step_index_shape < m_descLength);
				assert(shapeToColorStride + (desc_index + 1) * (m_params.colorBins + 1) + step_index_color >= 0 && shapeToColorStride + (desc_index + 1) * (m_params.colorBins + 1) + step_index_color < m_descLength);
				shot[(desc_index + 1) * (m_params.shapeBins + 1) + step_index_shape] -= inclinationDistance;
				shot[shapeToColorStride + (desc_index + 1) * (m_params.colorBins + 1) + step_index_color] -= inclinationDistance;
			}

		}
		else {

			double inclinationDistance = (inclination - DEG_45_TO_RAD) / DEG_90_TO_RAD;
			if (inclination < DEG_45_TO_RAD)
			{
				intWeightShape += 1 + inclinationDistance;
				intWeightColor += 1 + inclinationDistance;
			}
			else {
				intWeightShape += 1 - inclinationDistance;
				intWeightColor += 1 - inclinationDistance;
				assert((desc_index - 1) * (m_params.shapeBins + 1) + step_index_shape >= 0 && (desc_index - 1) * (m_params.shapeBins + 1) + step_index_shape < m_descLength);
				assert(shapeToColorStride + (desc_index - 1) * (m_params.colorBins + 1) + step_index_color >= 0 && shapeToColorStride + (desc_index - 1) * (m_params.colorBins + 1) + step_index_color < m_descLength);
				shot[(desc_index - 1) * (m_params.shapeBins + 1) + step_index_shape] += inclinationDistance;
				shot[shapeToColorStride + (desc_index - 1) * (m_params.colorBins + 1) + step_index_color] += inclinationDistance;
			}
		}

		if (yInFeatRef != 0.0 || xInFeatRef != 0.0)
		{
			//Interpolation on the azimuth (adjacent horizontal volumes)
			double azimuth = atan2(yInFeatRef, xInFeatRef);

			int sel = desc_index >> 2;
			double angularSectorSpan = DEG_45_TO_RAD;
			double angularSectorStart = -DEG_168_TO_RAD;

			double azimuthDistance = (azimuth - (angularSectorStart + angularSectorSpan*sel)) / angularSectorSpan;

			azimuthDistance = std::max(-0.5, std::min(azimuthDistance, 0.5));

			assert((azimuthDistance < 0.5 || areEquals(azimuthDistance, 0.5)) && (azimuthDistance > -0.5 || areEquals(azimuthDistance, -0.5)));

			if (azimuthDistance > 0) {
				intWeightShape += 1 - azimuthDistance;
				intWeightColor += 1 - azimuthDistance;
				int interp_index = (desc_index + 4) % m_maxAngularSectors;
				assert(interp_index * (m_params.shapeBins + 1) + step_index_shape >= 0 && interp_index * (m_params.shapeBins + 1) + step_index_shape < m_descLength);
				assert(shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color >= 0 && shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color < m_descLength);
				shot[interp_index * (m_params.shapeBins + 1) + step_index_shape] += azimuthDistance;
				shot[shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color] += azimuthDistance;
			}
			else {
				int interp_index = (desc_index - 4 + m_maxAngularSectors) % m_maxAngularSectors;
				intWeightShape += 1 + azimuthDistance;
				intWeightColor += 1 + azimuthDistance;
				assert(interp_index * (m_params.shapeBins + 1) + step_index_shape >= 0 && interp_index * (m_params.shapeBins + 1) + step_index_shape < m_descLength);
				assert(shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color >= 0 && shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color < m_descLength);
				shot[interp_index * (m_params.shapeBins + 1) + step_index_shape] -= azimuthDistance;
				shot[shapeToColorStride + interp_index * (m_params.colorBins + 1) + step_index_color] -= azimuthDistance;
			}

		}

		assert(volume_index_shape + step_index_shape >= 0 && volume_index_shape + step_index_shape < m_descLength);
		assert(volume_index_color + step_index_color >= 0 && volume_index_color + step_index_color < m_descLength);
		shot[volume_index_shape + step_index_shape] += intWeightShape;
		shot[volume_index_color + step_index_color] += intWeightColor;
	}
}

void SHOTDescriptor::computeDesc(SHOTParams params, Mesh *cloud, KDtree* index, Feature3D *feat, double **desc, int startIndex, int endIndex) {


	double centralPoint[3];
	double point_f[3];
	double point[3], x, y, z;

	std::cout << "compute desc1 startIndex = " << startIndex << " endIndex = " << endIndex << std::endl;
	for (int i = startIndex; i<endIndex; i++) {

		std::vector<double> binDistanceShape;
		std::vector<double> binDistanceColor;
		std::vector<double> distances;

		centralPoint[0] = cloud->point(cloud->vertex_handle(feat[i].index))[0];
		centralPoint[1] = cloud->point(cloud->vertex_handle(feat[i].index))[1];
		centralPoint[2] = cloud->point(cloud->vertex_handle(feat[i].index))[2];

		point_f[0] = double(centralPoint[0]);
		point_f[1] = double(centralPoint[1]);
		point_f[2] = double(centralPoint[2]);
		std::vector<std::pair<int, double>> NNpoints;
		if (areEquals(m_params.localRFradius, m_params.radius)) {
			findpointsInradius(centralPoint, cloud, m_params.radius, NNpoints);
			getSHOTLocalRF(cloud, NNpoints, m_params.localRFradius, feat[i].index, feat[i].rf);
		}
		else {
			getSHOTLocalRF(cloud, index, m_params.localRFradius, feat[i].index, feat[i].rf);
			findpointsInradius(centralPoint, cloud, m_params.radius, NNpoints);
		}
		int nNeighbors = NNpoints.size();

		if (nNeighbors <  m_params.minNeighbors)
		{
			printf("Warning! Neighborhood has less than 5 vertexes. Aborting description of feature point %d, index %d\n", i, feat[i].index);
			continue;
		}
		//Compute all distances once; useful in case both color and shape are described, so they are computed only once instead than twice
		distances.resize(nNeighbors);
		for (int j = 0; j<nNeighbors; j++) {

			int id = NNpoints[j].first;
			point[0] = cloud->point(cloud->vertex_handle(id))[0];
			point[1] = cloud->point(cloud->vertex_handle(id))[1];
			point[2] = cloud->point(cloud->vertex_handle(id))[2];

			x = point[0] - centralPoint[0];
			y = point[1] - centralPoint[1];
			z = point[2] - centralPoint[2];
			distances[j] = x*x + y*y + z*z;
		}

		//Compute current shape bin values
		if (params.describeShape)
		{
			binDistanceShape.resize(nNeighbors);


			for (int j = 0; j<nNeighbors; j++) {

				if (areEquals(distances[j], 0.0))
					continue;

				OpenMesh::VertexHandle idh = cloud->vertex_handle(NNpoints[j].first);

				double cosineDesc = feat[i].rf[6] * cloud->normal(idh)[0] + feat[i].rf[7] * cloud->normal(idh)[1] + feat[i].rf[8] * cloud->normal(idh)[2];

				if (cosineDesc > 1.0) cosineDesc = 1.0;
				if (cosineDesc < -1.0) cosineDesc = -1.0;

				binDistanceShape[j] = ((1.0 + cosineDesc) * params.shapeBins) / 2;
			}
		}

		//Color Description
		if (params.describeColor)
		{

			if (!cloud->has_vertex_colors()) {
				printf("Error: no color array found in the processed 3D data..\n");
				throw std::runtime_error("no color array found in the processed 3D data");
			}

			binDistanceColor.resize(nNeighbors);

			unsigned char r_RGB = cloud->color(cloud->vertex_handle(feat[i].index))[0];
			unsigned char g_RGB = cloud->color(cloud->vertex_handle(feat[i].index))[1];
			unsigned char b_RGB = cloud->color(cloud->vertex_handle(feat[i].index))[2];	//componenti RGB del punto centrale

			double l_LAB, a_LAB, b_LAB, l_LAB_norm, a_LAB_norm, b_LAB_norm;

			SHOTDescriptor::RGB2CIELAB(r_RGB, g_RGB, b_RGB, l_LAB, a_LAB, b_LAB);

			l_LAB_norm = l_LAB / 100.0;
			a_LAB_norm = a_LAB / 120.0;
			b_LAB_norm = b_LAB / 120.0;		//normalized LAB components (0<L<1, -1<a<1, -1<b<1)


			for (int j = 0; j<nNeighbors; j++) {	//per ogni vicino

				OpenMesh::VertexHandle idh = cloud->vertex_handle(NNpoints[j].first);

				if (areEquals(distances[j], 0.0))	//same point as the central one
					continue;

				unsigned char r_RGB_cur = cloud->point(idh)[0];
				unsigned char g_RGB_cur = cloud->point(idh)[1];
				unsigned char b_RGB_cur = cloud->point(idh)[2];

				double l_LAB_cur, a_LAB_cur, b_LAB_cur, l_LAB_norm_cur, a_LAB_norm_cur, b_LAB_norm_cur;

				SHOTDescriptor::RGB2CIELAB(r_RGB_cur, g_RGB_cur, b_RGB_cur, l_LAB_cur, a_LAB_cur, b_LAB_cur);

				l_LAB_norm_cur = l_LAB_cur / 100.0;
				a_LAB_norm_cur = a_LAB_cur / 120.0;
				b_LAB_norm_cur = b_LAB_cur / 120.0;		//normalized LAB components (0<L<1, -1<a<1, -1<b<1)


				double colorDesc = (std::abs(l_LAB_norm - l_LAB_norm_cur) + ((std::abs(a_LAB_norm - a_LAB_norm_cur) + std::abs(b_LAB_norm - b_LAB_norm_cur)) / 2)) / 3;


				if (colorDesc > 1.0) colorDesc = 1.0;
				if (colorDesc < 0.0) colorDesc = 0.0;

				binDistanceColor[j] = colorDesc * params.colorBins;

			}
		}

		if (params.describeShape && params.describeColor)
		{
			interpolateDoubleChannel(cloud, NNpoints, distances, centralPoint, feat[i].rf, binDistanceShape, binDistanceColor, desc[i]);
		}
		else if (params.describeColor)
		{
			interpolateSingleChannel(cloud, NNpoints, distances, centralPoint, feat[i].rf, binDistanceColor, m_params.colorBins, desc[i]);
		}
		else
		{
			interpolateSingleChannel(cloud, NNpoints, distances, centralPoint, feat[i].rf, binDistanceShape, m_params.shapeBins, desc[i]);
		}

		double accNorm = 0;
		for (int j = 0; j< m_descLength; j++)
			accNorm += desc[i][j] * desc[i][j];
		accNorm = sqrt(accNorm);
		for (int j = 0; j< m_descLength; j++)
			desc[i][j] /= accNorm;

	}
}

void SHOTDescriptor::describe(Mesh *cloud, KDtree* index, Feature3D *&feat, double **&desc, int nPoints)
{
	if (m_params.describeShape && !cloud->has_vertex_normals())
		cloud->update_vertex_normals();

	if (m_params.nThreads == 0)
		m_params.nThreads = 1;//omp_get_max_threads();

	setParams(m_params, nPoints);
	int step = nPoints / m_params.nThreads;
	computeDesc(m_params, cloud, index, feat, m_desc, (m_params.nThreads - 1)*step, nPoints);

	desc = m_desc;

}
