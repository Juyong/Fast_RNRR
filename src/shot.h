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

#ifndef SHOT_H
#define SHOT_H

#include "MeshDefinition.h"
#include "nanoflann.h"


#define EQUAL_ZERO_TH 1e-15f

#define DEG_45_TO_RAD 0.78539816339744830961566084581988f
#define DEG_90_TO_RAD 1.5707963267948966192313216916398f
#define DEG_135_TO_RAD 2.3561944901923449288469825374596f
#define DEG_168_TO_RAD 2.7488935718910690836548129603696f

#define DEG_180_TO_RAD 3.1415926535897932384626433832795f
#define DEG_360_TO_RAD 6.283185307179586476925286766558f

/*!
* \brief Check if val1 and val2 are equals.
*
* \param val1 first number to check.
* \param val2 second number to check.
* \return true if val1 is equal to val2, false otherwise.
*/
inline bool areEquals(double val1, double val2, double zeroDoubleEps = EQUAL_ZERO_TH){
    return ( (val1 - val2)<zeroDoubleEps && (val1 - val2)>-zeroDoubleEps );
};



struct Feature3D
{
	double x,y,z;
	double scale;
    int index;
	double rf[9];
	double score;
};

//class MeshPointsIndex
//{
//protected:
//    vtkPointLocator*		m_index;
//    vtkIdList*				m_pointsList;
//    vtkPolyData*			m_polyData;
//
//public:
//    MeshPointsIndex(Mesh* polyData);
//    ~MeshPointsIndex();
//
//    void SetPolyData(vtkPolyData* polyData);
//
//    vtkPointLocator* getMeshPointsIndex(){ return m_index; }
//    vtkIdList* getFoundPoints(){ return m_pointsList; }
//
//    vtkIdList* findPointsWithinRadius(const double* point, double radius);
//
//
//    int findNearestPoint(double* point);
//    int findNearestPointWithinRadius(double* point, double radius, double & dist);
//    int findNearestPointWithinRadius(float* point, float radius, double & dist);
//
//    double findNearestPoint(double* point, int &nearestPointId);
//    vtkIdList* findNearestNPoints(int n, const double* point);
//};

void getSHOTLocalRF(Mesh *cloud, std::vector<std::pair<int, double>> &NNpoints, double radius, int index, double *rfc);
void getSHOTLocalRF(Mesh *cloud, KDtree* pointsIndex, double radius, int index, double *rfc);
/*!
* \brief Round val to nearest integer.
*
* \param val double to round.
*/
//inline int round(double val){ return (int)floor(val + 0.5);};

/*!
* \brief Round val to nearest integer.
*
* \param val float to round.
*/
//inline int round(float val){ return (int)floor(val + 0.5);};

/*!
* \brief Compute minimum.
*
* \param a,b instances to compare.
*/
template <typename Type> inline Type Min(Type a, Type b) { return (a <= b)? a : b; };
/*!
* \brief Compute maximum.
*
* \param a,b instances to compare.
*/
template <typename Type> inline Type Max(Type a, Type b) { return (a >= b)? a : b; };



class  Random3DDetector
{
public:
    Random3DDetector(int nPoints, bool random = false, double radius = -1.0,  int minNeigh = 0, double borderDistance = -1.0); //computes "nPoints" randomly that have at least "minNeigh" neighboring points in a sphere of radius "radius" and that are not closer to any border than "borderDistance"

    virtual ~Random3DDetector();

    virtual int extract(Mesh* cloud, KDtree* index, Feature3D* & feat);

private:
    int m_minNeigh;
    double m_radius;
    int m_requestedNumFeat;

    Feature3D* m_feat;
    int m_numFeat;
    int m_featCapacity;
    double m_borderDistance;

    static void updateFeatCapacity(Feature3D * & feat, int & capacity, const int nFeat);
};



struct SHOTParams
{
	double radius;			///< radius of the sphere that defines the local neighborhood
	double localRFradius;	///< radius of the support to be used for the estimation of the local RF
    int shapeBins;				///< quantization bins for shape description
    int colorBins;				///< quantization bins for color description
    int minNeighbors;

    bool describeShape;		///< describe the shape channel
    bool describeColor;		///< describe the color channel

    int nThreads;			///< nThreads for parallel execution (0 for hardware parallelism)

    SHOTParams()	// default values
    {
        radius = 15;
        localRFradius = radius;
        shapeBins = 10;
        colorBins = 30;
        minNeighbors = 10;

        describeShape = true;
        describeColor = true;

        nThreads = 0;

    };

};


int findpointsInradius(double* query_point, Mesh* cloud, double radius, std::vector<std::pair<int, double>>& indices);

class SHOTDescriptor
{
public:
    SHOTDescriptor(SHOTParams params);
    virtual ~SHOTDescriptor();

    void setParams(SHOTParams params, int nPoints = 1000);

    static void RGB2CIELAB(unsigned char R, unsigned char G, unsigned char B, double &L, double &A, double &B2);

    virtual void describe(Mesh* cloud, KDtree* index, Feature3D* & feat, double** & desc, int nPoints);

    virtual int getDescriptorLength() { return m_descLength; };

    int getNumDescriptors() { return m_numDesc; };

private:
    double** m_desc;
    int m_numDesc;
    int m_descCapacity;
    int m_descLength;

    SHOTParams m_params;
    int m_k; //number of onion husks
    const int m_maxAngularSectors;
    const int m_nSect;

    //double m_sqradius;
    //double m_sqradius4;
    double m_radius3_4;
    double m_radius1_4;
    double m_radius1_2;

    static double sRGB_LUT[256];
    static double sXYZ_LUT[4000];

    void computeDesc (SHOTParams params, Mesh *cloud, KDtree* index, Feature3D *feat, double **desc, int startIndex, int endIndex);

    void interpolateDoubleChannel(Mesh* cloud, std::vector<std::pair<int, double>>& nnpoints, const std::vector<double> & distances, const double centralPoint[3], double rf[9], std::vector<double> & binDistanceShape, std::vector<double> & binDistanceColor, double* shot );

    void interpolateSingleChannel(Mesh* cloud,std::vector<std::pair<int, double>>& NNpoints, const std::vector<double> & distances, const double centralPoint[3], double rf[9], std::vector<double> & binDistance, int nr_bins, double* shot );

    void updateDescCapacity(double ** & desc, int & capacity, int nDesc, int descLength);
};

#endif
