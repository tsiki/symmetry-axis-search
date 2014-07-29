/*
---------------------------------------------------------------------------------------------------

Searches the symmetry axises of a set of 2d points.


---------------------------------------------------------------------------------------------------
*/
#define _USE_MATH_DEFINES

#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <complex>
#include <unordered_set>

#include <boost/functional/hash.hpp>
using namespace std;

// Thresholds
const double DISTANCE_FROM_AXIS_TH = 0.0000001;
const double FUZZINESS_TH = 0.0000001;
const double EIGENVALUE_TH = 0.0000001;
const double COVARIANCE_TH = 0.0000001;


#include <boost/geometry.hpp>
namespace bg = boost::geometry;

typedef bg::model::point
<
  int, 2, bg::cs::cartesian
> integer_point;

typedef bg::model::point
<
  double, 2, bg::cs::cartesian
> double_point;


struct Hash {
   size_t operator() (const std::pair<double, double> &doublePair) const {
		std::size_t seed = 0;
		boost::hash_combine(seed, std::get<0>(doublePair));
		boost::hash_combine(seed, std::get<1>(doublePair));
		return seed;
  }
};




double_point getCenterPoint(vector<double_point> pointVector) {
	  double xTotal = 0, yTotal = 0, total = 0;
		vector<double_point>::iterator iter;
    for (iter = pointVector.begin(); iter != pointVector.end(); iter++) {
        xTotal += (*iter).get<0>();
        yTotal += (*iter).get<1>();
        total++;
    }
		double_point center;
    center.set<0>(xTotal/total);
    center.set<1>(yTotal/total);
		return center;
}



double* getCovarianceMatrix(vector<double_point> points) {
	double *matrix = new double [4] (); // TODO: delete // TODO: 2d vector? 2d array?

	for (int i = 0; i != points.size(); i++) {
		matrix[0] += pow(points.at(i).get<0>(), 2);
		matrix[1] += points.at(i).get<0>()*points.at(i).get<1>();
		matrix[3] += pow(points.at(i).get<1>(), 2); 
	}
	matrix[2] = matrix[1];
	cout << matrix[0] << matrix[3];
	return matrix;
};


std::pair<double, double> getEigenValues(vector<double_point> points, double *covMatrix) {
	// Calculate eigenvalues, ie. solve (a-x)*(d-x)-b*c = 0 = ad - ax - dx + x^2 - bc
	// = x^2 - (a+d)x - bc + ad
	// TODO: ambiguous a,b,c
	double a = 1;
	double b = -(covMatrix[0]+covMatrix[3]);
	double c = -covMatrix[1]*covMatrix[2] + covMatrix[0]*covMatrix[3];
	double x1, x2;
	double delta = b*b-4*a*c;
	x1=(-b+sqrt(delta))/(2*a);
	x2=(-b-sqrt(delta))/(2*a);
	
	if (x1 > x2) {
		return std::make_pair(x1, x2);
	}
	return std::make_pair(x2, x1);
};

// Return the eigenvectors for centered points.
std::pair<double_point, double_point> getEigenVectors(double *covMatrix, std::pair<double, double> eigenValues) {

	double eig1 = eigenValues.first;
	double eig2 = eigenValues.second;

	//cout << "cov matrix: " << covMatrix[0] << " " << covMatrix[1] << " " << covMatrix[2] << " " << covMatrix[3] << "\n"; 

	double tempMatrix1[2][2];
	tempMatrix1[0][0] = covMatrix[0] - eig1;
	tempMatrix1[0][1] = covMatrix[1];
	tempMatrix1[1][0] = covMatrix[2];
	tempMatrix1[1][1] = covMatrix[3] - eig1;

	double tempMatrix2[2][2];
	tempMatrix2[0][0] = covMatrix[0] - eig2;
	tempMatrix2[0][1] = covMatrix[1];
	tempMatrix2[1][0] = covMatrix[2];
	tempMatrix2[1][1] = covMatrix[3] - eig2;

	// We have:
	// tempMatrix[0][0]*v11 + tempMatrix[0][1]*v12 = 0
	// tempMatrix[1][0]*v11 + tempMatrix[1][1]*v12 = 0
	// -> v11 = (tempMatrix[0][1]/tempMatrix[0][0])*v12
	// -> (tempMatrix[1][0]*tempMatrix[0][1]/tempMatrix[0][0] + tempMatrix[1][1])*v12 = 0
	// set v12 = 1

	double_point eigenVector1;
	double_point eigenVector2;

	// http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
	if (abs(covMatrix[2]) > COVARIANCE_TH) {
		eigenVector1 = double_point(eig1 - covMatrix[3], covMatrix[2]);
		eigenVector2 = double_point(eig2 - covMatrix[3], covMatrix[2]);
	} else if (abs(covMatrix[1]) > COVARIANCE_TH) {
		eigenVector1 = double_point(covMatrix[2], eig1 - covMatrix[0]);
		eigenVector2 = double_point(covMatrix[2], eig2 - covMatrix[0]);
	} else {
		eigenVector1 = double_point(1,0);
		eigenVector2 = double_point(0,1);
	}

	// TODO: check for zero division and zero eigenvalues and zero length eigenvectors

	return std::make_pair(eigenVector1, eigenVector2);
};


bool isSymmetryAxis(double_point axisDirection, vector<double_point> centeredPoints,
				double epsilon = std::numeric_limits<double>::epsilon()) {
	
	double eigx = axisDirection.get<0>();
	double eigy = axisDirection.get<1>();

	cout << "tryin with " << eigx << " " << eigy << endl;

	// Create vector orthogonal to axisDirection // TODO: make sure both dimensions are positive
	double x = cos(M_PI/2)*eigx - sin(M_PI/2)*eigy;
	double y = sin(M_PI/2)*eigx + cos(M_PI/2)*eigy;
	
	std::unordered_set<std::pair<double, double>, Hash> set;
	std::vector<std::pair<double, double>> positiveDotProductPoints;
	for (auto point = centeredPoints.begin(); point != centeredPoints.end(); ++point) {
		double dotProduct1 = (*point).get<0>() * x + (*point).get<1>() * y;
		double dotProduct2 = (*point).get<0>() * eigx + (*point).get<1>() * eigy;
		if (dotProduct1 < -1*DISTANCE_FROM_AXIS_TH) {

			cout << "adding " << dotProduct1 << " and " << dotProduct2 << endl;

			// Due to floating point inaccuracies, insert both the ceiling and floor
			double floorVal1 = floor(-1*dotProduct1/FUZZINESS_TH)*FUZZINESS_TH;
			double ceilVal1 = ceil(-1*dotProduct1/FUZZINESS_TH)*FUZZINESS_TH;
			double ceilVal2 = ceil(dotProduct2/FUZZINESS_TH)*FUZZINESS_TH;
			double floorVal2 = floor(dotProduct2/FUZZINESS_TH)*FUZZINESS_TH;
			set.insert(std::make_pair(floorVal1, floorVal2));
			set.insert(std::make_pair(floorVal1, ceilVal2));
			set.insert(std::make_pair(ceilVal1, floorVal2));
			set.insert(std::make_pair(ceilVal1, ceilVal2));
		} else if (dotProduct1 > DISTANCE_FROM_AXIS_TH) {
			positiveDotProductPoints.push_back(std::make_pair(dotProduct1, dotProduct2));
		}
	}

	cout << "positives: " << positiveDotProductPoints.size() << endl;
	cout << "set size: " << set.size() << endl;
	
	for (auto pair = positiveDotProductPoints.begin(); pair != positiveDotProductPoints.end(); ++pair) {
		double dotProduct1 = (*pair).first;
		double dotProduct2 = (*pair).second;

		double floorVal1 = floor(dotProduct1/FUZZINESS_TH)*FUZZINESS_TH;
		double ceilVal1 = ceil(dotProduct1/FUZZINESS_TH)*FUZZINESS_TH;
		double ceilVal2 = ceil(dotProduct2/FUZZINESS_TH)*FUZZINESS_TH;
		double floorVal2 = floor(dotProduct2/FUZZINESS_TH)*FUZZINESS_TH;

		cout << "testing " << dotProduct1 << " and " << dotProduct2 << endl;

		if (!set.count(std::make_pair(floorVal1, floorVal2)) && !set.count(std::make_pair(floorVal1, ceilVal2)) &&
				!set.count(std::make_pair(ceilVal1, floorVal2)) && !set.count(std::make_pair(ceilVal1, ceilVal2))) {
			cout << "not there!" << endl;
			return false;
		}
	}
	
	return true;
};


vector<double_point> centerPoints(vector<double_point> points) {
	std::vector<double_point> centeredPoints;
	double_point center = getCenterPoint(points);
	for (int i = 0; i != points.size(); i++) {
		double x = points.at(i).get<0>() - center.get<0>();
		double y = points.at(i).get<1>() - center.get<1>();
		double_point cp(x,y);
		centeredPoints.push_back(cp);
	}
	return centeredPoints;
};


vector<double_point> getSymmetryAxes(vector<double_point> points) {

	// Center given points.
	std::vector<double_point> centeredPoints = centerPoints(points);

	// Get eigenvalues ang eigenvectors.
	double *covMatrix = getCovarianceMatrix(centeredPoints);
	std::pair<double, double> eigenValues = getEigenValues(centeredPoints, covMatrix);
	std::pair<double_point, double_point> eigenVectors = getEigenVectors(covMatrix, eigenValues);

	cout << "eigen1: " << eigenVectors.first.get<0>() << " " << eigenVectors.first.get<1>() << "\n";
	cout << "eigen2: " << eigenVectors.second.get<0>() << " " << eigenVectors.second.get<1>() << "\n";

	// TODO: zero length checks etc. be here
	// TODO: if eigenvalues are the ~same, default to nlogn search
	// TODO: both eigenvectors need to be normalized to length of 1

	vector<double_point> returnVector;
	if (abs(eigenValues.first - eigenValues.second) < EIGENVALUE_TH) {
		std::cout << "Falling back to nlogn search... someday" << endl;
		// TODO: fall back to nlogn check here
	} else {
		if (isSymmetryAxis(eigenVectors.first, centeredPoints)) {
			returnVector.push_back(eigenVectors.first);
		}
		if (isSymmetryAxis(eigenVectors.second, centeredPoints)) {
			returnVector.push_back(eigenVectors.second);
		}
	}
	
	return returnVector;
}






int main() {

	/*
	double_point point1(2,1);
	double_point point2(-1,1);
	double_point point3(-1,-1);
	double_point point4(2,-1);
	double_point point5(5,0);

	vector<double_point> asd;
	asd.push_back(point1);
	asd.push_back(point2);
	asd.push_back(point3);
	asd.push_back(point4);
	asd.push_back(point5);
	*/

	double_point point1(1,0);
	double_point point2(0,1);
	double_point point3(4,2);
	double_point point4(2,4);

	vector<double_point> asd;
	asd.push_back(point1);
	asd.push_back(point2);
	asd.push_back(point3);
	asd.push_back(point4);


	vector<double_point> axes = getSymmetryAxes(asd);

	cout << "\nthe number of axes is " << axes.size();


	return 0;
};