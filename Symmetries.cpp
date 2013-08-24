/*
---------------------------------------------------------------------------------------------------

Searches the symmetry axises of a set of 2d points.


---------------------------------------------------------------------------------------------------
*/
#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <complex>
#include <boost/functional/hash.hpp> // TODO: why does removing this break everything?
using namespace std;


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


// Simple strutue which represents 2d line on infinite length. Contains two
// different points on the line.
struct Line {
	double_point p1;
	double_point p2;
	Line(double_point point1, double_point point2) {
		point1 = p1;
		point2 = p2;
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




// TODO: separate eigenvalue handling from eigenvectors (& decide if sensible)
// Returns two eigenvectors and corresponding eigenvalues
std::tuple<std::tuple<double, double>, std::tuple<double_point, double_point>>
					getEigenvectorsAndValues(vector<double_point> points) {
	double matrix[2][2] = { 0 };
	for (int i = 0; i != points.size(); i++) {
		matrix[0][0] += pow(points.at(i).get<0>(), 2);
		matrix[0][1] += points.at(i).get<0>()*points.at(i).get<1>();
		matrix[1][1] += pow(points.at(i).get<1>(), 2); 
	}
	matrix[1][0] = matrix[0][1];

	// Calculate eigenvalues, ie. solve (a-x)*(d-x)-b*c = 0 = ad - ax - dx + x^2 - bc
	// = x^2 - (a+d)x - bc + ad
	double a = 1;
	double b = -(matrix[0][0]+matrix[1][1]);
	double c = -matrix[0][1]*matrix[1][0] + matrix[0][0]*matrix[1][1];

	double x1, x2;
	double delta = b*b-4*a*c;
	if (delta<0) cout << "The solution is complex =(";
	else {
		x1=(-b+sqrt(delta))/(2*a);
		x2=(-b-sqrt(delta))/(2*a);
	}

	/*
	std::tuple<double,double> eigenvalues = getEigenvalues(points);
	double x1 = std::get<0>(eigenvalues);
	double x2 = std::get<1>(eigenvalues);
	*/

	double tempMatrix1[2][2];
	tempMatrix1[0][0] = matrix[0][0] - x1;
	tempMatrix1[0][1] = matrix[0][1];
	tempMatrix1[1][0] = matrix[1][0];
	tempMatrix1[1][1] = matrix[1][1] - x1;

	double tempMatrix2[2][2];
	tempMatrix2[0][0] = matrix[0][0] - x2;
	tempMatrix2[0][1] = matrix[0][1];
	tempMatrix2[1][0] = matrix[1][0];
	tempMatrix2[1][1] = matrix[1][1] - x2;

	// tempMatrix[0][0]*v11 + tempMatrix[0][1]*v12 = 0
	// tempMatrix[1][0]*v11 + tempMatrix[1][1]*v12 = 0
	// -> v11 = (tempMatrix[0][1]/tempMatrix[0][0])*v12
	// set v12 = 1
	
	double v11 = (tempMatrix1[0][1]/tempMatrix1[0][0]);
	double v12 = 1;
	double_point eigenvector1 = double_point(v11, v12);

	double v21 = (tempMatrix2[0][1]/tempMatrix2[0][0]);
	double v22 = 1;
	double_point eigenvector2 = double_point(v21, v22);

	if (abs(v11*v21 + v12*v22) > 0.0001) // TODO: get rid of magic constant
		cout << "oops, eigenvectors not orthogonal";

	//vector<double_point> eigenvectors(eigenvector1, eigenvector2);
	std::tuple<double_point, double_point> eigenvectors(eigenvector1, eigenvector2);
	std::tuple<double, double> eigenvalues(x1, x2);
	std::tuple<std::tuple<double, double>, std::tuple<double_point, double_point>> eigens(eigenvalues, eigenvectors);

	return eigens;
};



// Rotates points (vectors). Angle in radians.
vector<double_point> rotatePoints(vector<double_point> points, double angle) {

	double rotationMatrix[2][2] = {
		{cos(angle), -sin(angle)},
		{sin(angle), cos(angle)}
	};

	vector<double_point> rotatedPoints;
	for (int i = 0; i != points.size(); i++) {
		double x = points.at(i).get<0>() * rotationMatrix[0][0] + points.at(i).get<1>() * rotationMatrix[1][0];
		double y = points.at(i).get<0>() * rotationMatrix[0][1] + points.at(i).get<1>() * rotationMatrix[1][1];
		rotatedPoints.push_back(double_point(x, y));
	}

	return rotatedPoints;
};



// Assume line passes through center
bool isSymmetryaxis(Line line, vector<double_point> points, double_point center,
				double epsilon = std::numeric_limits<double>::epsilon()) {
	
	// TODO: center that shit?

	// Calculate angle of line.
	double xdiff = line.p1.get<0>() - line.p1.get<1>();
	double ydiff = line.p2.get<0>() - line.p2.get<1>();

	double angle;
	if (abs(xdiff) < std::numeric_limits<double>::epsilon()) {
		angle = (ydiff < 0 ? -1 : 1) * boost::math::constants::pi<double>()/2;
	} else {
		double angle = atan(ydiff/xdiff); // TODO: xdiff/ydiff? ydiff/xdiff? Is this shit correct?
	}

	// Rotate point so the given line is the x-axis.
	vector<double_point> rotatedPoints = rotatePoints(points, angle);

	// Now the vector formed by x-components if "points" should be orthogonal to
	// that of y-component only if "line" is a symmetry axis.
	double dotProduct = 0;
	for (int i = 0; i != rotatedPoints.size(); i++) {
		dotProduct += rotatedPoints.at(i).get<0>() * rotatedPoints.at(i).get<1>();
	}

	return dotProduct <= epsilon;
};




vector<Line> getSymmetryAxes(vector<double_point> points) {

	// Center given points.
	double_point center = getCenterPoint(points);
	for (int i = 0; i != points.size(); i++) {
		double x = points.at(i).get<0>() - center.get<0>();
		double y = points.at(i).get<1>() - center.get<1>();
		points.at(i).set<0>(x);
		points.at(i).set<1>(y);
	}

	// Get eigenvalues ang eigenvectors.
	std::tuple<std::tuple<double, double>, std::tuple<double_point, double_point>> eigens;
	eigens = getEigenvectorsAndValues(points);
	
	std::tuple<double, double> eigenvalues = std::get<0>(eigens);
	std::tuple<double_point, double_point> eigenvectors = std::get<1>(eigens);

	double eig1 = std::get<0>(eigenvalues);
	double eig2 = std::get<1>(eigenvalues);

	Line vector1 = Line(double_point(0,0), std::get<0>(eigenvectors));
	Line vector2 = Line(double_point(0,0), std::get<1>(eigenvectors));

	// Test if eigenvectors are symmetry axes
	bool isVector1Symmetric = isSymmetryaxis(vector1, points, double_point(0,0));
	bool isVector2Symmetric = isSymmetryaxis(vector2, points, double_point(0,0));


	vector<Line> symmetryAxes;

	if (isVector1Symmetric != isVector2Symmetric) {
		if (isVector1Symmetric) {
			// TODO: add center
			symmetryAxes.push_back(vector1);
		} else {
			// TODO: add center
			symmetryAxes.push_back(vector2);
		}
	}

	if (abs(eig1 - eig2) < std::numeric_limits<double>::epsilon()) {
		// TODO: handle >= 2 symmetry axes
	}
	return symmetryAxes;
}


bool isSymmetryaxis(Line line, vector<double_point> points,
				double epsilon = std::numeric_limits<double>::epsilon()) {
	double_point center = getCenterPoint(points);
	return isSymmetryaxis(line, points, center, epsilon);
};





int main() {
	std::cout << "QWEQWEQWE";
	double_point point1(1,1);
	double_point point2(-1,1);
	double_point point3(-1,-1);
	double_point point4(1,-1);

	vector<double_point> asd; //(point1, point2, point3, point4);
	asd.push_back(point1);
	asd.push_back(point2);
	asd.push_back(point3);
	asd.push_back(point4);

	vector<Line> axes = getSymmetryAxes(asd);

	return 0;
};