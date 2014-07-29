/*
Implementation of the symmetry axis finding algorithm presented by Highnam:
http://www.ri.cmu.edu/pub_files/pub3/highnam_p_t_1985_1/highnam_p_t_1985_1.pdf
*/


#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <complex>

#include <boost/geometry.hpp>
namespace bg = boost::geometry;

const double ANGLE_TH = 0.000001;
const double DISTANCE_TH = 0.0000001;

// TODO: remove redundancy
typedef bg::model::point
<
  double, 2, bg::cs::cartesian
> double_point;

bool polarPointComparator (std::complex<double> i, std::complex<double> j) {
	double angleDiff = std::arg(i) - std::arg(j);
	if (angleDiff < -ANGLE_TH) {
		return true;
	} else if (angleDiff > ANGLE_TH) {
		return false;
	}
	return std::abs(i) < std::abs(j);
};

// Assumes points are centered.
std::vector<double_point> getSymmetryAxesHighnam(std::vector<double_point> points) {

	std::vector<std::complex<double>> polarPoints;
	for (auto point = points.begin(); point != points.end(); ++point) {
		double x = (*point).get<0>();
		double y = (*point).get<1>();
		double distance = sqrt(pow(x,2) + pow(y,2));
		if (distance > DISTANCE_TH) {
			polarPoints.push_back(std::polar(distance, atan(y/x)));
		}
	}
	// TODO: check that the angles are fine

	std::sort(polarPoints.begin(), polarPoints.end(), polarPointComparator);

};
