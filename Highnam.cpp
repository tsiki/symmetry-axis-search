/*
Implementation of the symmetry axis finding algorithm presented by Highnam:
http://www.ri.cmu.edu/pub_files/pub3/highnam_p_t_1985_1/highnam_p_t_1985_1.pdf
*/

#define _USE_MATH_DEFINES

#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <complex>
#include <unordered_set>
#include <algorithm>
#include <regex>

#include <boost/geometry.hpp>
namespace bg = boost::geometry;

const double ANGLE_TH = 0.000001;
const double DISTANCE_TH = 0.0000001;
const double COMPARISON_TH = 0.0000001;
const double PI = 3.141592653589793238463; // TODO: standardize PI usage

// TODO: remove redundancy
typedef bg::model::point
<
  double, 2, bg::cs::cartesian
> double_point;

bool polarPointComparator (std::complex<double> i, std::complex<double> j) {
	double angleDiff = std::arg(i) - std::arg(j);
	if (angleDiff < 0) {
		return true;
	} else if (angleDiff > 0) {
		return false;
	}
	return std::abs(i) < std::abs(j);
};


// Returns a vector which contains, in alternate steps, angles and vectors of distances.
std::vector<std::vector<double>> transformList(std::vector<std::complex<double>> points) {
	std::vector<std::vector<double>> retval; // TODO: preallocate size of 2 x points
	// TODO: check only one point special case
	double prevAngle;
	for (int i = 0; i != points.size(); i++) {
		double angle = getAngle(points[i]);
		if (i > 0) {
			double angleDiff = prevAngle - angle;
			if (angleDiff < ANGLE_TH) {
				retval.back().push_back(getDistance(points[i]));
			} else {
				std::sort(retval.back().begin(), retval.back().end());
				std::vector<double> newAngle;
				newAngle.push_back(angle);
				retval.push_back(newAngle);
				std::vector<double> newDistance;
				newDistance.push_back(getDistance(points[i]));
				retval.push_back(newDistance);
			}
		} else {
			std::vector<double> newDistance;
			newDistance.push_back(getDistance(points[i]));
			retval.push_back(newDistance);
		}
	}
	double firstAndLastPointDiff = 2*PI - getAngle(points.back()) + getAngle(points.front());
	std::vector<double> newAngle;
	newAngle.push_back(firstAndLastPointDiff);
	retval.push_back(newAngle);
	return retval;
};



bool areEqual(std::vector<double> a, std::vector<double> b) {
	if (a.size() != b.size()) {
		return false;
	}
	for (int i = 0; i != a.size(); i++) {
		if (abs(a[i] - b[i]) > COMPARISON_TH) {
			return false;
		}
	}
	return true;
}

// Uses Manacher’s algorithm to find the lengths of the palindroms in O(n) 
// http://leetcode.com/2011/11/longest-palindromic-substring-part-ii.html
std::vector<int> getPalindromeLengths(std::vector<std::vector<double>> anglesAndDistances) {
	std::vector<int> lp;
	lp.reserve(anglesAndDistances.size());
	int C = 0;
	int R = 0;
	for (int i = 0; i != anglesAndDistances.size(); i++) {
		int i_mirror = 2*C - i;

		lp[i] = (R > i) ? std::min(R - i, lp[i_mirror]) : 0;

		while (areEqual(anglesAndDistances[i + 1 + lp[i]], anglesAndDistances[i - 1 - lp[i]])) {
			lp[i]++;
		}

		if (i + lp[i] > R) {
			C = i;
			R = i + lp[i];
		}
	}

	return lp;
}



// Assumes the set of points is centered.
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

	std::vector<std::vector<double>> anglesAndDistances = transformList(polarPoints);
	int n = anglesAndDistances.size();
	for (int i = 0; i != n; i++) {
		anglesAndDistances.push_back(anglesAndDistances[i]);
	}

	std::vector<int> palindromeLengths = getPalindromeLengths(anglesAndDistances);

	// TODO: now when "palindromeLengths" is over the threshold value, we have a symmetry axis
	
	return (std::vector<double_point>)NULL;

};





double getAngle(std::complex<double> point) {
	return std::arg(point);
}

double getDistance(std::complex<double> point) {
	return std::abs(point);
}

