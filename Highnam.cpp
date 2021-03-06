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
#include <iostream>

namespace symmetryAxisSearch {

	// TODO: fix comments
	// Points with radian value under this are considered to have the same angle as measured from the centroid.
	constexpr double ANGLE_TH = 1e-7;
	// Points whose distance from the centroid is plus minus this are considered to be at the same distance.
	constexpr double DISTANCE_TH = 1e-7;
	constexpr double PI = 3.141592653589793238463;

	using namespace std;

	struct AngleAndDistances {
		double angle;
		vector<double> distances;

		AngleAndDistances(double angle, double distance) {
			this->angle = angle;
			this->distances.push_back(distance);
		}

		AngleAndDistances(double angle) {
			this->angle = angle;
		}

		AngleAndDistances() {}

		bool operator==(const AngleAndDistances & rhs) {
			if (abs(this->angle - rhs.angle) >= ANGLE_TH || this->distances.size() != rhs.distances.size()) {
				return false;
			}
			for (int i = 0; i < this->distances.size(); i++) {
				if (abs(this->distances[i] - rhs.distances[i]) >= DISTANCE_TH) {
					return false;
				}
			}
			return true;
		}
	};


	template <class T> std::pair<double,double> getCenterPoint(std::vector<pair<T,T>> points) {
		double xTotal = 0, yTotal = 0, total = 0;
		for (auto p : points) {
			xTotal += p.first;
			yTotal += p.second;
		}
		auto s = points.size();
		return {xTotal/s, yTotal/s};
	}

	template <class T> std::vector<std::pair<double,double>> centerPoints(std::vector<pair<T,T>> points, std::pair<double,double> center) {
		std::vector<std::pair<double,double>> centeredPoints;
		for (int i = 0; i != points.size(); i++) {
			centeredPoints.push_back({points[i].first - center.first, points[i].second - center.second});
		}
		return centeredPoints;
	};

	std::vector<AngleAndDistances> transformList(std::vector<std::pair<double, double>> points) {
		std::vector<AngleAndDistances> ans;
		for (int i = 0; i != points.size(); i++) {
			double angle, distance;
			tie(angle, distance) = points[i];
			if (ans.size() > 0) {
				double prevAngle = points[i-1].first;
				double angleDiff = angle - prevAngle;
				if (angleDiff < ANGLE_TH) {
					ans.back().distances.push_back(distance);
					continue;
				}
				ans.push_back(AngleAndDistances(angleDiff));
			}
			ans.push_back(AngleAndDistances(0, distance));
		}
		double firstLastDiff = points.front().first + (2*PI - points.back().first);
		ans.push_back(AngleAndDistances(firstLastDiff));
		return ans;
	};


	std::vector<int> getPalindromeLengths(std::vector<AngleAndDistances> points) {
		int N = points.size();
		std::vector<int> L(N); // LPS Length vector
		if (N == 0)
		    return L;
		int C = 1; // centerPosition
		int R = 2; // centerRightPosition
		int i = 0; // currentRightPosition
		int iMirror; // currentLeftPosition
		int maxLPSLength = 0;
		int maxLPSCenterPosition = 0;
		int start = -1;
		int end = -1;
		int diff = -1;
		for (i = 1; i < N; i++) {
		    //get currentLeftPosition iMirror for currentRightPosition i
		    iMirror  = 2*C-i;
		    L[i] = 0;
		    diff = R - i;
		    //If currentRightPosition i is within centerRightPosition R
		    if(diff > 0)
		        L[i] = min(L[iMirror], diff);
	 
		    //Attempt to expand palindrome centered at currentRightPosition i
		    //Here for odd positions, we compare characters and 
		    //if match then increment LPS Length by ONE
		    //If even position, we just increment LPS by ONE without 
		    //any character comparison
			while (i + L[i] < N && i - L[i] > 0 && points[(i + L[i] + 1)] == points[(i - L[i] - 1)]) {
		        L[i]++;
			}

		    //If palindrome centered at currentRightPosition i 
		    //expand beyond centerRightPosition R,
		    //adjust centerPosition C based on expanded palindrome.
		    if (i + L[i] > R) {
		        C = i;
		        R = i + L[i];
		    }
		}
		return L;
	}


	std::vector<std::pair<double, double>> getPolarPoints(std::vector<std::pair<double,double>> & points) {
		std::vector<std::pair<double, double>> polarPoints;
		for (auto & p : points) {
			double x,y;
			tie(x,y) = p;
			double distance = sqrt(pow(x,2) + pow(y,2));
			if (distance > DISTANCE_TH) {
				double angle = atan2(y,x);
				if (angle < 0) {
					angle += 2*PI;
				}
				polarPoints.push_back({angle, distance});
			}
		}
		return polarPoints;
	}

	std::vector<std::pair<double, double>> getSymmetryAxes(std::vector<AngleAndDistances> & ad,
			std::vector<int> & palindromeLengths, double startAngle) {
		int pfloor = ad.size()/4 - 1;
		double angleAcc = startAngle;
		vector<pair<double, double>> ans;
		int start = palindromeLengths.size()/4;
		int end = palindromeLengths.size()/2;

		for (int i = 0; i < start; i++) angleAcc += ad[i].angle;
		for (int i = start; i < end; i++) {
//cout << ad[i].angle << " zzz " << ad[i].distances.size() << endl;
//if (ad[i].distances.size()) cout << ad[i].distances[0] << endl;
			if (palindromeLengths[i] >= pfloor) {
				double angle = angleAcc + ad[i].angle/2;
//cout << angleAcc << " " << ad[i].angle << " " << angle << " QQQ " << palindromeLengths[i] << " " << startAngle << endl;
				ans.push_back({cos(angle), sin(angle)});
			}
			angleAcc += ad[i].angle;
		}
		return ans;
	}

	// Returns the symmetry axes of a vector of 2D points.
	// The return value is a vector of points where the first point is the centroid of the given
	// set of points, and the other points are directions of the symmetry axes from that centroid
	// (as all symmetry axes pass through the centroid of the group, it can be used for each).
	template <class T> vector<pair<double, double>> getSymmetryAxesHighnam(std::vector<pair<T, T>> points) {
		auto center = getCenterPoint(points);
		auto centeredPoints = centerPoints(points, center);
		std::vector<std::pair<double, double>> polarPoints = getPolarPoints(centeredPoints);
		std::sort(polarPoints.begin(), polarPoints.end());
		double startAngle = polarPoints.front().first;

		std::vector<AngleAndDistances> anglesAndDistances = transformList(polarPoints); // TODO check that contains at least 2
		if (anglesAndDistances.size() < 2) {
			return vector<pair<double, double>> {center};
		}

		// append vector to itself for calculating the palindrome lengths
		auto old_count = anglesAndDistances.size();
		anglesAndDistances.resize(2 * old_count);
		std::copy_n(anglesAndDistances.begin(), old_count, anglesAndDistances.begin() + old_count);

		std::vector<int> palindromeLengths = getPalindromeLengths(anglesAndDistances);
//for (auto i : palindromeLengths) cout << i << " "; cout << endl;
		auto axes = getSymmetryAxes(anglesAndDistances, palindromeLengths, startAngle);

		for (auto & axis : axes) {
			axis.first += center.first;
			axis.second += center.second;
		}

		axes.insert(axes.begin(), center);

		return axes;
	};

}




int main() {
	using namespace std;
	vector<pair<double, double>> pts {{0,0}, {1,1}, {0,1}, {1,0}};
	//vector<pair<double, double>> pts {{0,0}, {1,1}, {0,1}, {2,0}};
	//vector<pair<double, double>> pts {{0,0}, {2,1}, {0,1}, {2,0}};
	//vector<pair<double, double>> pts {{0,0}, {0,1}, {1,0}};
	auto v = symmetryAxisSearch::getSymmetryAxesHighnam(pts);
	for (auto p : v) cout << p.first << "x" << p.second << " "; cout << endl;

	return 0;
};





