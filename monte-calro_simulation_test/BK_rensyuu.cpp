//practice practice practice 
#include <boost/math/special_functions/bessel.hpp>
#include <random>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <iomanip>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include<complex>
//std::bind
#include <functional>
#include "Parameters.h"
#include "Spline.h"
#include "INTDE2_functional.h"
#include <stdlib.h>
#include <stdio.h>  

using namespace std;

//parameter set
	const double alpha_s_bar = 0.2;
	const double anomalus_dimA1 = 11.0 / 12.0;
	const double Color_number = 3;
	const double Flavor_number = 3;
	const double Constant_Ca = 2.5;
	const double Constant_Q0 = 0.4;
	const double Constant_Rp = 0.7;
	const double Constant_p = 2.3;
	const double Lambda_Nf = 1.0;

double Integration_Kernel_DLA(double rho) {
	double result = _j1(2.0*sqrt(alpha_s_bar*rho*rho)) / sqrt(alpha_s_bar*rho*rho);
	if (std::abs(rho) < 1.0e-5) {
		result = 1.0 - alpha_s_bar*rho*rho / 2.0 + alpha_s_bar*rho*rho*alpha_s_bar*rho*rho / 12.0;
	}
	return result;
}

double Integration_Kernel_pm(const double x_y_square, const double x_z_square, const double y_z_square) {
	double result = 0;
	double min_x_z_y_z = (std::min)(abs(x_z_square), abs(y_z_square));
	double index_power = alpha_s_bar * anomalus_dimA1;
	if (abs(x_y_square) < min_x_z_y_z) {
		result = pow(x_y_square*x_y_square / min_x_z_y_z / min_x_z_y_z, index_power);
	}
	else {
		result = pow(x_y_square*x_y_square / min_x_z_y_z / min_x_z_y_z, - index_power);
	}
	return result;
}

double Integration_Kernel_BK(const double x_y_square, const double x_z_square, const double y_z_square) {
	return x_y_square*x_y_square / x_z_square / x_z_square / y_z_square / y_z_square;
}

double Function_L(const double x_y_square, const double x_z_square) {
	return log(x_z_square*x_z_square / x_y_square / x_y_square);
}

double alpha_bar_running(const double relative_distance) {
	double b_Nf = (11.0 * Color_number - 2.0*Flavor_number) / 12.0 / Pi;
	return Color_number / (b_Nf*log(4.0*Constant_Ca*Constant_Ca / (relative_distance*relative_distance*Lambda_Nf*Lambda_Nf)))/Pi;
}

double running_alpha_fastest(const double x_y_square, const double x_z_square, const double y_z_square) {
	double factor = (x_z_square*x_z_square - y_z_square*y_z_square) / (x_y_square * x_y_square);
	return 1.0/(1.0/alpha_bar_running(x_y_square) + factor*(alpha_bar_running(x_z_square) alpha_bar_running() ))
}