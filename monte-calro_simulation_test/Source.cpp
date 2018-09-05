﻿//From https://qiita.com/Ushio/items/0040b3c74a480c46c80c


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

using namespace std;

#define K_0_4_K_2_4_FUNC

//#define MAKE_FUNCTION

#define POLAR_COORDINATE

#define MULTI_PARTON

int nfunc;

/**
* The evolution step size.
*/
#define DELTA_T         0.01
#define OUTPUT_DELTA_T  0.2
#define END_T           1.0


//number of ansamble?
const int N = 1100000;

const int NPHI = 40;

const double h_phi_intmeasure = 2.0*Pi / NPHI;

/**
* coupling constant
*/
const double alpha_s = 0.2;

/**
* The size of nuclear included in the solution of the BK equation
*/
const double Radius_Nuc = 1.0;

/**
* The parameter in Gauss function
*/
const double Gauss_param = 2.0;

/**
* The parameter in Gauss function
*/
const double Gauss_param_B = 2.0;

//Upper limoit momentum of the produced particles
const double Upper_momk = 5.0;


/**
* The parameter \Lamda_QCD
*/
const double Lamda_QCD = 1.0;

/**
* The parameter \Lamda_QCD for confinement scale of the  nucleus
*/
const double Lamda_QCD_nucleus = 1.0/6.0;

//d^2(x',y')=d^2(x,-x) -> x
vector<double> vector_x_val;

//Solution of BK equation
vector<double> sol_BK;

//The spline interpolation of the solution of BK equation using  Tino Kluge code.
tk::spline sol_BK_sp;

vector<double> c_2_2(200, 0);
vector<double> error_c_2_2(200, 0);
vector<double> s_2_4(200, 0);
vector<double> error_s_2_4(200, 0);
vector<double> c_2_4(200, 0);
vector<double> error_c_2_4(200, 0);


const int number_impact = 1000;

const double stepsize_impact = 1.0 / number_impact;

vector<double> K_0_4(number_impact - 2, 0);
vector<double> error_K_0_4(number_impact - 2, 0);
vector<double> K_2_4(number_impact - 2, 0);
vector<double> error_K_2_4(number_impact - 1, 0);
vector<double> K_vector(number_impact - 2, 0);

tk::spline K_0_4_sp, K_2_4_sp;


//caclulate d^2(x,y)=d^2(r,b,co_theta(b-r))
double dist2_rad_imp(const double relat, const double impact_param, const double cos_theta){
	double ret_d2 = Radius_Nuc*Radius_Nuc*relat*relat
		/ ((Radius_Nuc*Radius_Nuc + impact_param*impact_param + relat*relat / 4.0)
			*(Radius_Nuc*Radius_Nuc + impact_param*impact_param + relat*relat / 4.0)
			- impact_param*impact_param*relat*relat*cos_theta*cos_theta);
	double ret_d2_v2 = Radius_Nuc*Radius_Nuc*relat*relat
		/ ((Radius_Nuc*Radius_Nuc + impact_param*impact_param + relat*relat / 4.0 - impact_param*relat*cos_theta)
			*(Radius_Nuc*Radius_Nuc + impact_param*impact_param + relat*relat / 4.0 + impact_param*relat*cos_theta)
			);
	if (abs(ret_d2_v2) > 1.0) {
		cerr << "d2>1 !!!??? \n";
		ret_d2_v2 = 1.0; 
	}
	return ret_d2_v2;
}

//caclulate d^2(x,y)=d^2(bx,by,rx,ry)
double dist2_rad_imp_ncoor(const double b_x, const double b_y, const double r_x, const double r_y){


	return Radius_Nuc*Radius_Nuc*(r_x*r_x+r_y*r_y)
		/ ((Radius_Nuc*Radius_Nuc + (b_x*b_x + b_y*b_y) + (r_x*r_x + r_y*r_y) / 4.0)
		*(Radius_Nuc*Radius_Nuc + (b_x*b_x + b_y*b_y) + (r_x*r_x + r_y*r_y) / 4.0)
		- (b_x*r_x + b_y*r_y)*(b_x*r_x + b_y*r_y));
}


//calculate d^2(x,-x) to x. (d^2 = 0 -> x=0)
double gain_x_val(double dist_2){
	if (dist_2 == 0){
		return 0;
	}
	else if (abs(dist_2)>1.0) {
		cerr << "The distance d^2 is over the 1.0" << "\n";
		return Radius_Nuc;
	}
	else{
		return Radius_Nuc*sqrt(-1.0 + 2.0 / dist_2 - 2.0 / dist_2*sqrt(1.0 - dist_2));
	}
}

void initial(double tau){


	ostringstream ifilename;
	ifilename << "res_501_1001_0.001_R_1_t_" << tau << "_hipre.txt";
	//imput and output file
	std::ifstream ifs(ifilename.str().c_str());


	double in_vector_x_val, in_sol_BK;
	int i = 0, max;

	char str[256];
	if (ifs.fail())
	{
		std::cerr << "失敗" << std::endl;
	}
	ifs.getline(str, 256 - 1);
	//      std::cout << str << std::endl;

	vector_x_val.erase(vector_x_val.begin(), vector_x_val.end());
	sol_BK.erase(sol_BK.begin(), sol_BK.end());


	while ((ifs >> in_vector_x_val >> in_sol_BK)){
		vector_x_val.push_back(in_vector_x_val);
		sol_BK.push_back(in_sol_BK);
		i++;
	}
	max = i;
	double size_vec = vector_x_val.size();


	sol_BK_sp.set_points(vector_x_val, sol_BK);

}

inline double f(double x) {
	return std::exp(-x*x / (0.5));
}

inline double p(double x) {
	return 1.0 - 0.5 * x;
}
inline double inverse_of_cdf(double x) {
	return 2.0 - 2.0 * std::sqrt(1.0 - x);
}

//Bessel
double bessel_ink_0_4(const double x_val)
{
	double bessel_in=0;
		bessel_in = _j0(x_val);

	return bessel_in;

}


//Bessel
double bessel_ink_2_4(const double x_val)
{
	double bessel_in_2 = 0;
	
	if(0<=x_val && x_val <1e-3){
		bessel_in_2 = x_val*x_val*x_val / 32.0 - x_val*x_val*x_val*x_val*x_val / 576.0;
	}
	else{

		bessel_in_2 = (2.0 - 2.0*_j0(x_val) - x_val*_j1(x_val)) / x_val;
	}

	return bessel_in_2;

}

//multi particle distribution : original
double multi_particle_dist(const double b_1, const double b_2, const double b_3, const double b_4,
	const double phi_1, const double phi_2, const double phi_3, const double phi_4)
{
	double f_A = exp(-(b_1*b_1 + b_2*b_2 + b_3*b_3 + b_4*b_4) / Gauss_param_B / Gauss_param_B
		- Lamda_QCD*Lamda_QCD*(3.0*(b_1*b_1 + b_2*b_2 + b_3*b_3 + b_4*b_4) -
		2.0*(cos(phi_1 - phi_2)*b_1*b_2 + cos(phi_1 - phi_3)*b_1*b_3 + cos(phi_1 - phi_4)*b_1*b_4 + cos(phi_3 - phi_2)*b_2*b_3 + cos(phi_4 - phi_2)*b_2*b_4 + cos(phi_3 - phi_4)*b_3*b_4
		)));
				
		return f_A;

}

//multi particle distribution : original
double multi_particle_dist_ov_normal_dist(const double b_1, const double b_2, const double b_3, const double b_4,
	const double phi_1, const double phi_2, const double phi_3, const double phi_4)
{
	double f_A = Gauss_param_B*sqrt(Pi)* b_1*b_2*b_3*b_4*
		exp(
		- Lamda_QCD*Lamda_QCD*(3.0*(b_1*b_1 + b_2*b_2 + b_3*b_3 + b_4*b_4) -
		2.0*(cos(phi_1 - phi_2)*b_1*b_2 + cos(phi_1 - phi_3)*b_1*b_3 + cos(phi_1 - phi_4)*b_1*b_4 + cos(phi_3 - phi_2)*b_2*b_3 + cos(phi_4 - phi_2)*b_2*b_4 + cos(phi_3 - phi_4)*b_3*b_4
		)));

	if (b_1 < 0 || b_2 < 0 || b_3 < 0 || b_4 < 0){
		f_A = 0;
	}

	return f_A;

}

//multi particle distribution : original
double multi_particle_dist_ov_normal_dist_ncoor(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{
	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));

	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8){ phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8){ phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8){ phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8){ phi_4 = 0; }

	double factor_exp = (b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
		+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4) + (b_3 - b_2)*(b_3 - b_2) + (by_3 - by_2)*(by_3 - by_2)
		+ (b_4 - b_2)*(b_4 - b_2) + (by_4 - by_2)*(by_4 - by_2) + (b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4);

	double f_A = Gauss_param_B*sqrt(Pi)* sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_2*b_2 + by_2*by_2)*sqrt(b_3*b_3 + by_3*by_3)*sqrt(b_4*b_4 + by_4*by_4)*
		exp(
		-Lamda_QCD*Lamda_QCD*(3.0*((b_1-b_2)*(b_1-b_2) + (by_1-by_2)*(by_1-by_2) + b_2*b_2 + by_2*by_2 + b_3*b_3 + by_3*by_3 + b_4*b_4 + by_4*by_4) -
		2.0*(cos(phi_1 - phi_2)*sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_2*b_2 + by_2*by_2) + cos(phi_1 - phi_3)*sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_3*b_3 + by_3*by_3)
		+ cos(phi_1 - phi_4)*sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_4*b_4 + by_4*by_4) + cos(phi_3 - phi_2)*sqrt(b_2*b_2 + by_2*by_2)*sqrt(b_3*b_3 + by_3*by_3)
		+ cos(phi_4 - phi_2)*sqrt(b_2*b_2 + by_2*by_2)*sqrt(b_4*b_4 + by_4*by_4) + cos(phi_3 - phi_4)*sqrt(b_3*b_3 + by_3*by_3)*sqrt(b_4*b_4 + by_4*by_4)
		)));

	return f_A;

}


double k_0_4_fuc(const double p_mom, const double r_dis, const double b_imp, const double cos_theta)
{
	double rtn_function = p_mom*exp(-r_dis*r_dis / Gauss_param / Gauss_param)*bessel_ink_0_4(p_mom*r_dis)*(1.0 - sol_BK_sp(gain_x_val(dist2_rad_imp(r_dis, b_imp, cos(cos_theta)))));
	if (isinf(rtn_function) || isnan(rtn_function)) { 
		cerr << setprecision(std::numeric_limits<double>::max_digits10) 
			<< "error nan p_mom " << p_mom << "\t r_dis " << r_dis << "\t b_imp "<<b_imp << "\t cos(cos_theta) " << cos(cos_theta) << "\n"; 
	}
	return rtn_function;
}

double k_0_4_factfunc(const double p_mom, const double b_imp, double *error)
{
	double integration = 0;

	double max_error = 0;

	for (int phi = 0; phi <= NPHI; ++phi)
	{
		double c_1 = h_phi_intmeasure;

		if (phi == 0 || phi == NPHI){ c_1 = c_1*0.5; }

		int lenaw;
		double tiny, aw[8000], i, err;

		lenaw = 8000;
		tiny = 1.0e-307;


		intdeiini(lenaw, tiny, 1.0e-15, aw);
		nfunc = 0;
		intdei(std::bind(k_0_4_fuc, p_mom, placeholders::_1, b_imp, phi*h_phi_intmeasure), 0.0, aw, &i, &err);
		//printf("I_5=int_0^infty sin(x)/x dx\n");
		//printf(" I_5= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);

		if (isinf(i) || isnan(i)) { cout << setprecision(std::numeric_limits<double>::max_digits10)
			<< "k_0_4_isinf isnan b_imp " << b_imp << "\t phi" << phi*h_phi_intmeasure << "\t error" << err << "\n"; }
		
		if (abs(err) > max_error) max_error = err;

		integration += c_1*i;

	}

	*error = max_error;


	return integration;
}


double calculate_k_0_4_ncvmp(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	 double& temp2_0_4, double& temp_0_2, double& temp2_0_2)
{
	double error_temp = 0;

	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));

	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8) { phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8) { phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8) { phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8) { phi_4 = 0; }


	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) >1.0) { phi_1 = 0; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) >1.0) { phi_2 = 0; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) >1.0) { phi_3 = 0; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) >1.0) { phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) < -1.0) { phi_1 = Pi; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) < -1.0) { phi_2 = Pi; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) < -1.0) { phi_3 = Pi; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) < -1.0) { phi_4 = Pi; }


	//double factor_exp = 3.0/4.0*((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
	//	+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4))
	//	-1.0/2.0*((b_1 - b_2)*(b_1 - b_3)+ (by_1 - by_2)*(by_1 - by_3) + (b_1 - b_4)*(b_1 - b_3) + (by_1 - by_4)*(by_1 - by_3) 
	//		+ (b_1 - b_2)*(b_1 - b_4) + (by_1 - by_2)*(by_1 - by_4)); 
	double factor_exp = -1.0 / 2.0*(b_1*(b_2 + b_3 + b_4) + by_1*(by_2 + by_3 + by_4)
				+ b_2*b_3 + by_2*by_3 + b_4*b_3 + by_4*by_3 + b_2*b_4 + by_2*by_4);

	double factor_exp_k_0_2 = 1.0/4.0 *(b_1*b_1+by_1*by_1+b_2*b_2+by_2*by_2)-2.0*(b_1*b_2+ by_1*by_2);

	double k_func_b_1 = k_0_4_factfunc(Upper_momk, sqrt(b_1*b_1 + by_1*by_1), &error_temp);
	double k_func_b_2 = k_0_4_factfunc(Upper_momk, sqrt(b_2*b_2 + by_2*by_2), &error_temp);

	double nuclear_confinement_factor_2 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus *(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2));
	double nuclear_confinement_factor_4 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus
		*(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2 + b_3*b_3 + by_3*by_3 + b_4*b_4 + by_4*by_4));

	double f_A =  exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*nuclear_confinement_factor_4
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_0_4_factfunc(Upper_momk, sqrt(b_3*b_3 + by_3*by_3), &error_temp)
		*k_0_4_factfunc(Upper_momk, sqrt(b_4*b_4 + by_4*by_4), &error_temp);



	double f_A_k_0_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*nuclear_confinement_factor_2
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;

	double temp_test = (phi_1 + phi_2 - phi_3 - phi_4);

	if (isinf(temp_test) || isnan(temp_test)) {
		cout << "temp_test isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}
	if (isinf(f_A) || isnan(f_A)) {
		cout << "isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}

	temp_0_2 = f_A_k_0_2;
	temp2_0_2 = f_A_k_0_2*f_A_k_0_2;


	temp2_0_4 = f_A*f_A;

	return f_A;
}


double calculate_k_0_4_ncvmp_sp(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	double& temp2_0_4, double& temp_0_2, double& temp2_0_2)
{
	double error_temp = 0;

	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));

	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8) { phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8) { phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8) { phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8) { phi_4 = 0; }


	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) >1.0) { phi_1 = 0; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) >1.0) { phi_2 = 0; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) >1.0) { phi_3 = 0; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) >1.0) { phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) < -1.0) { phi_1 = Pi; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) < -1.0) { phi_2 = Pi; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) < -1.0) { phi_3 = Pi; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) < -1.0) { phi_4 = Pi; }

	double factor_exp = -1.0 / 2.0*(b_1*(b_2 + b_3 + b_4) + by_1*(by_2 + by_3 + by_4)
		+ b_2*b_3 + by_2*by_3 + b_4*b_3 + by_4*by_3 + b_2*b_4 + by_2*by_4);

	double factor_exp_k_0_2 = 1.0 / 4.0 *(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2) - 2.0*(b_1*b_2 + by_1*by_2);

	double constant_integration_1 = Pi * sqrt(Pi)*exp(-Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0)
		* boost::math::cyl_bessel_i(0, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0) * sqrt(Gauss_param*Gauss_param*Upper_momk*Upper_momk);

	double epsilon_over_p2 = Gauss_param*Gauss_param*Upper_momk*Upper_momk;

	double constant_integration_2 = Pi*sqrt(Pi) / epsilon_over_p2 / sqrt(epsilon_over_p2)
		*(8.0*epsilon_over_p2 - exp(-1.0 / epsilon_over_p2 / 8.0)*(1.0 + 8.0*epsilon_over_p2) * boost::math::cyl_bessel_i(0, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0)
			+ boost::math::cyl_bessel_i(1, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0));


	double k_func_b_1 = constant_integration_1 - K_0_4_sp( exp(-1.0/sqrt(b_1*b_1 + by_1*by_1)));
	double k_func_b_2 = constant_integration_1 - K_0_4_sp(exp(-1.0/sqrt(b_2*b_2 + by_2*by_2)));
	double k_func_b_3 = constant_integration_1 - K_0_4_sp(exp(-1.0 / sqrt(b_3*b_3 + by_3*by_3)));
	double k_func_b_4 = constant_integration_1 - K_0_4_sp(exp(-1.0 / sqrt(b_4*b_4 + by_4*by_4)));


	double nuclear_confinement_factor_2 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus *(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2));
	double nuclear_confinement_factor_4 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus
		*(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2 + b_3*b_3 + by_3*by_3 + b_4*b_4 + by_4*by_4));

	double f_A = exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*nuclear_confinement_factor_4
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_func_b_3
		*k_func_b_4;

	double f_A_k_0_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*nuclear_confinement_factor_2
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;

	double temp_test = (phi_1 + phi_2 - phi_3 - phi_4);

	if (isinf(temp_test) || isnan(temp_test)) {
		cout << "temp_test isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}
	if (isinf(f_A) || isnan(f_A)) {
		cout << "isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}

	temp_0_2 = f_A_k_0_2;
	temp2_0_2 = f_A_k_0_2*f_A_k_0_2;


	temp2_0_4 = f_A*f_A;

	return f_A;
}

double calculate_k_0_4_radisalcvmp(const double b_1, const double bp_1, const double b_2, const double bp_2,
	const double b_3, const double bp_3, const double b_4, const double bp_4,
	double& temp2_0_4, double& temp_0_2, double& temp2_0_2)
{
	double error_temp = 0;



	//double factor_exp = 3.0/4.0*((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
	//	+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4))
	//	-1.0/2.0*((b_1 - b_2)*(b_1 - b_3)+ (by_1 - by_2)*(by_1 - by_3) + (b_1 - b_4)*(b_1 - b_3) + (by_1 - by_4)*(by_1 - by_3) 
	//		+ (b_1 - b_2)*(b_1 - b_4) + (by_1 - by_2)*(by_1 - by_4)); 
	double factor_exp = -1.0 / 2.0*(b_1*b_2*cos(bp_1 - bp_2) + b_1*b_3*cos(bp_1 - bp_3) + b_1*b_4*cos(bp_1 - bp_4)
		+ b_2*b_3*cos(bp_3 - bp_2) + b_4*b_3*cos(bp_3 - bp_4) + b_2*b_4*cos(bp_4 - bp_2));

	double factor_exp_k_0_2 = 1.0 / 4.0 *(b_1*b_1 + b_2*b_2 ) - 2.0*(b_1*b_2*cos(bp_1-bp_2));

	double k_func_b_1 = k_0_4_factfunc(Upper_momk, b_1, &error_temp);
	double k_func_b_2 = k_0_4_factfunc(Upper_momk, b_2, &error_temp);

	double f_A = b_1*b_2*b_3*b_4*exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_0_4_factfunc(Upper_momk, b_3, &error_temp)
		*k_0_4_factfunc(Upper_momk, b_4, &error_temp);

	double f_A_k_0_2 = b_1*b_2*exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;


	temp_0_2 = f_A_k_0_2;
	temp2_0_2 = f_A_k_0_2*f_A_k_0_2;


	temp2_0_4 = f_A*f_A;

	return f_A;
}


double calculate_k_0_4_ncoor(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	const double r_1, const double ry_1, const double r_2, const double ry_2,
	const double r_3, const double ry_3, const double r_4, const double ry_4, double &temp_0_4)
{
	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));

	double phi_r1 = acos(ry_1 / sqrt(r_1*r_1 + ry_1*ry_1));
	double phi_r2 = acos(ry_2 / sqrt(r_2*r_2 + ry_2*ry_2));
	double phi_r3 = acos(ry_3 / sqrt(r_3*r_3 + ry_3*ry_3));
	double phi_r4 = acos(ry_4 / sqrt(r_4*r_4 + ry_4*ry_4));

	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8){ phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8){ phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8){ phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8){ phi_4 = 0; }

	if (abs(r_1) < 1e-8 && abs(ry_1) <1e-8){ phi_r1 = 0; }
	if (abs(r_2) < 1e-8 && abs(ry_2) <1e-8){ phi_r2 = 0; }
	if (abs(r_3) < 1e-8 && abs(ry_3) <1e-8){ phi_r3 = 0; }
	if (abs(r_4) < 1e-8 && abs(ry_4) <1e-8){ phi_r4 = 0; }


	double factor_exp = (b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
		+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4) + (b_3 - b_2)*(b_3 - b_2) + (by_3 - by_2)*(by_3 - by_2)
		+ (b_4 - b_2)*(b_4 - b_2) + (by_4 - by_2)*(by_4 - by_2) + (b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4);

	double f_A = Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Pi*Pi*Pi*Pi
		* sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_2*b_2 + by_2*by_2)*sqrt(b_3*b_3 + by_3*by_3)*sqrt(b_4*b_4 + by_4*by_4)
		*exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Pi*Pi*Pi*Pi
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*bessel_ink_0_4(Upper_momk*sqrt(r_1*r_1 + ry_1*ry_1))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_1, by_1, r_1, ry_1)))
		*bessel_ink_0_4(Upper_momk*sqrt(r_2*r_2 + ry_2*ry_2))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_2, by_2, r_2, ry_2)))
		*bessel_ink_0_4(Upper_momk*sqrt(r_3*r_3 + ry_3*ry_3))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_3, by_3, r_3, ry_3)))
		*bessel_ink_0_4(Upper_momk*sqrt(r_4*r_4 + ry_4*ry_4))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_4, by_4, r_4, ry_4)));

	temp_0_4 = 2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi*f_A
		*2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi*f_A;

	return f_A*2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi;
}


double k_2_4_fuc(const double p_mom, const double r_dis, const double b_imp, const double cos_theta)
{
	return p_mom*exp(-r_dis*r_dis / Gauss_param / Gauss_param)*bessel_ink_2_4(p_mom*r_dis)*(1.0 - sol_BK_sp(gain_x_val(dist2_rad_imp(r_dis, b_imp, cos(cos_theta)))))*cos(2.0*cos_theta);
}

double k_2_4_factfunc(const double p_mom, const double b_imp, double *error)
{
	double integration = 0;
	double max_error = 0;


	for (int phi = 0; phi <= NPHI; ++phi)
	{
		double c_1 = h_phi_intmeasure;

		if (phi == 0 || phi == NPHI){ c_1 = c_1*0.5; }

		int lenaw;
		double tiny, aw[8000], i, err;

		lenaw = 8000;
		tiny = 1.0e-307;


		intdeiini(lenaw, tiny, 1.0e-15, aw);
		nfunc = 0;
		intdei(std::bind(k_2_4_fuc, p_mom, placeholders::_1, b_imp, phi*h_phi_intmeasure), 0.0, aw, &i, &err);
		//printf("I_5=int_0^infty sin(x)/x dx\n");
		//printf(" I_5= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);

		if (isinf(i) || isnan(i)) { cout<< setprecision(std::numeric_limits<double>::max_digits10)
			<< "k_2_4_isinf isnan b_imp " << b_imp << "\t phi" << phi*h_phi_intmeasure << "\t error" << err << "\n"; }

		if (abs(err) > max_error) max_error = err;


		integration += c_1*i;

	}

	*error = max_error;


	return integration;
}


double calculate_k_2_4_ncvmp(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	 double& temp2_2_4, double& temp_2_2, double& temp2_2_2)
{
	double error_temp = 0;

	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));


	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8){ phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8){ phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8){ phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8){ phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) >1.0) { phi_1 = 0; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) >1.0) { phi_2 = 0; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) >1.0) { phi_3 = 0; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) >1.0) { phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) < -1.0) { phi_1 = Pi; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) < -1.0) { phi_2 = Pi; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) < -1.0) { phi_3 = Pi; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) < -1.0) { phi_4 = Pi; }




	//double factor_exp = 3.0 / 4.0*((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
	//	+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4))
	//	- 1.0 / 2.0*((b_1 - b_2)*(b_1 - b_3) + (by_1 - by_2)*(by_1 - by_3) + (b_1 - b_4)*(b_1 - b_3) + (by_1 - by_4)*(by_1 - by_3)
	//	+ (b_1 - b_2)*(b_1 - b_4) + (by_1 - by_2)*(by_1 - by_4));
	double factor_exp =  -1.0 / 2.0*(b_1*(b_2+b_3+b_4)+ by_1*(by_2 + by_3 + by_4) 
		+b_2*b_3+ by_2*by_3+ b_4*b_3 + by_4*by_3+ b_2*b_4 + by_2*by_4);
	double factor_exp_k_0_2 = 1.0/4.0 *(b_1*b_1+by_1*by_1+b_2*b_2+by_2*by_2)-2.0*(b_1*b_2 + by_1*by_2);

	double k_func_b_1 = k_2_4_factfunc(Upper_momk, sqrt(b_1*b_1 + by_1*by_1), &error_temp);
	double k_func_b_2 = k_2_4_factfunc(Upper_momk, sqrt(b_2*b_2 + by_2*by_2), &error_temp);

	double f_A =  exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*cos(2.0*(phi_1 + phi_2 - phi_3 - phi_4))
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_2_4_factfunc(Upper_momk, sqrt(b_3*b_3 + by_3*by_3), &error_temp)
		*k_2_4_factfunc(Upper_momk, sqrt(b_4*b_4 + by_4*by_4), &error_temp);

	double f_A_k_2_2 =  exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(phi_1 - phi_2))
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;

	double temp_test = (phi_1 + phi_2 - phi_3 - phi_4);

	if (isinf(temp_test) || isnan(temp_test)) {
		cout << "temp_test isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}
	if (isinf(f_A) || isnan(f_A)) {
		cout << "isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}

	temp_2_2 = f_A_k_2_2;
	temp2_2_2 = f_A_k_2_2*f_A_k_2_2;

	temp2_2_4 = f_A*f_A;

	return f_A;
}


double calculate_k_2_4_ncvmp_sp(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	double& temp2_2_4, double& temp_2_2, double& temp2_2_2)
{
	double error_temp = 0;

	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));


	if (abs(b_1) < 1e-8 && abs(by_1) < 1e-8) { phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) < 1e-8) { phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) < 1e-8) { phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) < 1e-8) { phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) > 1.0) { phi_1 = 0; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) > 1.0) { phi_2 = 0; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) > 1.0) { phi_3 = 0; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) > 1.0) { phi_4 = 0; }

	if (by_1 / sqrt(b_1*b_1 + by_1*by_1) < -1.0) { phi_1 = Pi; }
	if (by_2 / sqrt(b_2*b_2 + by_2*by_2) < -1.0) { phi_2 = Pi; }
	if (by_3 / sqrt(b_3*b_3 + by_3*by_3) < -1.0) { phi_3 = Pi; }
	if (by_4 / sqrt(b_4*b_4 + by_4*by_4) < -1.0) { phi_4 = Pi; }

	double factor_exp = -1.0 / 2.0*(b_1*(b_2 + b_3 + b_4) + by_1*(by_2 + by_3 + by_4)
		+ b_2*b_3 + by_2*by_3 + b_4*b_3 + by_4*by_3 + b_2*b_4 + by_2*by_4);
	double factor_exp_k_0_2 = 1.0 / 4.0 *(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2) - 2.0*(b_1*b_2 + by_1*by_2);

	double constant_integration_1 = Pi * sqrt(Pi)*exp(-Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0)
		* boost::math::cyl_bessel_i(0, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0) * sqrt(Gauss_param*Gauss_param*Upper_momk*Upper_momk);

	double epsilon_over_p2 = Gauss_param*Gauss_param*Upper_momk*Upper_momk;

	double constant_integration_2 = Pi*sqrt(Pi) / epsilon_over_p2 / sqrt(epsilon_over_p2)
		*(8.0*epsilon_over_p2 - exp(-1.0 / epsilon_over_p2 / 8.0)*(1.0 + 8.0*epsilon_over_p2) * boost::math::cyl_bessel_i(0, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0)
			+ boost::math::cyl_bessel_i(1, Gauss_param*Gauss_param*Upper_momk*Upper_momk / 8.0));


	double k_func_b_1 =- K_2_4_sp(exp(-1.0 / sqrt(b_1*b_1 + by_1*by_1)));
	double k_func_b_2 =- K_2_4_sp(exp(-1.0 / sqrt(b_2*b_2 + by_2*by_2)));
	double k_func_b_3 =- K_2_4_sp(exp(-1.0 / sqrt(b_3*b_3 + by_3*by_3)));
	double k_func_b_4 =- K_2_4_sp(exp(-1.0 / sqrt(b_4*b_4 + by_4*by_4)));

	double nuclear_confinement_factor_2 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus *(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2));
	double nuclear_confinement_factor_4 = exp(-Lamda_QCD_nucleus *Lamda_QCD_nucleus
		*(b_1*b_1 + by_1*by_1 + b_2*b_2 + by_2*by_2 + b_3*b_3 + by_3*by_3 + b_4*b_4 + by_4*by_4));

	double f_A = exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*cos(2.0*(phi_1 + phi_2 - phi_3 - phi_4))
		*nuclear_confinement_factor_4
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_func_b_3
		*k_func_b_4;

	double f_A_k_2_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(phi_1 - phi_2))
		*nuclear_confinement_factor_2
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;

	double temp_test = (phi_1 + phi_2 - phi_3 - phi_4);

	if (isinf(temp_test) || isnan(temp_test)) {
		cout << "temp_test isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}
	if (isinf(f_A) || isnan(f_A)) {
		cout << "isinf isnan b_1 " << sqrt(b_1*b_1 + by_1*by_1) << "\t b_2 " << sqrt(b_2*b_2 + by_2*by_2) << "\t b_3 " << sqrt(b_3*b_3 + by_3*by_3) << "\t b_4 " << sqrt(b_4*b_4 + by_4*by_4) << "\n";
	}

	temp_2_2 = f_A_k_2_2;
	temp2_2_2 = f_A_k_2_2*f_A_k_2_2;

	temp2_2_4 = f_A*f_A;

	return f_A;
}


double calculate_k_2_4_radisalcvmp(const double b_1, const double bp_1, const double b_2, const double bp_2,
	const double b_3, const double bp_3, const double b_4, const double bp_4,
	double& temp2_2_4, double& temp_2_2, double& temp2_2_2)
{
	double error_temp = 0;



	//double factor_exp = 3.0/4.0*((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
	//	+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4))
	//	-1.0/2.0*((b_1 - b_2)*(b_1 - b_3)+ (by_1 - by_2)*(by_1 - by_3) + (b_1 - b_4)*(b_1 - b_3) + (by_1 - by_4)*(by_1 - by_3) 
	//		+ (b_1 - b_2)*(b_1 - b_4) + (by_1 - by_2)*(by_1 - by_4)); 
	double factor_exp = -1.0 / 2.0*(b_1*b_2*cos(bp_1 - bp_2) + b_1*b_3*cos(bp_1 - bp_3) + b_1*b_4*cos(bp_1 - bp_4)
		+ b_2*b_3*cos(bp_3 - bp_2) + b_4*b_3*cos(bp_3 - bp_4) + b_2*b_4*cos(bp_4 - bp_2));

	double factor_exp_k_0_2 = 1.0 / 4.0 *(b_1*b_1 + b_2*b_2) - 2.0*(b_1*b_2*cos(bp_1 - bp_2));

	double k_func_b_1 = k_2_4_factfunc(Upper_momk, b_1, &error_temp);
	double k_func_b_2 = k_2_4_factfunc(Upper_momk, b_2, &error_temp);

	double f_A = b_1*b_2*b_3*b_4*exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*cos(2.0*(bp_1 + bp_2 - bp_3 - bp_4))
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2
		*k_2_4_factfunc(Upper_momk, b_3, &error_temp)
		*k_2_4_factfunc(Upper_momk, b_4, &error_temp);

	double f_A_k_2_2 = b_1*b_2* exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(bp_1 - bp_2))
		*Upper_momk*Upper_momk
		*k_func_b_1
		*k_func_b_2;


	temp_2_2 = f_A_k_2_2;
	temp2_2_2 = f_A_k_2_2*f_A_k_2_2;

	temp2_2_4 = f_A*f_A;

	return f_A;
}


double calculate_k_2_4_ncoor(const double b_1, const double by_1, const double b_2, const double by_2,
	const double b_3, const double by_3, const double b_4, const double by_4,
	const double r_1, const double ry_1, const double r_2, const double ry_2,
	const double r_3, const double ry_3, const double r_4, const double ry_4,double &temp_2_4)
{
	double phi_1 = acos(by_1 / sqrt(b_1*b_1 + by_1*by_1));
	double phi_2 = acos(by_2 / sqrt(b_2*b_2 + by_2*by_2));
	double phi_3 = acos(by_3 / sqrt(b_3*b_3 + by_3*by_3));
	double phi_4 = acos(by_4 / sqrt(b_4*b_4 + by_4*by_4));

	double phi_r1 = acos(ry_1 / sqrt(r_1*r_1 + ry_1*ry_1));
	double phi_r2 = acos(ry_2 / sqrt(r_2*r_2 + ry_2*ry_2));
	double phi_r3 = acos(ry_3 / sqrt(r_3*r_3 + ry_3*ry_3));
	double phi_r4 = acos(ry_4 / sqrt(r_4*r_4 + ry_4*ry_4));

	if (abs(b_1) < 1e-8 && abs(by_1) <1e-8){ phi_1 = 0; }
	if (abs(b_2) < 1e-8 && abs(by_2) <1e-8){ phi_2 = 0; }
	if (abs(b_3) < 1e-8 && abs(by_3) <1e-8){ phi_3 = 0; }
	if (abs(b_4) < 1e-8 && abs(by_4) <1e-8){ phi_4 = 0; }

	if (abs(r_1) < 1e-8 && abs(ry_1) <1e-8){ phi_r1 = 0; }
	if (abs(r_2) < 1e-8 && abs(ry_2) <1e-8){ phi_r2 = 0; }
	if (abs(r_3) < 1e-8 && abs(ry_3) <1e-8){ phi_r3 = 0; }
	if (abs(r_4) < 1e-8 && abs(ry_4) <1e-8){ phi_r4 = 0; }


	double factor_exp = (b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2) + (b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)
		+ (b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4) + (b_3 - b_2)*(b_3 - b_2) + (by_3 - by_2)*(by_3 - by_2)
		+ (b_4 - b_2)*(b_4 - b_2) + (by_4 - by_2)*(by_4 - by_2) + (b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4);

	double f_A = Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Gauss_param_B*Pi*Pi*Pi*Pi
		* sqrt(b_1*b_1 + by_1*by_1)*sqrt(b_2*b_2 + by_2*by_2)*sqrt(b_3*b_3 + by_3*by_3)*sqrt(b_4*b_4 + by_4*by_4)
		*exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*cos(2.0*(phi_1 + phi_2 - phi_3 - phi_4))
		*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Gauss_param*Pi*Pi*Pi*Pi
		*Upper_momk*Upper_momk*Upper_momk*Upper_momk
		*bessel_ink_2_4(Upper_momk*sqrt(r_1*r_1 + ry_1*ry_1))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_1, by_1, r_1, ry_1)))*cos(2.0*phi_r1)
		*bessel_ink_2_4(Upper_momk*sqrt(r_2*r_2 + ry_2*ry_2))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_2, by_2, r_2, ry_2)))*cos(2.0*phi_r2)
		*bessel_ink_2_4(Upper_momk*sqrt(r_3*r_3 + ry_3*ry_3))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_3, by_3, r_3, ry_3)))*cos(2.0*phi_r3)
		*bessel_ink_2_4(Upper_momk*sqrt(r_4*r_4 + ry_4*ry_4))*sol_BK_sp(gain_x_val(dist2_rad_imp_ncoor(b_4, by_4, r_4, ry_4)))*cos(2.0*phi_r4);

	temp_2_4 = 2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi*f_A
		*2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi*f_A;

	return f_A*2.0*Pi*2.0*Pi*2.0*Pi*2.0*Pi;
}

int print(){
#ifdef K_0_4_K_2_4_FUNC

	ostringstream ofilename;
	ofilename << "c_2_4_k_spline_initc_1_Nphi_" << NPHI << "_N_" << N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());

#else
	ostringstream ofilename;
	ofilename << "c_2_4_initc_1_Nphi_"<< NPHI <<"_N_"<< N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());
#endif

	ofs_res  << "#tau \t c_2_4 \t error \t s_2_4 \t error \t c_2_2 \t error \n";

	for (int i = 0; i  <= END_T / OUTPUT_DELTA_T + 1e-6; ++i){
		ofs_res << i*OUTPUT_DELTA_T << "\t" << c_2_4[i] <<"\t"<<error_c_2_4[i] << "\t" << s_2_4[i] << "\t" << error_s_2_4[i] << "\t" << c_2_2[i] << "\t" << error_c_2_2[i] << "\n";
		}
	return 0;
}


int print_impact(const double tau) {
	ostringstream ofilename_i;
	ofilename_i << "K_2_4_initc_1_tau_" << tau << "_Nphi_" << NPHI << ".txt";
	ofstream ofs_res_i(ofilename_i.str().c_str());


	ofs_res_i << "#tilde_b = exp(-1/b) \t K_0_4 \t error \t K_2_4 \t error \n";

	for (int i = 0; i < number_impact-1 ; ++i) {
		ofs_res_i << (i+1)*stepsize_impact << "\t" << K_0_4[i] << "\t" << error_K_0_4[i] << "\t" << K_2_4[i] << "\t" << error_K_2_4[i]  << "\n";
	}
	return 0;
}


double test_2pow_over_normal_function(const double x_value)
{
	return x_value*x_value;
}


int main()
{



	for (int int_tau = 0; int_tau <= END_T / OUTPUT_DELTA_T+1e-6;++int_tau){
		initial(int_tau*OUTPUT_DELTA_T);

#ifdef K_0_4_K_2_4_FUNC


#pragma omp parallel for num_threads(2)
		for (int n =0; n < number_impact-1; ++n)
		{
			double impact_parameter = -1.0 / log((n+1)*stepsize_impact);
			K_vector[n] = (n + 1)*stepsize_impact;
			double error = 0;

			K_0_4[n] = k_0_4_factfunc(Upper_momk, impact_parameter, &error_K_0_4[n]);
			K_2_4[n] = k_2_4_factfunc(Upper_momk, impact_parameter, &error_K_2_4[n]);

		}

		print_impact(int_tau*OUTPUT_DELTA_T);


		K_0_4_sp.set_points(K_vector, K_0_4);
		K_2_4_sp.set_points(K_vector, K_2_4);


#else

		//if not defined MAKE_FUNCTION
#ifndef MAKE_FUNCTION


#ifdef  POLAR_COORDINATE

		std::mt19937 mt_1, mt_2, mt_3, mt_4, mt_r1, mt_r2, mt_r3, mt_r4;
		std::uniform_real_distribution<double> uniform_01(0.0, 1.0), uniform_02(0.0, 1.0), uniform_03(0.0, 1.0), uniform_04(0.0, 1.0), 
			uniform_0r1(0.0, 1.0), uniform_0r2(0.0, 1.0), uniform_0r3(0.0, 1.0), uniform_0r4(0.0, 1.0);
		std::random_device seed_gen_1, seed_gen_2, seed_gen_3, seed_gen_4, seed_gen_r1, seed_gen_r2, seed_gen_r3, seed_gen_r4,
			seed_gen_p1, seed_gen_p2, seed_gen_p3, seed_gen_p4, seed_gen_pr1, seed_gen_pr2, seed_gen_pr3, seed_gen_pr4;
		std::default_random_engine engine_1(seed_gen_1()), engine_2(seed_gen_2()), engine_3(seed_gen_3()), engine_4(seed_gen_4()),
			engine_p1(seed_gen_p1()), engine_p2(seed_gen_p2()), engine_p3(seed_gen_p3()), engine_p4(seed_gen_p4());


		// 平均0.0、標準偏差Gauss_param_Bで分布させる
		std::normal_distribution<> dist_1(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_2(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_3(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
			dist_4(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0));



		double sum_k_0_4 = 0.0;
		double sum_k_2_4 = 0.0;
		double sum2_k_0_4 = 0.0;
		double sum2_k_2_4 = 0.0;
		double sum_k_0_2 = 0.0;
		double sum_k_2_2 = 0.0;
		double sum2_k_0_2 = 0.0;
		double sum2_k_2_2 = 0.0;

#pragma omp parallel for reduction(+:sum_k_0_4,sum_k_2_4,sum2_k_0_4,sum2_k_2_4,sum_k_0_2,sum_k_2_2,sum2_k_0_2,sum2_k_2_2) num_threads(4)
		for (int i = 0; i < N; ++i)
		{


			double bp1 = 2.0*Pi*uniform_01(engine_p1);
			double bp2 = 2.0*Pi*uniform_02(engine_p2);
			double bp3 = 2.0*Pi*uniform_01(engine_p3);
			double bp4 = 2.0*Pi*uniform_01(engine_p4);

			double temp2_0_4 = 0.0;
			double temp2_2_4 = 0.0;
			double temp_0_2 = 0.0;
			double temp_2_2 = 0.0;
			double temp2_0_2 = 0.0;
			double temp2_2_2 = 0.0;
			double temp_0_4 = calculate_k_0_4_radisalcvmp(abs(dist_1(engine_1)), bp1, abs(dist_2(engine_2)), bp2,
				abs(dist_3(engine_3)), bp3,abs( dist_4(engine_4)), bp4, temp2_0_4, temp_0_2, temp2_0_2);
			double temp_2_4 = calculate_k_2_4_radisalcvmp(abs(dist_1(engine_1)), bp1, abs(dist_2(engine_2)), bp2,
				abs(dist_3(engine_3)), bp3, abs(dist_4(engine_4)), bp4, temp2_2_4, temp_2_2, temp2_2_2);


			sum_k_0_4 += temp_0_4;

			sum_k_2_4 += temp_2_4;

			sum2_k_0_4 += temp2_0_4;
			sum2_k_2_4 += temp2_2_4;

			sum_k_0_2 += temp_0_2;

			sum_k_2_2 += temp_2_2;

			sum2_k_0_2 += temp2_0_2;
			sum2_k_2_2 += temp2_2_2;

		}


#else 

#ifndef MULTI_PARTON

		std::random_device seed_gen_1, seed_gen_2, seed_gen_3, seed_gen_4, seed_gen_r1, seed_gen_r2, seed_gen_r3, seed_gen_r4,
			seed_gen_p1, seed_gen_p2, seed_gen_p3, seed_gen_p4, seed_gen_pr1, seed_gen_pr2, seed_gen_pr3, seed_gen_pr4;
		std::default_random_engine engine_1(seed_gen_1()), engine_2(seed_gen_2()), engine_3(seed_gen_3()), engine_4(seed_gen_4()),
			engine_r1(seed_gen_r1()), engine_r2(seed_gen_r2()), engine_r3(seed_gen_r3()), engine_r4(seed_gen_r4()),
			engine_p1(seed_gen_p1()), engine_p2(seed_gen_p2()), engine_p3(seed_gen_p3()), engine_p4(seed_gen_p4()),
			engine_pr1(seed_gen_pr1()), engine_pr2(seed_gen_pr2()), engine_pr3(seed_gen_pr3()), engine_pr4(seed_gen_pr4());

		// 平均0.0、標準偏差Gauss_param_Bで分布させる
		std::normal_distribution<> dist_1(0.0, Gauss_param_B / sqrt(2.0)), dist_2(0.0, Gauss_param_B / sqrt(2.0)), dist_3(0.0, Gauss_param_B / sqrt(2.0)), dist_4(0.0, Gauss_param_B / sqrt(2.0)),
			dist_p1(0.0, Gauss_param_B / sqrt(2.0)), dist_p2(0.0, Gauss_param_B / sqrt(2.0)), dist_p3(0.0, Gauss_param_B / sqrt(2.0)), dist_p4(0.0, Gauss_param_B / sqrt(2.0));

		std::normal_distribution<> dist_r1(0.0, Gauss_param / sqrt(2.0)), dist_r2(0.0, Gauss_param / sqrt(2.0)), dist_r3(0.0, Gauss_param / sqrt(2.0)), dist_r4(0.0, Gauss_param / sqrt(2.0)),
			dist_pr1(0.0, Gauss_param / sqrt(2.0)), dist_pr2(0.0, Gauss_param / sqrt(2.0)), dist_pr3(0.0, Gauss_param / sqrt(2.0)), dist_pr4(0.0, Gauss_param / sqrt(2.0));

		double sum_k_0_4 = 0.0;
		double sum_k_2_4 = 0.0;
		double sum2_k_0_4 = 0.0;
		double sum2_k_2_4 = 0.0;
		double sim_pi = 0.0;
		double sim_pi_2 = 0;
		double sim_pi_3 = 0;

		for (int i = 0; i < N; ++i)
		{
			// 正規分布で乱数を生成する
			//double result = dist_1(engine_1);

			double temp2_0_4 = 0;
			double temp2_2_4 = 0;


			sum_k_0_4 += calculate_k_0_4_ncoor(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
				dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4),
				dist_r1(engine_r1), dist_pr1(engine_pr1), dist_r2(engine_r2), dist_pr2(engine_pr2),
				dist_r3(engine_r3), dist_pr3(engine_pr3), dist_r4(engine_r4), dist_pr4(engine_pr4),temp2_0_4);

			sum_k_2_4 += calculate_k_0_4_ncoor(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
				dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4),
				dist_r1(engine_r1), dist_pr1(engine_pr1), dist_r2(engine_r2), dist_pr2(engine_pr2),
				dist_r3(engine_r3), dist_pr3(engine_pr3), dist_r4(engine_r4), dist_pr4(engine_pr4),temp2_2_4);

			sum2_k_0_4 += temp2_0_4;
			sum2_k_2_4 += temp2_2_4;

			if (i != 0 && i % 10 == 0) {
				double ans = sum_k_0_4 / i;
				//	printf("%.10f\n", std::fabs(ans - 0.626617));
			}


		}

		cout << sum_k_2_4 / sum_k_0_4 << "\t" << sum2_k_2_4 / N - sum_k_2_4 / N*sum_k_2_4 / N << "\t" << sum2_k_0_4 / N - sum_k_0_4 / N*sum_k_0_4 / N << "\n";
		c_2_4[int_tau] = sum_k_2_4 / sum_k_0_4;
		double error_k_0_4 = sqrt(sum2_k_0_4 / N  - sum_k_0_4 / N*sum_k_0_4 / N ) / sqrt(N);
		double error_k_2_4 = sqrt(sum2_k_2_4 / N  - sum_k_2_4 / N*sum_k_2_4 / N )/sqrt(N);
		error_c_2_4[int_tau] = N / abs(sum_k_0_4)*sqrt(error_k_2_4*error_k_2_4 + c_2_4[int_tau] * c_2_4[int_tau] * error_k_0_4*error_k_0_4);

#else

std::random_device seed_gen_1, seed_gen_2, seed_gen_3, seed_gen_4, 
seed_gen_p1, seed_gen_p2, seed_gen_p3, seed_gen_p4;
std::default_random_engine engine_1(seed_gen_1()), engine_2(seed_gen_2()), engine_3(seed_gen_3()), engine_4(seed_gen_4()),
engine_p1(seed_gen_p1()), engine_p2(seed_gen_p2()), engine_p3(seed_gen_p3()), engine_p4(seed_gen_p4());

// 平均0.0、標準偏差Gauss_param_Bで分布させる generate 1/sigma/sqrt(2pi)exp(-(x-mu)^2/2/sigma^2)
//std::normal_distribution<> dist_1(0.0, Gauss_param_B / sqrt(2.0)), dist_2(0.0, Gauss_param_B / sqrt(2.0)), dist_3(0.0, Gauss_param_B / sqrt(2.0)), dist_4(0.0, Gauss_param_B / sqrt(2.0)),
//dist_p1(0.0, Gauss_param_B / sqrt(2.0)), dist_p2(0.0, Gauss_param_B / sqrt(2.0)), dist_p3(0.0, Gauss_param_B / sqrt(2.0)), dist_p4(0.0, Gauss_param_B / sqrt(2.0));

std::normal_distribution<> dist_1(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_2(0.0, 1.0/Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_3(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_4(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_p1(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_p2(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_p3(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_p4(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0));


double sum_k_0_4 = 0.0;
double sum_k_2_4 = 0.0;
double sum2_k_0_4 = 0.0;
double sum2_k_2_4 = 0.0;
double sum_k_0_2 = 0.0;
double sum_k_2_2 = 0.0;
double sum2_k_0_2 = 0.0;
double sum2_k_2_2 = 0.0;
//#pragma loop(hint_parallel(4))
#pragma omp parallel for reduction(+:sum_k_0_4,sum_k_2_4,sum2_k_0_4,sum2_k_2_4,sum_k_0_2,sum_k_2_2,sum2_k_0_2,sum2_k_2_2) num_threads(4)
for (int i = 0; i < N; ++i)
{
	// 正規分布で乱数を生成する
	//double result = dist_1(engine_1);

	double temp2_0_4 = 0.0;
	double temp2_2_4 = 0.0;
	double temp_0_2 = 0.0;
	double temp_2_2 = 0.0;
	double temp2_0_2 = 0.0;
	double temp2_2_2 = 0.0;
	double temp_0_4 = calculate_k_0_4_ncvmp(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_0_4, temp_0_2, temp2_0_2);
	double temp_2_4 = calculate_k_2_4_ncvmp(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_2_4, temp_2_2, temp2_2_2);


	sum_k_0_4 += temp_0_4;

	sum_k_2_4 += temp_2_4;

	sum2_k_0_4 += temp2_0_4;
	sum2_k_2_4 += temp2_2_4;

	sum_k_0_2 += temp_0_2;

	sum_k_2_2 += temp_2_2;

	sum2_k_0_2 += temp2_0_2;
	sum2_k_2_2 += temp2_2_2;

	//if (isinf(temp_2_4) || isnan(temp_2_4)) {	cout << "N " <<i << "\n";}
	//cout << "tau " << int_tau*OUTPUT_DELTA_T << " Number " << i << "\n";



}


#endif
			 
#endif


#else



#endif


cout << sum_k_2_4 / sum_k_0_4 << "\t" << sum2_k_2_4 / N - sum_k_2_4 / N*sum_k_2_4 / N << "\t" << sum2_k_0_4 / N - sum_k_0_4 / N*sum_k_0_4 / N << "\n";

c_2_4[int_tau] = sum_k_2_4 / sum_k_0_4 - 2.0*(sum_k_2_2 / sum_k_0_2)*(sum_k_2_2 / sum_k_0_2);
double c_2_4_4 = sum_k_2_4 / sum_k_0_4;
double c_2_4_2 = sum_k_2_2 / sum_k_0_2;
double error_k_0_4 = sqrt(sum2_k_0_4 / N - sum_k_0_4 / N*sum_k_0_4 / N) / sqrt(N-1);
double error_k_2_4 = sqrt(sum2_k_2_4 / N - sum_k_2_4 / N*sum_k_2_4 / N) / sqrt(N - 1);
double error_k_0_2 = sqrt(sum2_k_0_2 / N - sum_k_0_2 / N*sum_k_0_2 / N) / sqrt(N - 1);
double error_k_2_2 = sqrt(sum2_k_2_2 / N - sum_k_2_2 / N*sum_k_2_2 / N) / sqrt(N - 1);
double error_4 = N / abs(sum_k_0_4)*sqrt(error_k_2_4*error_k_2_4 + c_2_4_4 * c_2_4_4 * error_k_0_4*error_k_0_4);
double error_2 = N / abs(sum_k_0_2)*sqrt(error_k_2_2*error_k_2_2 + c_2_4_2 * c_2_4_2 * error_k_0_2*error_k_0_2);

c_2_2[int_tau] = c_2_4_2;
error_c_2_2[int_tau] = error_2;
s_2_4[int_tau] = c_2_4_4;
error_s_2_4[int_tau] = error_4;
error_c_2_4[int_tau] = sqrt(error_4*error_4 + 4.0*2.0*error_2*error_2*c_2_4_2*c_2_4_2);
	
#endif

std::random_device seed_gen_1, seed_gen_2, seed_gen_3, seed_gen_4,
seed_gen_p1, seed_gen_p2, seed_gen_p3, seed_gen_p4;
std::default_random_engine engine_1(seed_gen_1()), engine_2(seed_gen_2()), engine_3(seed_gen_3()), engine_4(seed_gen_4()),
engine_p1(seed_gen_p1()), engine_p2(seed_gen_p2()), engine_p3(seed_gen_p3()), engine_p4(seed_gen_p4());

// 平均0.0、標準偏差Gauss_param_Bで分布させる generate 1/sigma/sqrt(2pi)exp(-(x-mu)^2/2/sigma^2)
//std::normal_distribution<> dist_1(0.0, Gauss_param_B / sqrt(2.0)), dist_2(0.0, Gauss_param_B / sqrt(2.0)), dist_3(0.0, Gauss_param_B / sqrt(2.0)), dist_4(0.0, Gauss_param_B / sqrt(2.0)),
//dist_p1(0.0, Gauss_param_B / sqrt(2.0)), dist_p2(0.0, Gauss_param_B / sqrt(2.0)), dist_p3(0.0, Gauss_param_B / sqrt(2.0)), dist_p4(0.0, Gauss_param_B / sqrt(2.0));

std::normal_distribution<> dist_1(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_2(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_3(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_4(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_p1(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_p2(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)), dist_p3(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0)),
dist_p4(0.0, 1.0 / Lamda_QCD*sqrt(2.0) / sqrt(3.0));


double sum_k_0_4 = 0.0;
double sum_k_2_4 = 0.0;
double sum2_k_0_4 = 0.0;
double sum2_k_2_4 = 0.0;
double sum_k_0_2 = 0.0;
double sum_k_2_2 = 0.0;
double sum2_k_0_2 = 0.0;
double sum2_k_2_2 = 0.0;
//#pragma loop(hint_parallel(4))
#pragma omp parallel for reduction(+:sum_k_0_4,sum_k_2_4,sum2_k_0_4,sum2_k_2_4,sum_k_0_2,sum_k_2_2,sum2_k_0_2,sum2_k_2_2) num_threads(2)
for (int i = 0; i < N; ++i)
{
	// 正規分布で乱数を生成する
	//double result = dist_1(engine_1);

	double temp2_0_4 = 0.0;
	double temp2_2_4 = 0.0;
	double temp_0_2 = 0.0;
	double temp_2_2 = 0.0;
	double temp2_0_2 = 0.0;
	double temp2_2_2 = 0.0;
	double temp_0_4 = calculate_k_0_4_ncvmp_sp(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_0_4, temp_0_2, temp2_0_2);
	double temp_2_4 = calculate_k_2_4_ncvmp_sp(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_2_4, temp_2_2, temp2_2_2);


	sum_k_0_4 += temp_0_4;

	sum_k_2_4 += temp_2_4;

	sum2_k_0_4 += temp2_0_4;
	sum2_k_2_4 += temp2_2_4;

	sum_k_0_2 += temp_0_2;

	sum_k_2_2 += temp_2_2;

	sum2_k_0_2 += temp2_0_2;
	sum2_k_2_2 += temp2_2_2;

	//if (isinf(temp_2_4) || isnan(temp_2_4)) {	cout << "N " <<i << "\n";}
	//cout << "tau " << int_tau*OUTPUT_DELTA_T << " Number " << i << "\n";



}


cout << sum_k_2_4 / sum_k_0_4 << "\t" << sum2_k_2_4 / N - sum_k_2_4 / N*sum_k_2_4 / N << "\t" << sum2_k_0_4 / N - sum_k_0_4 / N*sum_k_0_4 / N << "\n";

c_2_4[int_tau] = sum_k_2_4 / sum_k_0_4 - 2.0*(sum_k_2_2 / sum_k_0_2)*(sum_k_2_2 / sum_k_0_2);
double c_2_4_4 = sum_k_2_4 / sum_k_0_4;
double c_2_4_2 = sum_k_2_2 / sum_k_0_2;
double error_k_0_4 = sqrt(sum2_k_0_4 / N - sum_k_0_4 / N*sum_k_0_4 / N) / sqrt(N - 1);
double error_k_2_4 = sqrt(sum2_k_2_4 / N - sum_k_2_4 / N*sum_k_2_4 / N) / sqrt(N - 1);
double error_k_0_2 = sqrt(sum2_k_0_2 / N - sum_k_0_2 / N*sum_k_0_2 / N) / sqrt(N - 1);
double error_k_2_2 = sqrt(sum2_k_2_2 / N - sum_k_2_2 / N*sum_k_2_2 / N) / sqrt(N - 1);
double error_4 = N / abs(sum_k_0_4)*sqrt(error_k_2_4*error_k_2_4 + c_2_4_4 * c_2_4_4 * error_k_0_4*error_k_0_4);
double error_2 = N / abs(sum_k_0_2)*sqrt(error_k_2_2*error_k_2_2 + c_2_4_2 * c_2_4_2 * error_k_0_2*error_k_0_2);

c_2_2[int_tau] = c_2_4_2;
error_c_2_2[int_tau] = error_2;
s_2_4[int_tau] = c_2_4_4;
error_s_2_4[int_tau] = error_4;
error_c_2_4[int_tau] = sqrt(error_4*error_4 + 4.0*2.0*error_2*error_2*c_2_4_2*c_2_4_2);

	}

	print();

	return 0;
}