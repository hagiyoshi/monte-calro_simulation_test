//From https://qiita.com/Ushio/items/0040b3c74a480c46c80c


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
int N = 30000;

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

vector<double> s_k_0_2(200, 0);
vector<double> error_s_k_0_2(200, 0);
vector<double> s_k_2_2(200, 0);
vector<double> error_s_k_2_2(200, 0);
vector<double> s_k_0_4(200, 0);
vector<double> error_s_k_0_4(200, 0);
vector<double> s_k_2_4(200, 0);
vector<double> error_s_k_2_4(200, 0);

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

	double f_A =  exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*k_func_b_1
		*k_func_b_2
		*k_0_4_factfunc(Upper_momk, sqrt(b_3*b_3 + by_3*by_3), &error_temp)
		*k_0_4_factfunc(Upper_momk, sqrt(b_4*b_4 + by_4*by_4), &error_temp);

	double f_A_k_0_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
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

	double k_func_b_1 = K_0_4_sp( exp(-1.0/sqrt(b_1*b_1 + by_1*by_1)));
	double k_func_b_2 = K_0_4_sp(exp(-1.0/sqrt(b_2*b_2 + by_2*by_2)));
	double k_func_b_3 = K_0_4_sp(exp(-1.0 / sqrt(b_3*b_3 + by_3*by_3)));
	double k_func_b_4 = K_0_4_sp(exp(-1.0 / sqrt(b_4*b_4 + by_4*by_4)));

	double f_A = exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*k_func_b_1
		*k_func_b_2
		*k_func_b_3
		*k_func_b_4;

	double f_A_k_0_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
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
		*k_func_b_1
		*k_func_b_2
		*k_0_4_factfunc(Upper_momk, b_3, &error_temp)
		*k_0_4_factfunc(Upper_momk, b_4, &error_temp);

	double f_A_k_0_2 = b_1*b_2*exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*k_func_b_1
		*k_func_b_2;


	temp_0_2 = f_A_k_0_2;
	temp2_0_2 = f_A_k_0_2*f_A_k_0_2;


	temp2_0_4 = f_A*f_A;

	return f_A;
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
		*k_func_b_1
		*k_func_b_2
		*k_2_4_factfunc(Upper_momk, sqrt(b_3*b_3 + by_3*by_3), &error_temp)
		*k_2_4_factfunc(Upper_momk, sqrt(b_4*b_4 + by_4*by_4), &error_temp);

	double f_A_k_2_2 =  exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(phi_1 - phi_2))
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

	double k_func_b_1 = K_2_4_sp(exp(-1.0 / sqrt(b_1*b_1 + by_1*by_1)));
	double k_func_b_2 = K_2_4_sp(exp(-1.0 / sqrt(b_2*b_2 + by_2*by_2)));
	double k_func_b_3 = K_2_4_sp(exp(-1.0 / sqrt(b_3*b_3 + by_3*by_3)));
	double k_func_b_4 = K_2_4_sp(exp(-1.0 / sqrt(b_4*b_4 + by_4*by_4)));

	double f_A = exp(-Lamda_QCD*Lamda_QCD*factor_exp)
		*cos(2.0*(phi_1 + phi_2 - phi_3 - phi_4))
		*k_func_b_1
		*k_func_b_2
		*k_func_b_3
		*k_func_b_4;

	double f_A_k_2_2 = exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(phi_1 - phi_2))
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
		*k_func_b_1
		*k_func_b_2
		*k_2_4_factfunc(Upper_momk, b_3, &error_temp)
		*k_2_4_factfunc(Upper_momk, b_4, &error_temp);

	double f_A_k_2_2 = b_1*b_2* exp(-Lamda_QCD*Lamda_QCD*factor_exp_k_0_2)
		*cos(2.0*(bp_1 - bp_2))
		*k_func_b_1
		*k_func_b_2;


	temp_2_2 = f_A_k_2_2;
	temp2_2_2 = f_A_k_2_2*f_A_k_2_2;

	temp2_2_4 = f_A*f_A;

	return f_A;
}


int print(){
#ifdef K_0_4_K_2_4_FUNC

	ostringstream ofilename;
	ofilename << "k_spline_initc_1_Nphi_" << NPHI << "_N_" << N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());

	ostringstream ofilename2;
	ofilename2 << "c_2_4_k_spline_initc_1_Nphi_" << NPHI << "_N_" << N << ".txt";
	ofstream ofs_res2(ofilename2.str().c_str());

#else
	ostringstream ofilename;
	ofilename << "c_2_4_initc_1_Nphi_"<< NPHI <<"_N_"<< N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());
#endif

	ofs_res  << "#tau \t s_k_0_2 \t error \t s_k_2_2 \t error \t s_k_0_4 \t error \t s_k_2_4 \t error \n";
	ofs_res2 << "#tau \t c_2_4 \t error \t s_2_4 \t error \t c_2_2 \t error \n";

	for (int i = 1; i  < 100 + 1e-6; ++i){
		ofs_res << i << "\t" << s_k_0_2[i] <<"\t"<<error_s_k_0_2[i] << "\t" << s_k_2_2[i] << "\t" << error_s_k_2_2[i] 
			<< "\t" << s_k_0_4[i] << "\t" << error_s_k_0_4[i] << "\t" << s_k_2_4[i] << "\t" << error_s_k_2_4[i] <<"\n";
		ofs_res2 << i << "\t" << c_2_4[i] << "\t" << error_c_2_4[i] << "\t" << s_2_4[i] << "\t" << error_s_2_4[i] << "\t" << c_2_2[i] << "\t" << error_c_2_2[i] << "\n";
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

	int int_tau = 2;
	initial(int_tau*OUTPUT_DELTA_T);

	ostringstream ofilename_i;
	ofilename_i << "k_0_4_impact_initc_1_Nphi_" << NPHI << ".txt";
	ofstream ofs_res_i(ofilename_i.str().c_str());

	ofs_res_i << "#inpact param \t k_0_4 \t error \t k_2_4 \t error \n";

	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	std::normal_distribution<> dist(0.0, 1.0);
	double x1_Gauss_integral = 0;
	double x2_Gauss_integral = 0;
	double vvvvv = 0;
	double errprprprprprr = 0;
	for (int i = 1; i < 200; ++i) {

		double impact_parameter = i*100.0;
		ofs_res_i << impact_parameter << "\t" << k_0_4_factfunc(Upper_momk, impact_parameter, &errprprprprprr) << "\t" << errprprprprprr
			<< "\t" << k_2_4_factfunc(Upper_momk, impact_parameter, &errprprprprprr) << "\t" << errprprprprprr << "\n";

	}

	


//		for (int int_N = 1; int_N <= 100+1e-6;++int_N){
//			N = int_N * 100000;
#ifdef K_0_4_K_2_4_FUNC


//#pragma omp parallel for num_threads(2)
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



#else 

#endif


#else



#endif

	
#endif


ostringstream ofilename_c;
ofilename_c << "c_2_4_initc_1.txt";
ofstream ofs_res_c(ofilename_c.str().c_str());
ostringstream ofilename_er;
ofilename_er << "c_2_4_er_initc_1.txt";
ofstream ofs_res_er(ofilename_er.str().c_str());


ofs_res_c << "#tau \t sum_k_0_4 \t sum2_k_0_4 \t sum_k_2_4 \t sum2_k_2_4 \t sum_k_0_2 \t sum2_k_0_2 \t sum_k_2_2 \t sum2_k_2_2 \n";

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

double factor_sigma_dist = 1.0 / sqrt(2.0*Pi)* Lamda_QCD*sqrt(3.0);

double sum_k_0_4 = 0.0;
double sum_k_2_4 = 0.0;
double sum2_k_0_4 = 0.0;
double sum2_k_2_4 = 0.0;
double sum_k_0_2 = 0.0;
double sum_k_2_2 = 0.0;
double sum2_k_0_2 = 0.0;
double sum2_k_2_2 = 0.0;
double subs_num = 0;
//#pragma loop(hint_parallel(4))
//#pragma omp parallel for reduction(+:sum_k_0_4,sum_k_2_4,sum2_k_0_4,sum2_k_2_4,sum_k_0_2,sum_k_2_2,sum2_k_0_2,sum2_k_2_2) num_threads(2)
for (int i = 1; i < N; ++i)
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

	//if (abs(temp_0_2)>4.0
	//	|| abs(temp_0_4)>1.0){ 
	//	subs_num += 1; 
	//	ofs_res_er << fixed << setprecision(std::numeric_limits<double>::max_digits10)
	//		<< "invalid integrand \t t_0_2 \t" << temp_0_2 << "\t t_0_4 \t" << temp_0_4 << "\n";
	//	continue;
	//}


	sum_k_0_4 += temp_0_4;

	sum_k_2_4 += temp_2_4;

	sum2_k_0_4 += temp2_0_4;
	sum2_k_2_4 += temp2_2_4;

	sum_k_0_2 += temp_0_2;

	sum_k_2_2 += temp_2_2;

	sum2_k_0_2 += temp2_0_2;
	sum2_k_2_2 += temp2_2_2;

	ofs_res_c << (i - subs_num) << "\t" << sum_k_0_4 / (i - subs_num) << "\t" << sum2_k_0_4 / (i - subs_num) << "\t" << sum_k_2_4 / (i - subs_num) << "\t" << sum2_k_2_4 / (i - subs_num) << "\t"
		<< sum_k_0_2 / (i - subs_num) << "\t" << sum2_k_0_2 / (i - subs_num) << "\t" << sum_k_2_2 / (i - subs_num) << "\t" << sum2_k_2_2 / (i - subs_num) << "\n";
}


	print();

	return 0;
}