//From the paper 1705.03051

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
#define OUTPUT_DELTA_T  0.1
#define END_T           2.0


//number of ansamble?
const int N = 10000000;

const int NPHI = 40;

const double h_phi_intmeasure = 2.0*Pi / NPHI;

#define Nc 3

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
const double Upper_momk = 2.0;

/**
* The parameter \Lamda_QCD
*/
const double bar_Lamda_QCD = 0.241;

/**
* The parameter \Lamda_QCD
*/
const double Lamda_QCD = 1.0/2.0;

double momentum_Q = 0.1;

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



void initial(double tau){

	momentum_Q = tau;
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

//gammma function in the paper 1705.03051
double gamma_kannsuu(const double absolute_x)
{
	return 1.0 / 4.0 / Pi*absolute_x*absolute_x*log(1.0 / absolute_x / bar_Lamda_QCD + exp(1.0));
}



double k_pm_0(double absolute_x)
{
	return Upper_momk / 2.0 / Pi / absolute_x*_j1(absolute_x*Upper_momk);
}


double k_pm_2(double absolute_x)
{
	return Upper_momk*Upper_momk / 2.0 / Pi *bessel_ink_2_4(absolute_x*Upper_momk);
}

double Delta_E_21(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{
	double gam_14 = gamma_kannsuu(sqrt((b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4)));
	double gam_23 = gamma_kannsuu(sqrt((b_2 - b_3)*(b_2 - b_3) + (by_2 - by_3)*(by_2 - by_3)));
	double gam_24 = gamma_kannsuu(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double gam_13 = gamma_kannsuu(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));

	return momentum_Q * momentum_Q*(gam_14 + gam_23 - gam_24 - gam_13);
}

double V_0_12(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{

	double gam_14 = gamma_kannsuu(sqrt((b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4)));
	double gam_23 = gamma_kannsuu(sqrt((b_2 - b_3)*(b_2 - b_3) + (by_2 - by_3)*(by_2 - by_3)));
	double gam_24 = gamma_kannsuu(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double gam_13 = gamma_kannsuu(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double gam_12 = gamma_kannsuu(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));
	double gam_34 = gamma_kannsuu(sqrt((b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4)));

	return momentum_Q * momentum_Q*Nc / (Nc*Nc - 1.0)*(gam_13 + gam_24 - gam_12 - gam_34);
}


double V_12_0(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{

	double gam_14 = gamma_kannsuu(sqrt((b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4)));
	double gam_23 = gamma_kannsuu(sqrt((b_2 - b_3)*(b_2 - b_3) + (by_2 - by_3)*(by_2 - by_3)));
	double gam_24 = gamma_kannsuu(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double gam_13 = gamma_kannsuu(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double gam_12 = gamma_kannsuu(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));
	double gam_34 = gamma_kannsuu(sqrt((b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4)));

	return momentum_Q * momentum_Q*Nc / (Nc*Nc - 1.0)*(gam_14 + gam_23 - gam_12 - gam_34);
}


double D_pm_0(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{
	double integrand = 0;

	double gam_14 = gamma_kannsuu(sqrt((b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4)));
	double gam_23 = gamma_kannsuu(sqrt((b_2 - b_3)*(b_2 - b_3) + (by_2 - by_3)*(by_2 - by_3)));
	double gam_24 = gamma_kannsuu(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double gam_13 = gamma_kannsuu(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double gam_12 = gamma_kannsuu(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));
	double gam_34 = gamma_kannsuu(sqrt((b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4)));

	double factor_1 = exp(-momentum_Q*momentum_Q*gam_13)*k_pm_0(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double factor_2 = exp(-momentum_Q*momentum_Q*gam_24)*k_pm_0(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double factor_3 = (1.0 - (1.0 - exp(-Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4))) / Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4));
	double factor_4 = V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Nc
		+ V_0_12(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)*V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);

	

	if (abs(Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)) < 1.0e-16)
	{
		factor_3 = 1.0 / 2.0;
		factor_4 = V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Nc*Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)
			+ V_0_12(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)*V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);
	}
	integrand = factor_1*factor_2*factor_3*factor_4;

	return integrand;


}

double D_pm_2(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4)
{
	double integrand = 0;


	double phi_1 = acos((by_1- by_3 )/ sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double phi_2 = acos((by_2- by_4 )/ sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));

	if (abs(by_1 - by_3) < 1e-8 && abs(b_1 - b_3) <1e-8) { phi_1 = 0; }
	if (abs(by_2 - by_4) < 1e-8 && abs(b_2 - b_4) <1e-8) { phi_2 = 0; }


	if ((by_1 - by_3) / sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)) >1.0) { phi_1 = 0; }
	if ((by_2 - by_4) / sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)) >1.0) { phi_2 = 0; }

	if ((by_1 - by_3) / sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)) < -1.0) { phi_1 = Pi; }
	if ((by_2 - by_4) / sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)) < -1.0) { phi_2 = Pi; }

	double gam_14 = gamma_kannsuu(sqrt((b_1 - b_4)*(b_1 - b_4) + (by_1 - by_4)*(by_1 - by_4)));
	double gam_23 = gamma_kannsuu(sqrt((b_2 - b_3)*(b_2 - b_3) + (by_2 - by_3)*(by_2 - by_3)));
	double gam_24 = gamma_kannsuu(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double gam_13 = gamma_kannsuu(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double gam_12 = gamma_kannsuu(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));
	double gam_34 = gamma_kannsuu(sqrt((b_3 - b_4)*(b_3 - b_4) + (by_3 - by_4)*(by_3 - by_4)));

	double factor_1 = exp(-momentum_Q*momentum_Q*gam_13)*k_pm_2(sqrt((b_1 - b_3)*(b_1 - b_3) + (by_1 - by_3)*(by_1 - by_3)));
	double factor_2 = exp(-momentum_Q*momentum_Q*gam_24)*k_pm_2(sqrt((b_2 - b_4)*(b_2 - b_4) + (by_2 - by_4)*(by_2 - by_4)));
	double factor_3 = (1.0 - (1.0 - exp(-Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4))) / Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4));
	double factor_4 = V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Nc
		+ V_0_12(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)*V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);


	if (abs(Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)) < 1.0e-16)
	{
		factor_3 = 1.0 / 2.0;
		factor_4 = V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4) / Nc*Delta_E_21(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)
			+ V_0_12(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4)*V_12_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);
	}

	integrand = factor_1*factor_2*factor_3*factor_4*cos(2.0*(phi_1 - phi_2));

	return integrand;


}

double D_0(const double b_1, const double b_2,
	const double by_1, const double by_2)
{
	double integrand = 0;

	double gam_12 = gamma_kannsuu(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));

	double factor_1 = exp(-momentum_Q*momentum_Q*gam_12)*k_pm_0(sqrt((b_1 - b_2)*(b_1 - b_2) + (by_1 - by_2)*(by_1 - by_2)));

	integrand = factor_1;

	return integrand;
}


double kappa_0_2(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4,double& sum2_0_2)
{
	double k_0_2 = 0;

	k_0_2 = D_0(b_1, by_1, b_2, by_2)*D_0(b_1, by_1, b_2, by_2) + D_pm_0(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);

	sum2_0_2 = k_0_2*k_0_2;

	return k_0_2;

}


double kappa_2_2(const double b_1, const double b_2, const double b_3, const double b_4,
	const double by_1, const double by_2, const double by_3, const double by_4, double& sum2_2_2)
{
	double k_2_2 = 0;

	k_2_2 =  D_pm_2(b_1, b_2, b_3, b_4, by_1, by_2, by_3, by_4);

	sum2_2_2 = k_2_2*k_2_2;

	return k_2_2;

}


int print_FH() {
#ifdef K_0_4_K_2_4_FUNC

	ostringstream ofilename;
	ofilename << "FH_result_maxQ_" << END_T <<"_N_" << N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());

#else
	ostringstream ofilename;
	ofilename << "c_2_4_initc_1_Nphi_" << NPHI << "_N_" << N << ".txt";
	ofstream ofs_res(ofilename.str().c_str());
#endif

	ofs_res << "#mom_Q \t bar_mom_Q \t c_2_2 \t error \n";

	for (int i = 1; i <= END_T / OUTPUT_DELTA_T + 1e-6; ++i) {
		ofs_res << i*OUTPUT_DELTA_T << "\t" << sqrt(1.0/2.0/gamma_kannsuu(abs(sqrt(2.0)/(i*OUTPUT_DELTA_T)))) << "\t" << c_2_2[i] << "\t" << error_c_2_2[i] << "\n";
	}
	return 0;
}

int print_impact(const double tau) {
	ostringstream ofilename_i;
	ofilename_i << "K_2_4_initc_1_tau_" << tau << "_Nphi_" << NPHI << ".txt";
	ofstream ofs_res_i(ofilename_i.str().c_str());


	ofs_res_i << "#tilde_b = exp(-1/b) \t K_0_4 \t error \t K_2_4 \t error \n";

	for (int i = 1; i < number_impact-1 ; ++i) {
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

std::random_device seed_gen_1, seed_gen_2, seed_gen_3, seed_gen_4,
seed_gen_p1, seed_gen_p2, seed_gen_p3, seed_gen_p4;
std::default_random_engine engine_1(seed_gen_1()), engine_2(seed_gen_2()), engine_3(seed_gen_3()), engine_4(seed_gen_4()),
engine_p1(seed_gen_p1()), engine_p2(seed_gen_p2()), engine_p3(seed_gen_p3()), engine_p4(seed_gen_p4());

// 平均0.0、標準偏差Gauss_param_Bで分布させる generate 1/sigma/sqrt(2pi)exp(-(x-mu)^2/2/sigma^2)
//std::normal_distribution<> dist_1(0.0, Gauss_param_B / sqrt(2.0)), dist_2(0.0, Gauss_param_B / sqrt(2.0)), dist_3(0.0, Gauss_param_B / sqrt(2.0)), dist_4(0.0, Gauss_param_B / sqrt(2.0)),
//dist_p1(0.0, Gauss_param_B / sqrt(2.0)), dist_p2(0.0, Gauss_param_B / sqrt(2.0)), dist_p3(0.0, Gauss_param_B / sqrt(2.0)), dist_p4(0.0, Gauss_param_B / sqrt(2.0));

std::normal_distribution<> dist_1(0.0, sqrt(Gauss_param_B)), dist_2(0.0, sqrt(Gauss_param_B)), dist_3(0.0, sqrt(Gauss_param_B)),
dist_4(0.0, sqrt(Gauss_param_B)),
dist_p1(0.0, sqrt(Gauss_param_B)), dist_p2(0.0, sqrt(Gauss_param_B)), dist_p3(0.0, sqrt(Gauss_param_B)),
dist_p4(0.0, sqrt(Gauss_param_B));


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
	double temp2_0_2 = 0.0;
	double temp2_2_2 = 0.0;
	double temp_0_2 = kappa_0_2(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_0_2);
	double temp_2_2 = kappa_2_2(dist_1(engine_1), dist_p1(engine_p1), dist_2(engine_2), dist_p2(engine_p2),
		dist_3(engine_3), dist_p3(engine_p3), dist_4(engine_4), dist_p4(engine_p4), temp2_2_2);



	sum_k_0_2 += temp_0_2;

	sum_k_2_2 += temp_2_2;

	sum2_k_0_2 += temp2_0_2;
	sum2_k_2_2 += temp2_2_2;

	//if (isinf(temp_2_4) || isnan(temp_2_4)) {	cout << "N " <<i << "\n";}
	//cout << "tau " << int_tau*OUTPUT_DELTA_T << " Number " << i << "\n";



}



double c_2_4_2 = sum_k_2_2 / sum_k_0_2;
double error_k_0_2 = sqrt(sum2_k_0_2 / N - sum_k_0_2 / N*sum_k_0_2 / N) / sqrt(N - 1);
double error_k_2_2 = sqrt(sum2_k_2_2 / N - sum_k_2_2 / N*sum_k_2_2 / N) / sqrt(N - 1);
double error_2 = N / abs(sum_k_0_2)*sqrt(error_k_2_2*error_k_2_2 + c_2_4_2 * c_2_4_2 * error_k_0_2*error_k_0_2);

c_2_2[int_tau] = c_2_4_2;
error_c_2_2[int_tau] = error_2;
	}

	print_FH();

	return 0;
}