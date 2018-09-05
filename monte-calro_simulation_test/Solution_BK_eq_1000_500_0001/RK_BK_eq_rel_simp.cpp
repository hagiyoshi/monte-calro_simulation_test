// This program is intended to solve the impact parameter dependent BK equation
// g_(Y+dY)(d^2(x,-x)) = g_Y(d^2(x,-x)) + dY (\int .... )
 
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
 
#include "Parameters.h"
#include "Spline.h"
 
using namespace std;
 
//#define DEBUG
 
//#define HIGH_PRECISION
 
//not to compare Berger Stasto
#define SOL_IN_D
 
/**
 * The size of nuclear
 */
const double Radius_Nuc = 1.0;
 
/**
 * The impact parameter
 */
const double IMPACTP_B = 1.0;
 
/**
 * The number of bins for calculating region
 * because of using Simpson rule in integral, we must take this number odd
 */
 
#define NREGION         1001
 
/**
 * The bins of integral for phi. 
 * because of using Simpson rule in integral, we must take this number odd
 */
 
#define NPHI            501
 
/**
 * The evolution step size.
 */
#define DELTA_T         0.001
#define OUTPUT_DELTA_T  0.02
#define END_T           2.0
 
//avoiding divergence
//#define epsilon 0.0007
 
/**
 * The solution of the BK equation
 */
vector<double> sol_BK_g(NREGION,0);
 
vector<double> vector_X_value(NREGION,0);
 
//For Runge Kutta calculation
vector<double> K1_g(NREGION, 0);
vector<double> K2_g(NREGION, 0);
vector<double> K3_g(NREGION, 0);
vector<double> K4_g(NREGION, 0);
 
/**
 * The evolution rapidity. tau = bar{alpha}_s * rapidity
 */
double tau = 0;
 
/**
 * the arbitrary parameter for initial conditon in the solution of the BK equation.
 */
const double init_c = 1.0;
 
// The grid size for 0 <= x <= R
const double h_leng_grid = Radius_Nuc / (NREGION-1.0);
 
// The grid size for phi
const double h_phi = 2.0 * Pi / (NPHI-1.0);
 
//cos_phi table
double cos_phi_table[NPHI];
double phi_table[NPHI];
 
// Initialization.
 
static struct TableInitializer {
  TableInitializer();
} table_initializer_;
 
 
TableInitializer::TableInitializer(){
  // 0 <= phi < 2 pi.
  cos_phi_table[0] = 1.0;
  cos_phi_table[NPHI-1] = 1.0;
  phi_table[0] = 0;
  phi_table[NPHI-1] = 2.0 * Pi;
  for (int i = 1; i < NPHI-1; i++) {
    double x = h_phi * i;
    phi_table[i] = x;
    cos_phi_table[i] = cos(x);
  }
}
 
//calculate x,y to d^2(x,y). (x cdot y = xy cos(phi) )
double dist_2(double x, double y, double cos_phi){
  return Radius_Nuc*Radius_Nuc*(x*x + y*y - 2.0*x*y*cos_phi) / (Radius_Nuc*Radius_Nuc + x*x) / (Radius_Nuc*Radius_Nuc + y*y);
}
 
//calculate d^2(x,-x) to x. (d^2 /= 0)
double gain_x_val(double dist_2){
  return sqrt(Radius_Nuc*Radius_Nuc*(-1.0 + 2.0 / dist_2 - 2.0 / dist_2*sqrt(1.0 - dist_2)));
}
 
 
void init(){
  //tau = 0;
 
  // Value_reg = x
  double Value_reg = 0;
  for (int i = 0; i < NREGION; i++){
 
    //initial condition g_(Y=0) = exp(-c*4*R^2*x^2/(R^2+x^2)/(R^2+x^2))
    sol_BK_g[i] = exp(-init_c * 4.0 * Radius_Nuc * Radius_Nuc * Value_reg * Value_reg 
		      / (Radius_Nuc * Radius_Nuc + Value_reg * Value_reg) / (Radius_Nuc * Radius_Nuc + Value_reg * Value_reg));
 
    vector_X_value[i] = Value_reg;
 
    Value_reg = Value_reg + h_leng_grid;
 
  }
#ifdef DEBUG
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK_g);
  ofstream ofs("debug.txt");
  for (int i = 0; i < NREGION; i++){
    ofs << vector_X_value[i] << "   " << sol_BK_g[i] << endl;
  }
 
  ofstream ofs2("debug_spline.txt");
  for (int j = 0; j < 2 * NREGION; ++j){
    double x_value = j*Radius_Nuc / (2.0 * NREGION -1.0);
    ofs2 << x_value << "    " << s_BK_sol(x_value) << endl;
  }
#endif
}
 
//for precise calculation
vector<double> Iterate_func(const vector<double> &OLD_BK, const vector<double> &NEW_BK , const double dtau){
  vector<double> N_NEW_BK(NREGION,0);
 
  tk::spline s_BK_OLD;
  s_BK_OLD.set_points(vector_X_value, OLD_BK);
 
  tk::spline s_BK_NEW;
  s_BK_NEW.set_points(vector_X_value, NEW_BK);
 
  for (int REGION = 0; REGION < NREGION; REGION++){
 
    double x_value = (REGION)*h_leng_grid;
 
    if (REGION < 1){
      N_NEW_BK[REGION] = 1.0;
      continue;
    }
 
    N_NEW_BK[REGION] = OLD_BK[REGION];
 
    // the integral
 
    for (int iphi = 0; iphi < NPHI - 1; iphi++){
      for (int iregion = 0; iregion < NREGION; iregion++){
	//measure of the integral * dtau
	double c_1 = dtau*h_phi / 2 / Pi*h_leng_grid;
 
	double z_value = (iregion)*h_leng_grid;
 
	if (iregion == 0 || iregion == NREGION - 1) {
	  c_1 *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iregion % 2 == 0){
	  c_1 *= 2.0 / 3.0;
	}
	else{
	  c_1 *= 4.0 / 3.0;
	}
 
	double simpson = 1.0;
 
	if (iphi == 0 || iphi == NPHI - 1) {
	  simpson *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iphi % 2 == 0){
	  simpson *= 2.0 / 3.0;
	}
	else{
	  simpson *= 4.0 / 3.0;
	}
 
	if (REGION == iregion && iphi == 0 || REGION == iregion && iphi == (NPHI - 1) / 2){
	  //avoiding divergence
	  continue;
	}
 
	//4x^2/(x+z)^2/(x-z)^2
                double factor_1 = z_value * 4.0 * x_value *x_value
		  / (x_value * x_value + z_value * z_value - 2.0*x_value*z_value*cos_phi_table[iphi])
		  / (x_value * x_value + z_value * z_value + 2.0*x_value*z_value*cos_phi_table[iphi]);
                //4(R^2/x)^2/(R^2/x - z)^2/(z + R^2/x)^2
                double factor_2 = z_value * 4.0 * Radius_Nuc*Radius_Nuc / x_value *Radius_Nuc*Radius_Nuc / x_value
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     - 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi])
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     + 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi]);
 
                //d(x,z) -> d(x',x') -> x'. 
                double d2_1x_z = gain_x_val(dist_2(vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x = gain_x_val(dist_2(-vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                // change cos_phi -> -cos_phi.
                double d2_1x2_z = gain_x_val(dist_2(Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x2 = gain_x_val(dist_2(-Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
 
#ifdef DEBUG
                double a = simpson * c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					  + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
                a += 1.0;
 
                ofs4 << d2_1x_z << endl << d2_z_2x << endl << d2_1x2_z << endl << d2_z_2x2 << endl;
#endif
 
                N_NEW_BK[REGION] += 1.0 /2.0 * ( c_1*(factor_1*(s_BK_OLD(d2_1x_z)*s_BK_OLD(d2_z_2x) - OLD_BK[REGION])
						      + factor_2*(s_BK_OLD(d2_1x2_z)*s_BK_OLD(d2_z_2x2) - OLD_BK[REGION])) )
		  + 1.0 / 2.0 * (c_1*(factor_1*(s_BK_NEW(d2_1x_z)*s_BK_NEW(d2_z_2x) - NEW_BK[REGION])
				      + factor_2*(s_BK_NEW(d2_1x2_z)*s_BK_NEW(d2_z_2x2) - NEW_BK[REGION])))
		  - 1.0 / 6.0*(c_1*(factor_1*(s_BK_NEW(d2_1x_z) - s_BK_OLD(d2_1x_z))*(s_BK_NEW(d2_z_2x) - s_BK_OLD(d2_z_2x))
				    + factor_2*(s_BK_NEW(d2_1x2_z) - s_BK_OLD(d2_1x2_z))*(s_BK_NEW(d2_z_2x2) - s_BK_OLD(d2_z_2x2))));
 
      }
    }
 
 
#ifdef HIGH_PRECISION
 
    double b = (N_NEW_BK[REGION] - OLD_BK[REGION]) / OLD_BK[REGION];
    if (b > 0.01){
      //b += 1.0;
      cout << tau << "    " << REGION << "    " << b << endl;
    }
#endif
 
  }
 
  return N_NEW_BK;
}
 
vector<double> f_one_step_BK(const vector<double> &sol_BK_g, const double dtau){
  vector<double> n_BK_g(NREGION,0);
 
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK_g);
 
#ifdef DEBUG
  ofstream ofs3("debug_ratio.txt");
  ofs3 << dtau << endl;
 
  ofstream ofs4("debug_x_value.txt");
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif 
  for (int REGION = 0; REGION < NREGION; REGION++){
 
    double x_value = (REGION )*h_leng_grid;
    if (REGION < 1)
      {
	n_BK_g[REGION] = 0;
	continue;
      }
 
 
 
    // the integral
 
    for (int iphi = 0; iphi < NPHI-1; iphi++){
      for (int iregion = 0; iregion < NREGION; iregion++){
	//measure of the integral * dtau
	double c_1 = h_phi / 2 / Pi*h_leng_grid;
 
	double z_value = (iregion )*h_leng_grid;
 
	if (iregion == 0 || iregion == NREGION - 1) {
	  c_1 *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iregion % 2 == 0){
	  c_1 *= 2.0 / 3.0;
	}
	else{
	  c_1 *= 4.0 / 3.0;
	}
 
	double simpson = 1.0;
 
	if (iphi == 0 || iphi == NPHI - 1) {
	  simpson *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iphi % 2 == 0){
	  simpson *= 2.0 / 3.0;
	}
	else{
	  simpson *= 4.0 / 3.0;
	}
 
	if (REGION == iregion && iphi == 0 || REGION == iregion && iphi == (NPHI -1)/2){
	  //avoiding divergence
	  continue;
	}
 
	//4x^2/(x+z)^2/(x-z)^2
                    double factor_1 = z_value * 4.0 * x_value *x_value
		      / (x_value * x_value + z_value * z_value - 2.0*x_value*z_value*cos_phi_table[iphi])
		      / (x_value * x_value + z_value * z_value + 2.0*x_value*z_value*cos_phi_table[iphi]);
                    //4(R^2/x)^2/(R^2/x - z)^2/(z + R^2/x)^2
                    double factor_2 = z_value * 4.0 * Radius_Nuc*Radius_Nuc / x_value *Radius_Nuc*Radius_Nuc / x_value
		      / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
			 - 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi])
		      / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
			 + 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi]);
 
 
                    //d(x,z) -> d(x',x') -> x'. 
                    double d2_1x_z = gain_x_val(dist_2(vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                    double d2_z_2x = gain_x_val(dist_2(-vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                    // change cos_phi -> -cos_phi.
                    double d2_1x2_z = gain_x_val(dist_2(Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
                    double d2_z_2x2 = gain_x_val(dist_2(-Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
 
#ifdef DEBUG
                    double a = simpson * c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					      + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
                    a += 1.0;
 
                    ofs4 << d2_1x_z << endl << d2_z_2x << endl << d2_1x2_z << endl << d2_z_2x2 << endl;
#endif
                     
                    n_BK_g[REGION] += c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					   + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
 
      }
    }
 
 
#ifdef HIGH_PRECISION
 
    double b =( n_BK_g[REGION] - sol_BK_g[REGION] )/ sol_BK_g[REGION];
    //if (b > 0.01){
    //b += 1.0;
    //  cout << REGION << " " << b << endl;
    n_BK_g = Iterate_func(sol_BK_g, n_BK_g, dtau);
    //}
#endif
 
  }
 
 
  return n_BK_g;
}
 
 
vector<double> s_one_step_BK(const vector<double> K1_g, const vector<double> &sol_BK_g, const double dtau){
  vector<double> n_BK_g(NREGION, 0);
 
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK_g);
 
  tk::spline s_K1_g;
  s_K1_g.set_points(vector_X_value, K1_g);
 
#ifdef DEBUG
  ofstream ofs3("debug_ratio.txt");
  ofs3 << dtau << endl;
 
  ofstream ofs4("debug_x_value.txt");
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif 
  for (int REGION = 0; REGION < NREGION; REGION++){
 
    double x_value = (REGION)*h_leng_grid;
    if (REGION < 1)
      {
	n_BK_g[REGION] = 0;
	continue;
      }
 
    // the integral
 
    for (int iphi = 0; iphi < NPHI - 1; iphi++){
      for (int iregion = 0; iregion < NREGION; iregion++){
	//measure of the integral * dtau
	double c_1 = h_phi / 2 / Pi*h_leng_grid;
 
	double z_value = (iregion)*h_leng_grid;
 
	if (iregion == 0 || iregion == NREGION - 1) {
	  c_1 *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iregion % 2 == 0){
	  c_1 *= 2.0 / 3.0;
	}
	else{
	  c_1 *= 4.0 / 3.0;
	}
 
	double simpson = 1.0;
 
	if (iphi == 0 || iphi == NPHI - 1) {
	  simpson *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iphi % 2 == 0){
	  simpson *= 2.0 / 3.0;
	}
	else{
	  simpson *= 4.0 / 3.0;
	}
 
	if (REGION == iregion && iphi == 0 || REGION == iregion && iphi == (NPHI - 1) / 2){
	  //avoiding divergence
	  continue;
	}
 
	//4x^2/(x+z)^2/(x-z)^2
                double factor_1 = z_value * 4.0 * x_value *x_value
		  / (x_value * x_value + z_value * z_value - 2.0*x_value*z_value*cos_phi_table[iphi])
		  / (x_value * x_value + z_value * z_value + 2.0*x_value*z_value*cos_phi_table[iphi]);
                //4(R^2/x)^2/(R^2/x - z)^2/(z + R^2/x)^2
                double factor_2 = z_value * 4.0 * Radius_Nuc*Radius_Nuc / x_value *Radius_Nuc*Radius_Nuc / x_value
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     - 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi])
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     + 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi]);
 
                //d(x,z) -> d(x',x') -> x'. 
                double d2_1x_z = gain_x_val(dist_2(vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x = gain_x_val(dist_2(-vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                // change cos_phi -> -cos_phi.
                double d2_1x2_z = gain_x_val(dist_2(Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x2 = gain_x_val(dist_2(-Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
 
#ifdef DEBUG
                double a = simpson * c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					  + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
                a += 1.0;
 
                ofs4 << d2_1x_z << endl << d2_z_2x << endl << d2_1x2_z << endl << d2_z_2x2 << endl;
#endif
 
                n_BK_g[REGION] += c_1*(factor_1*((s_BK_sol(d2_1x_z) + dtau / 2.0*s_K1_g(d2_1x_z))*(s_BK_sol(d2_z_2x) + dtau / 2.0*s_K1_g(d2_z_2x)) - (sol_BK_g[REGION]+dtau/2.0*K1_g[REGION]))
				       + factor_2*((s_BK_sol(d2_1x2_z) + dtau / 2.0*s_K1_g(d2_1x2_z))*(s_BK_sol(d2_z_2x2) + dtau / 2.0*s_K1_g(d2_z_2x2)) - (sol_BK_g[REGION] + dtau / 2.0*K1_g[REGION])));
 
 
      }
    }
 
 
#ifdef HIGH_PRECISION
 
    double b = (n_BK_g[REGION] - sol_BK_g[REGION]) / sol_BK_g[REGION];
    //if (b > 0.01){
    //b += 1.0;
    //  cout << REGION << " " << b << endl;
    n_BK_g = Iterate_func(sol_BK_g, n_BK_g, dtau);
    //}
#endif
 
  }
 
 
  return n_BK_g;
}
 
vector<double> t_one_step_BK(const vector<double> K2_g, const vector<double> &sol_BK_g, const double dtau){
  vector<double> n_BK_g(NREGION, 0);
 
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK_g);
 
  tk::spline s_K2_g;
  s_K2_g.set_points(vector_X_value, K2_g);
 
#ifdef DEBUG
  ofstream ofs3("debug_ratio.txt");
  ofs3 << dtau << endl;
 
  ofstream ofs4("debug_x_value.txt");
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif 
  for (int REGION = 0; REGION < NREGION; REGION++){
 
    double x_value = (REGION)*h_leng_grid;
    if (REGION < 1)
      {
	n_BK_g[REGION] = 0;
	continue;
      }
 
    // the integral
 
    for (int iphi = 0; iphi < NPHI - 1; iphi++){
      for (int iregion = 0; iregion < NREGION; iregion++){
	//measure of the integral * dtau
	double c_1 = h_phi / 2 / Pi*h_leng_grid;
 
	double z_value = (iregion)*h_leng_grid;
 
	if (iregion == 0 || iregion == NREGION - 1) {
	  c_1 *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iregion % 2 == 0){
	  c_1 *= 2.0 / 3.0;
	}
	else{
	  c_1 *= 4.0 / 3.0;
	}
 
	double simpson = 1.0;
 
	if (iphi == 0 || iphi == NPHI - 1) {
	  simpson *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iphi % 2 == 0){
	  simpson *= 2.0 / 3.0;
	}
	else{
	  simpson *= 4.0 / 3.0;
	}
 
	if (REGION == iregion && iphi == 0 || REGION == iregion && iphi == (NPHI - 1) / 2){
	  //avoiding divergence
	  continue;
	}
 
	//4x^2/(x+z)^2/(x-z)^2
                double factor_1 = z_value * 4.0 * x_value *x_value
		  / (x_value * x_value + z_value * z_value - 2.0*x_value*z_value*cos_phi_table[iphi])
		  / (x_value * x_value + z_value * z_value + 2.0*x_value*z_value*cos_phi_table[iphi]);
                //4(R^2/x)^2/(R^2/x - z)^2/(z + R^2/x)^2
                double factor_2 = z_value * 4.0 * Radius_Nuc*Radius_Nuc / x_value *Radius_Nuc*Radius_Nuc / x_value
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     - 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi])
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     + 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi]);
 
                //d(x,z) -> d(x',x') -> x'. 
                double d2_1x_z = gain_x_val(dist_2(vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x = gain_x_val(dist_2(-vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                // change cos_phi -> -cos_phi.
                double d2_1x2_z = gain_x_val(dist_2(Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x2 = gain_x_val(dist_2(-Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
 
#ifdef DEBUG
                double a = simpson * c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					  + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
                a += 1.0;
 
                ofs4 << d2_1x_z << endl << d2_z_2x << endl << d2_1x2_z << endl << d2_z_2x2 << endl;
#endif
 
                n_BK_g[REGION] += c_1*(factor_1*((s_BK_sol(d2_1x_z) + dtau / 2.0*s_K2_g(d2_1x_z))*(s_BK_sol(d2_z_2x) + dtau / 2.0*s_K2_g(d2_z_2x)) - (sol_BK_g[REGION] + dtau / 2.0*K2_g[REGION]))
				       + factor_2*((s_BK_sol(d2_1x2_z) + dtau / 2.0*s_K2_g(d2_1x2_z))*(s_BK_sol(d2_z_2x2) + dtau / 2.0*s_K2_g(d2_z_2x2)) - (sol_BK_g[REGION] + dtau / 2.0*K2_g[REGION])));
 
 
      }
    }
 
 
#ifdef HIGH_PRECISION
 
    double b = (n_BK_g[REGION] - sol_BK_g[REGION]) / sol_BK_g[REGION];
    //if (b > 0.01){
    //b += 1.0;
    //  cout << REGION << " " << b << endl;
    n_BK_g = Iterate_func(sol_BK_g, n_BK_g, dtau);
    //}
#endif
 
  }
 
 
  return n_BK_g;
}
 
vector<double> fourth_one_step_BK(const vector<double> K3_g, const vector<double> &sol_BK_g, const double dtau){
  vector<double> n_BK_g(NREGION, 0);
 
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK_g);
 
  tk::spline s_K3_g;
  s_K3_g.set_points(vector_X_value, K3_g);
 
#ifdef DEBUG
  ofstream ofs3("debug_ratio.txt");
  ofs3 << dtau << endl;
 
  ofstream ofs4("debug_x_value.txt");
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif 
  for (int REGION = 0; REGION < NREGION; REGION++){
 
    double x_value = (REGION)*h_leng_grid;
    if (REGION < 1)
      {
	n_BK_g[REGION] = 0;
	continue;
      }
 
    // the integral
 
    for (int iphi = 0; iphi < NPHI - 1; iphi++){
      for (int iregion = 0; iregion < NREGION; iregion++){
	//measure of the integral * dtau
	double c_1 = h_phi / 2 / Pi*h_leng_grid;
 
	double z_value = (iregion)*h_leng_grid;
 
	if (iregion == 0 || iregion == NREGION - 1) {
	  c_1 *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iregion % 2 == 0){
	  c_1 *= 2.0 / 3.0;
	}
	else{
	  c_1 *= 4.0 / 3.0;
	}
 
	double simpson = 1.0;
 
	if (iphi == 0 || iphi == NPHI - 1) {
	  simpson *= 1.0 / 3.0;  // Simpson rule.
	}
	else if (iphi % 2 == 0){
	  simpson *= 2.0 / 3.0;
	}
	else{
	  simpson *= 4.0 / 3.0;
	}
 
	if (REGION == iregion && iphi == 0 || REGION == iregion && iphi == (NPHI - 1) / 2){
	  //avoiding divergence
	  continue;
	}
 
	//4x^2/(x+z)^2/(x-z)^2
                double factor_1 = z_value * 4.0 * x_value *x_value
		  / (x_value * x_value + z_value * z_value - 2.0*x_value*z_value*cos_phi_table[iphi])
		  / (x_value * x_value + z_value * z_value + 2.0*x_value*z_value*cos_phi_table[iphi]);
                //4(R^2/x)^2/(R^2/x - z)^2/(z + R^2/x)^2
                double factor_2 = z_value * 4.0 * Radius_Nuc*Radius_Nuc / x_value *Radius_Nuc*Radius_Nuc / x_value
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     - 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi])
		  / (Radius_Nuc*Radius_Nuc / x_value * Radius_Nuc*Radius_Nuc / x_value + z_value * z_value
		     + 2.0*Radius_Nuc*Radius_Nuc / x_value*z_value*cos_phi_table[iphi]);
 
                //d(x,z) -> d(x',x') -> x'. 
                double d2_1x_z = gain_x_val(dist_2(vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x = gain_x_val(dist_2(-vector_X_value[REGION], vector_X_value[iregion], cos_phi_table[iphi]));
                // change cos_phi -> -cos_phi.
                double d2_1x2_z = gain_x_val(dist_2(Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
                double d2_z_2x2 = gain_x_val(dist_2(-Radius_Nuc*Radius_Nuc / vector_X_value[REGION], -vector_X_value[iregion], cos_phi_table[iphi]));
 
#ifdef DEBUG
                double a = simpson * c_1*(factor_1*(s_BK_sol(d2_1x_z)*s_BK_sol(d2_z_2x) - sol_BK_g[REGION])
					  + factor_2*(s_BK_sol(d2_1x2_z)*s_BK_sol(d2_z_2x2) - sol_BK_g[REGION]));
 
                a += 1.0;
 
                ofs4 << d2_1x_z << endl << d2_z_2x << endl << d2_1x2_z << endl << d2_z_2x2 << endl;
#endif
 
                n_BK_g[REGION] += c_1*(factor_1*((s_BK_sol(d2_1x_z) + dtau *s_K3_g(d2_1x_z))*(s_BK_sol(d2_z_2x) + dtau *s_K3_g(d2_z_2x)) - (sol_BK_g[REGION] + dtau *K3_g[REGION]))
				       + factor_2*((s_BK_sol(d2_1x2_z) + dtau *s_K3_g(d2_1x2_z))*(s_BK_sol(d2_z_2x2) + dtau *s_K3_g(d2_z_2x2)) - (sol_BK_g[REGION] + dtau *K3_g[REGION])));
 
 
      }
    }
 
 
#ifdef HIGH_PRECISION
 
    double b = (n_BK_g[REGION] - sol_BK_g[REGION]) / sol_BK_g[REGION];
    //if (b > 0.01){
    //b += 1.0;
    //  cout << REGION << " " << b << endl;
    n_BK_g = Iterate_func(sol_BK_g, n_BK_g, dtau);
    //}
#endif
 
  }
 
 
  return n_BK_g;
}
 
 
//Using Runge Kutta method
vector<double> one_step_BK(const vector<double> &sol_BK_g, const vector<double> K1_g, const vector<double> K2_g,
			   const vector<double> K3_g, const vector<double> K4_g, const double dtau){
 
  vector<double> next_BK_g(NREGION);
 
  next_BK_g[0] = 1.0;

#ifdef _OPENMP
#pragma omp parallel for
#endif 
  for (int i = 1; i < NREGION; ++i){
        next_BK_g[i] = sol_BK_g[i]
	  + dtau / 6.0*(K1_g[i]+2.0*K2_g[i]+2.0*K3_g[i]+K4_g[i]);
  }
 
  return next_BK_g;
 
}
 
void one_step(double dtau){
 
  //evolve using the BK equation
  K1_g = f_one_step_BK(sol_BK_g, dtau);
  K2_g = s_one_step_BK(K1_g, sol_BK_g, dtau);
  K3_g = s_one_step_BK(K2_g, sol_BK_g, dtau); 
  K4_g = fourth_one_step_BK(K3_g, sol_BK_g, dtau);
 
     
 
  tau += dtau; 
  sol_BK_g = one_step_BK(sol_BK_g, K1_g, K2_g, K3_g, K4_g, dtau);
}
 
 
// to compare Berger Stasto 1010.0671 Fig.2
vector<double> scatt_impact(const vector<double> &sol_BK_g,vector<double> &vector_X_val){
  const vector<double> &sol_BK = sol_BK_g;
  //log10 scale in x
  vector<double> S_BK_b(NREGION, 0);
 
  //The spline interpolation of the solution of BK equation using  Tino Kluge code.
  tk::spline s_BK_sol;
  s_BK_sol.set_points(vector_X_value, sol_BK);
 
  //dipole size in log10 scale [10^-3,10^3]
  const double d_dipole = 6.0 / (NREGION - 1);
  double dipole_size_l = -3.0;
 
  //b=(1,0) -> x=(1-r/2) , y=(1+r/2,0)
  for (int i = 0; i < NREGION; ++i){
    double dipole_size = pow(10.0, i*d_dipole+dipole_size_l);
    double d_r_b = gain_x_val(dist_2(IMPACTP_B - dipole_size / 2.0, IMPACTP_B + dipole_size / 2.0, 1.0));
    S_BK_b[i] = s_BK_sol(d_r_b);
    vector_X_val[i] = i*d_dipole + dipole_size_l;
  }
 
  return S_BK_b;
}
 
/**
 * Prints g_Y.
 */
void print_g() {
  ostringstream ofilename;
  ofilename << "res_" << NPHI << "_" << NREGION << "_" << DELTA_T << "_R_" << Radius_Nuc << "_t_" << tau << "_hipre.txt";
  ofstream ofs_res(ofilename.str().c_str());
 
  //evaluate the number of time step
  static bool first = true;
  if (first) {
    const double EPS = 1e-12;
    first = false;
    int n = 0;
    double next_tau = 0;
    for (;;) {
      n++;
      if (next_tau >= END_T - EPS) {
	break;
      }
      next_tau = min(next_tau + OUTPUT_DELTA_T, END_T);
    }
  }
 
#ifndef SOL_IN_D
  sol_BK_g = scatt_impact(sol_BK_g,vector_X_value);
  //print the solution of the BK equation
  cout << "# t = " << tau << endl;
  for (int i = 0; i < NREGION; i++){
    ofs_res << vector_X_value[i] << "   " << sol_BK_g[i] << endl;
  }
#else
 
  //print the solution of the BK equation
  ofs_res << "# t = " << tau << endl;
  for (int i = 0; i < NREGION; i++){
    ofs_res << vector_X_value[i] << "   " << sol_BK_g[i] << endl;
  }
 
#endif
}
 
int main(){
  time_t t0 = time(NULL);
  init();
 
  const double EPS = 1e-12;
 
  double next_tau = 0;
  print_g();
 
  for (;;) {
    int reunit_count = 0;
    if (tau >= END_T - EPS) {
      break;
    }
    next_tau = min(next_tau + OUTPUT_DELTA_T, END_T);
    while (tau < next_tau - EPS) {
      one_step(min(DELTA_T, next_tau - tau));
    }
    print_g();
  }
 
  time_t t1 = time(NULL);
}
