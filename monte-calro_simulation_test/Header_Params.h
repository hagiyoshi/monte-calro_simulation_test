#pragma once

//#define OUTPUT_FT_NOISE
#define EVOLUTION

#define END_Y	10.0
#define OUTPUT_Y	0.5
#define DELTA_Y	0.01/3.0
#define EPS 1.0e-8


//GPU Block size
#define BSZ  32

//You should choose NX the power of 2.
#define NX  256
#define LATTICE_SIZE  8

//The number of the initial 
#define INITIAL_N  100

#define BATCH  10
#define M_PI  3.141592653589793238462643383
#define Nc 3
#define ADJNc 8
#define ALPHA_S	0.3
#define P_UPPER 5.0

#define NUMBER 256

#define MONTE_CARLO_NUMBER 500000

/**
* The parameter in Gauss function
*/
#define Gauss_param 2.0

/**
* The parameter \Lamda_QCD
*/
#define Lamda_QCD 1.0

/**
* The parameter \Lamda_QCD for confinement scale of the  nucleus
*/
//#define Lamda_QCD_nucleus  (Lamda_QCD / 1.0)
