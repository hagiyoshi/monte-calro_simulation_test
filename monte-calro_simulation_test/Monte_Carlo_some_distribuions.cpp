#include <random>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "Parameters.h"

//You should choose NX the power of 2.
#define NX  256
#define LATTICE_SIZE  8

void calculation_matrix(double Mat[2][2])
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			Mat[i][j] += 2;
		}
	}
}

/**
* Return the matrix exponential.
*/
inline Eigen::Matrix3cd exp_U(const Eigen::Matrix3cd m) {
	// Scaling and squaring + Taylor expansion.
	const double eps = std::numeric_limits<double>::epsilon() * 5;
	const int k = 1;
	// Find s such that |m/2^s| = |m|/2^s < 1/2^k.
	double norm = std::sqrt(m.squaredNorm());
	int s = std::max((int)(std::log(norm) / std::log(2.0)) + (k + 1), 0);
	// Scaling.
	const double scale = std::pow((double)2, -s);
	Eigen::Matrix3cd a;
	a = m;
	a = a.array()*scale;
	// Taylor expansion to get exp(m/2^s).
	Eigen::Matrix3cd sum = Eigen::MatrixXd::Identity(3,3);
	Eigen::Matrix3cd x=a;
	sum += x;
	for (int i = 2;; i++) {
		Eigen::Matrix3cd old=sum;
		x = x*a;
		double factor = 1.0/i;
		x = x.array() * factor;
		sum += x;
		Eigen::Matrix3cd osubs = old - sum;
		double osubs_norm2 = osubs.squaredNorm();
		if (osubs_norm2 < eps * eps) {
			break;
		}
	}
	// Squaring to get exp(m) = [exp(m/2^s)]^(2^s).
	for (int i = 0; i < s; i++) {
		sum *= sum;
	}
	return sum;
}

void noise_generation(double* noise)
{
	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());

#pragma omp parallel for num_threads(6)
	for (int j = 0; j < NX; j++) {
		for (int i = 0; i < NX; i++)
		{

			std::normal_distribution<double> dist(0.0, 1.0 );
			noise[NX*j + i] = dist(engine);
		}
	}
}

int main()
{
	double Mat[2][2];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			Mat[i][j] = 0;
		}
	}

	calculation_matrix(Mat);
	Eigen::Matrix3cd A,A_exp,A_exp_U;
	A << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, Pi*1.0), std::complex<double>(0.0, 0.0),
		std::complex<double>(0.0, Pi*1.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
		std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0);
	for (int i = 0; i < 1000; ++i) {
		A_exp = A.exp();
		A_exp_U = exp_U(A);
	}

	std::vector<double> noise_testa(NX*NX, 0);

	double noissss = 0;
	double noi = 0;
	noise_generation(noise_testa.data());

	for (int i = 0; i < NX*NX; ++i) {
		noissss += noise_testa[i] * noise_testa[i]/((double)NX*NX);
		noi += noise_testa[i] / ((double)NX*NX);
	}

	std::cout << noissss <<"\t"<<noi <<"\n" ;

	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());

	// 平均mu = 0.0、標準偏差sigma = 1.0で分布させる (1/(sqrt(2 pi sigma^2)))*exp( - (x - mu)^2/(2 sigma^2))
	std::normal_distribution<> dist(0.0, 1.0);
	int N = 10;
	double expectation_value = 0;
	double expectpatronam = 0;
	double expectpatronam2 = 0;
	for (int i = 0; i < N; ++i) {
		expectation_value += dist(engine)*dist(engine)/N;
		double gauss = dist(engine);
		double gauss2 = dist(engine);
		expectpatronam += gauss*gauss/N;
		expectpatronam2 += gauss*gauss2 / N;
	}

	double Radius_proton = 2.0;
	std::normal_distribution<double> dist2(0.0, Radius_proton / sqrt(2.0));
	double test_expectation = 0;
	for (int k = 0; k < N; k++) {
		double value = dist2(engine);
		test_expectation += value*value / N;
	}


	std::vector<double> test1(2, 0);
	for (int k = 0; k < N; k++) {

		std::random_device seed_gen2;
		std::default_random_engine engine2(seed_gen2());

		std::vector<double> v(2,0);

		std::generate(v.begin(), v.end(),
			[dist2,engine2]()mutable->double {
			return dist2(engine2);
		}
		);

		double value = v[0];
		double value2 = v[1];
		test1[0] += value*value / N;
		test1[1] += value*value2 / N;
	}

	std::ofstream file("normal_distribution.tsv");
	for (std::size_t n = 0; n < 10 * 10; ++n) {
		// 正規分布で乱数を生成する
		double result = dist(engine);
		file << result << "\t\n";
	}
}
