//one should use row = column.


#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
*This class aims to interpolate 2-dimentional function. 
*/
class interpolation_2dim
{
private:
	Eigen::MatrixXd M_func;
	vector<double> vec_x, vec_y;
	double x_max,x_min, y_max,y_min;

public:
	void set_points(const Eigen::MatrixXd& A,
		const std::vector<double>& y1, const std::vector<double>& y2);
	double operator() (double x,double y) const;
};


void interpolation_2dim::set_points(const Eigen::MatrixXd& A,
	const std::vector<double>& y1, const std::vector<double>& y2)
{
	int Acols = A.cols();
	int y1size = y1.size();
	int Arows = A.rows();
	int y2size = y2.size();
	//assert(Acols == y1size);
	//assert(Arows == y2size);
	//assert(A.cols() == y1.size())
	assert(A.cols()>1);
	assert(A.rows()>1);
	M_func = A;
	vec_x = y1;
	vec_y = y2;

	x_max = vec_x[y1.size() -1];
	y_max = vec_y[y2.size() - 1];
	x_min = vec_x[0];
	y_min = vec_y[0];

}

double interpolation_2dim::operator() (double x, double y) const
{
	// find the closest point vec_x[idx] < x, idx=0 even if x<vec_x[0]
	std::vector<double>::const_iterator it;
	it = std::lower_bound(vec_x.begin(), vec_x.end(), x);
	int idx = std::max(int(it - vec_x.begin()) - 1, 0);

	// find the closest point vec_y[idy] < y, idy=0 even if y<vec_y[0]
	std::vector<double>::const_iterator ity;
	ity = std::lower_bound(vec_y.begin(), vec_y.end(), y);
	int idy = std::max(int(ity - vec_y.begin()) - 1, 0);

	int it_d = int(it - vec_x.begin());
	int ity_d = int(ity - vec_y.begin());

	//cout <<"y_max "<< y_max <<" x "<<x <<" y " <<y << " it "<< it_d << " idx " << idx << " ity " << ity_d << " idy " << idy<< endl;

	double t, u;

	if (idx+1 <= vec_x.size() - 1 && idy+1  <= vec_y.size() - 1){
		//bilinear interpolation
		 t = (x - vec_x[idx]) / (vec_x[idx + 1] - vec_x[idx]);
		 u = (y - vec_y[idy]) / (vec_y[idy + 1] - vec_y[idy]);
	}
	else{

		//bilinear interpolation
		t = (x - vec_x[vec_x.size() - 1 - 1]) / (vec_x[vec_x.size() - 1] - vec_x[vec_x.size() - 1 - 1]);
		u = (y - vec_y[vec_y.size() - 1 - 1]) / (vec_y[vec_y.size() - 1] - vec_y[vec_y.size() - 1 - 1]);
	}

	double interpole;


	if (x > x_max || y > y_max || x < x_min || y < y_min)
	{
		cout << "out of range" << endl;
		cout << "x " << x << " y " << y << endl;
		interpole = 0;
		//assert(1);
	}
	else
	{
		//f(x,y) = (1-t)(1-u)f(x_i,y_k)+t(1-u)f(x_i+1,y_k)+tuf(x_i+1,y_k+1)+(1-t)uf(x_i,y_k+1)
		interpole = (1.0 - t)*(1.0 - u)*M_func(idx, idy) + t*(1.0 - u)*M_func(idx + 1, idy) 
			+ t*u*M_func(idx + 1, idy + 1) + (1.0 - t)*u*M_func(idx, idy + 1);
	}

	return interpole;
}