#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
  	rmse << 0,0,0,0;

  	// check the validity of the following inputs:
  	//  * the estimation vector size should not be zero
  	if (estimations.size() < 1)
  	{
  	    cout << "Invalid est vec" << endl;
  	    return rmse;
  	}

  	//  * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size())
  	{
  	    cout << "Size mismatch" << endl;
  	    return rmse;
  	}

  	//accumulate squared residuals
  	VectorXd acc_err(4);
  	acc_err<<0,0,0,0;
  	for(int i=0; i < estimations.size(); ++i){
        VectorXd err  = (estimations[i]-ground_truth[i]);
        VectorXd se = err.array()*err.array();
  	    acc_err = acc_err+se;

  	}

  	//calculate the mean
  	acc_err /= estimations.size();

  	//calculate the squared root
    rmse = acc_err.array().sqrt();
  	return rmse;
}

double Tools::AngleWrap(double& angle)
{
  return atan2(sin(angle),cos(angle));
}
