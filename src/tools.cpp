#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd residual(4);
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  for(uint i=0; i < estimations.size(); ++i){
    VectorXd err = estimations[i] - ground_truth[i];
    residual = err.array()*err.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();
  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  // initialize the Jacobian Matrix
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if (abs(px) < 0.001 && abs(py) < 0.001) {
	// if both px and py are 0, return a 0 matrix
	// in such a case, the state estimation will be the prediction, and the P (uncertainty) will stay the same
    cout << "NOTE - px and py are both 0! Returning a 0 matrix" << endl;
	Hj = MatrixXd(3,4);
  }
  else {
    //compute the Jacobian matrix
    float px2py2 = px*px+py*py;
	  float s_px2py2 = sqrt(px2py2);
	  float px2py2_32 = px2py2*s_px2py2;
	  float vp = vx*py-vy*px;

    Hj << px/s_px2py2, py/s_px2py2, 0, 0,
          -py/px2py2, px/px2py2, 0, 0,
          py*vp/px2py2_32, -px*vp/px2py2_32, px/s_px2py2, py/s_px2py2;
  }
  return Hj;
}
