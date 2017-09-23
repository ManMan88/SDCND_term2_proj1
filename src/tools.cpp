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
	// if both px and py are small, replace them with agreed small values
    cout << "NOTE - px and py are both small! Returning px = 0.001, py = 0.001" << endl;
    px = 0.001;
    py = 0.001;
  }

  //compute the Jacobian matrix
  float px2py2 = px*px+py*py;
  float s_px2py2 = sqrt(px2py2);
  float px2py2_32 = px2py2*s_px2py2;
  float vp = vx*py-vy*px;

  Hj << px/s_px2py2, py/s_px2py2, 0, 0,
        -py/px2py2, px/px2py2, 0, 0,
        py*vp/px2py2_32, -px*vp/px2py2_32, px/s_px2py2, py/s_px2py2;
  return Hj;
}
