#include "kalman_filter.h"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter(): window(0.2) {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &Rl_in, MatrixXd &Rr_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  Q_ = MatrixXd(4,4);
  H_ = H_in;
  Hj_ = MatrixXd(3,4);
  Rl_ = Rl_in;
  Rr_ = Rr_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_*x_;
  P_ = F_*P_*(F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  // perform the kalman filter measurement update
  VectorXd y = z - H_*x_;
  commonUpdate(P_, x_, H_, Rl_, y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float sq_x2y2 = sqrt(px*px + py*py);

  // calculate h(x)
  VectorXd h(3);
  h(0) = sq_x2y2;
  if (px != 0)
    h(1) = atan2(py,px);
  else if (py > 0)
    h(1) = M_PI/2.0;
  else if (py < 0)
    h(1) = -M_PI/2.0;
  else
    h(1) = 0;

  if (abs(px) < 0.0001 && abs(py) < 0.0001)
    h(2) = 0;
  else
    h(2) = (px*vx + py*vy)/sq_x2y2;


  // fix +-pi errors between measurement (z) and prediction (h)
  if (z(1) > M_PI - window && h(1) < 0 )
    h(1) += 2*M_PI;
  else if (z(1) < -M_PI + window && h(1) > 0 )
    h(1) -= 2*M_PI;

  // perform the kalman filter measurement update
  VectorXd y = z - h;
  commonUpdate(P_, x_, Hj_, Rr_, y);

}

void KalmanFilter::commonUpdate(MatrixXd &P, VectorXd &x, MatrixXd &H, MatrixXd &R, VectorXd &y) {

  MatrixXd Ht = H.transpose();
  MatrixXd I = MatrixXd::Identity(4,4);
  MatrixXd P_Ht = P*Ht;

  MatrixXd S = H*P_Ht + R;
  MatrixXd K = P_Ht*S.inverse();

  x = x + K*y;
  P = (I - K*H)*P;
}
