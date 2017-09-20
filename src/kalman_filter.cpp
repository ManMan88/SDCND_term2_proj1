#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht = H_.transpose;
  MatrixXd I = MatrixXd::Identity(3,3);

  // perform the kalman filter measurement update
  VectorXd y = z - H_*x_;
  MatrixXd S = H_*P*Ht + R_;
  MatrixXd K = P*Ht*S.inverse();

  x = x + K*y;
  P = (I - K*H)*P;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  MatrixXd Ht = H_.transpose;
  MatrixXd I = MatrixXd::Identity(3,3);
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
  // perform the kalman filter measurement update
  VectorXd y = z - h;
  MatrixXd S = H_*P*Ht + R_;
  MatrixXd K = P*Ht*S.inverse();

  x = x + K*y;
  P = (I - K*H)*P;
}
