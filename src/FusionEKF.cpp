#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() : noise_ax2(9*9), noise_ay2(9*9)  {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // set H_laser matrix
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    VectorXd x_i(4);
	  MatrixXd F_i = MatrixXd::Identity(4,4);
	  MatrixXd P_i(4,4);
	  P_i << 1,0,0,0,
           0,1,0,0,
           0,0,1000,0,
           0,0,1000,0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	    float r = measurement_pack.raw_measurements_(0);
	    float phi = measurement_pack.raw_measurements_(1);

	    x_i(0) = r*cos(phi);
	    x_i(1) = r*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	    x_i(0) = measurement_pack.raw_measurements_(0);
	    x_i(1) = measurement_pack.raw_measurements_(1);
    }
	  // Init Kalman Filter
	  ekf_.Init(x_i,P_i,F_i,H_laser_,R_laser_,R_radar_);

	  //update time stamp
	  previous_timestamp_ = measurement_pack.timestamp_;

      // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // calculate elapsed time since last measurement
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;
  // update F Matrix
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // update Q Matrix
  ekf_.Q_ << dt_4/4*noise_ax2, 0, dt_3/2*noise_ax2, 0,
             0, dt_4/4*noise_ay2, 0, dt_3/2*noise_ay2,
             dt_3/2*noise_ax2, 0, dt_2*noise_ax2, 0,
             0, dt_3/2*noise_ay2, 0, dt_2*noise_ay2;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

	  Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.Hj_ = Hj_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else {
    //Laser updates
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
