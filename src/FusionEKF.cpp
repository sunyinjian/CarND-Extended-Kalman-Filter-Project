#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    /**
    TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
    */
    MatrixXd F = MatrixXd::Identity(4, 4);

    VectorXd x = VectorXd::Zero(4);

    MatrixXd P = MatrixXd(4, 4);
    P << 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1000.0, 0.0,
            0.0, 0.0, 0.0, 1000.0;

    MatrixXd Q = MatrixXd::Zero(4, 4);

    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    noise_ax_ = 9.0;
    noise_ay_ = 9.0;

    ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
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
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        const float ro = measurement_pack.raw_measurements_[0];
        const float theta = measurement_pack.raw_measurements_[1];
        const float ro_dot = measurement_pack.raw_measurements_[2];

        ekf_.x_ << ro*cos(theta), ro*sin(theta), 0.0, 0.0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;
    }

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

    const double delta_T = (measurement_pack.timestamp_ - previous_timestamp_) * 1e-6;
    ekf_.F_(0, 2) = delta_T;
    ekf_.F_(1, 3) = delta_T;

    double t1 = delta_T;
    double t2 = t1 * t1;
    double t3 = t2 * t1;
    double t4 = t3 * t1;

    ekf_.Q_ << t4*noise_ax_/4.0, 0.0, t3*noise_ax_/2.0, 0.0,
            0.0, t4*noise_ay_/4.0, 0.0, t3*noise_ay_/2.0,
            t3*noise_ax_/2.0, 0.0, t2*noise_ax_, 0.0,
            0.0, t3*noise_ay_/2.0, 0.0, t2*noise_ay_;

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
        ekf_.R_ = R_radar_;
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

        // Radar updates
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;

        // Laser updates
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
