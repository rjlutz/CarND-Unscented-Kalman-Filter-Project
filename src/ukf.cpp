#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  use_laser_ = true;   // laser measurements will be ignored (except during init)
  use_radar_ = true;   // radar measurements will be ignored (except during init)
  x_ = VectorXd(5);    // initial state vector ->
                       // [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  P_ = MatrixXd(5, 5); // initial covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_lidar_ = MatrixXd(2, 2);
  std_laspx_ = 0.15;   // Laser measurement noise standard deviation position1 in m
  std_laspy_ = 0.15;   // Laser measurement noise standard deviation position2 in m
  std_radr_ = 0.3;     // Radar measurement noise standard deviation radius in m
  std_radphi_ = 0.03;  // Radar measurement noise standard deviation angle in rad
  std_radrd_ = 0.3;    // Radar measurement noise standard deviation radius change in m/s

  is_initialized_ = false; // initially set to false, set to true in first call of ProcessMeasurement
  std_a_ = 1.5;                             // Process noise std longitudinal acceleration in m/s^2
  std_yawdd_ = 1;                           // Process noise std deviation yaw acceleration in rad/s^2

  time_us_ = 0;                             // time when the state is true, in milliseconds
  n_x_ = 5;                                 // State dimension
  n_aug_ = 7;                               // Augmented state dimension
  n_sigma_pts = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_pts); // predicted sigma points matrix
  lambda_ = 3 - n_aug_;                     // Sigma point spreading parameter
  tools = Tools();

  weights_ = VectorXd(n_sigma_pts);
  weights_(0) = lambda_/(lambda_+n_aug_);   // Weights of sigma points
  for (int i=1; i<n_sigma_pts; i++)         // Weights of sigma points
    weights_(i) = 0.5/(n_aug_+lambda_);

  P_ <<   1, 0, 0, 0, 0,                    // state covariance matrix
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          1, 0, 0, 0, 1;

  R_radar_ << std_radr_ * std_radr_, 0,                         0,
              0,                     std_radphi_ * std_radphi_, 0,
              0,                     0,                         std_radrd_ * std_radrd_;


  R_lidar_ << std_laspx_*std_laspx_, 0,
              0,                     std_laspy_*std_laspy_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double x = rho * cos(theta);
      double y = rho * sin(theta);
      double v = sqrt(x * x + y * y);
      x_ << x, y, v, 0, 0; //polar to cartesian coords
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);


  // Update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    UpdateRadar(meas_package);
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    UpdateLidar(meas_package);

}

/**
 * Predicts sigma points, the state (x_), and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  pair<VectorXd, MatrixXd> result =
          tools.AugmentSigmaPoints(x_, P_, n_aug_-n_x_, lambda_, std_a_, std_yawdd_);
  VectorXd x_aug = result.first;
  MatrixXd P_aug = result.second;

  MatrixXd Xsig_aug = tools.GenerateSigmaPoints(x_aug, P_aug, lambda_);
  Xsig_pred_ = tools.SigmaPointPrediction(Xsig_aug, n_x_, delta_t);

  x_ = Xsig_pred_ * weights_;                                  // predict the state mean
  P_.fill(0.0);                                                // predicted state covariance matrix
  for (int i = 0; i < n_sigma_pts; i++) {                      // iterate over sigma points
    VectorXd residual = Xsig_pred_.col(i) - x_;                // state difference
    residual(3) = tools.ConstrainAngle(residual(3));           // angle normalization
    P_ += weights_(i) * residual * residual.transpose();
  }

}

/**
 * Updates the state (the belief about the object's
 * position) and the state covariance matrix using a laser measurement.
 *
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  int n_z = 2;
  MatrixXd Zsig = Xsig_pred_.topRows(n_z);

  VectorXd z_pred = VectorXd(n_z);                  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_pts; i++)
    z_pred +=  weights_(i) * Zsig.col(i);


  MatrixXd S = MatrixXd(n_z,n_z);                   // measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sigma_pts; i++) {           // iterate over sigma points
    VectorXd z_residual = Zsig.col(i) - z_pred;
    S += weights_(i) * z_residual * z_residual.transpose();
  }
  S += R_lidar_;                                    //innovation measurement covariance matrix S

  VectorXd z = meas_package.raw_measurements_;      // incoming lidar measurement

  MatrixXd Tc = MatrixXd(n_x_, n_z);                //create matrix for cross correlation Tc
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_pts; i++) {           // iterate over sigma points
    VectorXd z_residual = Zsig.col(i) - z_pred;     //residual
    VectorXd x_residual = Xsig_pred_.col(i) - x_;   // state difference
    Tc += weights_(i) * x_residual * z_residual.transpose();
  }

  MatrixXd K = Tc * S.inverse();                    //Kalman gain K
  VectorXd x_residual = z - z_pred;                 //residual

  x_ += K * x_residual;                             //update state mean and covariance matrix
  P_ -= K*S*K.transpose();
  NIS_laser_ = x_residual.transpose() * S.inverse() * x_residual;  // NIS Update

}

/**
 * Updates the state (the belief about the object's
 * position) and the state covariance matrix using a radar measurement.
 *
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {


  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, n_sigma_pts);   // Matrix for sigma points, in measurement space
  for (int i = 0; i < n_sigma_pts; i++) {       // Transform sigma points into measurement space
    VectorXd x = Xsig_pred_.col(i);
    double p_x = x(0);
    double p_y = x(1);
    double nu = x(2);
    double psi = x(3);

    Zsig.col(i) << sqrt(p_x*p_x + p_y*p_y),                                      // rho
            atan2(p_y, p_x),                                                     // phi
            (nu * (p_x * cos (psi) + p_y * sin(psi))) / sqrt(p_x*p_x + p_y*p_y); // rho_dot
  }

  VectorXd z_pred = VectorXd(n_z);                          //create vector for mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < n_sigma_pts; i++)
    z_pred += weights_(i) * Zsig.col(i);

  MatrixXd S = MatrixXd(n_z,n_z);                           // innovation measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sigma_pts; i++) {                   // iterate over sigma points
    VectorXd A = Zsig.col(i) - z_pred;
    A(1) = tools.ConstrainAngle(A(1));                      // angle normalization
    S+= weights_(i) * A * A.transpose();
  }
  S += R_radar_; // add on the measurement error

  VectorXd z = meas_package.raw_measurements_;              // incoming radar measurement
                                                            // rho in m, phi in radians, rho_dot in m/s

  MatrixXd Tc = MatrixXd(n_x_, n_z);                        // matrix for cross correlation Tc
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_pts; i++) {
    VectorXd x_residual = Xsig_pred_.col(i) - x_;           // x residual, state difference
    VectorXd z_residual = Zsig.col(i) - z_pred;             // z residual
    x_residual(3) = tools.ConstrainAngle(x_residual(3));    // angle normalization x
    z_residual(1) = tools.ConstrainAngle(z_residual(1));    // angle normalization z
    Tc += weights_(i) * x_residual * z_residual.transpose();
  }

  MatrixXd K = Tc * S.inverse();                            // Kalman gain K
  VectorXd z_residual = z - z_pred;

  x_ += K * z_residual;                                     // update state mean and covariance matrix
  P_ -= K * S * K.transpose();
  NIS_radar_ = z_residual.transpose() * S.inverse() * z_residual;  //NIS Update

}







