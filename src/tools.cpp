#include <iostream>
#include "tools.h"

#define EPSILON 0.0001

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
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

MatrixXd Tools::SigmaPointPrediction(MatrixXd Xsig, int n_x, double delta_t) {

  // code from Sigma Point Prediction Assignment

  int n_aug = Xsig.rows();

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  for (int i = 0; i < 2 * n_aug + 1; i++) {

    VectorXd x = Xsig.col(i);
    double p_x = x[0];
    double p_y = x[1];
    double nu = x[2];
    double psi = x[3];
    double psi_dot = x[4];
    double nu_a = x[5];
    double nu_psi_ddot = x[6];

    //std::cout << "x = " << std::endl << x << std::endl;

    // process model
    VectorXd model = VectorXd(5);
    model.fill(0);
    model(3) = psi_dot * delta_t;

    if (abs(psi_dot) > EPSILON) { //avoid division by zero
      model(0) = nu / psi_dot * ( sin(psi + psi_dot * delta_t) - sin(psi));
      model(1) = nu / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi));
    } else {
      model(0) = nu * cos(psi) * delta_t;
      model(1) = nu * sin(psi) * delta_t;
    }

    // process noise
    VectorXd noise = VectorXd(5);
    noise.fill(0);
    noise(0) = 0.5 * delta_t * delta_t * cos(psi) * nu_a;
    noise(1) = 0.5 * delta_t * delta_t * sin(psi) * nu_a;
    noise(2) = delta_t * nu_a;
    noise(3) = 0.5 * delta_t * delta_t * nu_psi_ddot;
    noise(4) = delta_t * nu_psi_ddot;

    // overall, write predicted sigma points into approriate column
    for (int j = 0; j < n_x; j++) {
      Xsig_pred(j, i) = x(j) + model(j) + noise(j);
    }
  }

  return Xsig_pred;

}

MatrixXd Tools::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda) {

  int n = x.size();

  MatrixXd Xsig = MatrixXd( n, 15);     // sigma point matrix
  MatrixXd A = P.llt().matrixL();       // calculate square root of P

  Xsig.col(0) = x;

  double radical_term = sqrt(lambda + n);
  for (int i = 0; i < n; i++){
    Xsig.col( i + 1 ) = x + radical_term * A.col(i);
    Xsig.col( i + 1 + n ) = x - radical_term * A.col(i);
  }
  return Xsig;

}

pair<VectorXd, MatrixXd> Tools::AugmentSigmaPoints(VectorXd x, MatrixXd P,
                                                   int n_incr, double lambda,
                                                   double std_a, double std_yawdd) {

// adapted from Augmentation Assignment

  int n_aug = x.rows() + n_incr;
  VectorXd x_aug = VectorXd(n_aug);                    //create augmented mean vector
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);             //create augmented state covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);  //create sigma point matrix

  x_aug.head(5) = x;                                   //create augmented mean state
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.setZero(7, 7);                                 //create augmented covariance matrix
  P_aug.topLeftCorner(5, 5) = P;
  P_aug(5, 5) = std_a * std_a;
  P_aug(6, 6) = std_yawdd * std_yawdd;

  auto result = make_pair(x_aug, P_aug);
  return result;

}

double Tools::ConstrainAngle(double x) {
//  const double TWO_M_PI = 2.0 * M_PI;
//  x = fmod(x + M_PI, TWO_M_PI);
//  if (x < 0) x += TWO_M_PI;
//  return x - M_PI;

  double alpha = fmod(x, 2.*M_PI);
  if (alpha > M_PI) alpha = alpha - 2.*M_PI;
  return alpha;


}




