//
// Created by Bob Lutz on 12/24/17.
//

// ###### Sigma Point  Assignment

#include <iostream>
#include "Dense"
#include <vector>
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main() {

  //Create a UKF instance
  UKF ukf;

/*******************************************************************************
* Programming assignment calls
*******************************************************************************/

  MatrixXd Xsig = MatrixXd(5, 11);
  ukf.GenerateSigmaPoints(&Xsig);

  //print result
  std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  return 0;
}

#include <iostream>
#include "ukf.h"

UKF::UKF() {
  //TODO Auto-generated constructor stub
  Init();
}

UKF::~UKF() {
  //TODO Auto-generated destructor stub
}

void UKF::Init() {

}

/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/


void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //define spreading parameter
  double lambda = 3 - n_x;

  //set example state
  VectorXd x = VectorXd(n_x);
  x <<   5.7441,
          1.3800,
          2.2049,
          0.5015,
          0.3528;

  //set example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
          0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //your code goes here
  Xsig.col(0) = x; // first column is the current state vector

//   MatrixXd next = MatrixXd(5,11);
//   next << 0,0,0,0,0,0,0,0,0,0,0,
//   0,0,0,0,0,0,0,0,0,0,0,
//   0,0,0,0,0,0,0,0,0,0,0,
//   0,0,0,0,0,0,0,0,0,0,0,
//   0,0,0,0,0,0,0,0,0,0,0;

  // TODO rewrite as for loop
  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig
  Xsig.col(1) = x + (sqrt(lambda + n_x) * A.col(0));
  Xsig.col(2) = x + (sqrt(lambda + n_x) * A.col(1));
  Xsig.col(3) = x + (sqrt(lambda + n_x) * A.col(2));
  Xsig.col(4) = x + (sqrt(lambda + n_x) * A.col(3));
  Xsig.col(5) = x + (sqrt(lambda + n_x) * A.col(4));
  Xsig.col(6) = x - (sqrt(lambda + n_x) * A.col(0));
  Xsig.col(7) = x - (sqrt(lambda + n_x) * A.col(1));
  Xsig.col(8) = x - (sqrt(lambda + n_x) * A.col(2));
  Xsig.col(9) = x - (sqrt(lambda + n_x) * A.col(3));
  Xsig.col(10) = x - (sqrt(lambda + n_x) * A.col(4));

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  //write result
  *Xsig_out = Xsig;

/* expected result:
   Xsig = n
    5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441
      1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38
    2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049
    0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015
    0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879
*/

}