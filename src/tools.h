#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    double ConstrainAngle(double x);
    MatrixXd GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda);
    MatrixXd SigmaPointPrediction(MatrixXd Xsig, int n_x, double delta_t);
    pair<VectorXd, MatrixXd> AugmentSigmaPoints(VectorXd x, MatrixXd P,
                                               int n_incr, double lambda,
                                               double std_a, double std_yawdd);

};

#endif /* TOOLS_H_ */