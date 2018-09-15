#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen/Core"
#include <cppad/cppad.hpp>

using namespace std;

class MPC {
 public:
  MPC();

    virtual ~MPC();

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuatotions.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
    void SetStateVariables(Eigen::VectorXd cur_state);
    void SetStateVariableBounds();
    void SetConstraintBounds(Eigen::VectorXd cur_state);
    double deg2rad(double x);
    double rad2deg(double x);
    size_t N_TIMESTEPS_;
    double dt_;
    size_t x_start_;
    size_t y_start_;
    size_t psi_start_;
    size_t vel_start_;
    size_t cte_start_;
    size_t psi_err_start_;
    size_t del_start_;
    size_t acc_start_;
  private:
    size_t N_VARS_;
    size_t N_CONSTRAINTS_;
    CppAD::vector<double> state_vars_;
    CppAD::vector<double> state_vars_lowerbound_;
    CppAD::vector<double> state_vars_upperbound_;
    CppAD::vector<double> constraints_lowerbound_;
    CppAD::vector<double> constraints_upperbound_;
};

#endif /* MPC_H */
