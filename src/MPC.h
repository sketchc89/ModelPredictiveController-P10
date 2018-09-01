#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  void SetStateVariables(CppAD::vector<double> mpc_state, Eigen::VectorXd cur_state);
  void SetStateVariableBounds(CppAD::vector<double> lowerbound, CppAD::vector<double> upperbound);
  void SetControlBounds(CppAD::vector<double> lowerbound, CppAD::vector<double> upperbound,
    Eigen::VectorXd cur_state);
};

#endif /* MPC_H */
