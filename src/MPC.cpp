#include "MPC.h"
#include "FG_eval.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen/Core"

//
// MPC class definition implementation.
//
MPC::MPC() {
  N_TIMESTEPS_ = 10;
  dt_ = 0.1;
  L_f_ = 2.67;

  N_VARS_ = 6*N_TIMESTEPS_ + 2*(N_TIMESTEPS_-1);
  N_CONSTRAINTS_ = 6*N_TIMESTEPS_;
  state_vars_.resize(N_VARS_);
  state_vars_lowerbound_.resize(N_VARS_);
  state_vars_upperbound_.resize(N_VARS_);
  constraints_lowerbound_.resize(N_CONSTRAINTS_);
  constraints_upperbound_.resize(N_CONSTRAINTS_);

  x_start_ = 0;
  y_start_ = N_TIMESTEPS_;
  psi_start_ = 2*N_TIMESTEPS_;
  vel_start_ = 3*N_TIMESTEPS_;
  cte_start_ = 4*N_TIMESTEPS_;
  psi_err_start_ = 5*N_TIMESTEPS_;
  del_start_ = 6*N_TIMESTEPS_;
  acc_start_ = 7*N_TIMESTEPS_ - 1;
}
MPC::~MPC() {}

std::vector<double> MPC::Solve(Eigen::VectorXd cur_state, Eigen::VectorXd poly_coeffs) {
  bool ok = true;

  SetStateVariables(cur_state);
  SetStateVariableBounds();
  SetConstraintBounds(cur_state);

  // object that computes objective and constraints
  FG_eval fg_eval(poly_coeffs, this);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<CppAD::vector<double>> solution;

  // solve the problem
  CppAD::ipopt::solve<CppAD::vector<double>, FG_eval>(
      options, state_vars_, state_vars_lowerbound_, state_vars_upperbound_,
      constraints_lowerbound_, constraints_upperbound_, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<CppAD::vector<double> >::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  std::vector<double> result_vector;
  for (size_t i=0; i < solution.x.size(); ++i) {
    result_vector.push_back(solution.x[i]);
  }
  return result_vector;
}

void MPC::SetStateVariables(Eigen::VectorXd cur_state) 
{
  for (size_t i = 0; i < state_vars_.size(); i++) {
    state_vars_[i] = 0;
  }
  state_vars_[x_start_] = cur_state[0];
  state_vars_[y_start_] = cur_state[1];
  state_vars_[psi_start_] = cur_state[2];
  state_vars_[vel_start_] = cur_state[3];
  state_vars_[cte_start_] = cur_state[4];
  state_vars_[psi_err_start_] = cur_state[5];
}

void MPC::SetStateVariableBounds() 
{
  if (state_vars_lowerbound_.size() != state_vars_upperbound_.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (size_t i = 0; i < del_start_; i++) {
    state_vars_lowerbound_[i] = std::numeric_limits<double>::min();
    state_vars_upperbound_[i] = std::numeric_limits<double>::max();
  }
  for (size_t i = del_start_; i < acc_start_; i++) {
    state_vars_lowerbound_[i] = deg2rad(-25);
    state_vars_upperbound_[i] = deg2rad(25);
  }
  for (size_t i = acc_start_; i < state_vars_lowerbound_.size(); i++) {
    state_vars_lowerbound_[i] = 0.0;
    state_vars_upperbound_[i] = 0.5;
  }
}

void MPC::SetConstraintBounds(Eigen::VectorXd cur_state) 
{
  if (constraints_lowerbound_.size() != constraints_upperbound_.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (size_t i = 0; i < constraints_upperbound_.size(); i++) {
    constraints_lowerbound_[i] = 0;
    constraints_upperbound_[i] = 0;
  }
  constraints_lowerbound_[x_start_] = cur_state[0];
  constraints_lowerbound_[y_start_] = cur_state[1];
  constraints_lowerbound_[psi_start_] = cur_state[2];
  constraints_lowerbound_[vel_start_] = cur_state[3];
  constraints_lowerbound_[cte_start_] = cur_state[4];
  constraints_lowerbound_[psi_err_start_] = cur_state[5];

  constraints_upperbound_[x_start_] = cur_state[0];
  constraints_upperbound_[y_start_] = cur_state[1];
  constraints_upperbound_[psi_start_] = cur_state[2];
  constraints_upperbound_[vel_start_] = cur_state[3];
  constraints_upperbound_[cte_start_] = cur_state[4];
  constraints_upperbound_[psi_err_start_] = cur_state[5];
}
double MPC::deg2rad(double x) { return x * M_PI / 180; }
double MPC::rad2deg(double x) { return x * 180 / M_PI; }