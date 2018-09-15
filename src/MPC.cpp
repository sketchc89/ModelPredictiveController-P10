#include "MPC.h"
#include "FG_eval.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#ifndef HELPER_H
  #include "Helper.h"
#endif

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd cur_state, Eigen::VectorXd poly_coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  Helper helper;
  size_t n_vars = 6*helper.N_TIMESTEPS + 2*(helper.N_TIMESTEPS-1);
  size_t n_constraints = 6*helper.N_TIMESTEPS;

  Dvector state_vars(n_vars);
  SetStateVariables(state_vars, cur_state);
  
  Dvector state_vars_lowerbound(n_vars);
  Dvector state_vars_upperbound(n_vars);
  SetStateVariableBounds(state_vars_lowerbound, state_vars_upperbound);

  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  SetControlBounds(constraints_lowerbound, constraints_upperbound, cur_state);

  // object that computes objective and constraints
  FG_eval fg_eval(poly_coeffs);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, state_vars, state_vars_lowerbound, state_vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  std::vector<double> step_1;
  step_1.push_back(solution.x[helper.x_start+1]);
  step_1.push_back(solution.x[helper.y_start+1]);
  step_1.push_back(solution.x[helper.psi_start+1]);
  step_1.push_back(solution.x[helper.vel_start+1]);
  step_1.push_back(solution.x[helper.cte_start+1]); 
  step_1.push_back(solution.x[helper.psi_err_start]);
  step_1.push_back(solution.x[helper.del_start]);
  step_1.push_back(solution.x[helper.acc_start]);
  return step_1;
}

void MPC::SetStateVariables(
  CppAD::vector<double> mpc_state, 
  Eigen::VectorXd cur_state) 
{
  Helper helper;
  for (size_t i = 0; i < mpc_state.size(); i++) {
    mpc_state[i] = 0;
  }
  mpc_state[helper.x_start] = cur_state[0];
  mpc_state[helper.y_start] = cur_state[1];
  mpc_state[helper.psi_start] = cur_state[2];
  mpc_state[helper.vel_start] = cur_state[3];
  mpc_state[helper.cte_start] = cur_state[4];
  mpc_state[helper.psi_err_start] = cur_state[5];
  mpc_state[helper.del_start] = cur_state[6];
  mpc_state[helper.acc_start] = cur_state[7];
}

void MPC::SetStateVariableBounds(
  CppAD::vector<double> lowerbound,
  CppAD::vector<double> upperbound) 
{
  Helper helper;
  if (lowerbound.size() != upperbound.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (size_t i = 0; i < helper.del_start; i++) {
    lowerbound[i] = std::numeric_limits<double>::min();
    upperbound[i] = std::numeric_limits<double>::max();
  }
  for (size_t i = helper.del_start; i < helper.acc_start; i++) {
    lowerbound[i] = helper.deg2rad(-25);
    upperbound[i] = helper.deg2rad(25);
  }
  for (size_t i = helper.acc_start; i < lowerbound.size(); i++) {
    lowerbound[i] = 0;
    upperbound[i] = 1;
  }
}

void MPC::SetControlBounds(
  CppAD::vector<double> lowerbound,
  CppAD::vector<double> upperbound,
  Eigen::VectorXd cur_state) 
{
  Helper helper;
  if (lowerbound.size() != upperbound.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (size_t i = 0; i < upperbound.size(); i++) {
    lowerbound[i] = 0;
    upperbound[i] = 0;
  }
  lowerbound[helper.x_start] = cur_state[0];
  lowerbound[helper.y_start] = cur_state[1];
  lowerbound[helper.psi_start] = cur_state[2];
  lowerbound[helper.vel_start] = cur_state[3];
  lowerbound[helper.cte_start] = cur_state[4];
  lowerbound[helper.psi_err_start] = cur_state[5];

  upperbound[helper.x_start] = cur_state[0];
  upperbound[helper.y_start] = cur_state[1];
  upperbound[helper.psi_start] = cur_state[2];
  upperbound[helper.vel_start] = cur_state[3];
  upperbound[helper.cte_start] = cur_state[4];
  upperbound[helper.psi_err_start] = cur_state[5];
}