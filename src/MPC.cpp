#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

size_t N_TIMESTEPS = 25;
double dt = 0.05;

const double Lf = 2.67; // length from front to center of gravity
int x_start = 0;
int y_start = N_TIMESTEPS;
int psi_start = 2*N_TIMESTEPS;
int vel_start = 3*N_TIMESTEPS;
int cte_start = 4*N_TIMESTEPS;
int psi_err_start = 5*N_TIMESTEPS;
int del_start = 6*N_TIMESTEPS;
int acc_start = 7*N_TIMESTEPS - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
  void operator()(ADvector& cost_vars, const ADvector& state_vars) {
    cost_vars[0] = 0; // cost value

    for (int t=0; t<N_TIMESTEPS; ++t) {
      cost_vars[0] += CppAD::pow(state_vars[vel_start + t], 2);
      cost_vars[0] += CppAD::pow(state_vars[cte_start + t], 2);
      cost_vars[0] += CppAD::pow(state_vars[psi_err_start + t], 2);
    }
    for (int t=0; t<N_TIMESTEPS - 1; ++t) {
      cost_vars[0] += CppAD::pow(state_vars[del_start + t], 2);
      cost_vars[0] += CppAD::pow(state_vars[acc_start + t], 2);
    }
    for (int t=0; t<N_TIMESTEPS - 2; ++t) {
      cost_vars[0] += CppAD::pow(state_vars[del_start + t + 1] - state_vars[del_start + t], 2);
      cost_vars[0] += CppAD::pow(state_vars[acc_start + t + 1] - state_vars[acc_start + t], 2);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6*N_TIMESTEPS + 2*(N_TIMESTEPS-1);
  // TODO: Set the number of constraints
  size_t n_constraints = 4*N_TIMESTEPS;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector state_vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    state_vars[i] = 0;
  }

  Dvector state_vars_lowerbound(n_vars);
  Dvector state_vars_upperbound(n_vars);
  for (int i = del_start; i < acc_start; i++) {
    state_vars_lowerbound[i] = deg2rad(-25);
    state_vars_upperbound[i] = deg2rad(25);
  }
  for (int i = acc_start; i < n_vars; i++) {
    state_vars_lowerbound[i] = 0;
    state_vars_upperbound[i] = -1;
  }

  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  std::vector<double> step_1;
  step_1.push_back(solution.x[x_start+1]);
  step_1.push_back(solution.x[y_start+1]);
  step_1.push_back(solution.x[psi_start+1]);
  step_1.push_back(solution.x[vel_start+1]);
  step_1.push_back(solution.x[cte_start+1]); 
  step_1.push_back(solution.x[psi_err_start]);
  step_1.push_back(solution.x[del_start]);
  step_1.push_back(solution.x[acc_start]);
  return step_1;
}
