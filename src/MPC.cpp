#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

const double pi() { return M_PI;}
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

size_t N_TIMESTEPS = 25;
double dt = 0.05;

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
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }
  double ref_v = 50;
  const double L_f = 2.67; // length from front to center of gravity

  typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
  void operator()(ADvector& cost_vars, const ADvector& state_vars) {
    cost_vars[0] = 0; // cost value

    for (int t=0; t<N_TIMESTEPS; ++t) {
      CppAD::AD<double> x_0 = state_vars[x_start + t - 1];
      CppAD::AD<double> x_1 = state_vars[x_start + t];
      CppAD::AD<double> y_0 = state_vars[y_start + t - 1];
      CppAD::AD<double> y_1 = state_vars[y_start + t];
      CppAD::AD<double> psi_0 = state_vars[psi_start + t - 1];
      CppAD::AD<double> psi_1 = state_vars[psi_start + t];
      CppAD::AD<double> vel_0 = state_vars[vel_start + t - 1];
      CppAD::AD<double> vel_1 = state_vars[vel_start + t];
      CppAD::AD<double> cte_0 = state_vars[cte_start + t - 1];
      CppAD::AD<double> cte_1 = state_vars[cte_start + t];
      CppAD::AD<double> psi_err_0 = state_vars[psi_err_start + t - 1];
      CppAD::AD<double> psi_err_1 = state_vars[psi_err_start + t];
      CppAD::AD<double> del_0 = state_vars[del_start + t - 1];
      CppAD::AD<double> del_1 = state_vars[del_start + t];
      CppAD::AD<double> del_2 = state_vars[del_start + t + 1];
      CppAD::AD<double> acc_0 = state_vars[acc_start + t - 1];
      CppAD::AD<double> acc_1 = state_vars[acc_start + t];
      CppAD::AD<double> acc_2 = state_vars[acc_start + t + 1];
      CppAD::AD<double> psi_dest = CppAD::atan(coeffs[1]);
      CppAD::AD<double> cte_x_0 = coeffs[0] + coeffs[1]*x_0;  
        
      cost_vars[0] += CppAD::pow(psi_1, 2);
      cost_vars[0] += CppAD::pow(vel_1 - ref_v, 2);
      cost_vars[0] += CppAD::pow(cte_1, 2);
      
      if (t > 0) {
        cost_vars[1 + x_start + t] = x_1 - x_0 - vel_0*dt*CppAD::cos(psi_0);
        cost_vars[1 + y_start + t] = y_1 - y_0 - vel_0*dt*CppAD::cos(psi_0);
        cost_vars[1 + psi_start + t] = psi_1 - psi_0 - vel_0*del_0*dt/L_f;
        cost_vars[1 + vel_start + t] = vel_1 - acc_0*dt;
        cost_vars[1 + cte_start + t] = psi_err_1 - cte_x_0 + y_0 - vel_0*dt*CppAD::sin(psi_err_0);
        cost_vars[1 + psi_err_start + t] = psi_err_1 - psi_0 + psi_dest -vel_0*del_0*dt/L_f; 
      }
      
      if (t < N_TIMESTEPS - 1){
        cost_vars[0] += CppAD::pow(del_1, 2);
        cost_vars[0] += CppAD::pow(acc_1, 2);
      }
      
      if (t < N_TIMESTEPS - 2) {
        cost_vars[0] += CppAD::pow(del_2 - del_1, 2);
        cost_vars[0] += CppAD::pow(acc_2 - acc_1, 2);
      }
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd cur_state, Eigen::VectorXd poly_coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  size_t n_vars = 6*N_TIMESTEPS + 2*(N_TIMESTEPS-1);
  size_t n_constraints = 6*N_TIMESTEPS;

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

void MPC::SetStateVariables(
  CppAD::vector<double> mpc_state, 
  Eigen::VectorXd cur_state) 
{
  for (int i = 0; i < mpc_state.size(); i++) {
    mpc_state[i] = 0;
  }
  mpc_state[x_start] = cur_state[0];
  mpc_state[y_start] = cur_state[1];
  mpc_state[psi_start] = cur_state[2];
  mpc_state[vel_start] = cur_state[3];
  mpc_state[cte_start] = cur_state[4];
  mpc_state[psi_err_start] = cur_state[5];
  mpc_state[del_start] = cur_state[6];
  mpc_state[acc_start] = cur_state[7];
}

void SetStateVariableBounds(
  CppAD::vector<double> lowerbound,
  CppAD::vector<double> upperbound) 
{
  if (lowerbound.size() != upperbound.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (int i = 0; i < del_start; i++) {
    lowerbound[i] = std::numeric_limits<double>::min();
    upperbound[i] = std::numeric_limits<double>::max();
  }
  for (int i = del_start; i < acc_start; i++) {
    lowerbound[i] = deg2rad(-25);
    upperbound[i] = deg2rad(25);
  }
  for (int i = acc_start; i < lowerbound.size(); i++) {
    lowerbound[i] = 0;
    upperbound[i] = 1;
  }
}

void MPC::SetControlBounds(
  CppAD::vector<double> lowerbound,
  CppAD::vector<double> upperbound,
  Eigen::VectorXd cur_state) 
{
  if (lowerbound.size() != upperbound.size()) {
    throw "\nLowerbound and upperbound size must match";
  }
  for (int i = 0; i < upperbound.size(); i++) {
    lowerbound[i] = 0;
    upperbound[i] = 0;
  }
  lowerbound[x_start] = cur_state[0];
  lowerbound[y_start] = cur_state[1];
  lowerbound[psi_start] = cur_state[2];
  lowerbound[vel_start] = cur_state[3];
  lowerbound[cte_start] = cur_state[4];
  lowerbound[psi_err_start] = cur_state[5];

  upperbound[x_start] = cur_state[0];
  upperbound[y_start] = cur_state[1];
  upperbound[psi_start] = cur_state[2];
  upperbound[vel_start] = cur_state[3];
  upperbound[cte_start] = cur_state[4];
  upperbound[psi_err_start] = cur_state[5];
}