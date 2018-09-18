#include "MPC.h"
#include "FG_eval.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen/Core"

//
// MPC class definition implementation.
//
MPC::MPC()
{
    N_TIMESTEPS_ = 10;
    dt_ = 0.1;
    L_f_ = 2.67;

    N_VARS_ = 2 * (N_TIMESTEPS_ - 1);
    N_CONSTRAINTS_ = 0;

    del_start_ = 0;
    acc_start_ = N_TIMESTEPS_ - 1;
}
MPC::~MPC() {}

std::vector<double> MPC::Solve(Eigen::VectorXd cur_state, Eigen::VectorXd poly_coeffs)
{
    bool ok = true;
    Dvector state_vars(N_VARS_);
    SetStateVariables(state_vars);

    Dvector state_vars_lowerbound(N_VARS_);
    Dvector state_vars_upperbound(N_VARS_);
    SetStateVariableBounds(state_vars_lowerbound, state_vars_upperbound);

    Dvector constraints_lowerbound(N_CONSTRAINTS_);
    Dvector constraints_upperbound(N_CONSTRAINTS_);
    SetConstraintBounds(constraints_lowerbound, constraints_upperbound);

    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<CppAD::vector<double>> solution;

    ADvector control(N_VARS_);

    // object that computes objective and constraints
    FG_eval fg_eval(cur_state, poly_coeffs, this);

    // solve the problem
    CppAD::ipopt::solve<CppAD::vector<double>, FG_eval>(
        options, state_vars, state_vars_lowerbound, state_vars_upperbound,
        constraints_lowerbound, constraints_upperbound, fg_eval, solution);
    for (size_t i = 0; i < N_VARS_; ++i)
    {
        control[i] = solution.x[i];
    }
    fg_eval.KinematicModel(control);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<CppAD::vector<double>>::success;

    // Cost
    auto cost = solution.obj_value;
    // std::cout << "Cost " << cost << std::endl;

    std::vector<double> result_vector;
    for (size_t i = 0; i < N_TIMESTEPS_; ++i)
    {
        result_vector.push_back(0.0);//CppAD::Value(fg_eval.x_[i]));
    }
    for (size_t i = N_TIMESTEPS_; i < 2*N_TIMESTEPS_; ++i)
    {
        result_vector.push_back(0.0);//CppAD::Value(fg_eval.y_[i]));
    }
    // for (size_t i = 0; i < N_VARS_; ++i) {
    //     std:: cout << solution.x[i] << "\t";
    // }
    std::cout << "\n";
    result_vector.push_back(solution.x[del_start_]);
    result_vector.push_back(solution.x[acc_start_]);
    return result_vector;
}



void MPC::SetStateVariables(Dvector &state_vars)
{
    for (size_t i = 0; i < state_vars.size(); i++)
    {
        state_vars[i] = 0;
    }
}
void MPC::SetStateVariableBounds(Dvector &state_vars_lowerbound, Dvector &state_vars_upperbound)
{
    if (state_vars_lowerbound.size() != state_vars_upperbound.size())
    {
        throw "\nLowerbound and upperbound size must match";
    }
    for (size_t i = del_start_; i < acc_start_; i++)
    {
        state_vars_lowerbound[i] = deg2rad(-25);
        state_vars_upperbound[i] = deg2rad(25);
    }
    for (size_t i = acc_start_; i < state_vars_lowerbound.size(); i++)
    {
        state_vars_lowerbound[i] = 0.0;
        state_vars_upperbound[i] = 1.0;
    }
}

void MPC::SetConstraintBounds(Dvector constraints_lowerbound, Dvector constraints_upperbound)
{
    if (constraints_lowerbound.size() != constraints_upperbound.size())
    {
        throw "\nLowerbound and upperbound size must match";
    }
    for (size_t i = 0; i < constraints_upperbound.size(); i++)
    {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
}
double MPC::deg2rad(double x) { return x * M_PI / 180; }
double MPC::rad2deg(double x) { return x * 180 / M_PI; }