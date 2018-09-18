#include "FG_eval.h"
#include "MPC.h"
#include <cppad/cppad.hpp>
#include "Eigen/Core"

FG_eval::FG_eval(Dvector &vel, Dvector &cte, Dvector &psi_err, MPC *mpc) { 
    this->vel_ = vel;
    this->cte_ = cte;
    this->psi_err_ = psi_err;
    this->N_TIMESTEPS_ = mpc->N_TIMESTEPS_;
    this->dt_ = mpc->dt_;
    this->L_f_ = mpc->L_f_;
    this->del_start_ = mpc->del_start_;
    this->acc_start_ = mpc->acc_start_;
    this->ref_v_ = 50;
    this->L_f_ = 2.67;
}
FG_eval::~FG_eval() {}

void FG_eval::operator()(ADvector &cost_vars, const ADvector &state_vars) {
    cost_vars[0] = 0; // cost value

    int vel_weight = 10;
    int cte_weight = 1;
    int psi_err_weight = 1;
    int del_weight = 1;
    int acc_weight = 0;
    int del_change_weight = 1;
    int acc_change_weight = 0;
    CppAD::AD<double> vel_cost = 0.0;
    CppAD::AD<double> cte_cost = 0.0;
    CppAD::AD<double> psi_err_cost = 0.0;
    CppAD::AD<double> del_cost = 0.0;
    CppAD::AD<double> acc_cost = 0.0;
    CppAD::AD<double> del_change_cost = 0.0;
    CppAD::AD<double> acc_change_cost = 0.0;

    for (size_t t=0; t < N_TIMESTEPS_; ++t) {
        cost_vars[0] += vel_weight * CppAD::pow(vel_[t] - ref_v_, 2);
        cost_vars[0] += cte_weight * CppAD::pow(cte_[t], 2);
        cost_vars[0] += psi_err_weight * CppAD::pow(psi_err_[t], 2);
    }
    for (size_t t=0; t < N_TIMESTEPS_ - 1; ++t) {
        cost_vars[0] += del_weight * CppAD::pow(state_vars[del_start_ + t], 2);
        cost_vars[0] += acc_weight * CppAD::pow(state_vars[acc_start_ + t], 2);
    }
    for (size_t t=0; t < N_TIMESTEPS_ - 2; ++t) {
        cost_vars[0] += del_change_weight * CppAD::pow(state_vars[del_start_ + t + 1] - 
                                                        state_vars[del_start_ + t], 2);
        cost_vars[0] += acc_change_weight * CppAD::pow(state_vars[acc_start_ + t + 1] - 
                                                        state_vars[acc_start_ + t], 2);
    }
    // cost_vars[0] = vel_cost + cte_cost + psi_err_cost + del_cost + acc_cost + del_change_cost + acc_change_cost;
    // std::cout << "Total cost: " << cost_vars[0] << "Velocity: " << vel_cost << "\tCTE: " << cte_cost << 
    //         "\tPsi Err: " << psi_err_cost << "\tDel: " << del_cost << "\tAcc: " << acc_cost << 
    //         "\tDel_change: " << del_change_cost << "\tAcc change: " << acc_change_cost << "\n";  
}