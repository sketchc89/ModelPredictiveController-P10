#include "FG_eval.h"
#include "MPC.h"
#include <cppad/cppad.hpp>
#include "Eigen/Core"

FG_eval::FG_eval(Eigen::VectorXd cur_state, Eigen::VectorXd coeffs, MPC *mpc)
{
    this->coeffs_ = coeffs;
    this->cur_x_ = cur_state[0];
    this->cur_y_ = cur_state[1];
    this->cur_psi_ = cur_state[2];
    this->cur_vel_ = cur_state[3];
    this->cur_cte_ = cur_state[4];
    this->cur_psi_err_ = cur_state[5];
    this->N_TIMESTEPS_ = mpc->N_TIMESTEPS_;
    x_.resize(N_TIMESTEPS_);
    y_.resize(N_TIMESTEPS_);
    psi_.resize(N_TIMESTEPS_);
    vel_.resize(N_TIMESTEPS_);
    cte_.resize(N_TIMESTEPS_);
    psi_err_.resize(N_TIMESTEPS_);

    this->dt_ = mpc->dt_;
    this->L_f_ = mpc->L_f_;
    this->del_start_ = mpc->del_start_;
    this->acc_start_ = mpc->acc_start_;
    this->ref_v_ = 70;
    this->L_f_ = 2.67;
    this->LATENCY_ = 0.1;
}
FG_eval::~FG_eval() {}

void FG_eval::operator()(ADvector &cost_vars, const ADvector &state_vars)
{
    cost_vars[0] = 0; // cost value
    KinematicModel(state_vars);

    int vel_weight = 1;
    int cte_weight = 5;
    int psi_err_weight = 100;
    int del_weight = 30;
    int fast_turn_weight = 5;
    int acc_weight = 1;
    int del_change_weight = 1000;
    int acc_change_weight = 1;
    
    for (size_t t = 0; t < N_TIMESTEPS_; ++t)
    {
        cost_vars[0] += vel_weight * CppAD::pow(vel_[t] - ref_v_, 2);
        cost_vars[0] += cte_weight * CppAD::pow(cte_[t], 2);
        cost_vars[0] += psi_err_weight * CppAD::pow(psi_err_[t], 2);
    }
    for (size_t t = 0; t < N_TIMESTEPS_ - 1; ++t)
    {
        cost_vars[0] += del_weight * CppAD::pow(state_vars[del_start_ + t], 2);
        cost_vars[0] += acc_weight * CppAD::pow(state_vars[acc_start_ + t], 2);
        cost_vars[0] += fast_turn_weight * CppAD::pow(vel_[t]*state_vars[del_start_ + t], 2);
        
    }
    for (size_t t = 0; t < N_TIMESTEPS_ - 2; ++t)
    {
        cost_vars[0] += del_change_weight * CppAD::pow(state_vars[del_start_ + t + 1] -
                                                           state_vars[del_start_ + t],
                                                       2);
        cost_vars[0] += acc_change_weight * CppAD::pow(state_vars[acc_start_ + t + 1] -
                                                           state_vars[acc_start_ + t],
                                                       2);
    }
}

void FG_eval::KinematicModel(const ADvector &state_vars)
{
    x_[0] = cur_x_;
    y_[0] = cur_y_;
    psi_[0] = cur_psi_;
    vel_[0] = cur_vel_;
    cte_[0] = cur_cte_;
    psi_err_[0] = cur_psi_err_;
    //latency
    x_[0] += vel_[0] * CppAD::cos(psi_[0])*LATENCY_;
    y_[0] += vel_[0] * CppAD::sin(psi_[0])*LATENCY_;
    psi_[0] += vel_[0] * state_vars[del_start_]*LATENCY_ / L_f_;
    vel_[0] += state_vars[acc_start_]*LATENCY_;
    cte_[0] = coeffs_[3] * CppAD::pow(x_[0], 3) + coeffs_[2] * CppAD::pow(x_[0], 2) +
              coeffs_[1] * x_[0] + coeffs_[0] - y_[0];
    psi_err_[0] = -CppAD::atan(coeffs_[1] + 2 * coeffs_[2] * x_[0] +
                            3 * coeffs_[3] * CppAD::pow(x_[0], 2)) + psi_[0];          

    for (size_t t = 1; t < N_TIMESTEPS_; ++t)
    {
        CppAD::AD<double> del = state_vars[del_start_ + t - 1];
        CppAD::AD<double> acc = state_vars[acc_start_ + t - 1];
        x_[t] = x_[t - 1] + vel_[t - 1] * CppAD::cos(psi_[t - 1]) * dt_;
        y_[t] = y_[t - 1] + vel_[t - 1] * CppAD::sin(psi_[t - 1]) * dt_;
        psi_[t] = psi_[t - 1] + vel_[t - 1] * del * dt_ / L_f_;
        vel_[t] = vel_[t - 1] + acc * dt_;
        cte_[t] = coeffs_[3] * CppAD::pow(x_[t], 3) + coeffs_[2] * CppAD::pow(x_[t], 2) +
                  coeffs_[1] * x_[t] + coeffs_[0] - y_[t];
        psi_err_[t] = psi_[t] - CppAD::atan(coeffs_[1] + 2 * coeffs_[2] * x_[t] +
                        3 * coeffs_[3] * CppAD::pow(x_[t], 2));
    }
}