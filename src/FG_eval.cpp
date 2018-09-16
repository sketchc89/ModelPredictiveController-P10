#include "FG_eval.h"
#include "MPC.h"
#include <cppad/cppad.hpp>
#include "Eigen/Core"

FG_eval::FG_eval(Eigen::VectorXd coeffs, MPC* mpc) { 
    this->coeffs = coeffs;
    this->N_TIMESTEPS_ = mpc->N_TIMESTEPS_;
    this->dt_ = mpc->dt_;
    this->L_f_ = mpc->L_f_;
    this->x_start_ = mpc->x_start_;
    this->y_start_ = mpc->y_start_;
    this->psi_start_ = mpc->psi_start_;
    this->vel_start_ = mpc->vel_start_;
    this->cte_start_ = mpc->cte_start_;
    this->psi_err_start_ = mpc->psi_err_start_;
    this->del_start_ = mpc->del_start_;
    this->acc_start_ = mpc->acc_start_;

    this->ref_v_ = 25;
    this->L_f_ = 2.67;
}
FG_eval::~FG_eval() {}

void FG_eval::operator()(ADvector& cost_vars, const ADvector& state_vars) {
    cost_vars[0] = 0; // cost value

    int vel_weight = 10;
    int cte_weight = 5000;
    int psi_err_weight = 5000;
    int del_weight = 10;
    int acc_weight = 10;
    int del_change_weight = 1000;
    int acc_change_weight = 10;
    CppAD::AD<double> vel_cost = 0.0;
    CppAD::AD<double> cte_cost = 0.0;
    CppAD::AD<double> psi_err_cost = 0.0;
    CppAD::AD<double> del_cost = 0.0;
    CppAD::AD<double> acc_cost = 0.0;
    CppAD::AD<double> del_change_cost = 0.0;
    CppAD::AD<double> acc_change_cost = 0.0; 

    for (size_t t=0; t < N_TIMESTEPS_; ++t) {
        vel_cost += vel_weight * CppAD::pow(state_vars[vel_start_ + t] - ref_v_, 2);
        cte_cost += cte_weight * CppAD::pow(state_vars[cte_start_ + t], 2);
        psi_err_cost += psi_err_weight * CppAD::pow(state_vars[psi_err_start_ + t], 2);
    }
    for (size_t t=0; t < N_TIMESTEPS_ - 1; ++t) {
        del_cost += del_weight * CppAD::pow(state_vars[del_start_ + t], 2);
        acc_cost += acc_weight * CppAD::pow(state_vars[acc_start_ + t], 2);
    }
    for (size_t t=0; t < N_TIMESTEPS_ - 2; ++t) {
        del_change_cost += del_change_weight * CppAD::pow(state_vars[del_start_ + t + 1] - 
                                                        state_vars[del_start_ + t], 2);
        acc_change_cost += acc_change_weight * CppAD::pow(state_vars[acc_start_ + t + 1] - 
                                                        state_vars[acc_start_ + t], 2);
    }
    cost_vars[0] = vel_cost + cte_cost + psi_err_cost + del_cost + acc_cost + del_change_cost + acc_change_cost;
    std::cout << "Total cost: " << cost_vars[0] << "Velocity: " << vel_cost << "\tCTE: " << cte_cost << 
            "\tPsi Err: " << psi_err_cost << "\tDel: " << del_cost << "\tAcc: " << acc_cost << 
            "\tDel_change: " << del_change_cost << "\tAcc change: " << acc_change_cost << "\n";  
    
    cost_vars[x_start_ + 1]       = state_vars[x_start_];
    cost_vars[y_start_ + 1]       = state_vars[y_start_];
    cost_vars[vel_start_ + 1]     = state_vars[vel_start_];
    cost_vars[cte_start_ + 1]     = state_vars[cte_start_];
    cost_vars[psi_err_start_ + 1] = state_vars[psi_err_start_];
    
    for (size_t t=1; t < N_TIMESTEPS_; ++t) {
        // std::cout << "T=\t" << t << "\n";
        CppAD::AD<double> x_0 = state_vars[x_start_ + t - 1];
        CppAD::AD<double> x_1 = state_vars[x_start_ + t];
        CppAD::AD<double> y_0 = state_vars[y_start_ + t - 1];
        CppAD::AD<double> y_1 = state_vars[y_start_ + t];
        CppAD::AD<double> psi_0 = state_vars[psi_start_ + t - 1];
        CppAD::AD<double> psi_1 = state_vars[psi_start_ + t];
        CppAD::AD<double> vel_0 = state_vars[vel_start_ + t - 1];
        CppAD::AD<double> vel_1 = state_vars[vel_start_ + t];
        CppAD::AD<double> cte_0 = state_vars[cte_start_ + t - 1];
        CppAD::AD<double> cte_1 = state_vars[cte_start_ + t];
        CppAD::AD<double> psi_err_0 = state_vars[psi_err_start_ + t - 1];
        CppAD::AD<double> psi_err_1 = state_vars[psi_err_start_ + t];        
        CppAD::AD<double> del_0 = state_vars[del_start_ + t - 1];
        CppAD::AD<double> acc_0 = state_vars[acc_start_ + t - 1];
        CppAD::AD<double> psi_dest = CppAD::atan(coeffs[1] + 2*coeffs[2]*x_0 + 
                                                3*coeffs[3]*CppAD::pow(x_0, 2));
        CppAD::AD<double> cte_x_0 = coeffs[0] + coeffs[1]*x_0 + 
                                    coeffs[2]*CppAD::pow(x_0, 2) + 
                                    coeffs[3]*CppAD::pow(x_0, 3);  
        
        cost_vars[1 + x_start_ + t] = x_1 - x_0 - vel_0*dt_*CppAD::cos(psi_0);
        cost_vars[1 + y_start_ + t] = y_1 - y_0 - vel_0*dt_*CppAD::sin(psi_0);
        cost_vars[1 + psi_start_ + t] = psi_1 - psi_0 + (vel_0*del_0*dt_) / L_f_;
        cost_vars[1 + vel_start_ + t] = vel_1 - vel_0 - acc_0*dt_;
        cost_vars[1 + cte_start_ + t] = cte_1 - cte_x_0 + y_0 - 
                                        vel_0*dt_*CppAD::sin(psi_err_0);
        cost_vars[1 + psi_err_start_ + t] = psi_err_1 - psi_0 + psi_dest +
                                            (vel_0*del_0*dt_) / L_f_; 

    }
    // for (size_t t=0; t < cost_vars.size(); ++t) {
    //     std::cout << cost_vars[t] << "\t";
    // }
    std::cout << "\n";
}