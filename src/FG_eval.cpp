#include "FG_eval.h"
#include "MPC.h"
#include <cppad/cppad.hpp>
#include "Eigen/Core"

FG_eval::FG_eval(Eigen::VectorXd coeffs, MPC* mpc) { 
    this->coeffs = coeffs;
    this->N_TIMESTEPS_ = mpc->N_TIMESTEPS_;
    this->dt_ = mpc->dt_;
    this->x_start_ = mpc->x_start_;
    this->y_start_ = mpc->y_start_;
    this->psi_start_ = mpc->psi_start_;
    this->vel_start_ = mpc->vel_start_;
    this->cte_start_ = mpc->cte_start_;
    this->psi_err_start_ = mpc->psi_err_start_;
    this->del_start_ = mpc->del_start_;
    this->acc_start_ = mpc->acc_start_;

    this->ref_v_ = 50;
    this->L_f_ = 2.67;
}
FG_eval::~FG_eval() {}

void FG_eval::operator()(ADvector& cost_vars, const ADvector& state_vars) {
    cost_vars[0] = 0; // cost value

    for (size_t t=0; t<N_TIMESTEPS_; ++t) {
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
        CppAD::AD<double> del_1 = state_vars[del_start_ + t];
        CppAD::AD<double> del_2 = state_vars[del_start_ + t + 1];
        CppAD::AD<double> acc_0 = state_vars[acc_start_ + t - 1];
        CppAD::AD<double> acc_1 = state_vars[acc_start_ + t];
        CppAD::AD<double> acc_2 = state_vars[acc_start_ + t + 1];
        CppAD::AD<double> psi_dest = CppAD::atan(coeffs[1]);
        CppAD::AD<double> cte_x_0 = coeffs[0] + coeffs[1]*x_0;  
        
        cost_vars[0] += CppAD::pow(psi_1, 2);
        cost_vars[0] += CppAD::pow(vel_1 - ref_v_, 2);
        cost_vars[0] += CppAD::pow(cte_1, 2);
        
        if (t > 0) {
        cost_vars[1 + x_start_ + t] = x_1 - x_0 - vel_0*dt_*CppAD::cos(psi_0);
        cost_vars[1 + y_start_ + t] = y_1 - y_0 - vel_0*dt_*CppAD::cos(psi_0);
        cost_vars[1 + psi_start_ + t] = psi_1 - psi_0 - vel_0*del_0*dt_/L_f;
        cost_vars[1 + vel_start_ + t] = vel_1 - acc_0*dt_;
        cost_vars[1 + cte_start_ + t] = psi_err_1 - cte_x_0 + y_0 - vel_0*dt_*CppAD::sin(psi_err_0);
        cost_vars[1 + psi_err_start_ + t] = psi_err_1 - psi_0 + psi_dest -vel_0*del_0*dt_/L_f; 
        }
        
        if (t < N_TIMESTEPS_ - 1){
        cost_vars[0] += CppAD::pow(del_1, 2);
        cost_vars[0] += CppAD::pow(acc_1, 2);
        }
        
        if (t < N_TIMESTEPS_ - 2) {
        cost_vars[0] += CppAD::pow(del_2 - del_1, 2);
        cost_vars[0] += CppAD::pow(acc_2 - acc_1, 2);
        }
    }
}