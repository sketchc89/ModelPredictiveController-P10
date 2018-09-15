#include "FG_eval.h"
#include <cppad/cppad.hpp>
#include "Helper.h"

FG_eval::FG_eval(Eigen::VectorXd coeffs) { 
    this->coeffs = coeffs; 
}
FG_eval::~FG_eval() {}

void FG_eval::operator()(ADvector& cost_vars, const ADvector& state_vars) {
    cost_vars[0] = 0; // cost value

    Helper helper;
    for (size_t t=0; t<helper.N_TIMESTEPS; ++t) {
        CppAD::AD<double> x_0 = state_vars[helper.x_start + t - 1];
        CppAD::AD<double> x_1 = state_vars[helper.x_start + t];
        CppAD::AD<double> y_0 = state_vars[helper.y_start + t - 1];
        CppAD::AD<double> y_1 = state_vars[helper.y_start + t];
        CppAD::AD<double> psi_0 = state_vars[helper.psi_start + t - 1];
        CppAD::AD<double> psi_1 = state_vars[helper.psi_start + t];
        CppAD::AD<double> vel_0 = state_vars[helper.vel_start + t - 1];
        CppAD::AD<double> vel_1 = state_vars[helper.vel_start + t];
        CppAD::AD<double> cte_0 = state_vars[helper.cte_start + t - 1];
        CppAD::AD<double> cte_1 = state_vars[helper.cte_start + t];
        CppAD::AD<double> psi_err_0 = state_vars[helper.psi_err_start + t - 1];
        CppAD::AD<double> psi_err_1 = state_vars[helper.psi_err_start + t];
        CppAD::AD<double> del_0 = state_vars[helper.del_start + t - 1];
        CppAD::AD<double> del_1 = state_vars[helper.del_start + t];
        CppAD::AD<double> del_2 = state_vars[helper.del_start + t + 1];
        CppAD::AD<double> acc_0 = state_vars[helper.acc_start + t - 1];
        CppAD::AD<double> acc_1 = state_vars[helper.acc_start + t];
        CppAD::AD<double> acc_2 = state_vars[helper.acc_start + t + 1];
        CppAD::AD<double> psi_dest = CppAD::atan(coeffs[1]);
        CppAD::AD<double> cte_x_0 = coeffs[0] + coeffs[1]*x_0;  
        
        cost_vars[0] += CppAD::pow(psi_1, 2);
        cost_vars[0] += CppAD::pow(vel_1 - ref_v, 2);
        cost_vars[0] += CppAD::pow(cte_1, 2);
        
        if (t > 0) {
        cost_vars[1 + helper.x_start + t] = x_1 - x_0 - vel_0*helper.dt*CppAD::cos(psi_0);
        cost_vars[1 + helper.y_start + t] = y_1 - y_0 - vel_0*helper.dt*CppAD::cos(psi_0);
        cost_vars[1 + helper.psi_start + t] = psi_1 - psi_0 - vel_0*del_0*helper.dt/L_f;
        cost_vars[1 + helper.vel_start + t] = vel_1 - acc_0*helper.dt;
        cost_vars[1 + helper.cte_start + t] = psi_err_1 - cte_x_0 + y_0 - vel_0*helper.dt*CppAD::sin(psi_err_0);
        cost_vars[1 + helper.psi_err_start + t] = psi_err_1 - psi_0 + psi_dest -vel_0*del_0*helper.dt/L_f; 
        }
        
        if (t < helper.N_TIMESTEPS - 1){
        cost_vars[0] += CppAD::pow(del_1, 2);
        cost_vars[0] += CppAD::pow(acc_1, 2);
        }
        
        if (t < helper.N_TIMESTEPS - 2) {
        cost_vars[0] += CppAD::pow(del_2 - del_1, 2);
        cost_vars[0] += CppAD::pow(acc_2 - acc_1, 2);
        }
    }
}