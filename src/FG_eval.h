#ifndef FG_EVAL_H
#define FG_EVAL_H

#include "MPC.h"
#include "Eigen/Core"
#include <cppad/cppad.hpp>

class FG_eval {
  public:
    Eigen::VectorXd coeffs;
    FG_eval(Eigen::VectorXd coeffs, MPC* mpc);
    virtual ~FG_eval();
    const double L_f = 2.67; // length from front to center of gravity
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    void operator()(ADvector& cost_vars, const ADvector& state_vars);
  private:
    size_t N_TIMESTEPS_;
    double dt_;
    size_t x_start_;
    size_t y_start_;
    size_t psi_start_;
    size_t vel_start_;
    size_t cte_start_;
    size_t psi_err_start_;
    size_t del_start_;
    size_t acc_start_;

    double ref_v_;
    double L_f_;
};

#endif