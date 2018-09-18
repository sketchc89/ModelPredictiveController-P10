#ifndef FG_EVAL_H
#define FG_EVAL_H

#include "MPC.h"
#include "Eigen/Core"
#include <cppad/cppad.hpp>

class FG_eval {
  public:
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    FG_eval(Eigen::VectorXd cur_state, Eigen::VectorXd coeffs, MPC* mpc);
    virtual ~FG_eval();
    void operator()(ADvector &cost_vars, const ADvector &state_vars);
    void KinematicModel(const ADvector &state_vars);
    ADvector x_;
    ADvector y_;
  private:
    Eigen::VectorXd coeffs_;
    size_t N_TIMESTEPS_;
    double dt_;
    ADvector psi_;
    ADvector vel_;
    ADvector cte_;
    ADvector psi_err_;
    size_t del_start_;
    size_t acc_start_;
    double cur_x_;
    double cur_y_;
    double cur_psi_;
    double cur_vel_;
    double cur_cte_;
    double cur_psi_err_;
    double ref_v_;
    double L_f_;
};

#endif