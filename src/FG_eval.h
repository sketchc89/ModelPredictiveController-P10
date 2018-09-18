#ifndef FG_EVAL_H
#define FG_EVAL_H

#include "MPC.h"
#include "Eigen/Core"
#include <cppad/cppad.hpp>

class FG_eval {
  public:
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    FG_eval(Dvector &vel, Dvector &cte, Dvector &psi_err, MPC* mpc);
    virtual ~FG_eval();
    void operator()(ADvector &cost_vars, const ADvector &state_vars);
    void KinematicModel(const ADvector &state_vars,  ADvector &x, 
                        ADvector &y, ADvector &psi, ADvector &vel,
                        ADvector &cte, ADvector &psi_err);
  private:
    Eigen::VectorXd coeffs_;
    size_t N_TIMESTEPS_;
    double dt_;
    Dvector vel_;
    Dvector cte_;
    Dvector psi_err_;
    size_t del_start_;
    size_t acc_start_;

    double ref_v_;
    double L_f_;
};

#endif