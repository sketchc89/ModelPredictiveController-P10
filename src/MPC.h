#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen/Core"
#include <cppad/cppad.hpp>

class MPC
{
  public:
    MPC();
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    virtual ~MPC();

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuatotions.
    std::vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
    void SetStateVariables(Dvector &state_vars);
    void SetStateVariableBounds(Dvector &state_vars_lowerbound, Dvector &state_vars_upperbound);
    void SetConstraintBounds(Dvector constraints_lowerbound, Dvector constraints_upperbound);
    double deg2rad(double x);
    double rad2deg(double x);
    size_t N_TIMESTEPS_;
    double dt_;
    size_t del_start_;
    size_t acc_start_;
    double L_f_;

  private:
    size_t N_VARS_;
    size_t N_CONSTRAINTS_;
};

#endif /* MPC_H */
