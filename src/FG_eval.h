#ifndef FG_EVAL_H
#define FG_EVAL_H

#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

class FG_eval {
 public:
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs);
  virtual ~FG_eval();
  double ref_v = 50;
  const double L_f = 2.67; // length from front to center of gravity
  typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
  void operator()(ADvector& cost_vars, const ADvector& state_vars);
};

#endif