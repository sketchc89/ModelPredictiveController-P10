#include <gmock/gmock.h>
#include "FG_eval.h"
#include "MPC.h"
#include "Eigen-3.3/Eigen/Core"

TEST(FG_eval, Initializes) {
    Eigen::VectorXd v(1);
    v << 0;
    MPC* mpc;
    FG_eval fg(v, mpc);

}