#include <gmock/gmock.h>
#include "FG_eval.h"
#include "Eigen-3.3/Eigen/Core"

TEST(FG_eval, Initializes) {
    Eigen::VectorXd v;
    v << 0;
    FG_eval fg(v);
}