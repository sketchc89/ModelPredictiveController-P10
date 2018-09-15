#include <gmock/gmock.h>
#include "MPC.h"
#include "Helper.h"

Helper helper;

TEST(MPC, Initializes) {
    MPC mpc;
}

TEST(MPC, StateVariables) {
    MPC mpc;
    size_t n_vars = 6*helper.N_TIMESTEPS + 2*(helper.N_TIMESTEPS-1);
    CppAD::vector<double> state_vars(n_vars), state_vars_expected(n_vars);
    for (size_t i=0; i < state_vars_expected.size(); ++i) {
        state_vars_expected[i] = 0;
    }
    state_vars_expected[helper.vel_start] = 50;
    Eigen::VectorXd cur_state(8);
    
    cur_state << 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0;
    mpc.SetStateVariables(state_vars, cur_state);
    std::vector<double> v0, v1;
    for (size_t i=0; i<state_vars.size(); ++i) {
        v0.push_back(state_vars[i]);
        v1.push_back(state_vars_expected[i]);
    }
    EXPECT_EQ(v0, v1);
}