#include "Helper.h"
#include <math.h>

Helper::Helper(){
    N_TIMESTEPS = 25;
    dt = 0.05;

    x_start = 0;
    y_start = N_TIMESTEPS;
    psi_start = 2*N_TIMESTEPS;
    vel_start = 3*N_TIMESTEPS;
    cte_start = 4*N_TIMESTEPS;
    psi_err_start = 5*N_TIMESTEPS;
    del_start = 6*N_TIMESTEPS;
    acc_start = 7*N_TIMESTEPS - 1;
}

Helper::~Helper() {}
const double Helper::pi() { return M_PI;}
double Helper::deg2rad(double x) { return x * pi() / 180; }
double Helper::rad2deg(double x) { return x * 180 / pi(); }