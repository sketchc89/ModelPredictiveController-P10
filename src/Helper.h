#ifndef HELPER_H
#define HELPER_H
#include <cstddef>
#include <math.h>
class Helper {
    public:
        Helper();
        virtual ~Helper();
        const double pi();
        double deg2rad(double x);
        double rad2deg(double x);
        
        size_t N_TIMESTEPS;
        double dt;
        size_t x_start;
        size_t y_start;
        size_t psi_start;
        size_t vel_start;
        size_t cte_start;
        size_t psi_err_start;
        size_t del_start;
        size_t acc_start;
};

#endif