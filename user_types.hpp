/*
 * user_types.hpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include "complex.hpp"

typedef struct domain_data {
    int n_r;
    int n_theta;
    int n_phi;
    int nt;
    double R;
} d_data;

typedef struct time_data {
    double to;
    double tf;
} t_data;

typedef struct physical_params {
    double h;
    double m;
} p_params;

typedef struct boundaries_and_psi_init {
    Complex psiInit;
    Complex psi0;
    Complex psiR;
} bound_and_psi_init;

typedef struct solver_data {
    Complex*** psi;
    Complex*** prob_density;
    Complex* psi_p_top;
    Complex* psi_p_bottom;
    Complex* prob_density_p_top;
    Complex* prob_density_p_bottom;
    double* r_p;
    double* theta_p;
    double* phi_p;
    double error_real;
    double error_im;
} s_data;

#endif /* USER_TYPES_HPP_ */
