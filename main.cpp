/*
 * main.cpp
 *
 *  Created on: Dec 29, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include "mem_ops.hpp"
#include "solver.hpp"
#include "user_types.hpp"

double V(double r, double t) {
    return -1/(r + 1e-20);
}

int main(int argc, char* argv[]) {
    d_data domain_data;
    t_data time_data;
    p_params physical_params;
    bound_and_psi_init boundaries_and_psi_init;
    s_data solver_data;

    /* Parameters */
    domain_data.n_r = 10;                         //Number of nodes in the radial direction
    domain_data.n_theta = 10;                     //Number of nodes in the theta direction
    domain_data.n_phi = 10;                       //Number of nodes in the phi direction
    domain_data.nt = 10;                          //Number of timesteps
    domain_data.R = 1.0;                          //Length of domain in the radial direction

    time_data.to = 0.0;                           //Initial time
    time_data.tf = 1.0;                           //Final time

    physical_params.h = 1.0;                      //Constant
    physical_params.m = 1.0;                      //Mass of particle

    boundaries_and_psi_init.psiInit.a = 0.5;      //Real value of psi at time t = 0
    boundaries_and_psi_init.psiInit.b = 0.5;      //Imaginary value of psi at time t = 0
    boundaries_and_psi_init.psi0.a = 0.0;         //Real value of psi at the center of the domain
    boundaries_and_psi_init.psi0.b = 0.0;         //Imaginary value of psi at the center of the domain
    boundaries_and_psi_init.psiR.a = 0.0;         //Real value of psi at R
    boundaries_and_psi_init.psiR.b = 0.0;         //Imaginary value of psi at R

    /* Allocate memory for solver data */
    solver_data.psi = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
    solver_data.prob_density = mat3D(domain_data.n_r, domain_data.n_theta, domain_data.n_phi);
    solver_data.psi_p_top = new Complex[domain_data.n_r];
    solver_data.psi_p_bottom = new Complex[domain_data.n_r];
    solver_data.prob_density_p_top = new Complex[domain_data.n_r];
    solver_data.prob_density_p_bottom = new Complex[domain_data.n_r];
    solver_data.r_p = new double[domain_data.n_r];
    solver_data.theta_p = new double[domain_data.n_theta];
    solver_data.phi_p = new double[domain_data.n_phi];

    /* Execute solver */
    solver(domain_data, time_data, physical_params, boundaries_and_psi_init, &solver_data);

    /* Print some results */
    for(int i = 0; i < domain_data.n_r; ++i) {
        printf("psi_p_top[%i] real: %f, im: %f, psi_p_bottom[%i] real: %f, im: %f\n",
        		i, solver_data.psi_p_top[i].a, solver_data.psi_p_top[i].b, i, solver_data.psi_p_bottom[i].a, solver_data.psi_p_bottom[i].b);
    }

    /* Free allocates data */
    free_mat3D(solver_data.psi, domain_data.n_r, domain_data.n_theta);
    free_mat3D(solver_data.prob_density, domain_data.n_r, domain_data.n_theta);
    delete [] solver_data.r_p;
    delete [] solver_data.theta_p;
    delete [] solver_data.phi_p;
    delete [] solver_data.psi_p_top;
    delete [] solver_data.psi_p_bottom;
    delete [] solver_data.prob_density_p_top;
    delete [] solver_data.prob_density_p_bottom;

    printf("done\n");

    return 0;
}
