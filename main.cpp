/*
 * main.cpp
 *
 *  Created on: Dec 29, 2020
 *      Author: d-w-h
 *
 *      Numerical solution of the Schrodinger equation for the hydrogen atom:
 *
 *      ih*d(psi)/dt = -h*h/2m*div(grad(psi)) + V*psi
 *
 */

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "mem_ops.hpp"
#include "solver.hpp"
#include "user_types.hpp"

double V(double r, double t) {
    /* Potential */
    double k = 1.0;
    return -k/(r + 1e-20);
}

double psi_init_real(double r, double theta, double phi) {
    /* Initial Re(psi) */
    return sin(2*M_PI*r);
}

double psi_init_im(double r, double theta, double phi) {
    /* Initial Im(psi) */
    return 0.0;
}

double psi_init_p_top_real(double r) {
    /* Initial Re(psi) for top pole nodes */
    return sin(2*M_PI*r);
}

double psi_init_p_top_im(double r) {
    /* Initial Im(psi) for top pole nodes */
    return 0.0;
}

double psi_init_p_bottom_real(double r) {
    /* Initial Re(psi) for bottom pole nodes */
    return sin(2*M_PI*r);
}

double psi_init_p_bottom_im(double r) {
    /* Initial Im(psi) for bottom pole nodes */
    return 0.0;
}

int main(int argc, char* argv[]) {
    d_data domain_data;
    t_data time_data;
    p_params physical_params;
    b_data boundaries;
    s_data solver_data;
    clock_t t_start = clock();

    /* Parameters */
    domain_data.n_r = 21;            //Number of nodes in the radial direction
    domain_data.n_theta = 11;        //Number of nodes in the theta direction
    domain_data.n_phi = 11;          //Number of nodes in the phi direction
    domain_data.nt = 20;             //Number of timesteps
    domain_data.R = 1.0;             //Length of domain in the radial direction

    time_data.to = 0.0;              //Initial time
    time_data.tf = 2.0;              //Final time

    physical_params.h = 1.0;         //Constant
    physical_params.m = 1.0;         //Mass of particle

    boundaries.psi0.a = 0.0;         //Real value of psi at the center of the domain
    boundaries.psi0.b = 0.0;         //Imaginary value of psi at the center of the domain
    boundaries.psiR.a = 0.0;         //Real value of psi at R
    boundaries.psiR.b = 0.0;         //Imaginary value of psi at R

    /* Allocate memory for solver data ans initial psi */
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
    solver(domain_data, time_data, physical_params, boundaries, &solver_data);

    /* Print some results */
    for(int i = 0; i < domain_data.n_r; ++i) {
        printf("r_p: %f, psi_p_top[%i] real: %f, im: %f, psi_p_bottom[%i] real: %f, im: %f, psi[%i][%i][%i] real: %f, im: %f\n",
                solver_data.r_p[i],
                i, solver_data.psi_p_top[i].a, solver_data.psi_p_top[i].b, i, solver_data.psi_p_bottom[i].a, solver_data.psi_p_bottom[i].b,
                i, domain_data.n_theta/2, domain_data.n_phi/2, solver_data.psi[i][domain_data.n_theta/2][domain_data.n_phi/2].a,
                solver_data.psi[i][domain_data.n_theta/2][domain_data.n_phi/2].b);
    }

    /* Free allocated data */
    free_mat3D(solver_data.psi, domain_data.n_r, domain_data.n_theta);
    free_mat3D(solver_data.prob_density, domain_data.n_r, domain_data.n_theta);
    delete [] solver_data.r_p;
    delete [] solver_data.theta_p;
    delete [] solver_data.phi_p;
    delete [] solver_data.psi_p_top;
    delete [] solver_data.psi_p_bottom;
    delete [] solver_data.prob_density_p_top;
    delete [] solver_data.prob_density_p_bottom;

    /* Print time taken for execution */
    clock_t t_end = clock();
    double sim_time = double (t_end - t_start)/CLOCKS_PER_SEC;
    printf("time taken: %f\n", sim_time);

    return 0;
}
