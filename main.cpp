/*
 * main.cpp
 *
 *  Created on: Dec 29, 2020
 *      Author: d-w-h
 */

#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <time.h>

class Complex {
public:
    double a;
    double b;
    Complex();
    Complex(double a, double b);
    Complex operator+(const Complex& m);
    Complex operator*(const Complex& m);
    Complex operator-(const Complex& m);
};

Complex::Complex() {}

Complex::Complex(double a, double b) {
    this->a = a;
    this->b = b;
}

Complex Complex::operator+(const Complex& m) {
    Complex result_addition(0,0);
    result_addition.a = this->a + m.a;
    result_addition.b = this->b + m.b;

    return result_addition;
}

Complex Complex::operator*(const Complex& m) {
    Complex result_multiplication(0,0);
    result_multiplication.a = this->a*m.a - this->b*m.b;
    result_multiplication.b = this->a*m.b + this->b*m.a;

    return result_multiplication;
}

Complex Complex::operator-(const Complex& m) {
    Complex result_subtraction(0,0);
    result_subtraction.a = this->a - m.a;
    result_subtraction.b = this->b - m.b;

    return result_subtraction;
}

Complex*** mat3D(int n_r, int n_theta, int n_phi) {

    Complex*** f = new Complex**[n_r];

    for(int i = 0; i < n_r; ++i) {
        f[i] = new Complex*[n_theta];
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            f[i][j] = new Complex[n_phi];
        }
    }

    return f;
}

double V(double r, double t) {
    return -1/(r + 1e-20);
}

int main(int argc, char* argv[]) {
    int n_r, n_theta, n_phi, nt;
    double R, to, tf, h, m;
    Complex psi0, psiR, psiInit;
    clock_t t_start = clock();

    /* Parameters */
    n_r = 10;
    n_theta = 10;
    n_phi = 10;

    nt = 10;
    to = 0.0;
    tf = 1.0;

    R = 1.0;

    h = 1.0;
    m = 1.0;

    psi0.a = 0.0;
    psi0.b = 0.0;
    psiR.a = 0.0;
    psiR.b = 0.0;
    psiInit.a = 0.5;
    psiInit.b = 0.5;

    /* Start calculations */
    double* r = new double[n_r];
    double* theta = new double[n_theta];
    double* phi = new double[n_phi+1];
    double* r_p = new double[n_r];
    double* theta_p = new double[n_theta];
    double* phi_p = new double[n_phi];

    Complex*** psi = mat3D(n_r, n_theta, n_phi);
    Complex*** psi_prev = mat3D(n_r, n_theta, n_phi);
    Complex*** psio = mat3D(n_r, n_theta, n_phi);
    Complex*** prob_density = mat3D(n_r, n_theta, n_phi);
    Complex* psi_p_top = new Complex[n_r];
    Complex* psi_p_bottom = new Complex[n_r];
    Complex* psi_p_topo = new Complex[n_r];
    Complex* psi_p_bottomo = new Complex[n_r];
    Complex* prob_density_p_top = new Complex[n_r];
    Complex* prob_density_p_bottom = new Complex[n_r];

    double dr = R/n_r;
    double dtheta = 2*M_PI/n_theta;
    double dphi = M_PI/(n_phi+1);
    double dt = (tf - to)/nt;
    double t = 0.0;
    double alpha = h/(2*m);

    /* Initialize r vector */
    for(int i = 0; i < n_r; ++i) {
        r[i] = i*dr + 0.5*dr;
    }

    /* Initialize r_p vector */
    for(int i = 0; i < n_r; ++i) {
        r_p[i] = i*dr;
    }

    /* Initialize theta vector */
    for(int i = 0; i < n_theta; ++i) {
        theta[i] = i*dtheta;
    }

    /* Initialize theta_p vector */
    for(int i = 0; i < n_theta; ++i) {
        theta_p[i] = i*dtheta + 0.5*dtheta;
    }

    /* Initialize phi vector */
    for(int i = 0; i < n_phi + 1; ++i) {
        phi[i] = i*dphi + 0.5*dphi;
    }

    /* Initialize phi_p vector */
    for(int i = 0; i < n_phi; ++i) {
        phi_p[i] = (i + 1)*dphi;
    }

    /* Initialize probability density */
    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                prob_density[i][j][k].a = 0.0;
                prob_density[i][j][k].b = 0.0;
            }
        }
    }

    /* Initialize probability density top and bottom pole nodes */
    for(int i = 0; i < n_r; ++i) {
        prob_density_p_top[i].a = 0.0;
        prob_density_p_top[i].b = 0.0;
        prob_density_p_bottom[i].a = 0.0;
        prob_density_p_bottom[i].b = 0.0;
    }

    /* Initialize psi */
    for(int i = 1; i < n_r - 1; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                psi[i][j][k].a = 0.0;
                psi[i][j][k].b = 0.0;
                psio[i][j][k].a = psiInit.a;
                psio[i][j][k].b = psiInit.b;
            }
        }
    }

    for(int j = 0; j < n_theta; ++j) {
        for(int k = 0; k < n_phi; ++k) {
            psi[0][j][k].a = psi0.a;
            psi[0][j][k].b = psi0.b;
            psio[0][j][k].a = psi0.a;
            psio[0][j][k].b = psi0.b;
        }
    }

    for(int i = 0; i < n_r; ++i) {
        for(int j = 0; j < n_theta; ++j) {
            for(int k = 0; k < n_phi; ++k) {
                psi_prev[i][j][k].a = psi[i][j][k].a;
                psi_prev[i][j][k].b = psi[i][j][k].b;
            }
        }
    }


    for(int j = 0; j < n_theta; ++j) {
        for(int k = 0; k < n_phi; ++k) {
            psi[n_r-1][j][k].a = psiR.a;
            psi[n_r-1][j][k].b = psiR.b;
            psio[n_r-1][j][k].a = psiR.a;
            psio[n_r-1][j][k].b = psiR.b;
        }
    }

    /* Initialize pole nodes */
    psi_p_top[0].a = psi0.a;
    psi_p_top[0].b = psi0.b;
    psi_p_bottom[0].a = psi0.a;
    psi_p_bottom[0].b = psi0.b;
    for(int i = 1; i < n_r - 1; ++i) {
        psi_p_top[i].a = 0.0;
        psi_p_top[i].b = 0.0;
        psi_p_bottom[i].a = 0.0;
        psi_p_bottom[i].b = 0.0;
        psi_p_topo[i].a = psiInit.a;
        psi_p_topo[i].b = psiInit.b;
        psi_p_bottomo[i].a = psiInit.a;
        psi_p_bottomo[i].b = psiInit.b;
    }
    psi_p_top[n_r-1].a = psiR.a;
    psi_p_top[n_r-1].b = psiR.b;
    psi_p_bottom[n_r-1].a = psiR.a;
    psi_p_bottom[n_r-1].b = psiR.b;

    /* Start simulation */
    int max_gs = 500;
    int timestep = 0;
    Complex I(0,1);
    double max_real = -1e+8;
    double min_real = 1e+8;
    double max_im = -1e+8;
    double min_im = 1e+8;
    double max_pd = -1e+8;

    while(t < tf) {
        /* Start Gauss-Seidel iterations */
        int gs_it = 0;
        while(gs_it < max_gs) {
            /* Calculate 'previous values' for convergence check */
            for(int i = 0; i < n_r; ++i) {
                for(int j = 0; j < n_theta; ++j) {
                    for(int k = 0; k < n_phi; ++k) {
                        psi_prev[i][j][k].a = psi[i][j][k].a;
                        psi_prev[i][j][k].b = psi[i][j][k].b;
                    }
                }
            }

            /* Central top pole nodes */
            for(int i = 1; i < n_r - 1; ++i) {
                double Ars = 2*M_PI*r[i-1]*r[i-1]*(1 - cos(dphi));
                double Arn = 2*M_PI*r[i]*r[i]*(1 - cos(dphi));
                double Aphit = r_p[i]*sin(phi[0])*dtheta*dr;
                double dV = 2*M_PI/3*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1])*(1 - cos(dphi));

                double aS = Ars*alpha/(dr*dV);
                double aN = Arn*alpha/(dr*dV);
                double aT = Aphit*alpha/(r_p[i]*dphi*dV);

                Complex psi_s(0,0);
                psi_s.a = -aS*psi_p_top[i-1].a;
                psi_s.b = -aS*psi_p_top[i-1].b;

                Complex psi_n(0,0);
                psi_n.a = -aN*psi_p_top[i+1].a;
                psi_n.b = -aN*psi_p_top[i+1].b;

                Complex sum_psi_t(0,0);
                for(int j = 0; j < n_theta; ++j) {
                    sum_psi_t.a = sum_psi_t.a + aT*psi[i][j][0].a;
                    sum_psi_t.b = sum_psi_t.b + aT*psi[i][j][0].b;
                }

                Complex psi_t(0,0);
                psi_t.a = -sum_psi_t.a;
                psi_t.b = -sum_psi_t.b;

                Complex ihpsipo_dt(0,0);
                ihpsipo_dt.a = psi_p_topo[i].a*h/dt;
                ihpsipo_dt.b = psi_p_topo[i].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Complex Fp(0,0);
                Fp = psi_s + psi_n + psi_t + ihpsipo_dt;

                double a = aS + aN + n_theta*aT + V(r_p[i], t);
                double b = h/dt;

                Complex min_Fp(-Fp.a, -Fp.b);
                Complex ib_plus_a(a,b);

                psi_p_top[i] = min_Fp * ib_plus_a;
                psi_p_top[i].a = psi_p_top[i].a/(a*a+b*b);
                psi_p_top[i].b = psi_p_top[i].b/(a*a+b*b);
            }

            /* Central bottom pole nodes */
            for(int i = 1; i < n_r - 1; ++i) {
                double Ars = 2*M_PI*r[i-1]*r[i-1]*(1 - cos(dphi));
                double Arn = 2*M_PI*r[i]*r[i]*(1 - cos(dphi));
                double Aphib = r_p[i]*sin(phi[n_phi])*dtheta*dr;
                double dV = 2*M_PI/3*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1])*(1 - cos(dphi));

                double aS = Ars*alpha/(dr*dV);
                double aN = Arn*alpha/(dr*dV);
                double aB = Aphib*alpha/(r_p[i]*dphi*dV);

                Complex psi_s(0,0);
                psi_s.a = -aS*psi_p_bottom[i-1].a;
                psi_s.b = -aS*psi_p_bottom[i-1].b;

                Complex psi_n(0,0);
                psi_n.a = -aN*psi_p_bottom[i+1].a;
                psi_n.b = -aN*psi_p_bottom[i+1].b;

                Complex sum_psi_t(0,0);
                for(int j = 0; j < n_theta; ++j) {
                    sum_psi_t.a = sum_psi_t.a + aB*psi[i][j][n_phi-1].a;
                    sum_psi_t.b = sum_psi_t.b + aB*psi[i][j][n_phi-1].b;
                }

                Complex psi_t(0,0);
                psi_t.a = -sum_psi_t.a;
                psi_t.b = -sum_psi_t.b;

                Complex ihpsipo_dt(0,0);
                ihpsipo_dt.a = psi_p_bottomo[i].a*h/dt;
                ihpsipo_dt.b = psi_p_bottomo[i].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Complex Fp(0,0);
                Fp = psi_s + psi_n + psi_t + ihpsipo_dt;

                double a = aS + aN + n_theta*aB + V(r_p[i], t);
                double b = h/dt;

                Complex min_Fp(-Fp.a, -Fp.b);
                Complex ib_plus_a(a,b);

                psi_p_bottom[i] = min_Fp * ib_plus_a;
                psi_p_bottom[i].a = psi_p_bottom[i].a/(a*a+b*b);
                psi_p_bottom[i].b = psi_p_bottom[i].b/(a*a+b*b);
            }

            /* Central top cone nodes */
            for(int i = 1; i < n_r - 1; ++i) {
                for(int j = 1; j < n_theta - 1; ++j) {
                    double dV1 = r_p[i]*r_p[i]*dr;
                    double dV2 = r_p[i]*r_p[i]*sin(phi_p[0])*sin(phi_p[0])*dtheta;
                    double dV3 = r_p[i]*r_p[i]*sin(phi_p[0])*dphi;

                    double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                    double aN = r[i]*r[i]*alpha/(dV1*dr);
                    double aW = alpha/(dV2*dtheta);
                    double aE = alpha/(dV2*dtheta);
                    double aB = sin(phi[0])*alpha/(dV3*dphi);
                    double aT = sin(phi[1])*alpha/(dV3*dphi);

                    Complex psi_s(0,0);
                    psi_s.a = -aS*psi[i-1][j][0].a;
                    psi_s.b = -aS*psi[i-1][j][0].b;

                    Complex psi_n(0,0);
                    psi_n.a = -aN*psi[i+1][j][0].a;
                    psi_n.b = -aN*psi[i+1][j][0].b;

                    Complex psi_w(0,0);
                    psi_w.a = -aW*psi[i][j-1][0].a;
                    psi_w.b = -aW*psi[i][j-1][0].b;

                    Complex psi_e(0,0);
                    psi_e.a = -aE*psi[i][j+1][0].a;
                    psi_e.b = -aE*psi[i][j+1][0].b;

                    Complex psi_b(0,0);
                    psi_b.a = -aB*psi_p_top[i].a;
                    psi_b.b = -aB*psi_p_top[i].b;

                    Complex psi_t(0,0);
                    psi_t.a = -aT*psi[i][j][1].a;
                    psi_t.b = -aT*psi[i][j][1].b;

                    Complex ihpsipo_dt(0,0);
                    ihpsipo_dt.a = psio[i][j][0].a*h/dt;
                    ihpsipo_dt.b = psio[i][j][0].b*h/dt;
                    ihpsipo_dt = ihpsipo_dt*I;

                    Complex Fp(0,0);
                    Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                    double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                    double b = h/dt;

                    Complex min_Fp(-Fp.a, -Fp.b);
                    Complex ib_plus_a(a,b);

                    psi[i][j][0] = min_Fp * ib_plus_a;
                    psi[i][j][0].a = psi[i][j][0].a/(a*a+b*b);
                    psi[i][j][0].b = psi[i][j][0].b/(a*a+b*b);
                }

                /* First node central top cone nodes */
                double dV1 = r_p[i]*r_p[i]*dr;
                double dV2 = r_p[i]*r_p[i]*sin(phi_p[0])*sin(phi_p[0])*dtheta;
                double dV3 = r_p[i]*r_p[i]*sin(phi_p[0])*dphi;

                double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                double aN = r[i]*r[i]*alpha/(dV1*dr);
                double aW = alpha/(dV2*dtheta);
                double aE = alpha/(dV2*dtheta);
                double aB = sin(phi[0])*alpha/(dV3*dphi);
                double aT = sin(phi[1])*alpha/(dV3*dphi);

                Complex psi_s(0,0);
                psi_s.a = -aS*psi[i-1][0][0].a;
                psi_s.b = -aS*psi[i-1][0][0].b;

                Complex psi_n(0,0);
                psi_n.a = -aN*psi[i+1][0][0].a;
                psi_n.b = -aN*psi[i+1][0][0].b;

                Complex psi_w(0,0);
                psi_w.a = -aW*psi[i][n_theta-1][0].a;
                psi_w.b = -aW*psi[i][n_theta-1][0].b;

                Complex psi_e(0,0);
                psi_e.a = -aE*psi[i][1][0].a;
                psi_e.b = -aE*psi[i][1][0].b;

                Complex psi_b(0,0);
                psi_b.a = -aB*psi_p_top[i].a;
                psi_b.b = -aB*psi_p_top[i].b;

                Complex psi_t(0,0);
                psi_t.a = -aT*psi[i][0][1].a;
                psi_t.b = -aT*psi[i][0][1].b;

                Complex ihpsipo_dt(0,0);
                ihpsipo_dt.a = psio[i][0][0].a*h/dt;
                ihpsipo_dt.b = psio[i][0][0].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Complex Fp(0,0);
                Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                double b = h/dt;

                Complex min_Fp(-Fp.a, -Fp.b);
                Complex ib_plus_a(a,b);

                psi[i][0][0] = min_Fp * ib_plus_a;
                psi[i][0][0].a = psi[i][0][0].a/(a*a+b*b);
                psi[i][0][0].b = psi[i][0][0].b/(a*a+b*b);

                /* Last node central top cone nodes */
                dV1 = r_p[i]*r_p[i]*dr;
                dV2 = r_p[i]*r_p[i]*sin(phi_p[0])*sin(phi_p[0])*dtheta;
                dV3 = r_p[i]*r_p[i]*sin(phi_p[0])*dphi;

                aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                aN = r[i]*r[i]*alpha/(dV1*dr);
                aW = alpha/(dV2*dtheta);
                aE = alpha/(dV2*dtheta);
                aB = sin(phi[0])*alpha/(dV3*dphi);
                aT = sin(phi[1])*alpha/(dV3*dphi);

                psi_s.a = -aS*psi[i-1][n_theta-1][0].a;
                psi_s.b = -aS*psi[i-1][n_theta-1][0].b;

                psi_n.a = -aN*psi[i+1][n_theta-1][0].a;
                psi_n.b = -aN*psi[i+1][n_theta-1][0].b;

                psi_w.a = -aW*psi[i][n_theta-1-1][0].a;
                psi_w.b = -aW*psi[i][n_theta-1-1][0].b;

                psi_e.a = -aE*psi[i][0][0].a;
                psi_e.b = -aE*psi[i][0][0].b;

                psi_b.a = -aB*psi_p_top[i].a;
                psi_b.b = -aB*psi_p_top[i].b;

                psi_t.a = -aT*psi[i][n_theta-1][1].a;
                psi_t.b = -aT*psi[i][n_theta-1][1].b;

                ihpsipo_dt.a = psio[i][n_theta-1][0].a*h/dt;
                ihpsipo_dt.b = psio[i][n_theta-1][0].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                b = h/dt;

                min_Fp.a = -Fp.a;
                min_Fp.b = -Fp.b;

                ib_plus_a.a = a;
                ib_plus_a.b = b;

                psi[i][n_theta-1][0] = min_Fp * ib_plus_a;
                psi[i][n_theta-1][0].a = psi[i][n_theta-1][0].a/(a*a+b*b);
                psi[i][n_theta-1][0].b = psi[i][n_theta-1][0].b/(a*a+b*b);
            }

            /* Central bottom cone nodes */
            for(int i = 1; i < n_r - 1; ++i) {
                for(int j = 1; j < n_theta - 1; ++j) {
                    double dV1 = r_p[i]*r_p[i]*dr;
                    double dV2 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*sin(phi_p[n_phi-1])*dtheta;
                    double dV3 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*dphi;

                    double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                    double aN = r[i]*r[i]*alpha/(dV1*dr);
                    double aW = alpha/(dV2*dtheta);
                    double aE = alpha/(dV2*dtheta);
                    double aB = sin(phi[n_phi-1])*alpha/(dV3*dphi);
                    double aT = sin(phi[n_phi])*alpha/(dV3*dphi);

                    Complex psi_s(0,0);
                    psi_s.a = -aS*psi[i-1][j][n_phi-1].a;
                    psi_s.b = -aS*psi[i-1][j][n_phi-1].b;

                    Complex psi_n(0,0);
                    psi_n.a = -aN*psi[i+1][j][n_phi-1].a;
                    psi_n.b = -aN*psi[i+1][j][n_phi-1].b;

                    Complex psi_w(0,0);
                    psi_w.a = -aW*psi[i][j-1][n_phi-1].a;
                    psi_w.b = -aW*psi[i][j-1][n_phi-1].b;

                    Complex psi_e(0,0);
                    psi_e.a = -aE*psi[i][j+1][n_phi-1].a;
                    psi_e.b = -aE*psi[i][j+1][n_phi-1].b;

                    Complex psi_b(0,0);
                    psi_b.a = -aB*psi[i][j][n_phi-2].a;
                    psi_b.b = -aB*psi[i][j][n_phi-2].b;

                    Complex psi_t(0,0);
                    psi_t.a = -aT*psi_p_bottom[i].a;
                    psi_t.b = -aT*psi_p_bottom[i].b;

                    Complex ihpsipo_dt(0,0);
                    ihpsipo_dt.a = psio[i][j][n_phi-1].a*h/dt;
                    ihpsipo_dt.b = psio[i][j][n_phi-1].b*h/dt;
                    ihpsipo_dt = ihpsipo_dt*I;

                    Complex Fp(0,0);
                    Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                    double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                    double b = h/dt;

                    Complex min_Fp(-Fp.a, -Fp.b);
                    Complex ib_plus_a(a,b);

                    psi[i][j][n_phi-1] = min_Fp * ib_plus_a;
                    psi[i][j][n_phi-1].a = psi[i][j][n_phi-1].a/(a*a+b*b);
                    psi[i][j][n_phi-1].b = psi[i][j][n_phi-1].b/(a*a+b*b);
                }

                /* First node central bottom cone nodes */
                double dV1 = r_p[i]*r_p[i]*dr;
                double dV2 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*sin(phi_p[n_phi-1])*dtheta;
                double dV3 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*dphi;

                double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                double aN = r[i]*r[i]*alpha/(dV1*dr);
                double aW = alpha/(dV2*dtheta);
                double aE = alpha/(dV2*dtheta);
                double aB = sin(phi[n_phi-1])*alpha/(dV3*dphi);
                double aT = sin(phi[n_phi])*alpha/(dV3*dphi);

                Complex psi_s(0,0);
                psi_s.a = -aS*psi[i-1][0][n_phi-1].a;
                psi_s.b = -aS*psi[i-1][0][n_phi-1].b;

                Complex psi_n(0,0);
                psi_n.a = -aN*psi[i+1][0][n_phi-1].a;
                psi_n.b = -aN*psi[i+1][0][n_phi-1].b;

                Complex psi_w(0,0);
                psi_w.a = -aW*psi[i][n_theta-1][n_phi-1].a;
                psi_w.b = -aW*psi[i][n_theta-1][n_phi-1].b;

                Complex psi_e(0,0);
                psi_e.a = -aE*psi[i][1][n_phi-1].a;
                psi_e.b = -aE*psi[i][1][n_phi-1].b;

                Complex psi_b(0,0);
                psi_b.a = -aB*psi[i][0][n_phi-2].a;
                psi_b.b = -aB*psi[i][0][n_phi-2].b;

                Complex psi_t(0,0);
                psi_t.a = -aT*psi_p_bottom[i].a;
                psi_t.b = -aT*psi_p_bottom[i].b;

                Complex ihpsipo_dt(0,0);
                ihpsipo_dt.a = psio[i][0][n_phi-1].a*h/dt;
                ihpsipo_dt.b = psio[i][0][n_phi-1].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Complex Fp(0,0);
                Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                double b = h/dt;

                Complex min_Fp(-Fp.a, -Fp.b);
                Complex ib_plus_a(a,b);

                psi[i][0][n_phi-1] = min_Fp * ib_plus_a;
                psi[i][0][n_phi-1].a = psi[i][0][n_phi-1].a/(a*a+b*b);
                psi[i][0][n_phi-1].b = psi[i][0][n_phi-1].b/(a*a+b*b);

                /* Last node central bottom cone nodes */
                dV1 = r_p[i]*r_p[i]*dr;
                dV2 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*sin(phi_p[n_phi-1])*dtheta;
                dV3 = r_p[i]*r_p[i]*sin(phi_p[n_phi-1])*dphi;

                aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                aN = r[i]*r[i]*alpha/(dV1*dr);
                aW = alpha/(dV2*dtheta);
                aE = alpha/(dV2*dtheta);
                aB = sin(phi[n_phi-1])*alpha/(dV3*dphi);
                aT = sin(phi[n_phi])*alpha/(dV3*dphi);

                psi_s.a = -aS*psi[i-1][n_theta-1][n_phi-1].a;
                psi_s.b = -aS*psi[i-1][n_theta-1][n_phi-1].b;

                psi_n.a = -aN*psi[i+1][n_theta-1][n_phi-1].a;
                psi_n.b = -aN*psi[i+1][n_theta-1][n_phi-1].b;

                psi_w.a = -aW*psi[i][n_theta-1-1][n_phi-1].a;
                psi_w.b = -aW*psi[i][n_theta-1-1][n_phi-1].b;

                psi_e.a = -aE*psi[i][0][n_phi-1].a;
                psi_e.b = -aE*psi[i][0][n_phi-1].b;

                psi_b.a = -aB*psi[i][n_theta-1][n_phi-2].a;
                psi_b.b = -aB*psi[i][n_theta-1][n_phi-2].b;

                psi_t.a = -aT*psi_p_bottom[i].a;
                psi_t.b = -aT*psi_p_bottom[i].b;

                ihpsipo_dt.a = psio[i][n_theta-1][n_phi-1].a*h/dt;
                ihpsipo_dt.b = psio[i][n_theta-1][n_phi-1].b*h/dt;
                ihpsipo_dt = ihpsipo_dt*I;

                Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                b = h/dt;

                min_Fp.a = -Fp.a;
                min_Fp.b = -Fp.b;

                ib_plus_a.a = a;
                ib_plus_a.b = b;

                psi[i][n_theta-1][n_phi-1] = min_Fp * ib_plus_a;
                psi[i][n_theta-1][n_phi-1].a = psi[i][n_theta-1][n_phi-1].a/(a*a+b*b);
                psi[i][n_theta-1][n_phi-1].b = psi[i][n_theta-1][n_phi-1].b/(a*a+b*b);
            }

            /* Central nodes */
            for(int i = 1; i < n_r - 1; ++i) {
                for(int k = 1; k < n_phi - 1; ++k) {
                    for(int j = 1; j < n_theta - 1; ++j) {
                        double dV1 = r_p[i]*r_p[i]*dr;
                        double dV2 = r_p[i]*r_p[i]*sin(phi_p[k])*sin(phi_p[k])*dtheta;
                        double dV3 = r_p[i]*r_p[i]*sin(phi_p[k])*dphi;

                        double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                        double aN = r[i]*r[i]*alpha/(dV1*dr);
                        double aW = alpha/(dV2*dtheta);
                        double aE = alpha/(dV2*dtheta);
                        double aB = sin(phi[k])*alpha/(dV3*dphi);
                        double aT = sin(phi[k+1])*alpha/(dV3*dphi);

                        Complex psi_s(0,0);
                        psi_s.a = -aS*psi[i-1][j][k].a;
                        psi_s.b = -aS*psi[i-1][j][k].b;

                        Complex psi_n(0,0);
                        psi_n.a = -aN*psi[i+1][j][k].a;
                        psi_n.b = -aN*psi[i+1][j][k].b;

                        Complex psi_w(0,0);
                        psi_w.a = -aW*psi[i][j-1][k].a;
                        psi_w.b = -aW*psi[i][j-1][k].b;

                        Complex psi_e(0,0);
                        psi_e.a = -aE*psi[i][j+1][k].a;
                        psi_e.b = -aE*psi[i][j+1][k].b;

                        Complex psi_b(0,0);
                        psi_b.a = -aB*psi[i][j][k-1].a;
                        psi_b.b = -aB*psi[i][j][k-1].b;

                        Complex psi_t(0,0);
                        psi_t.a = -aT*psi[i][j][k+1].a;
                        psi_t.b = -aT*psi[i][j][k+1].b;

                        Complex ihpsipo_dt(0,0);
                        ihpsipo_dt.a = psio[i][j][k].a*h/dt;
                        ihpsipo_dt.b = psio[i][j][k].b*h/dt;
                        ihpsipo_dt = ihpsipo_dt*I;

                        Complex Fp(0,0);
                        Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                        double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                        double b = h/dt;

                        Complex min_Fp(-Fp.a, -Fp.b);
                        Complex ib_plus_a(a,b);

                        psi[i][j][k] = min_Fp * ib_plus_a;
                        psi[i][j][k].a = psi[i][j][k].a/(a*a+b*b);
                        psi[i][j][k].b = psi[i][j][k].b/(a*a+b*b);
                    }

                    /* First node central nodes */
                    double dV1 = r_p[i]*r_p[i]*dr;
                    double dV2 = r_p[i]*r_p[i]*sin(phi_p[k])*sin(phi_p[k])*dtheta;
                    double dV3 = r_p[i]*r_p[i]*sin(phi_p[k])*dphi;

                    double aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                    double aN = r[i]*r[i]*alpha/(dV1*dr);
                    double aW = alpha/(dV2*dtheta);
                    double aE = alpha/(dV2*dtheta);
                    double aB = sin(phi[k])*alpha/(dV3*dphi);
                    double aT = sin(phi[k+1])*alpha/(dV3*dphi);

                    Complex psi_s(0,0);
                    psi_s.a = -aS*psi[i-1][0][k].a;
                    psi_s.b = -aS*psi[i-1][0][k].b;

                    Complex psi_n(0,0);
                    psi_n.a = -aN*psi[i+1][0][k].a;
                    psi_n.b = -aN*psi[i+1][0][k].b;

                    Complex psi_w(0,0);
                    psi_w.a = -aW*psi[i][n_theta-1][k].a;
                    psi_w.b = -aW*psi[i][n_theta-1][k].b;

                    Complex psi_e(0,0);
                    psi_e.a = -aE*psi[i][1][k].a;
                    psi_e.b = -aE*psi[i][1][k].b;

                    Complex psi_b(0,0);
                    psi_b.a = -aB*psi[i][0][k-1].a;
                    psi_b.b = -aB*psi[i][0][k-1].b;

                    Complex psi_t(0,0);
                    psi_t.a = -aT*psi[i][0][k+1].a;
                    psi_t.b = -aT*psi[i][0][k+1].b;

                    Complex ihpsipo_dt(0,0);
                    ihpsipo_dt.a = psio[i][0][k].a*h/dt;
                    ihpsipo_dt.b = psio[i][0][k].b*h/dt;
                    ihpsipo_dt = ihpsipo_dt*I;

                    Complex Fp(0,0);
                    Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                    double a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                    double b = h/dt;

                    Complex min_Fp(-Fp.a, -Fp.b);
                    Complex ib_plus_a(a,b);

                    psi[i][0][k] = min_Fp * ib_plus_a;
                    psi[i][0][k].a = psi[i][0][k].a/(a*a+b*b);
                    psi[i][0][k].b = psi[i][0][k].b/(a*a+b*b);

                    /* Last node central nodes */
                    dV1 = r_p[i]*r_p[i]*dr;
                    dV2 = r_p[i]*r_p[i]*sin(phi_p[k])*sin(phi_p[k])*dtheta;
                    dV3 = r_p[i]*r_p[i]*sin(phi_p[k])*dphi;

                    aS = r[i-1]*r[i-1]*alpha/(dV1*dr);
                    aN = r[i]*r[i]*alpha/(dV1*dr);
                    aW = alpha/(dV2*dtheta);
                    aE = alpha/(dV2*dtheta);
                    aB = sin(phi[k])*alpha/(dV3*dphi);
                    aT = sin(phi[k+1])*alpha/(dV3*dphi);

                    psi_s.a = -aS*psi[i-1][n_theta-1][k].a;
                    psi_s.b = -aS*psi[i-1][n_theta-1][k].b;

                    psi_n.a = -aN*psi[i+1][n_theta-1][k].a;
                    psi_n.b = -aN*psi[i+1][n_theta-1][k].b;

                    psi_w.a = -aW*psi[i][n_theta-1-1][k].a;
                    psi_w.b = -aW*psi[i][n_theta-1-1][k].b;

                    psi_e.a = -aE*psi[i][0][k].a;
                    psi_e.b = -aE*psi[i][0][k].b;

                    psi_b.a = -aB*psi[i][n_theta-1][k-1].a;
                    psi_b.b = -aB*psi[i][n_theta-1][k-1].b;

                    psi_t.a = -aT*psi[i][n_theta-1][k+1].a;
                    psi_t.b = -aT*psi[i][n_theta-1][k+1].b;

                    ihpsipo_dt.a = psio[i][n_theta-1][k].a*h/dt;
                    ihpsipo_dt.b = psio[i][n_theta-1][k].b*h/dt;
                    ihpsipo_dt = ihpsipo_dt*I;

                    Fp = psi_s + psi_n + psi_w + psi_e + psi_b + psi_t + ihpsipo_dt;

                    a = aS + aN + aW + aE + aB + aT + V(r_p[i], t);
                    b = h/dt;

                    min_Fp.a = -Fp.a;
                    min_Fp.b = -Fp.b;

                    ib_plus_a.a = a;
                    ib_plus_a.b = b;

                    psi[i][n_theta-1][k] = min_Fp * ib_plus_a;
                    psi[i][n_theta-1][k].a = psi[i][n_theta-1][k].a/(a*a+b*b);
                    psi[i][n_theta-1][k].b = psi[i][n_theta-1][k].b/(a*a+b*b);
                }
            }

            gs_it++;
        }

        /* Normalize psi */
        Complex integral(0,0);
        /* Compute integral of top pole nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            double dV = 2*M_PI/3*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1])*(1 - cos(dphi));
            Complex psi_conj(psi_p_top[i].a, -psi_p_top[i].b);
            Complex elem(0,0);
            elem = psi_p_top[i] * psi_conj;
            elem.a = elem.a*dV;
            elem.b = elem.b*dV;
            integral = integral + elem;
        }

        /* Compute integral of bottom pole nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            double dV = 2*M_PI/3*(r[i]*r[i]*r[i] - r[i-1]*r[i-1]*r[i-1])*(1 - cos(dphi));
            Complex psi_conj(psi_p_bottom[i].a, -psi_p_bottom[i].b);
            Complex elem(0,0);
            elem = psi_p_bottom[i] * psi_conj;
            elem.a = elem.a*dV;
            elem.b = elem.b*dV;
            integral = integral + elem;
        }

        /* Compute integral of central nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    double dV = r_p[i]*r_p[i]*sin(phi_p[k])*dr*dtheta*dphi;
                    Complex psi_conj(psi[i][j][k].a, -psi[i][j][k].b);
                    Complex elem(0,0);
                    elem = psi[i][j][k] * psi_conj;
                    elem.a = elem.a*dV;
                    elem.b = elem.b*dV;
                    integral = integral + elem;
                }
            }
        }

        /* Update psi */
        /* Top pole nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            psi_p_top[i].a = psi_p_top[i].a/sqrt(integral.a);
            psi_p_top[i].b = psi_p_top[i].b/sqrt(integral.a);
        }

        /* Bottom pole nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            psi_p_bottom[i].a = psi_p_bottom[i].a/sqrt(integral.a);
            psi_p_bottom[i].b = psi_p_bottom[i].b/sqrt(integral.a);
        }

        /* Central nodes */
        for(int i = 1; i < n_r - 1; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    psi[i][j][k].a = psi[i][j][k].a/sqrt(integral.a);
                    psi[i][j][k].b = psi[i][j][k].b/sqrt(integral.a);
                    psi_prev[i][j][k].a = psi_prev[i][j][k].a/sqrt(integral.a);
                    psi_prev[i][j][k].b = psi_prev[i][j][k].b/sqrt(integral.a);
                }
            }
        }

        /* Update psi previous timestep */
        for(int i = 1; i < n_r - 1; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    psio[i][j][k].a = psi[i][j][k].a;
                    psio[i][j][k].b = psi[i][j][k].b;
                }
            }
        }

        for(int i = 1; i < n_r - 1; ++i) {
            psi_p_topo[i].a = psi_p_top[i].a;
            psi_p_topo[i].b = psi_p_top[i].b;
            psi_p_bottomo[i].a = psi_p_bottom[i].a;
            psi_p_bottomo[i].b = psi_p_bottom[i].b;
        }

        /* Calculate min and max */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    if(psi[i][j][k].a > max_real) {
                        max_real = psi[i][j][k].a;
                    }
                    if(psi[i][j][k].b > max_im) {
                        max_im = psi[i][j][k].b;
                    }
                    if(psi[i][j][k].a < min_real) {
                        min_real = psi[i][j][k].a;
                    }
                    if(psi[i][j][k].b < min_im) {
                        min_im = psi[i][j][k].b;
                    }
                }
            }
        }

        /* Calculate probability density central nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    Complex psi_conj(psi[i][j][k].a, -psi[i][j][k].b);
                    prob_density[i][j][k] = psi[i][j][k]*psi_conj;
                }
            }
        }

        /* Calculate probability density top pole nodes */
        for(int i = 0; i < n_r; ++i) {
            Complex psi_conj(psi_p_top[i].a, -psi_p_top[i].b);
            prob_density_p_top[i] = psi_p_top[i]*psi_conj;
        }

        /* Calculate probability density bottom pole nodes */
        for(int i = 0; i < n_r; ++i) {
            Complex psi_conj(psi_p_bottom[i].a, -psi_p_bottom[i].b);
            prob_density_p_bottom[i] = psi_p_bottom[i]*psi_conj;
        }

        /* Calculate max probability density */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < 1; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    if(prob_density[i][j][k].a > max_pd) {
                        max_pd = prob_density[i][j][k].a;
                    }
                    if(prob_density_p_top[i].a > max_pd) {
                        max_pd = prob_density_p_top[i].a;
                    }
                    if(prob_density_p_bottom[i].a > max_pd) {
                        max_pd = prob_density_p_bottom[i].a;
                    }
                }
            }
        }

        /* Check convergence */
        double error_real = 0.0;
        double error_im = 0.0;
        for(int i = 1; i < n_r - 1; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    error_real = error_real + fabs(1.0 - psi_prev[i][j][k].a/(psi[i][j][k].a+(1.0e-20)));
                    error_im = error_im + fabs(1.0 - psi_prev[i][j][k].b/(psi[i][j][k].b+(1.0e-20)));
                }
            }
        }
        error_real = error_real / (n_r*n_theta*n_phi);
        error_im = error_im / (n_r*n_theta*n_phi);

        printf("error real: %E, error im: %E\n", error_real, error_im);


        /* Export psi data */
        std::ofstream myfile;
        std::string file_prefix = "psi_vs_t_";
        std::string time_step = std::to_string(timestep);
        std::string file_name = file_prefix + time_step + ".txt";
        myfile.open(file_name);

        /* Export central nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile << r_p[i] << " " << theta_p[j] << " " << phi_p[k] << " " << psi[i][j][k].a << " " << psi[i][j][k].b << "\n";
                }
            }
        }

        /* Export top pole nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile << r_p[i] << " " << theta_p[j] << " " << 0.0 << " " << psi_p_top[i].a << " " << psi_p_top[i].b << "\n";
                }
            }
        }

        /* Export bottom pole nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile << r_p[i] << " " << theta_p[j] << " " << M_PI << " " << psi_p_bottom[i].a << " " << psi_p_bottom[i].b << "\n";
                }
            }
        }

        myfile.close();

        /* Export probability density data */
        std::ofstream myfile_pd;
        file_prefix = "pd_vs_t_";
        time_step = std::to_string(timestep);
        file_name = file_prefix + time_step + ".txt";
        myfile_pd.open(file_name);

        /* Export central nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile_pd << r_p[i] << " " << theta_p[j] << " " << phi_p[k] << " " << prob_density[i][j][k].a << " " << prob_density[i][j][k].b << "\n";
                }
            }
        }

        /* Export top pole nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile_pd << r_p[i] << " " << theta_p[j] << " " << 0.0 << " " << prob_density_p_top[i].a << " " << prob_density_p_top[i].b << "\n";
                }
            }
        }

        /* Export bottom pole nodes */
        for(int i = 0; i < n_r; ++i) {
            for(int j = 0; j < n_theta; ++j) {
                for(int k = 0; k < n_phi; ++k) {
                    myfile_pd << r_p[i] << " " << theta_p[j] << " " << M_PI << " " << prob_density_p_bottom[i].a << " " << prob_density_p_bottom[i].b << "\n";
                }
            }
        }

        myfile_pd.close();

        timestep++;
        t = t + dt;
    }

    /* Export number of timesteps */
    std::ofstream file_timestep;
    std::string file_name = "number_of_timesteps.txt";
    file_timestep.open(file_name);
    file_timestep << nt;
    file_timestep.close();

    /* Export limits */
    std::ofstream file_limits;
    std::string file_name_limits = "limits.txt";
    file_limits.open(file_name_limits);
    file_limits << max_real << " "
                << min_real << " "
                << max_im << " "
                << min_im << " "
                << R << " "
                << max_pd;
    file_limits.close();

    /* Export dimensions */
    std::ofstream file_dims;
    std::string file_name_dims = "dims.txt";
    file_dims.open(file_name_dims);
    file_dims << n_r << " "
              << n_theta << " "
              << n_phi;
    file_dims.close();

    /* Print some results */
    for(int i = 0; i < n_r; ++i) {
        printf("psi_p_top[%i] real: %f, im: %f, psi_p_bottom[%i] real: %f, im: %f\n", i, psi_p_top[i].a, psi_p_top[i].b, i, psi_p_bottom[i].a, psi_p_bottom[i].b);
    }

    /* Print time taken for execution */
    clock_t t_end = clock();
    double sim_time = double (t_end - t_start)/CLOCKS_PER_SEC;
    printf("time taken: %f\n", sim_time);

    printf("done\n");
    return 0;
}
