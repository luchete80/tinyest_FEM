#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

#define m_dim 2
#define m_nodxelem 4
#define m_gp_count 1

double E = 206e9;  // Young's modulus in Pa
double nu = 0.3;   // Poisson's ratio
double rho = 7850.0;
double mat_G;
double K_mod;
int red_int = 0;
double element_length = 1.0; // Length of the element
int axi_symm = 0;            //FALSE: PLAIN STRAIN

double dt = 0.8e-5;
double tf = 0.8e-5;
//double tf = 1.0e-3;
double x[m_nodxelem][m_dim] = {{0.0,0.0},{0.0,0.1},{0.1,0.1},{0.0,0.1}};
double v[m_nodxelem][m_dim];
double a[m_nodxelem][m_dim];
double u[m_nodxelem][m_dim];

// double gauss_points[m_nodxelem][2]={{-0.577350269, -0.577350269},
                                    // {0.577350269, -0.577350269},
                                    // { 0.577350269,  0.577350269},
                                    // {-0.577350269,  0.577350269}};
                                    double gauss_points[1][2] = {{0,0}};
// double gauss_weights[m_gp_count] = {1.,1.,1.,1.};
double gauss_weights[m_gp_count] = {1.};
double detJ[m_gp_count];
double dNdX[m_gp_count][m_dim][m_nodxelem];
double str_rate[m_gp_count][3][3];
double rot_rate[m_gp_count][3][3];
double strain[m_gp_count][3][3];
double tau[m_gp_count][3][3];
double pres[m_gp_count];
double stress[m_gp_count][3][3];
double radius[m_gp_count];

void impose_bc(double vel[m_nodxelem][m_dim], double accel[m_nodxelem][m_dim]) {
    vel[2][1] = vel[3][1] = -1.0;
    vel[0][0] = vel[0][1] = vel[1][1] = 0.0;

    accel[2][1] = accel[3][1] = 0.0;
    accel[0][0] = accel[0][1] = accel[1][1] = 0.0;
}

void shape_functions(double xi, double eta, double N[m_nodxelem], double dNdX_[m_dim][m_nodxelem]) {
    N[0] = (1 - xi) * (1 - eta) / 4;
    N[1] = (1 + xi) * (1 - eta) / 4;
    N[2] = (1 + xi) * (1 + eta) / 4;
    N[3] = (1 - xi) * (1 + eta) / 4;

    dNdX_[0][0] = -(1 - eta) / 4;
    dNdX_[0][1] = (1 - eta) / 4;
    dNdX_[0][2] = (1 + eta) / 4;
    dNdX_[0][3] = -(1 + eta) / 4;

    dNdX_[1][0] = -(1 - xi) / 4;
    dNdX_[1][1] = -(1 + xi) / 4;
    dNdX_[1][2] = (1 + xi) / 4;
    dNdX_[1][3] = (1 - xi) / 4;
    
}

void calc_jacobian(double pos[m_nodxelem][m_dim], double J[m_gp_count][2][2]) {
    double N[m_nodxelem];
    double dNdX_[m_dim][m_nodxelem];
    double xi, eta;
    for (int gp = 0; gp < m_gp_count; gp++) {
        xi = gauss_points[gp][0];
        eta = gauss_points[gp][1];
        shape_functions(xi, eta, N, dNdX_);        
        // for (int i = 0; i < m_dim; i++) {
            // for (int j = 0; j < m_dim; j++) {
                // J[gp][i][j] = 0.0;
                // for (int k = 0; k < m_nodxelem; k++) {
                    // J[gp][i][j] += dNdX_[i][k] * pos[k][j];
                    for (int i = 0; i < m_dim; i++){
                      J[gp][0][i] = -pos[0][i]+pos[1][i]+pos[2][i]-pos[3][i];
                      J[gp][1][i] = -pos[0][i]-pos[1][i]+pos[2][i]+pos[3][i] ;                     
                    }
                      
        // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
        // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)                
                
  
                    //printf("pos %.6e", pos[k][j]);
                    // printf ("J %.6e", J[gp][i][j]);
                //}
            //}
        //}
         printf ("J %.6e %.6e \n %.6e %.6e\n", J[gp][0][0], J[gp][0][1], J[gp][1][0], J[gp][1][1] );
        detJ[gp] = J[gp][0][0] * J[gp][1][1] - J[gp][0][1] * J[gp][1][0];
        printf ("detJ %.6e\n", detJ[gp]);
    }
}

double calc_vol(double detJ[m_nodxelem]) {
    double vol = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        vol += detJ[gp] * gauss_weights[gp];
    }
    return vol;
}

void velocity_gradient_tensor(double dNdX[m_nodxelem][m_dim][m_nodxelem], double vel[m_nodxelem][m_dim], double grad_v[m_nodxelem][m_dim][m_dim]) {
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int I = 0; I < m_dim; I++) {
            for (int J = 0; J < m_dim; J++) {
                for (int k = 0; k < m_nodxelem; k++) {
                    grad_v[gp][I][J] += dNdX[gp][J][k] * vel[k][I];
                }
            }
        }
    }
}

void calc_str_rate(double dNdX[m_nodxelem][m_dim][m_nodxelem], double v[m_nodxelem][m_dim], double str_rate[m_nodxelem][3][3]) {
    double grad_v[m_nodxelem][m_dim][m_dim];
    velocity_gradient_tensor(dNdX, v, grad_v);
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_dim; j++) {
                str_rate[gp][i][j] = 0.5 * (grad_v[gp][i][j] + grad_v[gp][j][i]);
            }
        }
    }
}

void calc_strain(double str_rate[m_nodxelem][3][3], double dt, double strain[m_nodxelem][3][3]) {
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                strain[gp][i][j] = dt * str_rate[gp][i][j];
            }
        }
    }
}

void calc_pressure(double K_, double dstr[m_nodxelem][3][3], double stress[m_nodxelem][3][3], double pres[m_nodxelem]) {
    double pi_ = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            pi_ += dstr[gp][i][i];
        }
    }
    pi_ = -pi_ / m_nodxelem;
    for (int gp = 0; gp < m_gp_count; gp++) {
        pres[gp] = -1.0 / 3.0 * (stress[gp][0][0] + stress[gp][1][1] + stress[gp][2][2]) + K_ * pi_;
    }
}

double dev(double t[3][3]) {
    double d= 1.0 / 3.0 * (t[0][0] + t[1][1] + t[2][2]);
    
    return d;
}

void calc_stress2(double str_rate[m_nodxelem][3][3], double rot_rate[m_nodxelem][3][3], double tau[m_nodxelem][3][3], double p[m_nodxelem], double dt, double stress[m_nodxelem][3][3]) {
    double srt[m_nodxelem][3][3];
    double rs[m_nodxelem][3][3];

    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                srt[gp][i][j] = tau[gp][i][j] * rot_rate[gp][i][j];
                rs[gp][i][j] = rot_rate[gp][i][j] * tau[gp][i][j];
                tau[gp][i][j] += dt * (2.0 * mat_G * (dev(str_rate[gp])) + rs[gp][i][j] + srt[gp][i][j]);
                stress[gp][i][j] = tau[gp][i][j] - p[gp] * (i == j);
            }
        }
    }
}

void calc_forces(double stress[m_nodxelem][3][3], double dNdX[m_nodxelem][m_dim][m_nodxelem], double J[m_nodxelem][2][2], double forces[m_nodxelem][m_dim]) {
    double B[m_dim][m_nodxelem];

    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_nodxelem; j++) {
                B[i][j] = dNdX[gp][i][j];
            }
        }
        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                forces[gp][i] += B[j][i] * (stress[gp][0][0] + stress[gp][1][1]) * J[gp][0][0] * J[gp][1][1] * gauss_weights[gp];
            }
        }
    }
}

int main() {
    printf("Begin..\n");
    mat_G = E / (2.0 * (1 + nu));
    K_mod = E / (3.0 * (1.0 - 2.0 * nu));
    
    printf("Imposing bc..\n");
    impose_bc(v, a);
    printf("Done");
    double u_tot[m_nodxelem][m_dim];
    double prev_a[m_nodxelem][m_dim];

    double J[m_gp_count][2][2];
    calc_jacobian(x, J);

    double vol_0 = calc_vol(detJ);
    cout << "vol 0 "<<vol_0<<endl;
    double nod_mass = vol_0 * rho / m_nodxelem;

    double rho_b = 0.8182;
    double alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
    double beta = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
    double gamma = 1.5 - alpha;

    double t = 0.0;
    while (t < tf) {
      printf ("Time: %.6e\n", t);

        // PREDICTION PHASE
        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u[i][j] = dt * (v[i][j] + (0.5 - beta) * dt * prev_a[i][j]);
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                v[i][j] += (1.0 - gamma) * dt * prev_a[i][j];
            }
        }

        impose_bc(v, a);

        calc_jacobian(x, J);

        calc_str_rate(dNdX, v, str_rate);
        double str_inc[m_nodxelem][3][3];
        calc_strain(rot_rate, dt, str_inc);

        calc_pressure(K_mod, str_inc, stress, pres);

        calc_stress2(str_rate, rot_rate, tau, pres, dt, stress);

        calc_forces(stress, dNdX, J, a);

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                a[i][j] = -a[i][j] / nod_mass - alpha * prev_a[i][j];
                a[i][j] /= (1.0 - alpha);
                v[i][j] += gamma * dt * a[i][j];
            }
        }

        impose_bc(v, a);

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u[i][j] += beta * dt * dt * a[i][j];
                x[i][j] += u[i][j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                prev_a[i][j] = a[i][j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u_tot[i][j] += u[i][j];
            }
        }

        t += dt;
    }

    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", u_tot[i][j]);
        }
      printf("\n");
    }
    return 0;
}
