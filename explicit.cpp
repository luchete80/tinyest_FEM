#include <iostream>
#include <vector>
#include <cmath>

// Material properties
const double density = 1.0;
const double yield_stress = 0.1; // Yield stress
const double hourglass_coefficient = 0.1; // Hourglass coefficient

// Velocity gradient tensor computation using shape function derivatives
std::vector<std::vector<double>> velocity_gradient_tensor(const std::vector<std::vector<double>>& dNdX,
                                                          const std::vector<double>& nodal_velocities) {
    std::vector<std::vector<double>> grad_v(3, std::vector<double>(3, 0.0));
    for (int I = 0; I < 3; ++I) {
        for (int J = 0; J < 3; ++J) {
            for (int a = 0; a < 8; ++a) {
                grad_v[I][J] += dNdX[I][a] * nodal_velocities[J * 8 + a];
            }
        }
    }
    return grad_v;
}

// Strain rate tensor computation from velocity gradient tensor
std::vector<std::vector<double>> strain_rate_tensor(const std::vector<std::vector<double>>& grad_v) {
    // Strain rate tensor (symmetric)
    std::vector<std::vector<double>> D(3, std::vector<double>(3, 0.0));
    D[0][0] = grad_v[0][0]; // du/dx
    D[1][1] = grad_v[1][1]; // dv/dy
    D[2][2] = grad_v[2][2]; // dw/dz
    D[0][1] = 0.5 * (grad_v[0][1] + grad_v[1][0]); // (du/dy + dv/dx) / 2
    D[0][2] = 0.5 * (grad_v[0][2] + grad_v[2][0]); // (du/dz + dw/dx) / 2
    D[1][2] = 0.5 * (grad_v[1][2] + grad_v[2][1]); // (dv/dz + dw/dy) / 2
    D[1][0] = D[0][1];
    D[2][0] = D[0][2];
    D[2][1] = D[1][2];
    return D;
}

// Von Mises stress computation
double von_mises_stress(const std::vector<double>& sigma) {
    double s1 = sigma[0], s2 = sigma[1], s3 = sigma[2], s4 = sigma[3], s5 = sigma[4], s6 = sigma[5];
    return sqrt(0.5 * ((s1 - s2) * (s1 - s2) + (s2 - s3) * (s2 - s3) + (s3 - s1) * (s3 - s1) + 6 * (s4 * s4 + s5 * s5 + s6 * s6)));
}

// Plasticity update using von Mises criterion with Jaumann rate
void plasticity_update(std::vector<double>& sigma, std::vector<double>& epsilon,
                       const std::vector<double>& strain_rate, double dt) {
    // Compute von Mises stress
    double von_mises_stress_old = von_mises_stress(sigma);

    // Compute Jaumann rate of stress
    std::vector<double> sigma_increment(6, 0.0);
    for (int i = 0; i < 6; ++i) {
        sigma_increment[i] = 2 * yield_stress * strain_rate[i] * dt;
    }

    // Update stress using Jaumann rate
    for (int i = 0; i < 6; ++i) {
        sigma[i] += sigma_increment[i];
    }

    // Compute von Mises stress after update
    double von_mises_stress_new = von_mises_stress(sigma);

    // Check for yielding and adjust plastic strain if necessary
    if (von_mises_stress_new > yield_stress) {
        // Compute plastic strain increment
        double factor = yield_stress / von_mises_stress_old;
        for (int i = 0; i < 6; ++i) {
            epsilon[i] += factor * strain_rate[i] * dt;
        }
    }
}

// Hourglass control force
std::vector<double> hourglass_force(const std::vector<double>& u, const std::vector<double>& strain_rate) {
    std::vector<double> F(8, 0.0);
    for (int i = 0; i < 8; ++i) {
        F[i] = -hourglass_coefficient * strain_rate[i];
    }
    return F;
}
´
// von Mises yield criterion
double von_mises(double sigma[3]) {
    return sqrt(0.5 * (pow(sigma[0] - sigma[1], 2) + pow(sigma[1] - sigma[2], 2) + pow(sigma[2] - sigma[0], 2))
                 + 3 * pow(sigma[3], 2));
}

// Plasticity update using Jaumann rate
void plasticity_update2(std::vector<double>& sigma, std::vector<double>& epsilon, const std::vector<double>& stress_increment) {
    double s[3] = {sigma[0], sigma[1], sigma[2]};
    double s_dev[3] = {s[0] - (s[0] + s[1] + s[2]) / 3.0, s[1] - (s[0] + s[1] + s[2]) / 3.0, s[2] - (s[0] + s[1] + s[2]) / 3.0};
    double J2 = 0.5 * (pow(s_dev[0], 2) + pow(s_dev[1], 2) + pow(s_dev[2], 2));
    double eq_stress = sqrt(3 * J2);

    // Check if yielding
    if (eq_stress >= yield_stress) {
        // Plastic correction using Jaumann rate
        double factor = yield_stress / eq_stress;
        for (int i = 0; i < 3; ++i) {
            sigma[i] += factor * stress_increment[i];
            epsilon[i] += factor * stress_increment[i] / Young_modulus;
        }
    } else {
        // Elastic update
        for (int i = 0; i < 3; ++i) {
            sigma[i] += stress_increment[i];
            epsilon[i] += stress_increment[i] / Young_modulus;
        }
    }
}


// Explicit finite element solver with hourglass control and von Mises plasticity using Jaumann rate
void explicit_finite_element_solver(std::vector<double>& u, std::vector<double>& v,
                                    const std::vector<std::vector<double>>& M,
                                    const std::vector<std::vector<double>>& dNdX,
                                    const std::vector<double>& nodal_velocities) {
    // Initialize stress, strain, and rotation
    std::vector<double> sigma(6, 0.0); // Stress
    std::vector<double> epsilon(6, 0.0); // Strain

    // Time integration
    double dt = 0.01; // Time step (adjust as needed)
    int num_steps = 500; // Number of time steps (adjust as needed)
    for (int step = 0; step < num_steps; ++step) {
        // Compute velocity gradient tensor
        std::vector<std::vector<double>> grad_v = velocity_gradient_tensor(dNdX, nodal_velocities);

        // Compute strain rate tensor
        std::vector<std::vector<double>> D = strain_rate_tensor(grad_v);

      std::vector<double>sr(6,0.0);
        // Update stress, strain using plasticity with von Mises criterion and Jaumann rate
        plasticity_update(sigma, epsilon, D, dt);

        // Hourglass control force
        std::vector<double> F_hourglass = hourglass_force(u, epsilon);

        // Update velocity
        for (int i = 0; i < 8; ++i) {
            v[i] += dt * (F_hourglass[i] / M[i][i]);
        }

        // Update displacement
        for (int i = 0; i < 8; ++i) {
            u[i] += dt * v[i];
        }
    }
}

int main() {
    // Define nodal velocities (dummy data for demonstration)
    std::vector<double> nodal_velocities(24, 0.1); // Adjust as needed based on your problem

    // Define shape function derivatives (dummy data for demonstration)
    std::vector<std::vector<double>> dNdX(3, std::vector<double>(8, 0.0));
    // Assign shape function derivatives here

    // Solve the problem using explicit finite element method
    std::vector<double> u(8, 0.0); // Displacement vector
    std::vector<double> v(8, 0.0); // Velocity vector
    std::vector<std::vector<double>> M(8, std::vector<double>(8, 1.0)); // Mass matrix
    explicit_finite_element_solver(u, v, M, dNdX, nodal_velocities); // Provide appropriate parameters

    return 0;
}