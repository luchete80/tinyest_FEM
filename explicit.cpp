#include <iostream>
#include <vector>
#include <cmath>

// Material properties
const double density = 1.0;
const double Young_modulus = 1.0;
const double Poisson_ratio = 0.3;
const double yield_stress = 0.1; // Yield stress
const double hourglass_coefficient = 0.1; // Hourglass coefficient

// Mass matrix for 3D solid element (hexahedral)
std::vector<std::vector<double>> mass_matrix(double element_volume, double density) {
    std::vector<std::vector<double>> M(8, std::vector<double>(8, 0.0));
    double mass = density * element_volume / 8.0;
    for (int i = 0; i < 8; ++i)
        M[i][i] = mass;
    return M;
}

// von Mises yield criterion
double von_mises(const std::vector<double>& sigma) {
    double s1 = sigma[0], s2 = sigma[1], s3 = sigma[2], s4 = sigma[3], s5 = sigma[4], s6 = sigma[5];
    return sqrt(0.5 * ((s1 - s2) * (s1 - s2) + (s2 - s3) * (s2 - s3) + (s3 - s1) * (s3 - s1) + 6 * (s4 * s4 + s5 * s5 + s6 * s6)));
}

// Plasticity update using Jaumann rate with strain rate and rotation rate
void plasticity_update(std::vector<double>& sigma, std::vector<double>& epsilon, std::vector<double>& omega,
                       const std::vector<double>& strain_rate_increment, const std::vector<double>& rotation_rate_increment) {
    double d1 = strain_rate_increment[0], d2 = strain_rate_increment[1], d3 = strain_rate_increment[2],
           d4 = strain_rate_increment[3], d5 = strain_rate_increment[4], d6 = strain_rate_increment[5];
    double eq_strain_rate = sqrt(0.5 * (d1 * d1 + d2 * d2 + d3 * d3 + 2 * (d4 * d4 + d5 * d5 + d6 * d6)));

    double r1 = rotation_rate_increment[0], r2 = rotation_rate_increment[1], r3 = rotation_rate_increment[2];
    double eq_rotation_rate = sqrt(0.5 * (r1 * r1 + r2 * r2 + r3 * r3));

    // Check if yielding
    if (eq_strain_rate >= yield_stress || eq_rotation_rate >= yield_stress) {
        // Plastic correction using Jaumann rate
        double factor_strain = yield_stress / eq_strain_rate;
        double factor_rotation = yield_stress / eq_rotation_rate;
        for (int i = 0; i < 6; ++i) {
            epsilon[i] += factor_strain * strain_rate_increment[i];
        }
        for (int i = 0; i < 3; ++i) {
            omega[i] += factor_rotation * rotation_rate_increment[i];
        }
    } else {
        // Elastic update
        for (int i = 0; i < 6; ++i) {
            epsilon[i] += strain_rate_increment[i];
        }
        for (int i = 0; i < 3; ++i) {
            omega[i] += rotation_rate_increment[i];
        }
    }

    // Compute rotation tensor
    std::vector<std::vector<double>> omega_hat = {{0, -omega[2], omega[1]},
                                                  {omega[2], 0, -omega[0]},
                                                  {-omega[1], omega[0], 0}};

    // Compute strain rate tensor
    std::vector<std::vector<double>> D = {{d1, (d4 + d5) / 2, (d6 + d2) / 2},
                                          {(d4 + d5) / 2, d3, (d4 + d5) / 2},
                                          {(d6 + d2) / 2, (d4 + d5) / 2, d6}};

    // Compute updated stress
    std::vector<double> sigma_increment(6, 0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                sigma_increment[i] += 2 * Young_modulus * D[j][k] * omega_hat[j][i];
            }
        }
    }

    // Update stress
    for (int i = 0; i < 6; ++i) {
        sigma[i] += sigma_increment[i];
    }
}

// Hourglass control force
std::vector<double> hourglass_force(const std::vector<double>& u, const std::vector<double>& strain_rate,
                                    const std::vector<double>& rotation_rate) {
    std::vector<double> F(8, 0.0);
    for (int i = 0; i < 8; ++i) {
        F[i] = -hourglass_coefficient * (strain_rate[i] + rotation_rate[i]);
    }
    return F;
}

// Time integration using explicit central difference method with plasticity and hourglass control
void explicit_finite_element_solver(std::vector<double>& u, std::vector<double>& v, const std::vector<std::vector<double>>& M,
                                    const std::vector<double>& strain_rate, const std::vector<double>& rotation_rate) {
    // Initialize stress, strain, and rotation
    std::vector<double> sigma(6, 0.0); // Stress
    std::vector<double> epsilon(6, 0.0); // Strain
    std::vector<double> omega(3, 0.0); // Rotation

    // Time integration
    double dt = 0.01; // Time step (adjust as needed)
    int num_steps = 500; // Number of time steps (adjust as needed)
    for (int step = 0; step < num_steps; ++step) {
        // Compute strain rate and rotation rate increments (assume constant increments for simplicity)
        std::vector<double> strain_rate_increment(6, 0.0);
        std::vector<double> rotation_rate_increment(3, 0.0);
        for (int i = 0; i < 6; ++i)
            strain_rate_increment[i] = strain_rate[i] * dt;
        for (int i = 0; i < 3; ++i)
            rotation_rate_increment[i] = rotation_rate[i] * dt;

        // Update stress, strain, and rotation using plasticity
        plasticity_update(sigma, epsilon, omega, strain_rate_increment, rotation_rate_increment);

        // Hourglass control force
        std::vector<double> F_hourglass = hourglass_force(u, strain_rate, rotation_rate);

        // Update acceleration
        std::vector<double> a(8, 0.0);
        for (int i = 0; i < 8; ++i) {
            a[i] = (F_hourglass[i] / M[i][i]);
        }

        // Update velocity
        for (int i = 0; i < 8; ++i) {
            v[i] += dt * a[i];
        }

        // Update displacement
        for (int i = 0; i < 8; ++i) {
            u[i] += dt * v[i];
        }
    }
}

int main() {
    // Initial conditions
    std::vector<double> u(8, 0.0); // Displacement vector
    std::vector<double> v(8, 0.0); // Velocity vector

    // Compute element volume (for simplicity, assuming a unit volume for each element)
    double element_volume = 1.0;

    // Assemble global mass matrix (for simplicity, assuming uniform density)
    std::vector<std::vector<double>> M = mass_matrix(element_volume, density);

    // Define constant strain rate and rotation rate (for simplicity)
    std::vector<double> strain_rate = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; // Adjust as needed
    std::vector<double> rotation_rate = {0.1, 0.1, 0.1}; // Adjust as needed

    // Solve the problem using explicit finite element method
    explicit_finite_element_solver(u, v, M, strain_rate, rotation_rate);

    // Output the displacement
    for (int i = 0; i < 8; ++i) {
        std::cout << "Node " << i << ": " << u[i] << std::endl;
    }

    return 0;
}