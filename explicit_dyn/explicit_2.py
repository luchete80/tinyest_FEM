#include <iostream>
#include <vector>

// Material properties
const double E = 1.0;  // Young's modulus
const double nu = 0.3;  // Poisson's ratio
const double rho = 1.0;  // Density
const double yield_stress = 0.1;  // Yield stress

// Hexahedral element shape functions and their derivatives with respect to local coordinates
const double N[3][8] = {{1.0/8, -1.0/8, -1.0/8, 1.0/8, 1.0/8, -1.0/8, -1.0/8, 1.0/8},
                        {1.0/8, 1.0/8, -1.0/8, -1.0/8, 1.0/8, 1.0/8, -1.0/8, -1.0/8},
                        {1.0/8, 1.0/8, 1.0/8, 1.0/8, -1.0/8, -1.0/8, -1.0/8, -1.0/8}};

const double dN_dx[3][8] = {{-1.0/8, 1.0/8, 1.0/8, -1.0/8, -1.0/8, 1.0/8, 1.0/8, -1.0/8},
                             {-1.0/8, -1.0/8, 1.0/8, 1.0/8, -1.0/8, -1.0/8, 1.0/8, 1.0/8},
                             {-1.0/8, -1.0/8, -1.0/8, -1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8}};

// Integration point coordinates
const double xi = 0.5773502692;  // Gauss point coordinates
const double eta = 0.5773502692;
const double zeta = 0.5773502692;

// Initial conditions
const std::vector<double> initial_displacement(3, 0.0);
const std::vector<double> initial_velocity(3, 0.0);
const std::vector<double> initial_acceleration = {0.0, 0.0, -9.81};  // Gravitational acceleration
const std::vector<std::vector<double>> initial_stress(3, std::vector<double>(3, 0.0));  // Initial stress tensor

// Time parameters
const double total_time = 10.0;
const double dt = 0.01;
const int num_steps = total_time / dt;

// Compute Jacobian determinant at a single Gauss point
double computeJacobian(const std::vector<std::vector<double>>& coords) {
    double J = 0.0;
    for (int i = 0; i < 8; ++i) {
        double dN_dxi = dN_dx[0][i];
        double dN_deta = dN_dx[1][i];
        double dN_dzeta = dN_dx[2][i];

        J += dN_dxi * coords[0][i] + dN_deta * coords[1][i] + dN_dzeta * coords[2][i];
    }
    return J;
}

int main() {
    // Leapfrog time integration
    std::vector<double> displacement = initial_displacement;
    std::vector<double> velocity = initial_velocity;
    std::vector<std::vector<double>> stress = initial_stress;

    for (int step = 0; step < num_steps; ++step) {
        // Compute velocity gradients
        std::vector<std::vector<double>> grad_v(3, std::vector<double>(8, 0.0));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 8; ++j) {
                grad_v[i][j] = dN_dx[i][j] * velocity[0] + dN_dx[i][j] * velocity[1] + dN_dx[i][j] * velocity[2];
            }
        }

        // Compute strain rate tensor
        std::vector<std::vector<double>> strain_rate(3, std::vector<double>(3, 0.0));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 8; ++k) {
                    strain_rate[i][j] += 0.5 * (grad_v[i][k] + grad_v[j][k]);
                }
            }
        }

        // Compute nodal forces
        std::vector<double> nodal_forces(8, 0.0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 8; ++j) {
                for (int k = 0; k < 3; ++k) {
                    nodal_forces[j] += N[i][j] * stress[i][k] * dN_dx[i][j];
                }
            }
        }

        // Compute acceleration
        std::vector<double> acceleration(3, 0.0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 8; ++j) {
                acceleration[i] += nodal_forces[j] / element_mass;
            }
        }

        // Compute trial stress
        std::vector<std::vector<double>> trial_stress(3, std::vector<double>(3, 0.0));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                trial_stress[i][j] = stress[i][j] + 2 * dt * E / (1 + nu) * strain_rate[i][j];
            }
        }

        // Compute deviatoric trial stress and von Mises stress
        double dev_trial_stress_norm = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                trial_stress[i][j] -= (trial_stress[0][0] + trial_stress[1][1] + trial_stress[2][2]) / 3;
                dev_trial_stress_norm += trial_stress[i][j] * trial_stress[i][j];
            }
        }
        dev_trial_stress_norm = std::sqrt(dev_trial_stress_norm);

        // Check for yielding
        if (dev_trial_stress_norm > yield_stress) {
            // Plastic correction
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    trial_stress[i][j] *= yield_stress / dev_trial_stress_norm;
                }
            }
        }

        // Update stress using Jaumann rate
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                stress[i][j] = trial_stress[i][j] + 2 * dt * E / (1 + nu) * strain_rate[i][j];
            }
        }

        // Update velocity and displacement using leapfrog method
        for (int i = 0; i < 3; ++i) {
            velocity[i] += dt * acceleration[i];
            displacement[i] += dt * velocity[i];
        }

        // Hourglass correction
        double J = computeJacobian({{0.125, 0.125, 0.125, 0.125, -0.125, -0.125, -0.125, -0.125},
                                    {0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125},
                                    {0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125}});
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 3; ++j) {
                acceleration[j] -= (0.25 * (1.0 - J) * velocity[j]) / dt;
            }
        }

        // Output results
        std::cout << "Step: " << step << ", Time: " << step * dt << ", Displacement: [";
        for (int i = 0; i < 3; ++i) {
            std::cout << displacement[i];
            if (i < 2) std::cout << ", ";
        }
        std::cout << "], Velocity: [";
        for (int i = 0; i < 3; ++i) {
            std::cout << velocity[i];
            if (i < 2) std::cout << ", ";
        }
        std::cout << "], Stress: [";
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << stress[i][j];
                if (i < 2 || j < 2) std::cout << ", ";
            }
        }
        std::cout << "]\n";
    }

    return 0;
}
