import numpy as np

# Material properties
density = 1.0
yield_stress = 0.1
hourglass_coefficient = 0.1

# Velocity gradient tensor computation using shape function derivatives
def velocity_gradient_tensor(dNdX, nodal_velocities):
    grad_v = np.zeros((3, 3))
    for I in range(3):
        for J in range(3):
            for a in range(8):
                grad_v[I, J] += dNdX[I, a] * nodal_velocities[J * 8 + a]
    return grad_v

# Deformation gradient tensor computation from displacement gradient
def deformation_gradient_tensor(dNdx, nodal_displacements):
    F = np.zeros((3, 3))
    for I in range(3):
        for J in range(3):
            for a in range(8):
                F[I, J] += nodal_displacements[J * 3 + I] * dNdx[J, a]
    return F

# Strain rate tensor computation from deformation gradient tensor
def strain_rate_tensor(F):
    D = 0.5 * (F - F.T)
    return D

# Von Mises stress computation
def von_mises_stress(sigma):
    s1, s2, s3, s4, s5, s6 = sigma
    return np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2 + 6 * (s4 ** 2 + s5 ** 2 + s6 ** 2)))

# Plasticity update using von Mises criterion with Jaumann rate
def plasticity_update(sigma, epsilon, strain_rate, dt):
    von_mises_stress_old = von_mises_stress(sigma)
    sigma_increment = 2 * yield_stress * strain_rate * dt
    sigma += sigma_increment
    von_mises_stress_new = von_mises_stress(sigma)
    if von_mises_stress_new > yield_stress:
        factor = yield_stress / von_mises_stress_old
        epsilon += factor * strain_rate * dt

# Hourglass control force
def hourglass_force(u, strain_rate):
    F = np.zeros(8)
    for i in range(8):
        F[i] = -hourglass_coefficient * strain_rate[i]
    return F
    
def strain_rate_tensor_from_velocity_gradient(grad_v, dNdX):
    strain_rate = np.zeros((3, 3))
    for I in range(3):
        for J in range(3):
            for a in range(8):
                strain_rate[I, J] += grad_v[I, J] * dNdX[J, a]
    return strain_rate    

# Explicit finite element solver with hourglass control and von Mises plasticity using Jaumann rate
def explicit_finite_element_solver(u, v, M, dNdx, nodal_displacements, nodal_velocities):
    sigma = np.zeros(6)
    epsilon = np.zeros(6)
    dt = 0.01  # Time step
    num_steps = 500  # Number of time steps
    for step in range(num_steps):
        # Compute deformation gradient tensor
        F = deformation_gradient_tensor(dNdx, nodal_displacements)

        # Compute strain rate tensor
        D = strain_rate_tensor(F)

        # Update stress and strain using plasticity with von Mises criterion and Jaumann rate
        #plasticity_update(sigma, epsilon, D, dt)

        # Hourglass control force
        # F_hourglass = hourglass_force(u, epsilon)

        # Update velocity
        # v += dt * np.linalg.solve(M, F_hourglass)

        # Update displacement
        u += dt * v

        # Print or save nodal displacements
        if step % 10 == 0:  # Print every 10 steps
            print("Time step:", step)
            print("Nodal Displacements:")
            for i in range(len(nodal_displacements)):
                print("Node {}: {}".format(i, nodal_displacements[i]))
            print("\n")
            
# Main function
def main():
    # Define nodal velocities (dummy data for demonstration)
    nodal_velocities = np.full(24, 0.1)

    # Define shape function derivatives (dummy data for demonstration)
    dNdx = np.zeros((3, 8))
    # Assign shape function derivatives here

    # Solve the problem using explicit finite element method
    u = np.zeros(8)  # Displacement vector
    v = np.zeros(8)  # Velocity vector
    v[1] = -1.
    M = np.eye(8)  # Mass matrix (dummy identity matrix)
    nodal_displacements = np.zeros(24)  # Nodal displacements (dummy initial values)
    explicit_finite_element_solver(u, v, M, dNdx, nodal_displacements, nodal_velocities)

if __name__ == "__main__":
    main()
