import numpy as np

# Define material properties
E = 200e9  # Young's modulus in Pa
nu = 0.3   # Poisson's ratio
yield_stress = 250e6  # Yield stress in Pa
H = 5e9  # Hardening modulus in Pa

# Define element properties
num_nodes_element = 4  # Number of nodes per element
element_length = 1.0   # Length of the element

# Define shape functions and their derivatives for 2D quadrilateral element
def shape_functions(xi, eta):
    N = np.array([(1-xi)*(1-eta)/4,
                  (1+xi)*(1-eta)/4,
                  (1+xi)*(1+eta)/4,
                  (1-xi)*(1+eta)/4])
    dN_dxi = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
    dN_deta = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
    return N, dN_dxi, dN_deta

# Define material matrix for plane stress
def material_matrix():
    C = E / (1 - nu**2) * np.array([[1, nu, 0],
                                     [nu, 1, 0],
                                     [0, 0, (1 - nu) / 2]])
    return C

# Gauss quadrature points and weights
gauss_points = np.array([[-0.577350269, -0.577350269],
                         [ 0.577350269, -0.577350269],
                         [ 0.577350269,  0.577350269],
                         [-0.577350269,  0.577350269]])
gauss_weights = np.array([1, 1, 1, 1])


# Finite element strain rate calculation
def calculate_strain_rate(velocities):
    strain_rate = np.zeros((3, 3))
    for gp in range(len(gauss_points)):
        xi, eta = gauss_points[gp]
        weight = gauss_weights[gp]
        N, dN_dxi, dN_deta = shape_functions(xi, eta)
        J = np.dot(np.array([dN_dxi, dN_deta]).T, velocities)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        B = np.zeros((3, 8))
        for i in range(num_nodes_element):
            B[0, 2*i] = dN_dxi[i]
            B[1, 2*i+1] = dN_deta[i]
            B[2, 2*i] = dN_deta[i]
            B[2, 2*i+1] = dN_dxi[i]
        strain_rate += np.dot(B.T, np.dot(invJ.T, np.dot(B, velocities))) * detJ * weight
    strain_rate *= 0.5
    return strain_rate

# Leapfrog explicit integration to calculate strain from strain rates and stresses
def leapfrog_integration(strain_rate, stresses, dt):
    strain = np.zeros_like(strain_rate)
    strain_prev = np.zeros_like(strain_rate)
    strain_next = np.zeros_like(strain_rate)
    num_steps = 10  # Number of integration steps
    for i in range(num_steps):
        strain_next = strain_prev + 0.5 * dt * (3 * strain_rate - strain_prev)
        stresses_next = np.dot(material_matrix(), strain_next)
        stresses_next, _ = J2_plasticity(stresses_next, dt)
        strain = strain_prev + 0.5 * dt * (stresses + stresses_next)
        strain_prev = strain_next
        stresses = stresses_next
    return strain

# J2 plasticity function to update stresses
def j2_plasticity(strain_rate, stresses, dt):
    dev_strain_rate = strain_rate - np.trace(strain_rate) / 3 * np.eye(3)  # Deviatoric strain rate
    dev_stresses = stresses - np.trace(stresses) / 3 * np.eye(3)  # Deviatoric stresses
    dev_strain_rate_norm = np.sqrt(2 / 3 * np.sum(dev_strain_rate ** 2))  # Norm of deviatoric strain rate
    if dev_strain_rate_norm > 0:
        dev_stresses += 2 * E / (2 * E + 3 * yield_stress) * dev_strain_rate_norm * dev_strain_rate * dt
        dev_stresses /= np.sqrt(1 + 3 / (2 * E) * np.sum(dev_stresses ** 2))
    return dev_stresses + np.trace(stresses) / 3 * np.eye(3)  # Total stresses

# Calculate stresses from strain rate
def calculate_stresses(strain_rate, stresses, dt):
    updated_stresses = j2_plasticity(strain_rate, stresses, dt)
    return updated_stresses

# Calculate element forces from stresses
def calculate_element_forces(stresses):
    element_forces = np.dot(np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]]), stresses)
    return element_forces

# Calculate acceleration from element forces
def calculate_acceleration(element_forces):
    # Assuming constant mass for simplicity
    mass = 1.0
    acceleration = element_forces / mass
    return acceleration

# Calculate velocity from acceleration
def calculate_velocity(acceleration, dt):
    velocity = acceleration * dt
    return velocity

# Example usage
velocities = np.array([[1, 0], [2, 0], [2, 1], [1, 1]])  # Example velocities at nodes
strain_rate = calculate_strain_rate(velocities)
initial_stresses = np.zeros((3, 3))  # Initial stresses
dt = 0.01  # Time step
updated_stresses = calculate_stresses(strain_rate, initial_stresses, dt)
print("Updated Stresses:")
print(updated_stresses)

print("Strain:")
print(strain)
print("Element Forces:")
print(element_forces)
print("Acceleration:")
print(acceleration)
print("Velocity:")
print(velocity)
