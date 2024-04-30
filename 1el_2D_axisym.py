
import numpy as np

# Define material properties
E = 206e9  # Young's modulus in Pa
nu = 0.3   # Poisson's ratio
rho = 7850.0

# Define element properties
m_dim = 3  # Dimensionality for axisymmetric analysis
m_nodxelem = 4  # Number of nodes per element

# Define nodal velocity (dummy data for demonstration)
vel = np.full(m_dim * m_nodxelem, 0.1)
vel[5] = vel[7] = -1.0

dt = 0.1e-5
tf = 1.0e-3  # Final time
x = np.array([[0., 0.,0.0], [0.1, 0.,0.], [0.1, 0.1,0.], [0., 0.1,0.]])  # Initial nodal positions
v = np.array([[0, 0], [0, 0], [0, -1], [0, -1]])  # Example velocities at nodes

# Gauss quadrature points and weights
gauss_points = np.array([[0.3399810435848563, 0.8611363115940526],
                         [0.3399810435848563, -0.8611363115940526],
                         [0.8611363115940526, 0.3399810435848563],
                         [0.8611363115940526, -0.3399810435848563]])
gauss_weights = np.array([0.6521451548625461, 0.6521451548625461,
                          0.3478548451374538, 0.3478548451374538])
m_gp_count = 4

detJ = np.zeros((m_gp_count))
dNdX = np.zeros((m_gp_count, m_dim, m_nodxelem))
dNdrs = np.zeros((m_gp_count, m_dim, m_nodxelem))

strain = np.zeros((m_gp_count, m_dim, m_dim))


def impose_bc(vel, accel):
    vel[2, 1] = vel[3, 1] = -1.0
    vel[0, :] = vel[1, 1] = 0.0

    accel[2, 1] = accel[3, 1] = 0.0
    accel[0, :] = accel[1, 1] = 0.0


def shape_functions(xi, eta):
    dNdX_ = np.zeros((m_dim, m_nodxelem))
    N = np.array([(1 - xi) * (1 - eta) / 4,
                  (1 + xi) * (1 - eta) / 4,
                  (1 + xi) * (1 + eta) / 4,
                  (1 - xi) * (1 + eta) / 4])
    dNdX_[0, :] = np.array([-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4])
    dNdX_[1, :] = np.array([-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4])
    dNdX_[2, :] = np.array([0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), 0.25 * (1 + eta)])  # Update for axisymmetric
    return N, dNdX_


def calc_jacobian(pos):
    J = np.zeros((m_gp_count, m_dim, m_dim))
    detJ = np.zeros((m_gp_count))
    for gp in range(len(gauss_points)):
        xi, eta = gauss_points[gp]
        N, dNdrs[gp] = shape_functions(xi, eta)
        J[gp] = np.dot(dNdrs[gp], pos)
        detJ[gp] = np.linalg.det(J[gp])
        invJ = np.linalg.inv(J[gp])
        dNdX[gp] = np.dot(invJ, dNdrs[gp])
    return J, detJ, dNdX


def calc_vol(detJ):
    vol = 0.0
    for gp in range(len(gauss_points)):
        vol += detJ[gp] * gauss_weights[gp]
    return vol


J, detJ, dNdX = calc_jacobian(x)
vol_0 = calc_vol(detJ)
nod_mass = vol_0 * rho / m_nodxelem

# Main loop
t = 0.0
while t < tf:
    print("Time: ", t)
    u = dt * (v + 0.5 * dt * prev_a)  # Prediction phase
    v = v + dt * prev_a
    a = np.zeros((m_nodxelem, m_dim))
    impose_bc(v, a)

    J, detJ, dNdX = calc_jacobian(x)

    str_rate = np.zeros((m_gp_count, m_dim, m_dim))
    # Calculate strain rate
    for gp in range(m_gp_count):
        grad_v = np.zeros((m_dim, m_dim))
        for I in range(m_dim):
            for J in range(m_dim):
                for k in range(m_nodxelem):
                    grad_v[I, J] += dNdX[gp, J, k] * v[k, I]
        str_rate[gp] = 0.5 * (grad_v + grad_v.T)

    strain += dt * str_rate
    stress = np.zeros((m_gp_count, m_dim, m_dim))
    for gp in range(m_gp_count):
        c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))  # Axisymmetric material properties
        stress[gp, 0, 0] = c * ((1.0 - nu) * strain[gp, 0, 0] + nu * strain[gp, 1, 1])
        stress[gp, 1, 1] = c * ((1.0 - nu) * strain[gp, 1, 1] + nu * strain[gp, 0, 0])
        stress[gp, 0, 1] = stress[gp, 1, 0] = c * (1.0 - 2 * nu) * strain[gp, 0, 1]

    forces = np.zeros((m_nodxelem, m_dim))
    for gp in range(len(gauss_points)):
        B = np.zeros((m_dim, m_nodxelem))
        for i in range(m_nodxelem):
            B[0, i] = dNdX[gp, 0, i]
            B[1, i] = dNdX[gp, 1, i]
            B[2, i] = x[i, 0] * dNdX[gp, 2, i]  # Contribution from radial direction
        forces += np.dot(B.T, stress[gp]) * detJ[gp] * gauss_weights[gp]

    a = -forces / nod_mass
    v = v + dt * a
    u = u + dt * v
    x = x + dt * u
    prev_a = a
    t += dt

print("Final displacements:\n", u)