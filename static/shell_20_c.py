import numpy as np

# Material properties
E = 2.1e9  # Young's modulus (Pa)
nu = 0.3   # Poisson's ratio
t = 0.01   # Thickness (m)

# Geometry (node positions for a quadrilateral element)
nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
import numpy as np

def shape_functions(xi, eta):
    """Shape functions and their derivatives for MITC4 elements."""
    
    # Shape functions for a 4-node quadrilateral element
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)
    
    # First derivatives of shape functions w.r.t. xi and eta
    dN1_dxi = -0.25 * (1 - eta)
    dN1_deta = -0.25 * (1 - xi)
    
    dN2_dxi = 0.25 * (1 - eta)
    dN2_deta = -0.25 * (1 + xi)
    
    dN3_dxi = 0.25 * (1 + eta)
    dN3_deta = 0.25 * (1 + xi)
    
    dN4_dxi = -0.25 * (1 + eta)
    dN4_deta = 0.25 * (1 - xi)
    
    # Second derivatives of shape functions
    d2N1_dxi2 = 0  # ∂²N1/∂xi²
    d2N1_deta2 = 0  # ∂²N1/∂eta²
    d2N1_dxdy = 0.25  # ∂²N1/∂xi∂eta
    
    d2N2_dxi2 = 0
    d2N2_deta2 = 0
    d2N2_dxdy = -0.25
    
    d2N3_dxi2 = 0
    d2N3_deta2 = 0
    d2N3_dxdy = 0.25
    
    d2N4_dxi2 = 0
    d2N4_deta2 = 0
    d2N4_dxdy = -0.25

    # Returning the shape functions, first derivatives, and second derivatives
    N = np.array([N1, N2, N3, N4])
    dN_dxi = np.array([dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi])
    dN_deta = np.array([dN1_deta, dN2_deta, dN3_deta, dN4_deta])
    d2N_dxi2 = np.array([d2N1_dxi2, d2N2_dxi2, d2N3_dxi2, d2N4_dxi2])
    d2N_deta2 = np.array([d2N1_deta2, d2N2_deta2, d2N3_deta2, d2N4_deta2])
    d2N_dxdy = np.array([d2N1_dxdy, d2N2_dxdy, d2N3_dxdy, d2N4_dxdy])
    
    return N, dN_dxi, dN_deta, d2N_dxi2, d2N_deta2, d2N_dxdy

# Compute the membrane stiffness matrix
def compute_membrane_stiffness(N, dN_dxi, dN_deta, material_properties, thickness, detJ, weight):
    E, nu = material_properties
    # Plane stress assumption for shell, constitutive matrix D
    D_m = E / (1 - nu**2) * np.array([[1, nu, 0],
                                      [nu, 1, 0],
                                      [0, 0, (1 - nu) / 2]])

    # Create the strain-displacement matrix B for membrane
    B_m = np.zeros((3, 8))  # 4 nodes × 2 DOFs per node (u, v)
    for i in range(4):
        B_m[0, 2*i] = dN_dxi[i]
        B_m[1, 2*i+1] = dN_deta[i]
        B_m[2, 2*i] = dN_deta[i]
        B_m[2, 2*i+1] = dN_dxi[i]

    # Membrane stiffness contribution
    K_m = np.dot(B_m.T, np.dot(D_m, B_m)) * thickness * detJ * weight
    return K_m

def compute_bending_stiffness(N, d2N_dx2, d2N_dy2, d2N_dxdy, material_properties, thickness, detJ, weight):
    E, nu = material_properties
    
    # Bending stiffness matrix for thin shells
    D_b = E * thickness**3 / (12 * (1 - nu**2)) * np.array([[1, nu, 0],
                                                           [nu, 1, 0],
                                                           [0, 0, (1 - nu) / 2]])

    # Create the bending strain-displacement matrix B
    B_b = np.zeros((3, 8))  # 4 nodes × 2 DOFs per node (rotations)
    for i in range(4):
        B_b[0, 2*i] = d2N_dx2[i]    # Curvature κ_x
        B_b[1, 2*i+1] = d2N_dy2[i]  # Curvature κ_y
        B_b[2, 2*i] = d2N_dxdy[i]   # Twist curvature κ_xy
        B_b[2, 2*i+1] = d2N_dxdy[i] # Twist curvature κ_xy

    # Compute the bending stiffness matrix
    K_b = np.dot(B_b.T, np.dot(D_b, B_b)) * detJ * weight
    return K_b


# Compute the shear stiffness matrix
def compute_shear_stiffness(N, dN_dxi, dN_deta, material_properties, thickness, detJ, weight):
    E, nu = material_properties
    # Shear modulus (assuming isotropic material)
    G = E / (2 * (1 + nu))
    
    # Shear correction factor (approximated)
    k_shear = 5 / 6  # Typically used value for 4-node quadrilateral shells
    
    # Create the strain-displacement matrix S for shear
    S = np.zeros((2, 4))  # 4 nodes × 1 DOF per node (w and rotations)
    for i in range(4):
        S[0, i] = dN_dxi[i]
        S[1, i] = dN_deta[i]

    # Shear stiffness matrix with correction
    K_s = k_shear * np.dot(S.T, np.dot(G, S)) * thickness * detJ * weight
    return K_s

# Compute the Jacobian matrix
def compute_jacobian(dN_dxi, dN_deta, nodes):
    """Compute the Jacobian matrix for a 4-node quadrilateral element."""
    # Derivatives of shape functions w.r.t. xi and eta
    dN_dxi = dN_dxi.reshape(4, 1)
    dN_deta = dN_deta.reshape(4, 1)
    
    # Node coordinates (x, y)
    x = nodes[:, 0].reshape(4, 1)
    y = nodes[:, 1].reshape(4, 1)
    
    # Compute Jacobian matrix
    J = np.zeros((2, 2))
    J[0, 0] = np.dot(dN_dxi.T, x)  # dx/dxi
    J[0, 1] = np.dot(dN_dxi.T, y)  # dy/dxi
    J[1, 0] = np.dot(dN_deta.T, x)  # dx/deta
    J[1, 1] = np.dot(dN_deta.T, y)  # dy/deta
    
    return J
    


# Compute the global stiffness matrix
def compute_stiffness_matrix(nodes, material_properties, thickness):
    """Assembly of stiffness matrix for the MITC4 shell element"""
    # Initialize the global stiffness matrix for a shell element with 5 DOFs per node (20 DOFs)
    K_global = np.zeros((20, 20))  
    
    # Gauss quadrature points and weights for 2x2 integration
    gauss_points = [(-1/np.sqrt(3), -1/np.sqrt(3)), (1/np.sqrt(3), -1/np.sqrt(3)),
                    (1/np.sqrt(3), 1/np.sqrt(3)), (-1/np.sqrt(3), 1/np.sqrt(3))]
    gauss_weights = [1.0, 1.0, 1.0, 1.0]  # Weights for 2x2 Gauss quadrature

    # Loop over integration points
    for idx, (xi, eta) in enumerate(gauss_points):
        # Compute shape functions and derivatives
        N, dN_dxi, dN_deta, d2N_dxi2, d2N_deta2,d2N_dxdy = shape_functions(xi, eta)

        # Compute Jacobian matrix
        J = compute_jacobian(dN_dxi, dN_deta, nodes)
        detJ = np.linalg.det(J)  # Determinant of the Jacobian

        # Integration weight
        weight = gauss_weights[idx]

        # Compute membrane stiffness
        K_m = compute_membrane_stiffness(N, dN_dxi, dN_deta, material_properties, thickness, detJ, weight)

        # Compute bending stiffness
        
        K_b = compute_bending_stiffness(N,  d2N_dxi2, d2N_deta2,d2N_dxdy, material_properties, thickness, detJ, weight)
        print (K_b)
        
        # Compute shear stiffness
        K_s = compute_shear_stiffness(N, dN_dxi, dN_deta, material_properties, thickness, detJ, weight)

        # Assemble into the global stiffness matrix
        membrane_dofs = [0, 1, 5, 6, 10, 11, 15, 16]
        bending_dofs = [3, 4, 8, 9, 13, 14, 18, 19]
        shear_dofs = [2, 7, 12, 17]

        for i, dof_i in enumerate(membrane_dofs):
            for j, dof_j in enumerate(membrane_dofs):
                K_global[dof_i, dof_j] += K_m[i, j]

        for i, dof_i in enumerate(bending_dofs):
            for j, dof_j in enumerate(bending_dofs):
                K_global[dof_i, dof_j] += K_b[i, j]

        for i, dof_i in enumerate(shear_dofs):
            for j, dof_j in enumerate(shear_dofs):
                K_global[dof_i, dof_j] += K_s[i, j]
    
    return K_global

# Define the load vector (simple example, add forces to specific nodes)
F = np.zeros(20)  # Load vector for 20 DOFs
F[12] = -100  # Force at node 3 (e.g., applied in the z direction)

# Boundary conditions (fix some nodes)
fixed_dofs = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])  # Fix all DOFs at node 0
free_dofs = np.setdiff1d(np.arange(20), fixed_dofs)

# Compute the stiffness matrix
K_total = compute_stiffness_matrix(nodes, (E, nu), t)

# Solve the system of equations (for displacements)
K_total_reduced = K_total[np.ix_(free_dofs, free_dofs)]
F_reduced = F[free_dofs]

# Solve for displacements
U_reduced = np.linalg.solve(K_total_reduced, F_reduced)

# Reconstruct full displacement vector
U = np.zeros(20)
U[free_dofs] = U_reduced

# Output results (deformed shape, stresses, etc.)
print("Displacements:", U)
