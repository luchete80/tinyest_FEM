import numpy as np


# Material properties
E = 2.1e9  # Young's modulus (Pa)
nu = 0.3   # Poisson's ratio
t = 0.01   # Thickness (m)

# Geometry (node positions for a quadrilateral element)
nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])

# Nodal Displacements (u, v, w, rotations)
u = np.zeros(20)  # Placeholder for 5 DOFs per node, 4 nodes
# Fill in initial displacement values here

# Shape functions for MITC4 element (4-node quadrilateral)
def shape_functions(xi, eta):
    """Shape functions and their derivatives for MITC4 elements."""
    
    # Shape functions for a 4-node quadrilateral element
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)
    
    # Derivatives of shape functions w.r.t. xi and eta
    dN1_dxi = -0.25 * (1 - eta)
    dN1_deta = -0.25 * (1 - xi)
    
    dN2_dxi = 0.25 * (1 - eta)
    dN2_deta = -0.25 * (1 + xi)
    
    dN3_dxi = 0.25 * (1 + eta)
    dN3_deta = 0.25 * (1 + xi)
    
    dN4_dxi = -0.25 * (1 + eta)
    dN4_deta = 0.25 * (1 - xi)
    
    # Returning the shape functions and their derivatives
    N = np.array([N1, N2, N3, N4])
    dN_dxi = np.array([dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi])
    dN_deta = np.array([dN1_deta, dN2_deta, dN3_deta, dN4_deta])
    
    return N, dN_dxi, dN_deta

# Compute the stiffness matrix for membrane, shear, and bending components
def compute_membrane_stiffness(N, dN_dxi, dN_deta, material_properties, thickness):
    E, nu = material_properties
    # Plane stress assumption for shell, constitutive matrix D
    D_m = E / (1 - nu**2) * np.array([[1, nu, 0],
                                      [nu, 1, 0],
                                      [0, 0, (1 - nu) / 2]])

    # Create the MITC strain-displacement matrix B for membrane
    B_m = np.zeros((3, 8))  # Placeholder for the membrane B-matrix
    for i in range(4):
        B_m[0, 2*i] = dN_dxi[i]
        B_m[1, 2*i+1] = dN_deta[i]
        B_m[2, 2*i] = dN_deta[i]
        B_m[2, 2*i+1] = dN_dxi[i]

    # Membrane stiffness contribution using MITC formulation
    K_m = np.dot(B_m.T, np.dot(D_m, B_m)) * thickness
    return K_m
    
def compute_bending_stiffness(N, dN_dxi, dN_deta, material_properties, thickness):
    E, nu = material_properties
    # Bending stiffness matrix for thin shells (with MITC correction)
    D_b = E * thickness**3 / (12 * (1 - nu**2)) * np.array([[1, nu, 0],
                                                           [nu, 1, 0],
                                                           [0, 0, (1 - nu) / 2]])

    # Second derivatives of shape functions for bending
    d2N_dxi2 = np.zeros(4)  # Second derivatives wrt xi
    d2N_deta2 = np.zeros(4)  # Second derivatives wrt eta
    d2N_dxideta = np.zeros(4)  # Mixed second derivatives

    for i in range(4):
        d2N_dxi2[i] = dN_dxi[i] ** 2
        d2N_deta2[i] = dN_deta[i] ** 2
        d2N_dxideta[i] = dN_dxi[i] * dN_deta[i]

    # Bending strain-displacement matrix (second derivatives of N)
    B_b = np.zeros((3, 8))  # Placeholder for bending B-matrix
    for i in range(4):
        B_b[0, 2*i] = d2N_dxi2[i]
        B_b[1, 2*i+1] = d2N_deta2[i]
        B_b[2, 2*i] = d2N_dxideta[i]
        B_b[2, 2*i+1] = d2N_dxideta[i]

    # Bending stiffness contribution with MITC
    K_b = np.dot(B_b.T, np.dot(D_b, B_b)) * thickness
    return K_b
def compute_shear_stiffness(N, dN_dxi, dN_deta, material_properties, thickness):
    E, nu = material_properties
    # Shear modulus (assuming isotropic material)
    G = E / (2 * (1 + nu))
    
    # Shear correction factor (approximated)
    k_shear = 5 / 6  # Typically used value for 4-node quadrilateral shells
    
    # Compute shear strain-displacement matrix S (for 4 DOFs, not 8)
    S = np.zeros((4, 4))  # 4 DOFs for shear (transverse displacement)
    for i in range(4):
        S[0, i] = dN_dxi[i]
        S[1, i] = dN_deta[i]
        S[2, i] = dN_deta[i]
        S[3, i] = dN_dxi[i]

    # Shear stiffness matrix with correction
    K_s = k_shear * np.dot(S.T, np.dot(G, S)) * thickness
    return K_s


# Compute the global stiffness matrix
def compute_stiffness_matrix(nodes, material_properties):
    """Assembly of stiffness matrix for the MITC4 shell element"""
    # Initialize the global stiffness matrix for a shell element with 5 DOFs per node (20 DOFs)
    K_total = np.zeros((20, 20))  
    
    # Loop over integration points (typically Gauss quadrature for 2D integration)
    integration_points = [(-1/np.sqrt(3), -1/np.sqrt(3)), (1/np.sqrt(3), -1/np.sqrt(3)),
                          (1/np.sqrt(3), 1/np.sqrt(3)), (-1/np.sqrt(3), 1/np.sqrt(3))]  # 2x2 Gauss points

    for xi, eta in integration_points:
        N, dN_dxi, dN_deta = shape_functions(xi, eta)

        # Compute membrane stiffness (based on shape functions, material properties)
        K_m = compute_membrane_stiffness(N, dN_dxi, dN_deta, material_properties,t)

        # Compute bending stiffness (based on rotation and bending shape functions)
        K_b = compute_bending_stiffness(N, dN_dxi, dN_deta, material_properties,t)

        # Compute shear stiffness (modified shape functions for MITC shear correction)
        K_s = compute_shear_stiffness(N, dN_dxi, dN_deta, material_properties,t)

        # Assemble into the global stiffness matrix
        # Since this is a 4-node element with 5 DOFs per node, we need to add to the correct DOFs in the global matrix
        # (this step assumes you are assembling into the correct degrees of freedom)
        K_total[0:8, 0:8] += K_m
        
        # Bending DOFs (8 DOFs for bending)
        K_total[8:16, 8:16] += K_b  # Bending stiffness (8x8)
        
        # Shear DOFs (8 DOFs for shear)
        print (K_s)
        K_total[16:20, 16:20] += K_s  # Shear stiffness (8x8)
    
    return K_total

# Define the load vector (simple example, add forces to specific nodes)
F = np.zeros(20)  # Load vector for 20 DOFs
F[12] = -100  # Force at node 3 (e.g., applied in the z direction)

# Boundary conditions (fix some nodes)
fixed_dofs = np.array([0, 1, 2, 3, 4, 5,6,7,8,9])  # Fix all DOFs at node 0
free_dofs = np.setdiff1d(np.arange(20), fixed_dofs)

# Compute the stiffness matrix
K_total = compute_stiffness_matrix(nodes, (E, nu))

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
