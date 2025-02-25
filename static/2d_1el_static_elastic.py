import numpy as np

gauss_points = np.array([[-0.577350269, -0.577350269],
                         [ 0.577350269, -0.577350269],
                         [ 0.577350269,  0.577350269],
                         [-0.577350269,  0.577350269]])
gauss_weights = np.array([1, 1, 1, 1])

def compute_stiffness_matrix(node_coords, elasticity_matrix, thickness):
    m_dim = 2  # 2D problem
    m_nodxelem = 4  # Quadrilateral element
    
    K = np.zeros((m_nodxelem * m_dim, m_nodxelem * m_dim))
    
    for gp in range(4):
        xi, eta = gauss_points[gp]
        N, dNdrs = shape_functions(xi, eta)
        
        J = np.dot(dNdrs, node_coords)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        dNdX = np.dot(invJ, dNdrs)
        
        B = np.zeros((3, m_nodxelem * m_dim))
        for i in range(m_nodxelem):
            B[0, i*2]   = dNdX[0, i]
            B[1, i*2+1] = dNdX[1, i]
            B[2, i*2]   = dNdX[1, i]
            B[2, i*2+1] = dNdX[0, i]
        
        K += B.T @ elasticity_matrix @ B * detJ * gauss_weights[gp] * thickness
        #Equivalent to. K += np.matmul(np.matmul(B.T, elasticity_matrix), B) * detJ * gauss_weights[gp] * thickness
        # OR: K += np.dot(np.dot(B.T, elasticity_matrix), B) * detJ * gauss_weights[gp] * thickness
        # B.T @ elasticity_matrix → Multiplies the transposed strain-displacement matrix (B.T) with the elasticity matrix.
    return K

def shape_functions(xi, eta):
    dNdX_ = np.zeros((2, 4))
    N = np.array([(1-xi)*(1-eta)/4,
                  (1+xi)*(1-eta)/4,
                  (1+xi)*(1+eta)/4,
                  (1-xi)*(1+eta)/4])
    dNdX_[0,:] = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
    dNdX_[1,:] = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
    return N, dNdX_

def apply_boundary_condition(K, dof):
    """
    Apply a boundary condition to the stiffness matrix by modifying a specific degree of freedom (dof).
    """
    K[dof, :] = 0
    K[:, dof] = 0
    K[dof, dof] = 1
    
    return K

def elasticity_matrix_plane_strain(E, nu):
    """
    Compute the elasticity matrix for plane strain condition.
    """
    factor = E / ((1 + nu) * (1 - 2 * nu))
    C = factor * np.array([[1 - nu, nu, 0],
                            [nu, 1 - nu, 0],
                            [0, 0, (1 - 2 * nu) / 2]])
    return C

# The elasticity matrix for plane stress differs:
# C_plane_stress = (E / (1 - nu**2)) * np.array([[1, nu, 0],
#                                                 [nu, 1, 0],
#                                                 [0, 0, (1 - nu) / 2]])

# Plane strain elasticity matrix is 3x3, considering strain in z-direction as zero.
# Plane stress elasticity matrix is also 3x3 but assumes zero stress in z-direction.

def solve_displacement(K, F):
    """
    Solve the linear system K * X = F for the displacement vector X.
    """
    X = np.linalg.solve(K, F)
    return X
    
# Example usage
node_coords = np.array([[0, 0], [0.01, 0], [0.01, 0.01], [0, 0.01]])
E, nu = 200e9, 0.3  # Example values for steel
elasticity_matrix = elasticity_matrix_plane_strain(E, nu)
thickness = 1.0
K = compute_stiffness_matrix(node_coords, elasticity_matrix, thickness)

apply_boundary_condition(K,0)
apply_boundary_condition(K,1)
apply_boundary_condition(K,3)
print(K)

F = np.zeros(K.shape[0])
F[-2] = 1000  # Apply force at last node in x-direction
print(F)
# Solve for displacements
print("Displacements")
X = solve_displacement(K, F)
print(X)

