import numpy as np

def compute_stiffness_matrix(node_coords, elasticity_matrix, thickness):
    m_dim = 2  # 2D problem
    m_nodxelem = 4  # Quadrilateral element
    gauss_points = np.array([[-0.577350269, -0.577350269],
                             [ 0.577350269, -0.577350269],
                             [ 0.577350269,  0.577350269],
                             [-0.577350269,  0.577350269]])
    gauss_weights = np.array([1, 1, 1, 1])
    
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

def assemble_global_stiffness(nodal_coords, el_nodes, elasticity_matrix, thickness):
    num_nodes = nodal_coords.shape[0]
    global_dof = num_nodes * 2
    K_global = np.zeros((global_dof, global_dof))
    
    for elem in range(el_nodes.shape[0]):
        element_nodes = el_nodes[elem, :]
        node_coords = nodal_coords[element_nodes, :]
        K_local = compute_stiffness_matrix(node_coords, elasticity_matrix, thickness)
        
        for i in range(4):
            for j in range(4):
                K_global[element_nodes[i]*2, element_nodes[j]*2]     += K_local[i*2, j*2]
                K_global[element_nodes[i]*2+1, element_nodes[j]*2+1] += K_local[i*2+1, j*2+1]
                K_global[element_nodes[i]*2, element_nodes[j]*2+1]   += K_local[i*2, j*2+1]
                K_global[element_nodes[i]*2+1, element_nodes[j]*2]   += K_local[i*2+1, j*2]
    
    return K_global

# Example usage
node_coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
elasticity_matrix = np.eye(3)  # Replace with actual material matrix
thickness = 1.0
el_nodes = np.array([[0, 1, 2, 3]])  # Example element connectivity

global_stiffness_matrix = assemble_global_stiffness(node_coords, el_nodes, elasticity_matrix, thickness)

print(global_stiffness_matrix)
