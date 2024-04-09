import numpy as np

# Define material properties

E = 200e9  # Young's modulus in Pa
nu = 0.3   # Poisson's ratio
m_dim = 3
m_nodxelem = 8
# Define element properties

element_length = 1.0   # Length of the element
m_gp_count = 1

dNdX = np.zeros((m_dim, m_nodxelem)) 
dNdrs = np.zeros((m_dim, m_nodxelem)) 
# Define shape functions and their derivatives for 2D quadrilateral element
def shape_functions(xi, eta):
    dNdX_ = np.zeros((m_dim, m_nodxelem))
    N = np.array([0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125])
    # dNdX_[0,:] = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
    # dNdX_[1,:] = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
    return N, dNdX_
    print(dNdX)
# Gauss quadrature points and weights
gauss_points = np.array([[-0.577350269, -0.577350269],
                         [ 0.577350269, -0.577350269],
                         [ 0.577350269,  0.577350269],
                         [-0.577350269,  0.577350269]])

gauss_weights = np.array([1, 1, 1, 1])

gp_count = len(gauss_points)

# Finite element strain rate calculation
def calculate_jacobian(x2):
    J = np.zeros((m_dim, m_dim))
    # x2 = np.zeros((m_nodxelem,m_dim))
    # x2(i,:)=nod%x(x[:,i]
    # N, dNdrs[gp] = shap1_functions(xi, eta)
    #J[gp] = np.dot(dNdrs1gp], pos)
    print ("x2")
    print(x2)
    J[0,:] = -x2[0,:]+x2[1,:]+x2[2,:]-x2[3,:]-x2[4,:]+x2[5,:]+x2[6,:]-x2[7,:]
    J[1,:] = -x2[0,:]-x2[1,:]+x2[2,:]+x2[3,:]-x2[4,:]-x2[5,:]+x2[6,:]+x2[7,:]
    J[2,:] = -x2[0,:]-x2[1,:]-x2[2,:]-x2[3,:]+x2[4,:]+x2[5,:]+x2[6,:]+x2[7,:]
    
    detJ = np.linalg.det(J)
    print("det J\n", detJ)
    invJ = np.linalg.inv(J)
    print ("invJ", invJ)
    dNdX[:,0] = -invJ[:,0]-invJ[:,1]-invJ[:,2] #For each 3 rows of inv J and dHdxy
    dNdX[:,1] =  invJ[:,0]-invJ[:,1]-invJ[:,2]
    dNdX[:,2] =  invJ[:,0]+invJ[:,1]-invJ[:,2]
    dNdX[:,3] = -invJ[:,0]+invJ[:,1]-invJ[:,2]
    dNdX[:,4] = -invJ[:,0]-invJ[:,1]+invJ[:,2]
    dNdX[:,5] =  invJ[:,0]-invJ[:,1]+invJ[:,2]
    dNdX[:,6] =  invJ[:,0]+invJ[:,1]+invJ[:,2]
    dNdX[:,7] = -invJ[:,0]+invJ[:,1]+invJ[:,2]
        
    return J,dNdX 

def velocity_gradient_tensor(dNdX, vel):
    grad_v = np.zeros((m_gp_count,m_dim, m_dim))
    for gp in range (m_gp_count):
        for I in range(m_dim): 
            for J in range(m_dim):
                for k in range(m_nodxelem): 
                    #grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k * m_dim + I]
                    grad_v[gp,I, J] += dNdX[J, k] * vel[k, I]
    return grad_v

def calc_str_rate (dNdX,velocities):
    str_rate = np.zeros((m_gp_count,m_dim, m_dim))
    for gp in range (m_gp_count):
        grad_v = velocity_gradient_tensor(dNdX, velocities)
        print("Velocity gradients\n" ,grad_v[0])

        str_rate[gp] = 0.5*(grad_v[0]+grad_v[0].T)
    print("strain rate:\n" ,str_rate)
    return str_rate



# Define material matrix for plane stress
def material_matrix():
    C = E / (1 - nu**2) * np.array([[1, nu, 0],
                                     [nu, 1, 0],
                                     [0, 0, (1 - nu) / 2]])
    return C

def calculate_strain(str_rate,dt):
    strain = np.zeros((m_gp_count,m_dim, m_dim))
    strain = dt * str_rate
    return strain
    
def calculate_stress(eps,dNdX):
    stress = np.zeros((m_gp_count,m_dim, m_dim))
    # strain = np.zeros((m_gp_count,m_dim, m_dim))
    # eps[gp] +=  str_rate * dt
  # # PLAIN STRESS
    c = E / (1.0- nu*nu)
  
  # #!!!! PLAIN STRAIN
  # #c = dom%mat_E / ((1.0+dom%mat_nu)*(1.0-2.0*dom%mat_nu))
    for gp in range(len(gauss_points)):
        stress[gp,0,0] = c * ((1.0-nu)*eps[gp,0,0] + nu*eps[gp,1,1])
        stress[gp,1,1] = c * ((1.0-nu)*eps[gp,0,0] + nu*eps[gp,1,1])
        stress[gp,0,1] = stress[gp,1,0] = c * (1.0-2*nu)*eps[gp,0,1] 
    return stress
    

# Define nodal velocities (dummy data for demonstration)
vel = np.full(m_dim * m_nodxelem, 0.1)
vel[5] = vel[7] = -1.0
    
x    =  np.array([[0., 0.,0.], [0.1, 0.,0.], [0.1, 0.1,0.0], [0., 0.1,0.0],
                  [0., 0.,0.1], [0.1, 0.,0.1], [0.1, 0.1,0.1], [0., 0.1,0.1]])                  
velocities = np.array([[0, 0], [0, 0], [0, -1], [0, -1]])  # Example velocities at nodes
#strain_rate = calculate_strain_rate(velocities)
J = calculate_jacobian(x)
print ("Jacobian\n", J[0])

dt = 8.0e-5
str_rate = calc_str_rate (dNdX,velocities)

strain =  calculate_strain(str_rate,dt)
stress =  calculate_stress(strain,dt)

print (strain)
print("STRESS")
print (stress)
print("strain rate:\n" ,str_rate[0])
