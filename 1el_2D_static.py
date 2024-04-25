import numpy as np
# Define material properties
E   = 206e9  # Young's modulus in Pa
nu  = 0.3   # Poisson's ratio
rho = 7850.0
m_dim = 2
m_nodxelem = 4
# Define element properties

red_int = False
element_length = 1.0   # Length of the element
  

# Define nodal v (dummy data for demonstration)
vel = np.full(m_dim * m_nodxelem, 0.1)
vel[5] = vel[7] = -1.0

dt = 0.1e-5
# tf = dt
tf = 1.0e-3    
x      =  np.array([[0., 0.], [0.1, 0.], [0.1, 0.1], [0., 0.1]])
v      = np.array([[0, 0], [0, 0], [0, -1], [0, -1]])  # Example v at nodes

############################
# Gauss quadrature points and weights
a      = np.zeros((m_nodxelem,2)) 
u      = np.zeros((m_nodxelem,2)) 
u_tot  = np.zeros((m_nodxelem,2))
prev_a = np.zeros((m_nodxelem,2)) 

if red_int:
  gauss_points = np.array([[0.0, 0.0]])
  gauss_weights = np.array([4.0])
  m_gp_count = 1
else :
  gauss_points = np.array([[-0.577350269, -0.577350269],
                         [ 0.577350269, -0.577350269],
                         [ 0.577350269,  0.577350269],
                         [-0.577350269,  0.577350269]])

  gauss_weights = np.array([1, 1, 1, 1])
  m_gp_count = 4
  
detJ = np.zeros((m_gp_count))
dNdX = np.zeros((m_gp_count, m_dim, m_nodxelem)) 
dNdrs = np.zeros((m_gp_count, m_dim, m_nodxelem)) 

strain = np.zeros((m_gp_count,m_dim, m_dim))

def impose_bc(vel, accel):
  vel[2,1] = vel[3,1] = -1.0
  vel[0,:] = vel[1,1] = 0.0

  accel[2,1] = accel[3,1] = 0.0
  accel[0,:] = accel[1,1] = 0.0

  
# Define shape functions and their derivatives for 2D quadrilateral element
          # !!!!! J-1 = dr/dx
          # !!!! dHxy = J-1 x dHrs = [ ] x 0.25[-1 1 -1 1]
          # !!!!                               [-1 -1 1 1]
          # elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2) !For each 3 rows of inv J and dHdxy
          # elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)
          # elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)
          # elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)     
          
          # elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.25d0
def shape_functions(xi, eta):
    dNdX_ = np.zeros((m_dim, m_nodxelem))
    N = np.array([(1-xi)*(1-eta)/4,
                  (1+xi)*(1-eta)/4,
                  (1+xi)*(1+eta)/4,
                  (1-xi)*(1+eta)/4])
    dNdX_[0,:] = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
    dNdX_[1,:] = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
    return N, dNdX_
    # print(dNdX)


gp_count = len(gauss_points)

# Finite element JACOBIAN AND DERIVATIVES CALC
def calc_jacobian(pos):
    J = np.zeros((gp_count, 2, 2))
    detJ = np.zeros((gp_count))
    for gp in range(len(gauss_points)):
        xi, eta = gauss_points[gp]
        N, dNdrs[gp] = shape_functions(xi, eta)
        J[gp] = np.dot(dNdrs[gp], pos)
        detJ[gp] = np.linalg.det(J[gp])
        # print("det J\n", detJ)
        invJ = np.linalg.inv(J[gp])
        # print ("invJ", invJ)
        dNdX[gp] = np.dot(invJ,dNdrs[gp])
        # print ("test", -invJ[0,0]-invJ[0,1])
        # print ("deriv",dNdX[gp] )
    return J, detJ, dNdX

def calc_vol(detJ):
  vol = 0.0
  for gp in range(len(gauss_points)):
      vol += detJ[gp] * gauss_weights[gp]
      print ("vol " + str(vol))
  return vol

def velocity_gradient_tensor(dNdX, vel):
    grad_v = np.zeros((m_gp_count,m_dim, m_dim))
    for gp in range (m_gp_count):
        for I in range(m_dim): 
            for J in range(m_dim):
                for k in range(m_nodxelem): 
                    #grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k * m_dim + I]
                    grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k, I]
    return grad_v

def calc_str_rate (dNdX,v):
    str_rate = np.zeros((m_gp_count,m_dim, m_dim))
    for gp in range (m_gp_count):
        grad_v = velocity_gradient_tensor(dNdX, v)
        # print("Velocity gradients\n" ,grad_v[0])

        str_rate[gp] = 0.5*(grad_v[gp]+grad_v[gp].T)
    # print("strain rate:\n" ,str_rate)
    return str_rate


def calc_strain(str_rate,dt):
    strain = np.zeros((m_gp_count,m_dim, m_dim))
    strain = dt * str_rate
    return strain
    
def calc_stress(eps,dNdX):
    stress = np.zeros((m_gp_count,m_dim, m_dim))
    #c = E / (1.0- nu*nu)
    c = E / ((1.0+nu)*(1.0-2.0*nu)) # #!!!! PLAIN STRAIN
    for gp in range(len(gauss_points)):
        stress[gp,0,0] = c * ((1.0-nu)*eps[gp,0,0] + nu*eps[gp,1,1])
        stress[gp,1,1] = c * ((1.0-nu)*eps[gp,1,1] + nu*eps[gp,0,0])
        stress[gp,0,1] = stress[gp,1,0] = c * (1.0-2*nu)*eps[gp,0,1] 
    return stress

#We can calc with B matrix
#F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
#               = [dh2/dx dh2/dy ]   [ syx syy]
def calc_forces(stress,dNdX,J):
    forces = np.zeros((m_nodxelem,m_dim))
    B = np.zeros((m_dim, m_nodxelem))
    
    for gp in range(len(gauss_points)):
        for i in range(m_nodxelem):
            B[0, i] = dNdX[gp,0,i]
            B[1, i] = dNdX[gp,1,i]    
        forces +=  np.dot(B.T,stress[gp]) *  np.linalg.det(J[gp]) * gauss_weights[gp]
    # print ("forces")
    # print (forces)
    return forces


#strain_rate = calc_strain_rate(v)
J, detJ, dNdX = calc_jacobian(x)


str_rate = calc_str_rate (dNdX,v)
# print ("strain rate\n",str_rate)
# strain =  strain + calc_strain(str_rate,dt)
strain =  strain + calc_strain(str_rate,dt)
# print ("strain \n",strain)
stress =  calc_stress(strain,dt)
# print ("stress\n",stress)
forces =  calc_forces(stress,dNdX,J)
a = -forces/nod_mass




