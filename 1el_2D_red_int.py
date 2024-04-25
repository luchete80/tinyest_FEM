import numpy as np
# Define material properties
E   = 200e9  # Young's modulus in Pa
nu  = 0.3   # Poisson's ratio
rho = 7850.0
m_dim = 2
m_nodxelem = 4
# Define element properties

element_length = 1.0   # Length of the element
m_gp_count = 1

# Define nodal v (dummy data for demonstration)
vel = np.full(m_dim * m_nodxelem, 0.1)
vel[5] = vel[7] = -1.0

dt = 0.8e-5
    
x      =  np.array([[0., 0.], [0.1, 0.], [0.1, 0.1], [0., 0.1]])
v      = np.array([[0, 0], [0, 0], [0, -1], [0, -1]])  # Example v at nodes

############################
# Gauss quadrature points and weights
a      = np.zeros((m_nodxelem,2)) 
u      = np.zeros((m_nodxelem,2)) 
u_tot  = np.zeros((m_nodxelem,2))
prev_a = np.zeros((m_nodxelem,2)) 

gauss_points = np.array([[0.0, 0.0]])
gauss_weights = np.array([4.0])
detJ = np.zeros((m_gp_count))
dNdX = np.zeros((m_gp_count, m_dim, m_nodxelem)) 
dNdrs = np.zeros((m_gp_count, m_dim, m_nodxelem)) 

strain = np.zeros((m_gp_count,m_dim, m_dim))

def impose_bcv(vel):
  vel[2,1] = vel[2,1] = -1.0
  vel[0,:] = vel[1,0] = 0.0
  
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
        weight = gauss_weights[gp]
        N, dNdrs[gp] = shape_functions(xi, eta)
        J[gp] = np.dot(dNdrs[gp], pos)
        detJ[gp] = np.linalg.det(J[gp])
        print("det J\n", detJ)
        invJ = np.linalg.inv(J[gp])
        print ("invJ", invJ)
        dNdX[gp] = np.dot(invJ,dNdrs[gp])
    return J, detJ, dNdX

def calc_vol(detJ):
  vol = 0.0
  for gp in range(len(gauss_points)):
      vol += detJ[gp] * gauss_weights[gp]
      print ("vol " + str(vol))
  return vol
# !!!!! ASSUME VOLUME IS ALREADY CALCULATED
# subroutine calc_elem_density ()
  # implicit none
  # real(fp_kind), dimension(dim,dim) :: F !def gradient
  # real(fp_kind), dimension(nodxelem,dim) :: x !def gradient
  
  # integer :: e, n, gp
  # do e = 1, elem_count
    # do gp=1, elem%gausspc(e)
    # !if (elem%gausspc(e) .eq. 1) then
      # elem%rho(e,gp) = elem%rho_0(e,gp)*elem%vol_0(e)/elem%vol(e) !IS THE SAME
    # end do
  # end do
# end subroutine



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
        print("Velocity gradients\n" ,grad_v[0])

        str_rate[gp] = 0.5*(grad_v[0]+grad_v[0].T)
    print("strain rate:\n" ,str_rate)
    return str_rate



# Define material matrix for plane strain
def material_matrix():
    C = E / ((1+nu)*(1 - 2.0*nu)) * np.array([[1-nu, nu, 0],
                                     [nu, 1-nu, 0],
                                     [0, 0, (1 - 2.0*nu) ]])
    return C

def calc_strain(str_rate,dt):
    strain = np.zeros((m_gp_count,m_dim, m_dim))
    strain = dt * str_rate
    return strain
    
def calc_stress(eps,dNdX):
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
    print ("forces")
    print (forces)
    return forces


#strain_rate = calc_strain_rate(v)
J, detJ, dNdX = calc_jacobian(x)
print ("Jacobian\n", J[0])
print ("det J \n", detJ[0])
vol_0 = calc_vol(detJ)
nod_mass = vol_0 * rho / m_nodxelem 

# accel ()= forces/mass

rho_b = 0.8182  # DEFAULT SPECTRAL RADIUS

alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b)
beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b))
gamma = 1.5 - alpha;

print ("V", a)
################################# MAIN LOOP ###################################
for step in range(1):
    # !!! PREDICTION PHASE
    u = dt * (v + (0.5 - beta) * dt * prev_a)
    # !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
    v = v + (1.0-gamma)* dt * prev_a
    a[:,:] = 0.0
    impose_bcv(v)

    J, detJ, dNdX = calc_jacobian(x)
    # vol_0 = calc_vol(detJ)
    # nod_mass = vol_0 * rho / m_nodxelem 

    str_rate = calc_str_rate (dNdX,v)
    print ("strain rate\n",str_rate)
    # strain =  strain + calc_strain(str_rate,dt)
    strain =  calc_strain(str_rate,dt)
    stress =  calc_stress(strain,dt)
    print ("stress\n",stress)
    forces =  calc_forces(stress,dNdX,J)

    a = a - alpha * prev_a
    a = a / (1.0 - alpha)
    v = v + gamma * dt * a  

    impose_bcv (v) #REINFORCE VELOCITY BC

    u = u + beta * dt * dt * a   
    # nod%u = nod%u + u
    x = x + u
    
    # !call AverageData(elem%rho(:,1),nod%rho(:))  
    prev_a = a
    # time = time + dt
  
print ("DISPLACEMENTS\n",u)
print (strain)
print("STRESS")
print (stress)
print("strain rate:\n" ,str_rate[0])
