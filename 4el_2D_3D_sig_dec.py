### DECOMPOSED SIGMA

import numpy as np
# Define material properties
E   = 206e9  # Young's modulus in Pa
nu  = 0.3   # Poisson's ratio
rho = 7850.0
m_dim = 3

if (m_dim == 2):
  m_nodxelem = 4
else:
  m_nodxelem = 8
  
mat_G = E/(2.0*(1+nu))
print("mat G", mat_G)
K_mod = E / ( 3.0*(1.0 -2.0*nu) )
mat_cs = np.sqrt(K_mod/rho)
# Define element 

red_int = True
element_length = 1.0   # Length of the element
axi_symm = False #FALSE: PLAIN STRAIN 

# Define nodal v (dummy data for demonstration)
vel = np.full(m_dim * m_nodxelem, 0.1)
vel[5] = vel[7] = -1.0

dt = 0.8e-5
#tf = 2.0*dt
tf = 1.0e-3    
if (m_dim==2):
  x      =  np.array([[0.1, 0.], [0.2, 0.], [0.2, 0.1], [0.1, 0.1]])
  #x      =  np.array([[0., 0.], [0.1, 0.], [0.1, 0.1], [0.0, 0.1]]) #SHOULD BE RADIAL FIX
  v      = np.array([[0, 0], [0, 0], [0, -1], [0, -1]])  # Example v at nodes
else:
  x      =  np.array([[0.1, 0.,0.0], [0.2, 0.,0.0], [0.2, 0.1,0.0], [0.1, 0.1,0.0],
                    [0.1, 0.,0.1], [0.2, 0.,0.1], [0.2, 0.1,0.1], [0.1, 0.1,0.1]])
  v      = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],
                     [0, 0, -1], [0, 0, -1], [0, 0, -1], [0, 0, -1]])  # Example v at nodes

############################
# Gauss quadrature points and weights
a      = np.zeros((m_nodxelem,m_dim)) 
u      = np.zeros((m_nodxelem,m_dim)) 
u_tot  = np.zeros((m_nodxelem,m_dim))
prev_a = np.zeros((m_nodxelem,m_dim)) 

    
if red_int:
  gauss_points = np.array([[0.0, 0.0, 0.0]])
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
N    = np.zeros((m_gp_count, 1, m_nodxelem)) 
dNdrs = np.zeros((m_gp_count, m_dim, m_nodxelem)) 

rot_rate = np.zeros((m_gp_count,3, 3))
strain = np.zeros((m_gp_count,3, 3))
tau = np.zeros((m_gp_count,3, 3))
pres = np.zeros(m_gp_count)
stress = np.zeros((m_gp_count,3, 3))
radius = np.zeros(m_gp_count)

def impose_bc(vel, accel):
  vel[2,1] = vel[3,1] = -1.0
  vel[0,:] = vel[1,1] = 0.0

  accel[2,1] = accel[3,1] = 0.0
  accel[0,:] = accel[1,1] = 0.0

def calc_radius (N):
  for gp in range(len(gauss_points)):
    rv = np.dot(N[gp],x[:,0])
    radius[gp] = rv
  return radius
  
# Define shape functions and their derivatives for 2D quadrilateral element
          # !!!!! J-1 = dr/dx
          # !!!! dHxy = J-1 x dHrs = [ ] x 0.25[-1 1 -1 1]
          # !!!!                               [-1 -1 1 1]
          # elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2) !For each 3 rows of inv J and dHdxy
          # elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)
          # elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)
          # elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)     
          
          # elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.25d0
def shape_functions(gp):

    dNdX_ = np.zeros((m_dim, m_nodxelem))
    xi = gauss_points[gp][0]
    eta = gauss_points[gp][1]
    if (m_dim ==3):
      t = gauss_points[gp][2]
      
    N = np.array([(1-xi)*(1-eta)/4,
                  (1+xi)*(1-eta)/4,
                  (1+xi)*(1+eta)/4,
                  (1-xi)*(1+eta)/4])
    if (m_dim ==2):
      dNdX_[0,:] = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
      dNdX_[1,:] = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
    else:
      dNdX_[0,:] = np.array([-(1-eta)*(1-t)/8, (1-eta)*(1-t)/8, (1+eta)*(1-t)/8, -(1+eta)*(1-t)/8,
                             -(1-eta)*(1+t)/8, (1-eta)*(1+t)/8, (1+eta)*(1+t)/8, -(1+eta)*(1+t)/8])
      dNdX_[1,:] = np.array([-(1-eta)*(1-t)/8, (1-eta)*(1-t)/8, (1+eta)*(1-t)/8, -(1+eta)*(1-t)/8,
                             -(1-eta)*(1+t)/8, (1-eta)*(1+t)/8, (1+eta)*(1+t)/8, -(1+eta)*(1+t)/8])      
    
    return dNdX_
    # print(dNdX)


gp_count = len(gauss_points)

# Finite element JACOBIAN AND DERIVATIVES CALC
def calc_jacobian(pos):
    J = np.zeros((gp_count, m_dim, m_dim))
    detJ = np.zeros((gp_count))
    for gp in range(len(gauss_points)):
        dNdrs[gp] = shape_functions(gp)
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
      # print ("vol " + str(vol))
  return vol

def velocity_gradient_tensor(dNdX, vel):
    grad_v = np.zeros((m_gp_count,m_dim,m_dim))
    for gp in range (m_gp_count):
        for I in range(m_dim): 
            for J in range(m_dim):
                for k in range(m_nodxelem): 
                    #grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k * m_dim + I]
                    grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k, I]

    return grad_v

def calc_str_rate (dNdX,v):
    str_rate = np.zeros((m_gp_count,3,3))
    ss= np.zeros((m_gp_count,3,3))
    for gp in range (m_gp_count):
        grad_v = velocity_gradient_tensor(dNdX, v)
        # print("Velocity gradients\n" ,grad_v[0])

        str_rate[gp,0:m_dim,0:m_dim] = 0.5*(grad_v[gp]+grad_v[gp].T)
        rot_rate[gp,0:m_dim,0:m_dim] = 0.5*(grad_v[gp]-grad_v[gp].T)

        if (axi_symm):
            str_rate [gp,2,2] = np.dot(N[gp],v[:,0])/radius[gp]
    #print("strain rate:\n" ,str_rate)
    return str_rate, rot_rate


def calc_strain(str_rate,dt):
    strain = np.zeros((m_gp_count,3, 3))
    strain = dt * str_rate
    return strain

    
# def calc_stress(eps,dNdX):
    # stress = np.zeros((m_gp_count,m_dim, m_dim))
    # #c = E / (1.0- nu*nu)
    # c = E / ((1.0+nu)*(1.0-2.0*nu)) # #!!!! PLAIN STRAIN
    # for gp in range(len(gauss_points)):
        # stress[gp,0,0] = c * ((1.0-nu)*eps[gp,0,0] + nu*eps[gp,1,1])
        # stress[gp,1,1] = c * ((1.0-nu)*eps[gp,1,1] + nu*eps[gp,0,0])
        # stress[gp,0,1] = stress[gp,1,0] = c * (1.0-2*nu)*eps[gp,0,1] 
    # return stress
    
def calc_pressure(K_,dstr,stress):
  pr = np.zeros(m_gp_count)
  pi_= 0.0
  for gp in range(len(gauss_points)):
    pi_= pi_ + np.trace(dstr[gp])
  pi_ = -pi_/float(len(gauss_points))
  for gp in range(len(gauss_points)):
    pr[gp] = -1.0/3.0 *  np.trace(stress[gp]) + K_ * pi_
  return pr

def dev(t):
  d = np.zeros((3,3))
  d = t-1.0/3.0*np.trace(t)*np.identity(3) 
  return d

def calc_stress2(str_rate, rot_rate, tau, p, dt):

  for gp in range(len(gauss_points)):
    srt = np.dot(tau[gp],np.transpose(rot_rate[gp]))
    rs  = np.dot(rot_rate[gp],tau[gp])

    tau[gp] +=  dt * (2.0 * mat_G * (dev(str_rate[gp])) + rs + srt )
    stress[gp] =  tau[gp] - p[gp] * np.identity(3)

  return stress

#We can calc with B matrix
#F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
#               = [dh2/dx dh2/dy ]   [ syx syy]
def calc_forces(stress,dNdX,J):
  forces = np.zeros((m_nodxelem,m_dim))
  fax = np.zeros((m_nodxelem,m_dim))
  
  B = np.zeros((m_dim, m_nodxelem))
  f= 1.0

  for gp in range(len(gauss_points)):
    for i in range(m_nodxelem):
      for j in range (m_dim):
        B[j, i] = dNdX[gp,j,i]

    if (axi_symm):
        f = radius[gp]
    forces +=  np.dot(B.T,stress[gp,0:m_dim,0:m_dim]) *  np.linalg.det(J[gp]) * gauss_weights[gp] *f
  if (axi_symm):
    #(srr-stt)/r x r dr dz x dt
    for gp in range(len(gauss_points)):
      #RIGOUROUSLY SHOULD BE H[gp] * s(..)
      #is value taken to the nodes
      #print ("Ngp ", N[gp])
      #fax[:,0] += N[gp,:].T*(stress[gp,0,0]-stress[gp,2,2]) * gauss_weights[gp]  #SHAPE MAT
      #fax[:,1] += N[gp]* stress[gp,0,1] * gauss_weights[gp] #SHAPE MAT
      f = np.linalg.det(J[gp]) * gauss_weights[gp]
      fax[:,0] -= N[gp,0,:]*(stress[gp,0,0]-stress[gp,2,2]) * f
      #fax[:,0] += 0.25*(stress[gp,0,0]-stress[gp,2,2]) * f  #0,25 is from pass to nod, SHAPE MAT
      fax[:,1] -= N[gp,0,:]*stress[gp,0,1] * f
      #fax[:,1] += 0.25* stress[gp,0,1] * f #SHAPE MAT
    forces = (forces + fax)*2.0 * np.pi
  #print ("forces")
  #print (forces)
  return forces

def calc_hg_forces(rho, vol,cs):
  f_ = np.zeros((m_nodxelem,m_dim))
  Sig = np.array([[1.,-1.,1.,-1.],[1.,-1.,1.,-1.],[1.,-1.,1.,-1.],[1.,-1.,1.,-1.]])
  hmod = np.zeros((m_dim,4))
  jmax = 4
  for gp in range(len(gauss_points)):
    for j in range(jmax):
      for n in range(m_nodxelem):
        hmod[:,j] +=v[n,:]*Sig[j,n]
  for j in range(jmax):
    for n in range(m_nodxelem):
      f_[n,:] -=hmod[:,j]*Sig[j,n] 
  ch = 0.06 * pow (vol,0.66666) * rho * 0.25 * cs
  #print ("hg forces", f_)
  f_ *= ch
  return f_
#strain_rate = calc_strain_rate(v)
J, detJ, dNdX = calc_jacobian(x)
print ("Jacobian\n", J[0])
print ("det J \n", detJ[0])
vol_0 = calc_vol(detJ)
nod_mass = vol_0 * rho / m_nodxelem
radius =calc_radius(N)
#only for this example case
if (axi_symm):
    rr = 0.0
    for gp in range(m_gp_count):
        rr+=radius[gp]
    rr/=m_gp_count
    print ("avg radius", rr)
    nod_mass *=rr*2.0 *np.pi
print ("nodmass ",nod_mass)
# accel ()= forces/mass

rho_b = 0.8182  # DEFAULT SPECTRAL RADIUS

alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b)
beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b))
gamma = 1.5 - alpha;

print ("V", a)
################################# MAIN LOOP ###################################
t = 0.0
while (t < tf):
    print ("---Time: ", t)
    # !!! PREDICTION PHASE
    u = dt * (v + (0.5 - beta) * dt * prev_a)
    # !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
    v = v + (1.0-gamma)* dt * prev_a
    a[:,:] = 0.0
    impose_bc(v, a)

    J, detJ, dNdX = calc_jacobian(x)
    radius =calc_radius(N)
    #print ("radius", radius)
    str_rate,rot_rate = calc_str_rate (dNdX,v)

    str_inc = calc_strain(str_rate,dt)

    pres = calc_pressure(K_mod,str_inc, stress)
    strain =  strain + str_inc

    stress =  calc_stress2(str_rate,rot_rate,tau,pres,dt)


    forces =  calc_forces(stress,dNdX,J)
    if red_int:
      forces -= calc_hg_forces(rho,vol_0, mat_cs)
    a = -forces/nod_mass
    
    a = a - alpha * prev_a
    a = a / (1.0 - alpha)
    v = v + gamma * dt * a  

    impose_bc (v, a) #REINFORCE VELOCITY BC

    u = u + beta * dt * dt * a   
    # nod%u = nod%u + u
    x = x + u
    
    # !call AverageData(elem%rho(:,1),nod%rho(:))  
    prev_a = a
    # time = time + dt

    u_tot += u 
    t+= dt
str_rate = calc_str_rate (dNdX,v)
    

# print (strain)
# print("STRESS")
print ("veloc ", v)
print ("accel ", a)
print ("pressure", pres)
print ("stress",  stress)
print ("tau",  tau)
print ("Forces", forces)
print("strain rate:\n" ,str_rate[0])

print ("DISPLACEMENTS\n",u_tot)