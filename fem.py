import numpy as np
#  class Material:
#    # Define material properties
#    E   = 206e9  # Young's modulus in Pa
#    nu  = 0.3   # Poisson's ratio
#   rho = 7850.0
def dev(t):
  d = np.zeros((3,3))
  d = t-1.0/3.0*np.trace(t)*np.identity(3) 
  return d


class Domain:
    
  #gauss_points = np.zeros((4,2))
  #Booleans
# ############################
# # Gauss quadrature points and weights
# a      = np.zeros((m_nodxelem,2)) 
# u      = np.zeros((m_nodxelem,2)) 
# u_tot  = np.zeros((m_nodxelem,2))
# prev_a = np.zeros((m_nodxelem,2)) 
  def __init__(self):
    self.dim = 2
    self.nodxelem = 4
    self.red_int = False
    self.mat_E = 0.0;    self.mat_nu = 0.0;    self.mat_G = 0.0;
    self.mat_rho = 0.0
    self.tot_mass = 0.0
    print ("Domain!")
    self.dim = 2
    if self.red_int:
      self.gauss_points = np.array([[0.0, 0.0]])
      self.gauss_weights = np.array([4.0])
      self.gp_count = 1
    else :
      self.gauss_points = np.array([[-0.577350269, -0.577350269],
                             [ 0.577350269, -0.577350269],
                             [ 0.577350269,  0.577350269],
                             [-0.577350269,  0.577350269]])

      self.gauss_weights = np.array([1, 1, 1, 1])
      self.gp_count = 4
  def calcMatConst(self):
    self.mat_G = self.mat_E/(2.0*(1+self.mat_nu))
    print("mat G", self.mat_G)
    self.K_mod = self.mat_E / ( 3.0*(1.0 -2.0*self.mat_nu) )
    self.mat_cs = np.sqrt(self.K_mod/self.mat_rho)

  def allocateNodal(self,n):
    print ("Allocating nodal, dim ", self.dim)
    self.x = np.zeros((n,self.dim))
    self.u = np.zeros((n,self.dim))
    self.v = np.zeros((n,self.dim))
    self.a = np.zeros((n,self.dim))
    self.prev_a = np.zeros((n,self.dim))    
    self.forces = np.zeros((n,self.dim))
    self.mass   = np.zeros(n)
    self.is_bcv = np.zeros((n,self.dim),dtype='bool') 
    self.is_fix = np.zeros((n,self.dim),dtype='bool') 
  
  def allocateElemental(self,e):
    gp = self.gp_count
    self.sigma = np.zeros((e,gp,3,3))
    self.tau   = np.zeros((e,gp,3,3))
    self.p     = np.zeros((e,gp))
    
    self.str_rate = np.zeros((e,gp,3,3))
    self.rot_rate = np.zeros((e,gp,3,3))

    detJ = np.zeros((e,gp))
    dNdX = np.zeros((e,gp, self.dim, self.nodxelem)) 

    f_i = np.zeros((e,self.nodxelem,self.dim))
    # dNdrs = np.zeros((m_gp_count, m_dim, m_nodxelem)) 

    # strain = np.zeros((m_gp_count,m_dim, m_dim))


  def addBoxLength(self, lx, ly, le):
    ex = int(lx /le); ey = int(ly/le);
    ez = 1;
    r = le/2.0
    nel = np.array([ex,ey,ez])
    self.elem_count = ex*ey
    Xp = np.zeros(self.dim)
    self.node_count = (ex+1)*(ey+1)
    print ("Node count "  + str(self.node_count))
    self.allocateNodal(self.node_count)
    self.allocateElemental(self.elem_count)
    self.tot_mass = lx * ly * self.mat_rho
    print ("tot mass ", self.tot_mass)
    
    
    if (self.dim == 2):
      # !write(*,*) "Box Particle Count is ", node_count
      p = 0
    # !do while (Xp(3) <= (V(3)+Lz))
      # j = 1;         Xp(2) = V(2)
      for j in range(nel[1]+1):
        Xp[0] = 0.0
        for i in range(nel[1]+1):
          self.x[p,:] = Xp[:]
          # print *,"node ",p , "X: ",Xp(:)
          # p = p + 1
          Xp[0] += le
          p +=1
        Xp[1] += le
        # end do
        # Xp(2) = Xp(2) + 2 * r

      # end do 
      # Xp(3) = Xp(3) + 2 * r
      
      print (ex, ey, ez)
      print ("Generated " + str(p) + " nodes")
      print (self.x)
  
  def impose_bcv(self):
    for n in range(self.node_count):
      for d in range(self.dim):
        if (self.is_fix(n,d)):
          self.v[n,d] = 0.0

        # if (self.is_bcv(n,d)):
          # self.v[n,d] = 0.0

    
  # Define shape functions and their derivatives for 2D quadrilateral element
            # !!!!! J-1 = dr/dx
            # !!!! dHxy = J-1 x dHrs = [ ] x 0.25[-1 1 -1 1]
            # !!!!                               [-1 -1 1 1]
            # elem%dHxy_detJ(e,gp,:,1) = -invJ(:,1)-invJ(:,2) !For each 3 rows of inv J and dHdxy
            # elem%dHxy_detJ(e,gp,:,2) =  invJ(:,1)-invJ(:,2)
            # elem%dHxy_detJ(e,gp,:,3) =  invJ(:,1)+invJ(:,2)
            # elem%dHxy_detJ(e,gp,:,4) = -invJ(:,1)+invJ(:,2)     
            
            # elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.25d0
  def shape_functions(self,xi, eta):
      dNdX_ = np.zeros((self.dim, self.nodxelem))
      N = np.array([(1-xi)*(1-eta)/4,
                    (1+xi)*(1-eta)/4,
                    (1+xi)*(1+eta)/4,
                    (1-xi)*(1+eta)/4])
      dNdX_[0,:] = np.array([-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4])
      dNdX_[1,:] = np.array([-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4])
      return N, dNdX_
      # print(dNdX)


  #self.gp_count = len(self.gauss_points)

  # Finite element JACOBIAN AND DERIVATIVES CALC
  def calc_jacobian(pos):
    for e in range(self.elem_count):
      for gp in range(len(gauss_points)):
          xi, eta = gauss_points[gp]
          N, dNdrs[gp] = shape_functions(xi, eta)
          J[gp] = np.dot(dNdrs[gp], pos)
          self.detJ[e,gp] = np.linalg.det(J[gp])
          # print("det J\n", detJ)
          invJ = np.linalg.inv(J[gp])
          # print ("invJ", invJ)
          self.dNdX[e,gp] = np.dot(invJ,dNdrs[gp])
          # print ("test", -invJ[0,0]-invJ[0,1])
          # print ("deriv",dNdX[gp] )


  def calc_vol(detJ):
    vol = 0.0
    for gp in range(len(gauss_points)):
        vol += detJ[gp] * gauss_weights[gp]
        print ("vol " + str(vol))
    return vol

  def velocity_gradient_tensor():
      grad_v = np.zeros((self.gp_count,m_dim, m_dim))
      for gp in range (self.gp_count):
          for I in range(m_dim): 
              for J in range(m_dim):
                  for k in range(m_nodxelem): 
                      #grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k * m_dim + I]
                      grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k, I]
      return grad_v

  def calc_str_rate (self,dNdX, v):
    vel = np.zeros(self.nodxelem,self.dim)
    for e in range(self.elem_count):
      for gp in range (m_gp_count):
          grad_v = velocity_gradient_tensor(dNdX[e,gp], vel)
          # print("Velocity gradients\n" ,grad_v[0])

          str_rate[gp] = 0.5*(grad_v[gp]+grad_v[gp].T)
      # print("strain rate:\n" ,str_rate)


  def calc_strain(str_rate,dt):
      strain = np.zeros((m_gp_count,m_dim, m_dim))
      strain = dt * str_rate
      return strain

  def calc_stress(self, dt):
    for e in range(self.elem_count):
      for gp in range(self.gp_count):
        srt = np.dot(self.tau[e,gp],np.transpose(self.rot_rate[e,gp]))
        rs  = np.dot(self.rot_rate[e,gp],self.tau[e,gp])

        self.tau[e,gp] +=  dt * (2.0 * self.mat_G * (dev(self.str_rate[e,gp])) + rs + srt )
        self.sigma[e,gp] =  self.tau[e,gp] - self.p[e,gp] * np.identity(3)



  #We can calc with B matrix
  #F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
  #               = [dh2/dx dh2/dy ]   [ syx syy]
  def calc_forces(self):
    B = np.zeros((m_dim, m_nodxelem))
    self.fi = 0.0
    for e in range(self.elem_count):
      for gp in range(len(gauss_points)):
          for i in range(m_nodxelem):
              B[0, i] = dNdX[e,gp,0,i]
              B[1, i] = dNdX[e,gp,1,i]
          #TO OPTIMIZE MULTIPLY ONLY NON SPARSE ELEMENTS
          self.fi[e] +=  np.dot(B.T,self.stress[e,gp]) *  np.linalg.det(self.J[gp]) * self.gauss_weights[gp]
      # print ("forces")
      # print (forces)
      return forces


##strain_rate = calc_strain_rate(v)
#J, detJ, dNdX = calc_jacobian(x)
#print ("Jacobian\n", J[0])
#print ("det J \n", detJ[0])
#vol_0 = calc_vol(detJ)
#nod_mass = vol_0 * rho / m_nodxelem 
#print ("nodmass\n",nod_mass)
## accel ()= forces/mass



#print ("V", a)
#TO MODIFY


  def solve(self,tf, dt):
    self.mass[:] = self.tot_mass/self.node_count
    print ("NODAL  mass (FROM AVG): ", self.mass[0])

    rho_b = 0.8182  # DEFAULT SPECTRAL RADIUS
    alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b)
    beta  = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b))
    gamma = 1.5 - alpha;
################################# MAIN LOOP ###################################
    t = 0.0
    u_ = np.zeros((self.node_count,self.dim))
    while (t < tf):
        print ("Time: ", t)
        # !!! PREDICTION PHASE
        u_ = dt * (self.v + (0.5 - beta) * dt * self.prev_a)
        # !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
        self.v +=  (1.0-gamma)* dt * self.prev_a
        self.a[:,:] = 0.0
        self.impose_bc(self.v, self.a)

        self.calc_jacobian

        self.calc_str_rate 

        self.calc_stress(dt)
        #strain =  strain + calc_strain(str_rate,dt)
        # print ("strain \n",strain)

        # print ("stress\n",stress)
        print ("a ", self.a)
        self.calc_forces
        for d in range (self.dim):
          self.a[:,d] = -self.forces[:,d]/self.mass[:]

        print ("a ", self.a)
        self.a -= alpha * self.prev_a
        self.a = self.a / (1.0 - alpha)
        self.v = self.v + gamma * dt * self.a  

        self.impose_bc (self.v, self.a) #REINFORCE VELOCITY BC

        u_ += beta * dt * dt * self.a   
        # nod%u = nod%u + u
        self.x += self.u
        
        # !call AverageData(elem%rho(:,1),nod%rho(:))  
        self.prev_a = self.a
        # time = time + dt

        self.u += u_ 
        t += dt
    #str_rate = calc_str_rate (dNdX,v)
        
    print ("DISPLACEMENTS\n",self.u)
    # print (strain)
    # print("STRESS")
    # print (stress)
    # print("strain rate:\n" ,str_rate[0])


    