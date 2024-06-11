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
    
  def __init__(self, red_int):
    self.dim = 2
    self.nodxelem = 4
    self.red_int = red_int
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
    
  def AssemblyForces(self):
    self.f_i[:,:] = 0.0
    for e in range(self.elem_count):
      for ne in range(self.nodxelem):
        self.f_i[self.elnod[e,ne],:] += self.f_e[e,ne,:]
      

  def allocateNodal(self,n):
    print ("Allocating nodal, dim ", self.dim)
    self.x = np.zeros((n,self.dim))
    self.u = np.zeros((n,self.dim))
    self.v = np.zeros((n,self.dim))
    self.a = np.zeros((n,self.dim))
    self.prev_a = np.zeros((n,self.dim))    
    self.f_i = np.zeros((n,self.dim))
    self.mass   = np.zeros(n)
    self.is_bcv = np.zeros((n,self.dim),dtype='bool') 
    self.is_fix = np.zeros((n,self.dim),dtype='bool') 
    self.bcv    = np.zeros((n,self.dim)) 
  
  def allocateElemental(self,e):
    gp = self.gp_count
    self.sigma = np.zeros((e,gp,3,3))
    self.tau   = np.zeros((e,gp,3,3))
    self.p     = np.zeros((e,gp))
    
    self.str_rate = np.zeros((e,gp,3,3))
    self.rot_rate = np.zeros((e,gp,3,3))
    
    self.elnod = np.zeros((e,self.nodxelem), dtype=int)

    self.detJ = np.zeros((e,gp))
    self.J    = np.zeros((e,gp,self.dim, self.dim))
    self.dNdX = np.zeros((e,gp, self.dim, self.nodxelem)) 

    self.f_e = np.zeros((e,self.nodxelem,self.dim))
    # dNdrs = np.zeros((m_gp_count, dim, nodxelem)) 

    # strain = np.zeros((m_gp_count,dim, dim))


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
        for i in range(nel[0]+1):
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

      ei = 0
      for ey in range (nel[1]):
        for ex in range (nel[0]):
          nb1 = (nel[0]+1)* ey    + ex;    
          nb2 = (nel[0]+1)*(ey+1) + ex;
          self.elnod[ei,0] = nb1;                                
          self.elnod[ei,1] = nb1 + 1;                             
          self.elnod[ei,2] = nb2 + 1;                        
          self.elnod[ei,3] = (nel[0]+1)*(ey+1) + ex;           
			
      print ("ELNOD", self.elnod)
				 # for (int i=0;i<nodxelem;i++)cout << elnod_h[ei+i]<<", ";
					# cout << "Nel x : "<<nel[0]<<endl;
					# cout << "nodes "<<endl;
					# ei += nodxelem;
					 # }
       

      
      print (ex, ey, ez)
      print ("Generated " + str(p) + " nodes")
      print (self.x)
  
  def impose_bc(self):
    for n in range(self.node_count):
      for d in range(self.dim):
        if (self.is_fix[n,d]):
          self.v[n,d] = 0.0
          self.a[n,d] = 0.0
          
        if (self.is_bcv[n,d]):
          self.v[n,d] = self.bcv[n,d]
          self.a[n,d] = 0.0

    
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

  def getPos(self, e):
    pos = np.zeros((self.nodxelem,self.dim))
    for n in range (self.nodxelem):
      nglob = self.elnod[e,n]
      pos[n,:] = self.x[nglob,:]
    return pos
    
  def getVel(self, e):
    vel = np.zeros((self.nodxelem,self.dim))
    for n in range (self.nodxelem):
      nglob = self.elnod[e,n]
      vel[n,:] = self.v[nglob,:]
    return vel
    
  #self.gp_count = len(self.gauss_points)

  # Finite element JACOBIAN AND DERIVATIVES CALC
  def calc_jacobian(self):
    for e in range(self.elem_count):
      pos = self.getPos(e)
      for gp in range(self.gp_count):
          xi, eta = self.gauss_points[gp]
          N, dNdrs = self.shape_functions(xi, eta)
          self.J[gp] = np.dot(dNdrs, pos)
          self.detJ[e,gp] = np.linalg.det(self.J[gp])
          # print("det J\n", detJ)
          invJ = np.linalg.inv(self.J[gp])
          # print ("invJ", invJ)
          self.dNdX[e,gp] = np.dot(invJ,dNdrs)
          # print ("test", -invJ[0,0]-invJ[0,1])
          # print ("deriv",dNdX[gp] )
          # print ("dNdx ", self.dNdX)

  def calc_vol(detJ):
    vol = 0.0
    for gp in range(self.gp_count):
        vol += detJ[gp] * gauss_weights[gp]
        # print ("vol " + str(vol))
    return vol
  

    
  def velocity_gradient_tensor(self,e, gp):
    grad_v = np.zeros((self.gp_count,self.dim, self.dim))
    vel = np.zeros((self.nodxelem,self.dim))
    vel = self.getVel(e)
    # print ("vel", vel)
    for gp in range (self.gp_count):
      for I in range(self.dim): 
        for J in range(self.dim):
          for k in range(self.nodxelem): 
            #grad_v[gp,I, J] += dNdX[gp, J, k] * vel[k * dim + I]
            grad_v[gp,I, J] += self.dNdX[e,gp, J, k] * vel[k, I]
    return grad_v

  def calc_str_rate (self):
    for e in range(self.elem_count):
      for gp in range (self.gp_count):
        grad_v = self.velocity_gradient_tensor(e, gp)
        # print("Velocity gradients\n" ,grad_v[gp])
        self.str_rate[e,gp,0:self.dim,0:self.dim] = 0.5*(grad_v[gp]+grad_v[gp].T)
        self.rot_rate[e,gp,0:self.dim,0:self.dim] = 0.5*(grad_v[gp]-grad_v[gp].T)
    # print("strain rate:\n" ,self.str_rate)

  def calc_pressure(self, dt):
    
    for e in range(self.elem_count):
      pi_= 0.0
      for gp in range(self.gp_count):
        dstr = dt*self.str_rate[e,gp]
        pi_= pi_ + np.trace(dstr)
      pi_ = -pi_/float(self.gp_count)
      for gp in range(self.gp_count):
        self.p[e,gp] = -1.0/3.0 *  np.trace(self.sigma[e,gp]) + self.K_mod * pi_


  def calc_strain(str_rate,dt):
      strain = np.zeros((m_gp_count,dim, dim))
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
    B = np.zeros((self.dim, self.nodxelem))
    self.f_e[:,:] = 0.0
    for e in range(self.elem_count):
      for gp in range(self.gp_count):
          for i in range(self.nodxelem):
              B[0, i] = self.dNdX[e,gp,0,i]
              B[1, i] = self.dNdX[e,gp,1,i]
          #TO OPTIMIZE MULTIPLY ONLY NON SPARSE ELEMENTS
          self.f_e[e] +=  np.dot(B.T,self.sigma[e,gp,0:2,0:2]) *  self.detJ[e,gp] * self.gauss_weights[gp]
      # print ("forces")
      # print (self.f_e)
      # return forces


##strain_rate = calc_strain_rate(v)
#J, detJ, dNdX = calc_jacobian(x)
#print ("Jacobian\n", J[0])
#print ("det J \n", detJ[0])
#vol_0 = calc_vol(detJ)
#nod_mass = vol_0 * rho / nodxelem 
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
    
    self.impose_bc()
    
    print ("vel after bc", self.v)
    while (t < tf):
      print ("Time: ", t)
      # !!! PREDICTION PHASE
      u_ = dt * (self.v + (0.5 - beta) * dt * self.prev_a)
      # !!! CAN BE UNIFIED AT THE END OF STEP by v= (a(t+dt)+a(t))/2. but is not convenient for variable time step
      self.v +=  (1.0-gamma)* dt * self.prev_a
      self.a[:,:] = 0.0
      self.impose_bc()

      self.calc_jacobian()

      self.calc_str_rate()
      self.calc_pressure(dt)
      
      self.calc_stress(dt)

      #strain =  strain + calc_strain(str_rate,dt)
      # print ("strain \n",strain)

      # print ("stress\n",stress)
      #print ("a ", self.a)
      self.calc_forces()
      self.AssemblyForces()
      
      for d in range (self.dim):
        self.a[:,d] = -self.f_i[:,d]/self.mass[:]

      #print ("a ", self.a)
      self.a -= alpha * self.prev_a
      self.a = self.a / (1.0 - alpha)
      self.v += gamma * dt * self.a  

      self.impose_bc () #REINFORCE VELOCITY BC

      u_ += beta * dt * dt * self.a   
      # nod%u = nod%u + u
      self.x += u_
      
      # !call AverageData(elem%rho(:,1),nod%rho(:))  
      self.prev_a = self.a
      # time = time + dt

      self.u += u_ 
      t += dt
    #str_rate = calc_str_rate (dNdX,v)
    print ("END")
    print ("DISPLACEMENTS\n",self.u)
    print ("VEL\n",self.v)
    print ("ACC\n",self.a)
    # print (strain)
    print("STRESS")
    print (self.sigma)
    print("TAU")
    print (self.tau)
    
    # print("strain rate:\n" ,str_rate[0])


    