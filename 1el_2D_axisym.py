import numpy as np
import fem as fe
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

dom = fe.Domain() #Don't Forget parentheses
dom.addBoxLength(10,10,10)