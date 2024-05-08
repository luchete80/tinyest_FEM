import numpy as np
import fem as fe

dt = 0.8e-5

dom = fe.Domain()
dom.mat_E  = 206.e9 #PASS TO MAT CLASS.. 
dom.mat_nu = 0.3 #PASS TO MAT CLAS..
dom.mat_rho = 7850.0

dom.calcMatConst()
dom.addBoxLength(0.1,0.1,0.05)  #(self, ex, ey, ez):

dom.solve(dt,dt)
