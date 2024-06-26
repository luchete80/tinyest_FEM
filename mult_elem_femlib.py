    
import fem as fe
import numpy as np

dt = 0.8e-5

dom = fe.Domain(True) #red_int
dom.mat_E  = 206.e9 #PASS TO MAT CLASS.. 
dom.mat_nu = 0.3 #PASS TO MAT CLAS..
dom.mat_rho = 7850.0

dom.calcMatConst()
dx = 0.05
#def addBoxLength(self, lx, ly, le):
dom.addBoxLength(0.1,0.1,dx)  #(self, ex, ey, le):

for i in range (dom.node_count):

    if (dom.x[i,1]<dx/4.0):
        dom.is_fix[i,1] = True
        print("Node ", i, " fixed on y")
    if (dom.x[i,1]>0.1-dx/4.0): 
        dom.is_bcv[i,1] = True
        dom.bcv   [i,1] = -1.0
        print("Node ", i, " vel on y")




#AFTER CREATE
dom.is_fix[0,:] = True

print ("vel", dom.v)
tf = 1.0*dt
tf = 1.0e-3
dom.solve(tf, dt)