import numpy as np
import fem as fe

dt = 0.8e-5
dom = fe.Domain()
dom.addBoxLength(0.1,0.1,0.05)  #(self, ex, ey, ez):
dom.solve(dt,dt)
