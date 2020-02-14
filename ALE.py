# ==========================================
# Code created by Leandro Marques at 02/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to apply ALE scheme

import sys
import numpy as np


def rotate(_npoints, _x, _y, _dirichlet_pts, _neumann_edges):
 vx_mesh = np.zeros([_npoints,1], dtype = float)
 vy_mesh = np.zeros([_npoints,1], dtype = float)

 for i in range(0, _npoints): 
  vx_mesh[i] = _y[i]*0.001
  vy_mesh[i] = - _x[i]*0.001


 try: 
  for i in range(0, len(_dirichlet_pts)):
   v1 = _dirichlet_pts[i][1] - 1
   v2 = _dirichlet_pts[i][2] - 1
 

 except KeyError:
  pass

 try: 
  for i in range(0,len(_neumann_edges)):
   v1 = _neumann_edges[i][1] - 1
   v2 = _neumann_edges[i][2] - 1
 


 except KeyError:
  pass

 return vx_mesh, vy_mesh


def coord_mesh(_x, _y, _vxx, _vyy, _dirichlet_pts, _neumann_edges, _dt):
 x = _x - _vxx*_dt
 y = _y - _vyy*_dt

 try: 
  for i in range(0, len(_dirichlet_pts)):
   v1 = _dirichlet_pts[i][1] - 1
   v2 = _dirichlet_pts[i][2] - 1
 
   x[v1] = 0.0
   x[v2] = 0.0
   
   y[v1] = 0.0
   y[v2] = 0.0

 except KeyError:
  pass

 try: 
  for i in range(0,len(_neumann_edges)):
   v1 = _neumann_edges[i][1] - 1
   v2 = _neumann_edges[i][2] - 1
 
   x[v1] = 0.0
   x[v2] = 0.0


   y[v1] = 0.0
   y[v2] = 0.0


 except KeyError:
  pass


 return x, y
