# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used for to compute boundary condition


# ------------------------------------------------------------------
# Use:

# Applying vx condition
#condition_xvelocity = bc_apply.Linear(mesh.npoints,mesh.x,mesh.y,bc)
#condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
#condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
#condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)

# Applying psi condition
#condition_psi = bc_apply2.Poiseuille(mesh.npoints,mesh.x,mesh.y,bc)
#condition_psi.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
# ------------------------------------------------------------------


import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg


class Poiseuille:
 def __init__(_self, _npoints, _x, _y, _bc):
  _self.npoints = _npoints
  _self.x = _x
  _self.y = _y
  _self.bc = _bc

 
 def neumann_condition(_self, _neumann_edges):
  _self.bc_neumann = np.zeros([_self.npoints,1], dtype = float) 
  _self.neumann_edges = _neumann_edges 
 
  for i in range(0, len(_self.neumann_edges)):
   line = _self.neumann_edges[i][0] - 1
   v1 = _self.neumann_edges[i][1] - 1
   v2 = _self.neumann_edges[i][2] - 1

   x = _self.x[v1] - _self.x[v2]
   y = _self.y[v1] - _self.y[v2]
   length = np.sqrt(x**2 + y**2)
  
   _self.bc_neumann[v1] += (_self.bc[line]*length) / 2. 
   _self.bc_neumann[v2] += (_self.bc[line]*length) / 2. 


 def dirichlet_condition(_self, _dirichlet_pts):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  #_self.bc_1 = np.zeros([1,_self.npoints], dtype = float) #For numpy array solve
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.dirichlet_pts = _dirichlet_pts
 

  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1

   _self.bc_1[v1] = _self.bc[line]
   _self.bc_1[v2] = _self.bc[line]

   _self.bc_neumann[v1] = 0.0 #Dirichlet condition is preferential
   _self.bc_neumann[v2] = 0.0 #Dirichlet condition is preferential

   _self.ibc.append(v1)
   _self.ibc.append(v2)
   
  _self.ibc = np.unique(_self.ibc)


 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.neighbors_nodes = _neighbors_nodes

#  # Scipy sparse - Method 1
#  _self.LHS = _self.LHS.tolil() #For scipy array
#  for i in range(0,len(_self.ibc)):
#   mm = _self.ibc[i]
#   #_self.bc_dirichlet -= _self.LHS[:,mm]*_self.bc_1[mm] #For numpy array
#   _self.bc_dirichlet -= _self.LHS[:,mm].todense()*_self.bc_1[mm] #For scipy array
#   _self.LHS[:,mm] = 0.0
#   _self.LHS[mm,:] = 0.0
#   _self.LHS[mm,mm] = 1.0
#   #_self.bc_dirichlet[0][mm] = _self.bc_1[mm] #For numpy array solve
#   _self.bc_dirichlet[mm] = _self.bc_1[mm] #For scipy array solve
#   _self.bc_2[mm] = 0.0
#  #_self.bc_1 = np.transpose(_self.bc_1) #For numpy array solve

#  # Scipy sparse - Method 2
#  for i in range(0,len(_self.ibc)):
#   mm = _self.ibc[i]
# 
#   for nn in _self.neighbors_nodes[mm]:
#    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
#    _self.LHS[nn,mm] = 0.0
#    _self.LHS[mm,nn] = 0.0
#   
#   _self.LHS[mm,mm] = 1.0
#   _self.bc_dirichlet[mm] = _self.bc_1[mm]
#   _self.bc_2[mm] = 0.0
 
# # Scipy sparse - Method 3
  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 


 def streamfunction_condition(_self, _dirichlet_pts, _LHS0, _neighbors_nodes):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.dirichlet_pts = _dirichlet_pts
  _self.neighbors_nodes = _neighbors_nodes

  # Dirichlet condition
  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1

   if line == 8:
    _self.bc_1[v1] = 0.0
    _self.bc_1[v2] = 0.0
 
    _self.ibc.append(v1)
    _self.ibc.append(v2)

   elif line == 11:
    _self.bc_1[v1] = 1.0
    _self.bc_1[v2] = 1.0

    _self.ibc.append(v1)
    _self.ibc.append(v2)

   else:
    _self.bc_1[v1] = _self.y[v1]
    _self.bc_1[v2] = _self.y[v2]

    _self.ibc.append(v1)
    _self.ibc.append(v2)

  _self.ibc = np.unique(_self.ibc)


  # Gaussian elimination for psi
  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 
 



