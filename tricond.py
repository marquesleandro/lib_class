import numpy as np

def b_bc(_npoints, _x, _y, _bc, _neumann_edges,_dirichlet_pts, _LHS, _neighbors_nodes): #For Method 2
#def b_bc(_npoints, _x, _y, _bc, _neumann_edges,_dirichlet_pts, _LHS): #For Method 1
 bc_neumann = np.zeros([_npoints,1], dtype = float) 
 bc_dirichlet = np.zeros([_npoints,1], dtype = float) 
 ibc = [] 
 #bc_1 = np.zeros([1,_npoints], dtype = float) #For numpy array solve
 bc_1 = np.zeros([_npoints,1], dtype = float) #For scipy array solve
 bc_2 = np.ones([_npoints,1], dtype = float) 

 # Neumann condition 
 for i in range(0, len(_neumann_edges)):
  line = _neumann_edges[i][0] - 1
  v1 = _neumann_edges[i][1] - 1
  v2 = _neumann_edges[i][2] - 1

  x = _x[v1] - _x[v2]
  y = _y[v1] - _y[v2]
  length = np.sqrt(x**2 + y**2)
  
  bc_neumann[v1] += (_bc[line]*length) / 2. 
  bc_neumann[v2] += (_bc[line]*length) / 2. 

 # Dirichlet condition
 for i in range(0, len(_dirichlet_pts)):
  line = _dirichlet_pts[i][0] - 1
  v1 = _dirichlet_pts[i][1] - 1
  v2 = _dirichlet_pts[i][2] - 1

  bc_1[v1] = _bc[line]
  bc_1[v2] = _bc[line]

  bc_neumann[v1] = 0.0 #Dirichlet condition is preferential
  bc_neumann[v2] = 0.0 #Dirichlet condition is preferential

  ibc.append(v1)
  ibc.append(v2)
   
 ibc = np.unique(ibc)

 # Gaussian elimination
 # Scipy sparse - Method 1
# _LHS = _LHS.tolil() #For scipy array
# for i in range(0,len(ibc)):
#  mm = ibc[i]
#  #bc_1 -= _LHS[:,mm]*bc_dirichlet[mm] #For numpy array
#  bc_1 -= _LHS[:,mm].todense()*bc_dirichlet[mm] #For scipy array
#  _LHS[:,mm] = 0.0
#  _LHS[mm,:] = 0.0
#  _LHS[mm,mm] = 1.0
#  #bc_1[0][mm] = bc_dirichlet[mm] #For numpy array solve
#  bc_1[mm] = bc_dirichlet[mm] #For scipy array solve
#  bc_2[mm] = 0.0
# #bc_1 = np.transpose(bc_1) #For numpy array solve

# # Scipy sparse - Method 2
# for i in range(0,len(ibc)):
#  mm = ibc[i]
# 
#  for nn in _neighbors_nodes[mm]:
#   bc_1[nn] -= float(_LHS[nn,mm]*bc_dirichlet[mm])
#   _LHS[nn,mm] = 0.0
#   _LHS[mm,nn] = 0.0
#   
#  _LHS[mm,mm] = 1.0
#  bc_1[mm] = bc_dirichlet[mm]
#  bc_2[mm] = 0.0
 
# # Scipy sparse - Method 3
 for mm in ibc:
  for nn in _neighbors_nodes[mm]:
   bc_dirichlet[nn] -= float(_LHS[nn,mm]*bc_1[mm])
   _LHS[nn,mm] = 0.0
   _LHS[mm,nn] = 0.0
   
  _LHS[mm,mm] = 1.0
  bc_dirichlet[mm] = bc_1[mm]
  bc_2[mm] = 0.0
  
 
 return bc_neumann, bc_dirichlet, _LHS, bc_1, bc_2, ibc

