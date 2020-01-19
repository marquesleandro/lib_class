# ==========================================
# Code created by Leandro Marques at 01/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to solver governing equations



import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg
import semi_lagrangian 




# Diffusion implicit and convection explicit for 1D
class SemiImplicit_concentration_equation1D:
 def __init__(_self, _scheme):
  _self.scheme = _scheme

 def taylor_galerkin(_self, _c, _vx, _dt, _M, _K, _G, _LHS, _bc_dirichlet, _bc_2):
 
  _self.scheme_name = 'Taylor Galerkin' 
  
  c = _c
  vx = _vx
  dt = _dt
  
  M = _M
  K = _K
  G = _G
  
  LHS = _LHS
  bc_dirichlet = _bc_dirichlet
  bc_2 = _bc_2

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c) - np.multiply(vx,sps.lil_matrix.dot(G,c))\
                                - (dt/2.0)*np.multiply(vx,(np.multiply(vx,sps.lil_matrix.dot(K,c))))
 
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet
 
  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
  
  _self.c = c


 def semi_lagrangian_linear(_self, _npoints, _neighbors_elements, _IEN, _x, _vx, _dt, _c, _M, _LHS, _bc_dirichlet, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian' 
  
  npoints = _npoints
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  vx = _vx
  dt = _dt
  c = _c
  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_2 = _bc_2

  #c_d = semi_lagrangian.Linear1D_v2(npoints, nelem, IEN, x, vx, dt, c)
  c_d = semi_lagrangian.Linear1D(npoints, neighbors_elements, IEN, x, vx, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)
 
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c

 def semi_lagrangian_quad(_self, _npoints, _nelem, _neighbors_elements, _IEN, _x, _vx, _dt, _c, _M, _LHS, _bc_dirichlet, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian Quadratic' 
  
  npoints = _npoints
  nelem = _nelem
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  vx = _vx
  dt = _dt
  c = _c
  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_2 = _bc_2

  #c_d = semi_lagrangian.Quad1D_v2(npoints, nelem, IEN, x, vx, dt, c)
  c_d = semi_lagrangian.Quad1D(npoints, neighbors_elements, IEN, x, vx, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)
 
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c




# Diffusion implicit and convection explicit for 2D
class SemiImplicit_concentration_equation2D:
 def __init__(_self, _scheme):
  _self.scheme = _scheme

 def taylor_galerkin(_self, _c, _vx, _vy, _dt, _Re, _Sc, _M, _Kxx, _Kyx, _Kxy, _Kyy, _Gx, _Gy, _LHS, _bc_dirichlet, _bc_neumann, _bc_2):

  _self.scheme_name = 'Taylor Galerkin' 

  c = _c
  vx =  _vx
  vy =  _vy
  dt = _dt
  Re = _Re
  Sc = _Sc
  
  M = _M
  Kxx = _Kxx
  Kyx = _Kyx
  Kxy = _Kxy
  Kyy = _Kyy
  Gx = _Gx
  Gy = _Gy
  
  LHS = _LHS
  bc_dirichlet = _bc_dirichlet
  bc_neumann = _bc_neumann
  bc_2 = _bc_2

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c) - np.multiply(vx,sps.lil_matrix.dot(Gx,c))\
                                - np.multiply(vy,sps.lil_matrix.dot(Gy,c))\
                                - (dt/2.0)*np.multiply(vx,(np.multiply(vx,sps.lil_matrix.dot(Kxx,c))\
                                                         + np.multiply(vy,sps.lil_matrix.dot(Kyx,c))))\
                                - (dt/2.0)*np.multiply(vy,(np.multiply(vx,sps.lil_matrix.dot(Kxy,c))\
                                                         + np.multiply(vy,sps.lil_matrix.dot(Kyy,c))))

  RHS = RHS + (1.0/(Re*Sc))*bc_neumann
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet
  
  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))

  _self.c = c



 def semi_lagrangian_linear(_self, _npoints, _neighbors_nodes, _neighbors_elements, _IEN, _x, _y, _vx, _vy, _dt, _Re, _Sc, _c, _M, _LHS, _bc_dirichlet, _bc_neumann, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian' 
  
  npoints = _npoints
  neighbors_nodes = _neighbors_nodes
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  y = _y

  vx = _vx
  vy = _vy
  dt = _dt
  Re = _Re
  Sc = _Sc
  c = _c

  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_neumann = _bc_neumann
  bc_2 = _bc_2

  #c_d = semi_lagrangian.Linear2D_v2(npoints, nelem, IEN, x, y, vx, vy, dt, c)
  c_d = semi_lagrangian.Linear2D(npoints, neighbors_elements, IEN, x, y, vx, vy, dt, c)
  #c_d = semi_lagrangian.Linear2D_v3(npoints, neighbors_nodes, neighbors_elements, IEN, x, y, vx, vy, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)
 
  RHS = RHS + (1.0/(Re*Sc))*bc_neumann
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c


 def semi_lagrangian_mini(_self, _npoints, _nelem, _neighbors_elements, _IEN, _x, _y, _vx, _vy, _dt, _Re, _Sc, _c, _M, _LHS, _bc_dirichlet, _bc_neumann, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian Mini' 
  
  npoints = _npoints
  nelem = _nelem
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  y = _y

  vx = _vx
  vy = _vy
  dt = _dt
  c = _c
  Re = _Re
  Sc = _Sc
 
  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_neumann = _bc_neumann
  bc_2 = _bc_2

  c_d = semi_lagrangian.Mini2D(npoints, neighbors_elements, IEN, x, y, vx, vy, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)

  RHS = RHS + (1.0/(Re*Sc))*bc_neumann
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c








 def semi_lagrangian_quad(_self, _npoints, _nelem, _neighbors_elements, _IEN, _x, _y, _vx, _vy, _dt, _Re, _Sc, _c, _M, _LHS, _bc_dirichlet, _bc_neumann, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian Quadratic' 
  
  npoints = _npoints
  nelem = _nelem
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  y = _y

  vx = _vx
  vy = _vy
  dt = _dt
  c = _c
  Re = _Re
  Sc = _Sc
 
  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_neumann = _bc_neumann
  bc_2 = _bc_2

  #c_d = semi_lagrangian.Quad2D_v2(npoints, nelem, IEN, x, y, vx, vy, dt, c)
  c_d = semi_lagrangian.Quad2D(npoints, neighbors_elements, IEN, x, y, vx, vy, dt, c)
  #c_d = semi_lagrangian.Quad2D_v3(npoints, neighbors_elements, IEN, x, y, vx, vy, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)

  RHS = RHS + (1.0/(Re*Sc))*bc_neumann
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c


 def semi_lagrangian_cubic(_self, _npoints, _nelem, _neighbors_elements, _IEN, _x, _y, _vx, _vy, _dt, _Re, _Sc, _c, _M, _LHS, _bc_dirichlet, _bc_neumann, _bc_2):
  
  _self.scheme_name = 'Semi Lagrangian Cubic' 
  
  npoints = _npoints
  nelem = _nelem
  neighbors_elements = _neighbors_elements
  IEN = _IEN
  x = _x
  y = _y

  vx = _vx
  vy = _vy
  dt = _dt
  c = _c
  Re = _Re
  Sc = _Sc
 
  M = _M
  LHS = _LHS 
  bc_dirichlet = _bc_dirichlet
  bc_neumann = _bc_neumann
  bc_2 = _bc_2

  #c_d = semi_lagrangian.Cubic2D_v2(npoints, nelem, IEN, x, y, vx, vy, dt, c)
  c_d = semi_lagrangian.Cubic2D(npoints, neighbors_elements, IEN, x, y, vx, vy, dt, c)

  A = np.copy(M)/dt
  RHS = sps.lil_matrix.dot(A,c_d)

  RHS = RHS + (1.0/(Re*Sc))*bc_neumann
  RHS = np.multiply(RHS,bc_2)
  RHS = RHS - bc_dirichlet

  c = scipy.sparse.linalg.cg(LHS,RHS,c, maxiter=1.0e+05, tol=1.0e-05)
  c = c[0].reshape((len(c[0]),1))
 
  _self.c = c




