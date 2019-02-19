# ======================================
# Code created by Leandro Marques
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ======================================

# This code is used to assemble 
# the elementary arrays


# ------------------------------------------------------------------
# Use:
# linear = trielem.Linear(mesh.x, mesh.y, mesh.IEN)
# 
# for e in range(0,mesh.nelem):
#  linear.numerical(e)
#
#  for i in range(0,3):
#   ii = mesh.IEN[e][i]
#
#   for j in range(0,3):
#    jj = mesh.IEN[e][j]
#
#    K[ii][jj] += linear.kxx[i][j] + linear.kyy[i][j]
# ------------------------------------------------------------------


# Reference Analytic: Fundamentals of the Finite
#                     Element Method for Heat Transfer
#                     and Fluid Flow - Lewis, Nithiarasu,
#                     Seetharamu - pg. 196-200
#                     For Q_elem pg. 126
#                     For 1D pg. 193


import sys
import numpy as np

class Linear1D:
 def __init__(_self, _x, _IEN):
  _self.x = _x
  _self.IEN = _IEN
  _self.NUMNODE = 2  #Linear One-dimensional Element - 2 Nodes



 def GQNUM3(_self):
  _self.NUMGP = 3  #Number of Gauss Points


  #                                l1     
  _self.GQPoints = np.array([[-0.774596669], 
                             [ 0.000000000], 
                             [ 0.774596669]])


  #                                 w
  _self.GQWeights = np.array([[0.555555556], 
                              [0.888888889], 
                              [0.555555556]])



 
 def GQNUM4(_self):
  _self.NUMGP = 4  #Number of Gauss Points


  #                               l1     
  _self.GQPoints = np.array([[-0.861136], 
                             [-0.339981], 
                             [ 0.339981], 
                             [ 0.861136]])


  #                                w
  _self.GQWeights = np.array([[0.347855], 
                              [0.652145], 
                              [0.652145], 
                              [0.347855]])




 def GQNUM5(_self):
  _self.NUMGP = 5  #Number of Gauss Points


  #                               l1     
  _self.GQPoints = np.array([[-0.906180], 
                             [-0.538469], 
                             [ 0.000000], 
                             [ 0.538469], 
                             [ 0.906180]])


  #                                w
  _self.GQWeights = np.array([[0.236927], 
                              [0.478629], 
                              [0.568889], 
                              [0.478629], 
                              [0.236927]])





 def numerical(_self,_e):

  N = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGP,1], dtype = float)

  J = np.zeros([1,1], dtype = float)
  jacobian = np.zeros([_self.NUMGP,1], dtype = float)
  
  dNdx = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGP):
    
   # Area Coordinates
   L1 = (1.0/2.0)*(1.0 - _self.GQPoints[k][0])      #L1 = (1/2)*(1 - l1)
   L2 = (1.0/2.0)*(1.0 + _self.GQPoints[k][0])      #L2 = (1/2)*(1 + l1)

   # Shape Functions
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -0.5   #dN1/dl1
   dNdl1[k][1] =  0.5   #dN2/dl1


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...

   # Jacobian Matrix
   # Lewis pag. 64 Eq. 3.108
   J[0][0] = dxdl1[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x
   # Lewis pag. 63 Eq. 3.107
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*dNdl1[k][i]

 
  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  
  # Elementary Matrices 
  for k in range(0,_self.NUMGP): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.kx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.gx[i][j] += dNdx[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]
   



 def analytic(_self, _e):
  v1 = _self.IEN[_e][0]
  v2 = _self.IEN[_e][1]
  
  dx = _self.x[v2] - _self.x[v1]

  _self.kx = (1.0/dx)*np.array([[1,-1],[-1,1]])
  _self.mass = (dx/6.0)*np.array([[2,1],[1,2]])
  _self.gx = (1.0/2.0)*np.array([[-1,1],[-1,1]])








class Linear2D:
 def __init__(_self, _x, _y, _IEN, _GAUSSPOINTS):
  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN
  _self.NUMNODE = 3  #Linear Triangle Element - 3 Nodes

  if _GAUSSPOINTS == 4:
   _self.NUMGP = 4  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.33333333333333, 0.33333333333333], 
                              [0.60000000000000, 0.20000000000000], 
                              [0.20000000000000, 0.60000000000000], 
                              [0.20000000000000, 0.20000000000000]])


   #                                    w
   _self.GQWeights = np.array([[-0.56250000000000], 
                               [0.520833333333333], 
                               [0.520833333333333], 
                               [0.520833333333333]])

   _self.polynomial_order = 'Linear Element'
   _self.gausspoints = 4

  else:
   print ""
   print " Error: Gauss Points not found"
   print ""
   sys.exit()




 def numerical(_self,_e):

  N = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGP,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMGP,1], dtype = float)
  dydl1 = np.zeros([_self.NUMGP,1], dtype = float)
  dydl2 = np.zeros([_self.NUMGP,1], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGP,1], dtype = float)
  
  dNdx = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)
  dNdy = np.zeros([_self.NUMGP,_self.NUMNODE], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGP):
    
   # Area Coordinates
   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
   L2 = _self.GQPoints[k][0]                                #L2 = l1
   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2
   N[k][2] = L3  #N3 = L3

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -1.0   #dN1/dl1
   dNdl1[k][1] =  1.0   #dN2/dl1
   dNdl1[k][2] =  0.0   #dN3/dl1

   # Shape Functions Derivatives in respect to l2
   dNdl2[k][0] = -1.0   #dN1/dl2
   dNdl2[k][1] =  0.0   #dN2/dl2
   dNdl2[k][2] =  1.0   #dN3/dl2


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    dydl1[k] += _self.y[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dxdl2[k] += _self.x[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    dydl2[k] += _self.y[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...

   # Jacobian Matrix
   # Lewis pag. 64 Eq. 3.114
   J[0][0] = dxdl1[k]
   J[0][1] = dydl1[k]
   J[1][0] = dxdl2[k]
   J[1][1] = dydl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*dydl2[k] - dNdl2[k][i]*dydl1[k])
    dNdy[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dxdl2[k] + dNdl2[k][i]*dxdl1[k])

 
  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  
  # Elementary Matrices 
  for k in range(0,_self.NUMGP): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kxx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kxy[i][j] += dNdx[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyx[i][j] += dNdy[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyy[i][j] += dNdy[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gx[i][j] += dNdx[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gy[i][j] += dNdy[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
   

 
 def analytic(_self, _e):

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.y[j]-_self.y[k]
  bj = _self.y[k]-_self.y[i]
  bk = _self.y[i]-_self.y[j]
  ci = _self.x[k]-_self.x[j]
  cj = _self.x[i]-_self.x[k]
  ck = _self.x[j]-_self.x[i]


  A = 0.5*np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
 				  [1, _self.x[j], _self.y[j]],
				  [1, _self.x[k], _self.y[k]]]))


  _self.mass = (A/12.)*np.array([[2.,1.,1.],
                                 [1.,2.,1.],
                                 [1.,1.,2.]])

  _self.q = (A/3.)*np.ones([3,1], dtype = float)
  
  _self.gx = (1./6)*np.array([[bi,bj,bk],
                              [bi,bj,bk],
                              [bi,bj,bk]]) 
   
  _self.gy = (1./6)*np.array([[ci,cj,ck],
                              [ci,cj,ck],
                              [ci,cj,ck]])

  _self.kxx = (1./(4*A))*np.array([[bi*bi,bi*bj,bi*bk],
                                   [bj*bi,bj*bj,bj*bk],
                                   [bk*bi,bk*bj,bk*bk]])

  _self.kyy = (1./(4*A))*np.array([[ci*ci,ci*cj,ci*ck],
                                   [cj*ci,cj*cj,cj*ck],
                                   [ck*ci,ck*cj,ck*ck]])

  _self.kxy = (1./(4*A))*np.array([[bi*ci,bi*cj,bi*ck],
                                   [bj*ci,bj*cj,bj*ck],
                                   [bk*ci,bk*cj,bk*ck]])

  _self.kyx = (1./(4*A))*np.array([[ci*bi,ci*bj,ci*bk],
                                   [cj*bi,cj*bj,cj*bk],
                                   [ck*bi,ck*bj,ck*bk]])




 def axisymmetric(_self, _e):
  _self.r = _self.y
  _self.z = _self.x

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.z[j] - _self.z[k]
  bj = _self.z[k] - _self.z[i]
  bk = _self.z[i] - _self.z[j]
  ci = _self.r[k] - _self.r[j]
  cj = _self.r[i] - _self.r[k]
  ck = _self.r[j] - _self.r[i]

  A = 0.5*np.linalg.det(np.array([[1, _self.r[i], _self.z[i]],
 				  [1, _self.r[j], _self.z[j]],
				  [1, _self.r[k], _self.z[k]]]))

  r = (_self.r[i] + _self.r[j] + _self.r[k])/3.

  r_vec = np.array([[_self.r[i]],
                    [_self.r[j]],
                    [_self.r[k]]])

  _self.M_elem = (A/12.)*np.array([[2.,1.,1.],
				   [1.,2.,1.],
				   [1.,1.,2.]])

  _self.Q_elem = (2*np.pi)*np.dot(_self.M_elem,r_vec)
  
  _self.Gr_elem = (1./6)*np.array([[bi,bj,bk],
                                   [bi,bj,bk],
                                   [bi,bj,bk]]) 
   
  _self.Gz_elem = (1./6)*np.array([[ci,cj,ck],
                                   [ci,cj,ck],
                                   [ci,cj,ck]])

  _self.Kr_elem = ((2*np.pi*r)/(4*A))*np.array([[bi*bi,bj*bi,bk*bi],
                                                [bi*bj,bj*bj,bk*bj],
                                                [bi*bk,bj*bk,bk*bk]])

  _self.Kz_elem = ((2*np.pi*r)/(4*A))*np.array([[ci*ci,cj*ci,ck*ci],
                                                [ci*cj,cj*cj,ck*cj],
                                                [ci*ck,cj*ck,ck*ck]])








class Linear2D_v2:
 def __init__(_self, _x, _y, _IEN):

  _self.NUMRULE = 4
  _self.NUMGLEC = 3

  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN

  _self.gqPoints = np.array([[0.33333333333333, 0.33333333333333, 0.33333333333333], 
                             [0.60000000000000, 0.20000000000000, 0.20000000000000], 
                             [0.20000000000000, 0.60000000000000, 0.20000000000000], 
                             [0.20000000000000, 0.20000000000000, 0.60000000000000]])

  _self.gqWeights = np.array([[-0.562500000000000], 
                              [0.520833333333333], 
                              [0.520833333333333], 
                              [0.520833333333333]])

  _self.phiJ = np.array([[0.33333333333333, 0.33333333333333, 0.33333333333333], 
                         [0.60000000000000, 0.20000000000000, 0.20000000000000], 
                         [0.20000000000000, 0.60000000000000, 0.20000000000000], 
                         [0.20000000000000, 0.20000000000000, 0.60000000000000]])

  _self.dphiJdl1 = np.array([[1.0, 0.0,-1.0], 
                             [1.0, 0.0,-1.0], 
                             [1.0, 0.0,-1.0], 
                             [1.0, 0.0,-1.0]])

  _self.dphiJdl2 = np.array([[0.0, 1.0,-1.0], 
                             [0.0, 1.0,-1.0], 
                             [0.0, 1.0,-1.0], 
                             [0.0, 1.0,-1.0]])
 

 def numerical(_self,_e):

  _self.mass = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.kxx = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.kxy = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.kyx = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.kyy = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.gx = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.gy = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.dx = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  _self.dy = np.zeros([_self.NUMGLEC,_self.NUMGLEC], dtype = float)
  
  localx = np.zeros([_self.NUMRULE,1], dtype = float)
  localy = np.zeros([_self.NUMRULE,1], dtype = float)

  dxdl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl2 = np.zeros([_self.NUMRULE,1], dtype = float)

  dphiJdx = np.zeros([_self.NUMRULE,_self.NUMGLEC], dtype = float)
  dphiJdy = np.zeros([_self.NUMRULE,_self.NUMGLEC], dtype = float)
 

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  jacobian = np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
                                     [1, _self.x[j], _self.y[j]],
                                     [1, _self.x[k], _self.y[k]]]))


  for k in range(0,_self.NUMRULE):
   for i in range(0,_self.NUMGLEC):
    v = _self.IEN[_e][i]
    
    localx[k] += _self.x[v]*_self.phiJ[k][i]
    localy[k] += _self.y[v]*_self.phiJ[k][i]

    dxdl1[k] += _self.x[v]*_self.dphiJdl1[k][i]
    dxdl2[k] += _self.x[v]*_self.dphiJdl2[k][i]
    dydl1[k] += _self.y[v]*_self.dphiJdl1[k][i]
    dydl2[k] += _self.y[v]*_self.dphiJdl2[k][i]


  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEC):
    dphiJdx[k][i] = (1./jacobian)*(+_self.dphiJdl1[k][i]*dydl2[k]\
                                   -_self.dphiJdl2[k][i]*dydl1[k])

    dphiJdy[k][i] = (1./jacobian)*(-_self.dphiJdl1[k][i]*dxdl2[k]\
                                   +_self.dphiJdl2[k][i]*dxdl1[k])

  
  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEC):
    for j in range(0,_self.NUMGLEC):
     _self.mass[i][j] += 0.5*(_self.phiJ[k][i]*_self.phiJ[k][j]\
                            *jacobian*_self.gqWeights[k])
    
     _self.kxx[i][j] += 0.5*(dphiJdx[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kxy[i][j] += 0.5*(dphiJdx[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyx[i][j] += 0.5*(dphiJdy[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyy[i][j] += 0.5*(dphiJdy[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.gx[i][j] += 0.5*(_self.gqPoints[k][i]*dphiJdx[k][j]\
                          *jacobian*_self.gqWeights[k])

     _self.gy[i][j] += 0.5*(_self.gqPoints[k][i]*dphiJdy[k][j]\
                          *jacobian*_self.gqWeights[k])

     _self.dx[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dy[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])



 def analytic(_self, _e):

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.y[j]-_self.y[k]
  bj = _self.y[k]-_self.y[i]
  bk = _self.y[i]-_self.y[j]
  ci = _self.x[k]-_self.x[j]
  cj = _self.x[i]-_self.x[k]
  ck = _self.x[j]-_self.x[i]


  A = 0.5*np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
 				  [1, _self.x[j], _self.y[j]],
				  [1, _self.x[k], _self.y[k]]]))


  _self.mass = (A/12.)*np.array([[2.,1.,1.],
                                 [1.,2.,1.],
                                 [1.,1.,2.]])

  _self.q = (A/3.)*np.ones([3,1], dtype = float)
  
  _self.gx = (1./6)*np.array([[bi,bj,bk],
                              [bi,bj,bk],
                              [bi,bj,bk]]) 
   
  _self.gy = (1./6)*np.array([[ci,cj,ck],
                              [ci,cj,ck],
                              [ci,cj,ck]])

  _self.kxx = (1./(4*A))*np.array([[bi*bi,bi*bj,bi*bk],
                                   [bj*bi,bj*bj,bj*bk],
                                   [bk*bi,bk*bj,bk*bk]])

  _self.kyy = (1./(4*A))*np.array([[ci*ci,ci*cj,ci*ck],
                                   [cj*ci,cj*cj,cj*ck],
                                   [ck*ci,ck*cj,ck*ck]])

  _self.kxy = (1./(4*A))*np.array([[bi*ci,bi*cj,bi*ck],
                                   [bj*ci,bj*cj,bj*ck],
                                   [bk*ci,bk*cj,bk*ck]])

  _self.kyx = (1./(4*A))*np.array([[ci*bi,ci*bj,ci*bk],
                                   [cj*bi,cj*bj,cj*bk],
                                   [ck*bi,ck*bj,ck*bk]])


 def axisymmetric(_self, _e):
  _self.r = _self.y
  _self.z = _self.x

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.z[j] - _self.z[k]
  bj = _self.z[k] - _self.z[i]
  bk = _self.z[i] - _self.z[j]
  ci = _self.r[k] - _self.r[j]
  cj = _self.r[i] - _self.r[k]
  ck = _self.r[j] - _self.r[i]

  A = 0.5*np.linalg.det(np.array([[1, _self.r[i], _self.z[i]],
 				  [1, _self.r[j], _self.z[j]],
				  [1, _self.r[k], _self.z[k]]]))

  r = (_self.r[i] + _self.r[j] + _self.r[k])/3.

  r_vec = np.array([[_self.r[i]],
                    [_self.r[j]],
                    [_self.r[k]]])

  _self.M_elem = (A/12.)*np.array([[2.,1.,1.],
				   [1.,2.,1.],
				   [1.,1.,2.]])

  _self.Q_elem = (2*np.pi)*np.dot(_self.M_elem,r_vec)
  
  _self.Gr_elem = (1./6)*np.array([[bi,bj,bk],
                                   [bi,bj,bk],
                                   [bi,bj,bk]]) 
   
  _self.Gz_elem = (1./6)*np.array([[ci,cj,ck],
                                   [ci,cj,ck],
                                   [ci,cj,ck]])

  _self.Kr_elem = ((2*np.pi*r)/(4*A))*np.array([[bi*bi,bj*bi,bk*bi],
                                                [bi*bj,bj*bj,bk*bj],
                                                [bi*bk,bj*bk,bk*bk]])

  _self.Kz_elem = ((2*np.pi*r)/(4*A))*np.array([[ci*ci,cj*ci,ck*ci],
                                                [ci*cj,cj*cj,ck*cj],
                                                [ci*ck,cj*ck,ck*ck]])









class Mini:
 def __init__(_self, _x, _y, _IEN):

  _self.NUMRULE = 12
  _self.NUMGLEU = 4 
  _self.NUMGLEP = 3 

  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN

  _self.gqPoints= np.array([[0.873821971016996, 0.063089014491502, 0.063089014491502],
                            [0.063089014491502, 0.873821971016996, 0.063089014491502],
                            [0.063089014491502, 0.063089014491502, 0.873821971016996],
                            [0.501426509658179, 0.249286745170910, 0.249286745170911],
                            [0.249286745170910, 0.501426509658179, 0.249286745170911],
                            [0.249286745170910, 0.249286745170910, 0.50142650965818 ],
                            [0.636502499121399, 0.310352451033785, 0.053145049844816],
                            [0.636502499121399, 0.053145049844816, 0.310352451033785],
                            [0.310352451033785, 0.636502499121399, 0.053145049844816],
                            [0.310352451033785, 0.053145049844816, 0.636502499121399],
                            [0.053145049844816, 0.636502499121399, 0.310352451033785],
                            [0.053145049844816, 0.310352451033785, 0.636502499121399]])

  _self.gqWeights = np.array([[0.050844906370207],
                              [0.050844906370207],
                              [0.050844906370207],
                              [0.116786275726379],
                              [0.116786275726379],
                              [0.116786275726379],
                              [0.082851075618374],
                              [0.082851075618374],
                              [0.082851075618374],
                              [0.082851075618374],
                              [0.082851075618374],
                              [0.082851075618374]])

  _self.phiJ = np.array([[ 0.842519908360035, 0.031786951834541, 0.031786951834541, 0.093906187970883],
                         [ 0.031786951834541, 0.842519908360035, 0.031786951834541, 0.093906187970883],
                         [ 0.031786951834541, 0.031786951834541, 0.842519908360035, 0.093906187970883],
                         [ 0.220981204105529,-0.03115856038174 ,-0.031158560381739, 0.841335916657949],
                         [-0.03115856038174 , 0.220981204105529,-0.031158560381739, 0.841335916657949],
                         [-0.031158560381739,-0.031158560381739, 0.220981204105531, 0.841335916657947],
                         [ 0.542017987859968, 0.215867939772354,-0.041339461416615, 0.283453533784293],
                         [ 0.542017987859968,-0.041339461416615, 0.215867939772354, 0.283453533784293],
                         [ 0.215867939772354, 0.542017987859968,-0.041339461416615, 0.283453533784293],
                         [ 0.215867939772354,-0.041339461416615, 0.542017987859968, 0.283453533784293],
                         [-0.041339461416615, 0.542017987859968, 0.215867939772354, 0.283453533784293],
                         [-0.041339461416615, 0.215867939772354, 0.542017987859968, 0.283453533784293]])

  _self.dphiJdl1 = np.array([[ 1.460335089186776,  4.603350891867764e-01,-0.539664910813224,-1.381005267560329e+00],
                             [ 1.               ,  3.885780586188048e-16,-1.               ,-1.332267629550188e-15],
                             [ 0.539664910813224, -4.603350891867763e-01,-1.460335089186776, 1.381005267560329e+00],
                             [ 1.565695910954718,  5.656959109547176e-01,-0.434304089045282,-1.697087732864153e+00],
                             [ 0.999999999999996, -4.440892098500626e-15,-1.000000000000004, 1.332267629550188e-14],
                             [ 0.434304089045278, -5.656959109547223e-01,-1.565695910954722, 1.697087732864166e+00],
                             [ 2.62941772790624 ,  1.629417727906240e+00, 0.62941772790624 ,-4.888253183718719e+00],
                             [ 1.155999345062549,  1.559993450625485e-01,-0.844000654937452,-4.679980351876453e-01],
                             [ 2.473418382843691,  1.473418382843691e+00, 0.473418382843691,-4.420255148531075e+00],
                             [ 0.844000654937452, -1.559993450625484e-01,-1.155999345062548, 4.679980351876453e-01],
                             [-0.473418382843691, -1.473418382843691e+00,-2.473418382843691, 4.420255148531074e+00],
                             [-0.62941772790624 , -1.629417727906240e+00,-2.62941772790624 , 4.888253183718721e+00]])

  _self.dphiJdl2 = np.array([[ 4.440892098500626e-16, 1.               ,-1.               ,-1.332267629550188e-15],
                             [ 4.603350891867763e-01, 1.460335089186776,-0.539664910813224,-1.381005267560329e+00],
                             [-4.603350891867763e-01, 0.539664910813224,-1.460335089186776, 1.381005267560329e+00],
                             [-4.440892098500626e-15, 0.999999999999996,-1.000000000000004, 1.332267629550188e-14],
                             [ 5.656959109547176e-01, 1.565695910954718,-0.434304089045282,-1.697087732864153e+00],
                             [-5.656959109547223e-01, 0.434304089045278,-1.565695910954722, 1.697087732864166e+00],
                             [ 1.473418382843691e+00, 2.473418382843691, 0.473418382843691,-4.420255148531074e+00],
                             [-1.473418382843691e+00,-0.473418382843691,-2.473418382843691, 4.420255148531074e+00],
                             [ 1.629417727906240e+00, 2.62941772790624 , 0.62941772790624 ,-4.888253183718720e+00],
                             [-1.629417727906240e+00,-0.62941772790624 ,-2.62941772790624 , 4.888253183718721e+00],
                             [ 1.559993450625485e-01, 1.155999345062549,-0.844000654937452,-4.679980351876454e-01],
                             [-1.559993450625484e-01, 0.844000654937452,-1.155999345062548, 4.679980351876453e-01]])

  _self.dgqPointsdl1 = np.array([[1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1]])

  _self.dgqPointsdl2 = np.array([[0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1]])





 def numerical(_self,_e):

  _self.mass = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.gx = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.gy = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.dx = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  _self.dy = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  
  localx = np.zeros([_self.NUMRULE,1], dtype = float)
  localy = np.zeros([_self.NUMRULE,1], dtype = float)

  dxdl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl2 = np.zeros([_self.NUMRULE,1], dtype = float)

  dphiJdx = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
  dphiJdy = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
 

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  jacobian = np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
                                     [1, _self.x[j], _self.y[j]],
                                     [1, _self.x[k], _self.y[k]]]))


  for k in range(0,_self.NUMRULE):
   for i in range(0,_self.NUMGLEU):
    v = _self.IEN[_e][i]
    
    localx[k] += _self.x[v]*_self.phiJ[k][i]
    localy[k] += _self.y[v]*_self.phiJ[k][i]

    dxdl1[k] += _self.x[v]*_self.dphiJdl1[k][i]
    dxdl2[k] += _self.x[v]*_self.dphiJdl2[k][i]
    dydl1[k] += _self.y[v]*_self.dphiJdl1[k][i]
    dydl2[k] += _self.y[v]*_self.dphiJdl2[k][i]


  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    dphiJdx[k][i] = (1./jacobian)*(+_self.dphiJdl1[k][i]*dydl2[k]\
                                   -_self.dphiJdl2[k][i]*dydl1[k])

    dphiJdy[k][i] = (1./jacobian)*(-_self.dphiJdl1[k][i]*dxdl2[k]\
                                   +_self.dphiJdl2[k][i]*dxdl1[k])

  
  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    for j in range(0,_self.NUMGLEU):
     _self.mass[i][j] += 0.5*(_self.phiJ[k][i]*_self.phiJ[k][j]\
                            *jacobian*_self.gqWeights[k])
    
     _self.kxx[i][j] += 0.5*(dphiJdx[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])
     
     _self.kxy[i][j] += 0.5*(dphiJdx[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyx[i][j] += 0.5*(dphiJdy[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyy[i][j] += 0.5*(dphiJdy[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])
    
    for j in range(0,_self.NUMGLEP):
     _self.gx[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.gy[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dx[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dy[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])

class Quad:
 def __init__(_self, _x, _y, _IEN):

  _self.NUMRULE = 12
  _self.NUMGLEU = 6 
  _self.NUMGLEP = 3 

  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN

  _self.gqPoints = np.array([[0.249286745171, 0.249286745171, 0.501426509658 ],
                             [0.249286745171, 0.501426509658, 0.249286745171 ],
                             [0.501426509658, 0.249286745171, 0.249286745171 ],
                             [0.063089014492, 0.063089014492, 0.873821971017 ],
                             [0.063089014492, 0.873821971017, 0.063089014491 ],
                             [0.873821971017, 0.063089014492, 0.063089014491 ],
                             [0.310352451034, 0.636502499121, 0.053145049845 ],
                             [0.636502499121, 0.053145049845, 0.310352451034 ],
                             [0.053145049845, 0.310352451034, 0.636502499121 ],
                             [0.636502499121, 0.310352451034, 0.053145049845 ],
                             [0.310352451034, 0.053145049845, 0.636502499121 ],
                             [0.053145049845, 0.636502499121, 0.310352451034 ]])

  _self.gqWeights = np.array([[0.116786275726],
                              [0.116786275726],
                              [0.116786275726],
                              [0.050844906370],
                              [0.050844906370],
                              [0.050844906370],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618]])

  _self.phiJ = np.array([[-0.124998982535, -0.124998982535, 0.001430579518, 0.248575525272, 0.499995930140, 0.499995930140 ],
                         [-0.124998982535, 0.001430579518, -0.124998982535, 0.499995930140, 0.499995930140, 0.248575525272 ],
                         [0.001430579518, -0.124998982535, -0.124998982535, 0.499995930140, 0.248575525272, 0.499995930140 ],
                         [-0.055128566992, -0.055128566992, 0.653307703047, 0.015920894998, 0.220514267970, 0.220514267970 ],
                         [-0.055128566992, 0.653307703047, -0.055128566992, 0.220514267970, 0.220514267970, 0.015920894998 ],
                         [0.653307703047, -0.055128566992, -0.055128566992, 0.220514267970, 0.015920894998, 0.220514267970 ],
                         [-0.117715163308, 0.173768363654, -0.047496257199, 0.790160442766, 0.135307828169, 0.065974785919 ],
                         [0.173768363654, -0.047496257199, -0.117715163308, 0.135307828169, 0.065974785919, 0.790160442766 ],
                         [-0.047496257199, -0.117715163308, 0.173768363654, 0.065974785919, 0.790160442766, 0.135307828169 ],
                         [0.173768363654, -0.117715163308, -0.047496257199, 0.790160442766, 0.065974785919, 0.135307828169 ],
                         [-0.117715163308, -0.047496257199, 0.173768363654, 0.065974785919, 0.135307828169, 0.790160442766 ],
                         [-0.047496257199, 0.173768363654, -0.117715163308, 0.135307828169, 0.790160442766, 0.065974785919 ]])



  _self.dphiJdl1 = np.array([[-0.002853019316, 0.000000000000, -1.005706038632, 0.997146980684, -0.997146980684, 1.008559057948 ],
                             [-0.002853019316, 0.000000000000, 0.002853019316, 2.005706038632, -2.005706038632, 0.000000000000 ],
                             [1.005706038632, 0.000000000000, 0.002853019316, 0.997146980684,  -0.997146980684, -1.008559057948 ],
                             [-0.747643942034, 0.000000000000, -2.495287884068, 0.252356057966, -0.252356057966, 3.242931826102 ],
                             [-0.747643942034, 0.000000000000, 0.747643942034, 3.495287884068,  -3.495287884068, -0.000000000000 ],
                             [2.495287884068, 0.000000000000, 0.747643942034, 0.252356057966, -0.252356057966, -3.242931826102 ],
                             [0.241409804136, 0.000000000000, 0.787419800620, 2.546009996484, -2.546009996484, -1.028829604756 ],
                             [1.546009996484, 0.000000000000, -0.241409804137, 0.212580199379, -0.212580199379, -1.304600192347 ],
                             [-0.787419800621, 0.000000000000, -1.546009996485, 1.241409804136, -1.241409804136, 2.333429797106 ],
                             [1.546009996484, 0.000000000000, 0.787419800620, 1.241409804136, -1.241409804136, -2.333429797104 ],
                             [0.241409804136, 0.000000000000, -1.546009996485, 0.212580199379, -0.212580199379, 1.304600192349 ],
                             [-0.787419800621, 0.000000000000, -0.241409804137, 2.546009996484, -2.546009996484, 1.028829604758 ]])


  _self.dphiJdl2 = np.array([[0.000000000000, -0.002853019316, -1.005706038632, 0.997146980684, 1.008559057948, -0.997146980684 ],
                             [0.000000000000, 1.005706038632, 0.002853019316, 0.997146980684, -1.008559057948, -0.997146980684 ],
                             [0.000000000000, -0.002853019316, 0.002853019316, 2.005706038632, -0.000000000000, -2.005706038632 ],
                             [0.000000000000, -0.747643942034, -2.495287884068, 0.252356057966, 3.242931826102, -0.252356057966 ],
                             [0.000000000000, 2.495287884068, 0.747643942034, 0.252356057966, -3.242931826102, -0.252356057966 ],
                             [0.000000000000, -0.747643942034, 0.747643942034, 3.495287884068, -0.000000000000, -3.495287884068 ],
                             [0.000000000000, 1.546009996484, 0.787419800620, 1.241409804136, -2.333429797104, -1.241409804136 ],
                             [0.000000000000, -0.787419800621, -0.241409804137, 2.546009996484, 1.028829604758, -2.546009996484 ],
                             [0.000000000000, 0.241409804136, -1.546009996485, 0.212580199379, 1.304600192349, -0.212580199379 ],
                             [0.000000000000, 0.241409804136, 0.787419800620, 2.546009996484, -1.028829604756, -2.546009996484 ],
                             [0.000000000000, -0.787419800621, -1.546009996485, 1.241409804136, 2.333429797106, -1.241409804136 ],
                             [0.000000000000, 1.546009996484, -0.241409804137, 0.212580199379, -1.304600192347, -0.212580199379 ]])



  _self.dgqPointsdl1 = np.array([[1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1]])

  _self.dgqPointsdl2 = np.array([[0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1]])

 
 def numerical(_self,_e):

  _self.mass = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.gx = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.gy = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.dx = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  _self.dy = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  
  localx = np.zeros([_self.NUMRULE,1], dtype = float)
  localy = np.zeros([_self.NUMRULE,1], dtype = float)

  dxdl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl2 = np.zeros([_self.NUMRULE,1], dtype = float)

  dphiJdx = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
  dphiJdy = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
 

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  jacobian = np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
                                     [1, _self.x[j], _self.y[j]],
                                     [1, _self.x[k], _self.y[k]]]))


  for k in range(0,_self.NUMRULE):
   for i in range(0,_self.NUMGLEU):
    v = _self.IEN[_e][i]
    
    localx[k] += _self.x[v]*_self.phiJ[k][i]
    localy[k] += _self.y[v]*_self.phiJ[k][i]

    dxdl1[k] += _self.x[v]*_self.dphiJdl1[k][i]
    dxdl2[k] += _self.x[v]*_self.dphiJdl2[k][i]
    dydl1[k] += _self.y[v]*_self.dphiJdl1[k][i]
    dydl2[k] += _self.y[v]*_self.dphiJdl2[k][i]


  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    dphiJdx[k][i] = (1./jacobian)*(+_self.dphiJdl1[k][i]*dydl2[k]\
                                   -_self.dphiJdl2[k][i]*dydl1[k])

    dphiJdy[k][i] = (1./jacobian)*(-_self.dphiJdl1[k][i]*dxdl2[k]\
                                   +_self.dphiJdl2[k][i]*dxdl1[k])

  
  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    for j in range(0,_self.NUMGLEU):
     _self.mass[i][j] += 0.5*(_self.phiJ[k][i]*_self.phiJ[k][j]\
                            *jacobian*_self.gqWeights[k])
    
     _self.kxx[i][j] += 0.5*(dphiJdx[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])
     
     _self.kxy[i][j] += 0.5*(dphiJdx[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyx[i][j] += 0.5*(dphiJdy[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyy[i][j] += 0.5*(dphiJdy[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])
    
    for j in range(0,_self.NUMGLEP):
     _self.gx[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.gy[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dx[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dy[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])


class QuadBubble:
 def __init__(_self, _x, _y, _IEN):

  _self.NUMRULE = 12
  _self.NUMGLEU = 7 
  _self.NUMGLEP = 3 

  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN

  _self.gqPoints = np.array([[0.249286745171, 0.249286745171, 0.501426509658 ],
                             [0.249286745171, 0.501426509658, 0.249286745171 ],
                             [0.501426509658, 0.249286745171, 0.249286745171 ],
                             [0.063089014492, 0.063089014492, 0.873821971017 ],
                             [0.063089014492, 0.873821971017, 0.063089014491 ],
                             [0.873821971017, 0.063089014492, 0.063089014491 ],
                             [0.310352451034, 0.636502499121, 0.053145049845 ],
                             [0.636502499121, 0.053145049845, 0.310352451034 ],
                             [0.053145049845, 0.310352451034, 0.636502499121 ],
                             [0.636502499121, 0.310352451034, 0.053145049845 ],
                             [0.310352451034, 0.053145049845, 0.636502499121 ],
                             [0.053145049845, 0.636502499121, 0.310352451034 ]])

  _self.gqWeights = np.array([[0.116786275726],
                              [0.116786275726],
                              [0.116786275726],
                              [0.050844906370],
                              [0.050844906370],
                              [0.050844906370],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618]])

  _self.phiJ = np.array([[-0.031517214017514, -0.031517214017514, 0.094912348035192, -0.125351548798530, 0.126068856070057, 0.126068856070057, 0.841335916658253 ],
                         [-0.031517214017514, 0.094912348035192, -0.031517214017514, 0.126068856070057, 0.126068856070057, -0.125351548798530, 0.841335916658253 ],
                         [0.094912348035192, -0.031517214017514, -0.031517214017514, 0.126068856070057, -0.125351548798530, 0.126068856070057, 0.841335916658253 ],
                         [-0.044694546106830, -0.044694546106830, 0.663741723932723, -0.025815188544578, 0.178778184427318, 0.178778184427318, 0.093906187970878 ],
                         [-0.044694546106830, 0.663741723932723, -0.044694546106829, 0.178778184427318, 0.178778184427318, -0.025815188544578, 0.093906187970878 ],
                         [0.663741723932723, -0.044694546106830, -0.044694546106829, 0.178778184427318, -0.025815188544578, 0.178778184427318, 0.093906187970878 ],
                         [-0.086220326221122, 0.205263200740812, -0.016001420111690, 0.664181094416856, 0.009328479819991, -0.060004562430140, 0.283453533785293 ],
                         [0.205263200740714, -0.016001420111631, -0.086220326221172, 0.009328479819875, -0.060004562429953, 0.664181094417758, 0.283453533784409 ],
                         [-0.016001420111641, -0.086220326221231, 0.205263200741013, -0.060004562429953, 0.664181094417539, 0.009328479819959, 0.283453533784315 ],
                         [0.205263200740812, -0.086220326221122, -0.016001420111690, 0.664181094416856, -0.060004562430140, 0.009328479819991, 0.283453533785293 ],
                         [-0.086220326221231, -0.016001420111641, 0.205263200741013, -0.060004562429953, 0.009328479819959, 0.664181094417539, 0.283453533784315 ],
                         [-0.016001420111631, 0.205263200740714, -0.086220326221172, 0.009328479819875, 0.664181094417758, -0.060004562429953, 0.283453533784409 ]])


  _self.dphiJdl1 = np.array([[0.185712284335440, 0.188565303651440, -0.817140734980560, 0.242885766078239, -1.751408195289760, 0.254297843342239, 1.697087732862961    ],
                             [-0.002853019316000, 0.000000000000000, 0.002853019316000, 2.005706038632000, -2.005706038632000, -0.000000000000000, 0.000000000000000   ],
                             [0.817140734980560, -0.188565303651440, -0.185712284335440, 1.751408195289760, -0.242885766078240, -0.254297843342239, -1.697087732862961 ],
                             [-0.594198912305078, 0.153445029728922, -2.341842854339078, -0.361424060949687, -0.866136176881687, 2.629151707186313, 1.381005267560296  ],
                             [-0.747643942034000, -0.000000000000000, 0.747643942034000, 3.495287884068000, -3.495287884068000, 0.000000000000000, -0.000000000000001  ],
                             [2.341842854339078, -0.153445029728922, 0.594198912305078, 0.866136176881687, 0.361424060949687, -2.629151707186313, -1.381005267560296   ],
                             [-0.249729656811648, -0.491139460947648, 0.296280339672351, 4.510567840274594, -0.581452152693406, 0.935728239034594, -4.420255148528836  ],
                             [1.494010214796629, -0.051999781687371, -0.293409585824171, 0.420579326128683, -0.004581072629717, -1.096601065597717, -0.467998035186336 ],
                             [-0.244280557985181, 0.543139242635619, -1.002870753849181, -0.931147166406477, -3.413966774678477, 0.160872826563123, 4.888253183720574  ],
                             [1.002870753848753, -0.543139242635247, 0.244280557984753, 3.413966774676987, 0.931147166404987, -0.160872826563013, -4.888253183717221   ],
                             [0.293409585823434, 0.051999781687434, -1.494010214797366, 0.004581072629462, -0.420579326128938, 1.096601065599062, 0.467998035186910    ],
                             [-0.296280339672388, 0.491139460948412, 0.249729656811612, 0.581452152690350, -4.510567840277649, -0.935728239036049, 4.420255148535712   ]])

  _self.dphiJdl2 = np.array([[0.188565303651440, 0.185712284335440, -0.817140734980560, 0.242885766078239, 0.254297843342239, -1.751408195289760, 1.697087732862961    ],
                             [-0.188565303651440, 0.817140734980560, -0.185712284335440, 1.751408195289760, -0.254297843342239, -0.242885766078240, -1.697087732862961 ],
                             [-0.000000000000000, -0.002853019316000, 0.002853019316000, 2.005706038632000, -0.000000000000000, -2.005706038632000, -0.000000000000001 ],
                             [0.153445029728922, -0.594198912305078, -2.341842854339078, -0.361424060949687, 2.629151707186313, -0.866136176881687, 1.381005267560296  ],  
                             [-0.153445029728922, 2.341842854339078, 0.594198912305078, 0.866136176881687, -2.629151707186313, 0.361424060949687, -1.381005267560296   ],
                             [-0.000000000000000, -0.747643942034000, 0.747643942034000, 3.495287884068000, -0.000000000000001, -3.495287884068000, -0.000000000000001 ],
                             [-0.543139242635247, 1.002870753848753, 0.244280557984753, 3.413966774676987, -0.160872826563013, 0.931147166404987, -4.888253183717221   ],
                             [0.491139460948412, -0.296280339672388, 0.249729656811612, 0.581452152690350, -0.935728239036049, -4.510567840277649, 4.420255148535712   ],
                             [0.051999781687434, 0.293409585823434, -1.494010214797366, 0.004581072629462, 1.096601065599062, -0.420579326128938, 0.467998035186910    ],
                             [-0.491139460947648, -0.249729656811648, 0.296280339672351, 4.510567840274594, 0.935728239034594, -0.581452152693406, -4.420255148528836  ],
                             [0.543139242635619, -0.244280557985181, -1.002870753849181, -0.931147166406477, 0.160872826563123, -3.413966774678477, 4.888253183720574  ],
                             [-0.051999781687371, 1.494010214796629, -0.293409585824171, 0.420579326128683, -1.096601065597717, -0.004581072629717, -0.467998035186336 ]])


  _self.dgqPointsdl1 = np.array([[1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1]])

  _self.dgqPointsdl2 = np.array([[0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1]])

 def numerical(_self,_e):

  _self.mass = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.gx = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.gy = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.dx = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  _self.dy = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  
  localx = np.zeros([_self.NUMRULE,1], dtype = float)
  localy = np.zeros([_self.NUMRULE,1], dtype = float)

  dxdl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl2 = np.zeros([_self.NUMRULE,1], dtype = float)

  dphiJdx = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
  dphiJdy = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
 

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  jacobian = np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
                                     [1, _self.x[j], _self.y[j]],
                                     [1, _self.x[k], _self.y[k]]]))


  for k in range(0,_self.NUMRULE):
   for i in range(0,_self.NUMGLEU):
    v = _self.IEN[_e][i]
    
    localx[k] += _self.x[v]*_self.phiJ[k][i]
    localy[k] += _self.y[v]*_self.phiJ[k][i]

    dxdl1[k] += _self.x[v]*_self.dphiJdl1[k][i]
    dxdl2[k] += _self.x[v]*_self.dphiJdl2[k][i]
    dydl1[k] += _self.y[v]*_self.dphiJdl1[k][i]
    dydl2[k] += _self.y[v]*_self.dphiJdl2[k][i]


  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    dphiJdx[k][i] = (1./jacobian)*(+_self.dphiJdl1[k][i]*dydl2[k]\
                                   -_self.dphiJdl2[k][i]*dydl1[k])

    dphiJdy[k][i] = (1./jacobian)*(-_self.dphiJdl1[k][i]*dxdl2[k]\
                                   +_self.dphiJdl2[k][i]*dxdl1[k])

  
  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    for j in range(0,_self.NUMGLEU):
     _self.mass[i][j] += 0.5*(_self.phiJ[k][i]*_self.phiJ[k][j]\
                            *jacobian*_self.gqWeights[k])
    
     _self.kxx[i][j] += 0.5*(dphiJdx[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])
     
     _self.kxy[i][j] += 0.5*(dphiJdx[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyx[i][j] += 0.5*(dphiJdy[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyy[i][j] += 0.5*(dphiJdy[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])
    
    for j in range(0,_self.NUMGLEP):
     _self.gx[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.gy[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dx[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dy[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])


class Cubic:
 def __init__(_self, _x, _y, _IEN):

  _self.NUMRULE = 12
  _self.NUMGLEU = 10 
  _self.NUMGLEP = 3 

  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN

  _self.gqPoints = np.array([[0.249286745171, 0.249286745171, 0.501426509658 ],
                             [0.249286745171, 0.501426509658, 0.249286745171 ],
                             [0.501426509658, 0.249286745171, 0.249286745171 ],
                             [0.063089014492, 0.063089014492, 0.873821971017 ],
                             [0.063089014492, 0.873821971017, 0.063089014491 ],
                             [0.873821971017, 0.063089014492, 0.063089014491 ],
                             [0.310352451034, 0.636502499121, 0.053145049845 ],
                             [0.636502499121, 0.053145049845, 0.310352451034 ],
                             [0.053145049845, 0.310352451034, 0.636502499121 ],
                             [0.636502499121, 0.310352451034, 0.053145049845 ],
                             [0.310352451034, 0.053145049845, 0.636502499121 ],
                             [0.053145049845, 0.636502499121, 0.310352451034 ]])

  _self.gqWeights = np.array([[0.116786275726],
                              [0.116786275726],
                              [0.116786275726],
                              [0.050844906370],
                              [0.050844906370],
                              [0.050844906370],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618],
                              [0.082851075618]])

  _self.phiJ = np.array([[0.039351685817, 0.039351685817, -0.062673722052, -0.070510246199, -0.070510246199, -0.141827463079, 0.283654926158, 0.283654926158, -0.141827463079, 0.841335916658 ],
                         [0.039351685817, -0.062673722052, 0.039351685817, -0.141827463079, 0.283654926158, 0.283654926158, -0.141827463079, -0.070510246199, -0.070510246199, 0.841335916658 ],
                         [-0.062673722052, 0.039351685817, 0.039351685817, 0.283654926158, -0.141827463079, -0.070510246199, -0.070510246199, -0.141827463079, 0.283654926158, 0.841335916658 ],
                         [0.046307995391, 0.046307995391, 0.440268993399, -0.014521043556, -0.014521043556, -0.201125457481, 0.402250914961, 0.402250914961, -0.201125457481, 0.093906187971  ],
                         [0.046307995391, 0.440268993399, 0.046307995391, -0.201125457481, 0.402250914961, 0.402250914961, -0.201125457481, -0.014521043556, -0.014521043556, 0.093906187971  ],
                         [0.440268993399, 0.046307995391, 0.046307995391, 0.402250914961, -0.201125457481, -0.014521043556, -0.014521043556, -0.201125457481, 0.402250914961, 0.093906187971  ],
                         [0.011435826065, -0.026193226600, 0.041110728467, -0.061285221448, 0.808488952667, 0.138446419693, -0.127951879896, -0.062388096818, -0.005117035916, 0.283453533785 ],
                         [-0.026193226600, 0.041110728466, 0.011435826065, 0.138446419693, -0.127951879895, -0.062388096818, -0.005117035916, -0.061285221448, 0.808488952668, 0.283453533784 ],
                         [0.041110728466, 0.011435826065, -0.026193226600, -0.062388096818, -0.005117035916, -0.061285221448, 0.808488952668, 0.138446419693, -0.127951879895, 0.283453533784 ],
                         [-0.026193226600, 0.011435826065, 0.041110728467, 0.808488952667, -0.061285221448, -0.005117035916, -0.062388096818, -0.127951879896, 0.138446419693, 0.283453533785 ],
                         [0.011435826065, 0.041110728466, -0.026193226600, -0.005117035916, -0.062388096818, -0.127951879895, 0.138446419693, 0.808488952668, -0.061285221448, 0.283453533784 ],
                         [0.041110728466, -0.026193226600, 0.011435826065, -0.127951879895, 0.138446419693, 0.808488952668, -0.061285221448, -0.005117035916, -0.062388096818, 0.283453533784 ]])


  _self.dphiJdl1 = np.array([[-0.404638308747, 0.000000000000, 0.118553234987, 0.556094442315, -0.282847955477, 0.282847955477, -2.253182175178, -1.115316116704, 1.401401190464, 1.697087732863 ],
                             [-0.404638308747, 0.000000000000, 0.404638308747, 1.118553234987, 1.137866058474, -1.137866058474, -1.118553234987, -0.838942397792, 0.838942397792, 0.000000000000 ],
                             [-0.118553234987, 0.000000000000, 0.404638308747, 2.253182175178, -0.282847955477, 0.282847955477, -0.556094442315, -1.401401190464, 1.115316116704, -1.697087732863],
                             [0.485931890195, 0.000000000000, -3.443727560779, -0.176434523975, -0.230167544593, 0.230167544593, -1.204570743585, 5.171355686771, -2.213560016186, 1.381005267560],
                             [0.485931890195, 0.000000000000, -0.485931890195, -2.443727560779, 6.375926430356, -6.375926430356, 2.443727560779, -0.053733020618, 0.053733020618, -0.000000000000],
                             [3.443727560779, 0.000000000000, -0.485931890195, 1.204570743585, -0.230167544593, 0.230167544593, 0.176434523975, 2.213560016186, -5.171355686771, -1.381005267560 ],
                             [-0.492870367158, 0.000000000000, -0.559823901756, 2.469321742625, 2.605067077684, -2.605067077684, 1.950933405904, 0.750232850759, 0.302461418155, -4.420255148529 ],
                             [0.740805831639, 0.000000000000, 0.492870367158, 0.674175115836, -0.201023373941, 0.201023373941, -0.206177080649, -2.565606080133, 1.331929881335, -0.467998035186 ],
                             [0.559823901757, 0.000000000000, -0.740805831641, -0.951256224702, -0.096284337505, 0.096284337505, -3.936996959018, 1.930891961850, -1.749910031967, 4.888253183721],
                             [0.740805831639, 0.000000000000, -0.559823901756, 3.936996959017, -0.096284337505, 0.096284337505, 0.951256224701, 1.749910031962, -1.930891961845, -4.888253183717 ],
                             [-0.492870367158, 0.000000000000, -0.740805831641, 0.206177080649, -0.201023373941, 0.201023373941, -0.674175115836, -1.331929881332, 2.565606080131, 0.467998035187],
                             [0.559823901757, 0.000000000000, 0.492870367158, -1.950933405907, 2.605067077684, -2.605067077684, -2.469321742629, -0.302461418154, -0.750232850762, 4.420255148536]])

  _self.dphiJdl2 = np.array([[0.000000000000, -0.404638308747, 0.118553234987, -0.282847955477, 0.556094442315, 1.401401190464, -1.115316116704, -2.253182175178, 0.282847955477, 1.697087732863  ],
                             [0.000000000000, -0.118553234987, 0.404638308747, -0.282847955477, 2.253182175178, 1.115316116704, -1.401401190464, -0.556094442315, 0.282847955477, -1.697087732863 ],
                             [0.000000000000, -0.404638308747, 0.404638308747, 1.137866058474, 1.118553234987, 0.838942397792, -0.838942397792, -1.118553234987, -1.137866058474, -0.000000000000 ],
                             [0.000000000000, 0.485931890195, -3.443727560779, -0.230167544593, -0.176434523975, -2.213560016186, 5.171355686771, -1.204570743585, 0.230167544593, 1.381005267560 ],
                             [0.000000000000, 3.443727560779, -0.485931890195, -0.230167544593, 1.204570743585, -5.171355686771, 2.213560016186, 0.176434523975, 0.230167544593, -1.381005267560  ],
                             [0.000000000000, 0.485931890195, -0.485931890195, 6.375926430356, -2.443727560779, 0.053733020618, -0.053733020618, 2.443727560779, -6.375926430356, -0.000000000000 ],
                             [0.000000000000, 0.740805831639, -0.559823901756, -0.096284337505, 3.936996959017, -1.930891961845, 1.749910031962, 0.951256224701, 0.096284337505, -4.888253183717  ],
                             [0.000000000000, 0.559823901757, 0.492870367158, 2.605067077684, -1.950933405907, -0.750232850762, -0.302461418154, -2.469321742629, -2.605067077684, 4.420255148536 ],
                             [0.000000000000, -0.492870367158, -0.740805831641, -0.201023373941, 0.206177080649, 2.565606080131, -1.331929881332, -0.674175115836, 0.201023373941, 0.467998035187 ],
                             [0.000000000000, -0.492870367158, -0.559823901756, 2.605067077684, 2.469321742625, 0.302461418155, 0.750232850759, 1.950933405904, -2.605067077684, -4.420255148529  ],
                             [0.000000000000, 0.559823901757, -0.740805831641, -0.096284337505, -0.951256224702, -1.749910031967, 1.930891961850, -3.936996959018, 0.096284337505, 4.888253183721 ],
                             [0.000000000000, 0.740805831639, 0.492870367158, -0.201023373941, 0.674175115836, 1.331929881335, -2.565606080133, -0.206177080649, 0.201023373941, -0.467998035186  ]])



  _self.dgqPointsdl1 = np.array([[1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1],
                                 [1, 0, -1]])

  _self.dgqPointsdl2 = np.array([[0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1],
                                 [0, 1, -1]])

def numerical(_self,_e):

  _self.mass = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kxy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyx = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.kyy = np.zeros([_self.NUMGLEU,_self.NUMGLEU], dtype = float)
  _self.gx = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.gy = np.zeros([_self.NUMGLEU,_self.NUMGLEP], dtype = float)
  _self.dx = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  _self.dy = np.zeros([_self.NUMGLEP,_self.NUMGLEU], dtype = float)
  
  localx = np.zeros([_self.NUMRULE,1], dtype = float)
  localy = np.zeros([_self.NUMRULE,1], dtype = float)

  dxdl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl1 = np.zeros([_self.NUMRULE,1], dtype = float)
  dydl2 = np.zeros([_self.NUMRULE,1], dtype = float)

  dphiJdx = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
  dphiJdy = np.zeros([_self.NUMRULE,_self.NUMGLEU], dtype = float)
 

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  jacobian = np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
                                     [1, _self.x[j], _self.y[j]],
                                     [1, _self.x[k], _self.y[k]]]))


  for k in range(0,_self.NUMRULE):
   for i in range(0,_self.NUMGLEU):
    v = _self.IEN[_e][i]
    
    localx[k] += _self.x[v]*_self.phiJ[k][i]
    localy[k] += _self.y[v]*_self.phiJ[k][i]

    dxdl1[k] += _self.x[v]*_self.dphiJdl1[k][i]
    dxdl2[k] += _self.x[v]*_self.dphiJdl2[k][i]
    dydl1[k] += _self.y[v]*_self.dphiJdl1[k][i]
    dydl2[k] += _self.y[v]*_self.dphiJdl2[k][i]


  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    dphiJdx[k][i] = (1./jacobian)*(+_self.dphiJdl1[k][i]*dydl2[k]\
                                   -_self.dphiJdl2[k][i]*dydl1[k])

    dphiJdy[k][i] = (1./jacobian)*(-_self.dphiJdl1[k][i]*dxdl2[k]\
                                   +_self.dphiJdl2[k][i]*dxdl1[k])

  
  for k in range(0,_self.NUMRULE): 
   for i in range(0,_self.NUMGLEU):
    for j in range(0,_self.NUMGLEU):
     _self.mass[i][j] += 0.5*(_self.phiJ[k][i]*_self.phiJ[k][j]\
                            *jacobian*_self.gqWeights[k])
    
     _self.kxx[i][j] += 0.5*(dphiJdx[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])
     
     _self.kxy[i][j] += 0.5*(dphiJdx[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyx[i][j] += 0.5*(dphiJdy[k][i]*dphiJdx[k][j]\
                           *jacobian*_self.gqWeights[k])

     _self.kyy[i][j] += 0.5*(dphiJdy[k][i]*dphiJdy[k][j]\
                           *jacobian*_self.gqWeights[k])
    
    for j in range(0,_self.NUMGLEP):
     _self.gx[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.gy[i][j] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dx[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdx[k][i]\
                          *jacobian*_self.gqWeights[k])

     _self.dy[j][i] += 0.5*(-_self.gqPoints[k][j]*dphiJdy[k][i]\
                          *jacobian*_self.gqWeights[k])





