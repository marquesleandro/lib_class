# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used for to assembly global matrices

# ------------------------------------------------------------------
# Use:
# Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy = 
# assembly_linear(mesh.npoints, mesh.nelem, mesh.IEN, mesh.x, mesh.y)
# ------------------------------------------------------------------

import numpy as np
import gaussian_quadrature
import scipy.sparse as sps
from tqdm import tqdm



def Linear1D(_GL, _npoints, _nelem, _IEN, _x):
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 G = sps.lil_matrix((_npoints,_npoints), dtype = float)
 
 
 linear = gaussian_quadrature.Linear1D()

 for e in tqdm(range(0, _nelem)):
  dx = _x[e+1] - _x[e]

  for i in range(0,_GL): 
   ii = _IEN[e][i]
  
   for j in range(0,_GL):
    jj = _IEN[e][j]

    K[ii,jj] += (1.0/dx)*linear.K_elem[i][j]
    M[ii,jj] += (dx/6.0)*linear.M_elem[i][j]
    G[ii,jj] += (1.0/dx)*linear.G_elem[i][j]


 return K, M, G



def Linear2D(_GL, _npoints, _nelem, _IEN, _x, _y):
 
 Kxx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kxy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 MLump = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gy = sps.lil_matrix((_npoints,_npoints), dtype = float)


 linear = gaussian_quadrature.Linear2D(_x, _y, _IEN)

 for e in tqdm(range(0, _nelem)):
  linear.numerical(e)

  for i in range(0,_GL): 
   ii = _IEN[e][i]
  
   for j in range(0,_GL):
    jj = _IEN[e][j]

    Kxx[ii,jj] += linear.kxx[i][j]
    Kxy[ii,jj] += linear.kxy[i][j]
    Kyx[ii,jj] += linear.kyx[i][j]
    Kyy[ii,jj] += linear.kyy[i][j]
    K[ii,jj] += linear.kxx[i][j] + linear.kyy[i][j]
   
    M[ii,jj] += linear.mass[i][j]
    MLump[ii,ii] += linear.mass[i][j]

    Gx[ii,jj] += linear.gx[i][j]
    Gy[ii,jj] += linear.gy[i][j]


 return Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy



def Mini_NS2D(_GLV, _GLP, _NV, _NP, _nelem, _IEN, _x, _y):
 
 K = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 M = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 MLump = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 G = sps.lil_matrix((2*_NV,_NP), dtype = float)
 D = sps.lil_matrix((_NP,2*_NV), dtype = float)


 mini = gaussian_quadrature.Mini(_x, _y, _IEN)

 for e in tqdm(range(0, _nelem)):
  mini.numerical(e)

  for i in range(0, _GLV): 
   ii = _IEN[e][i]
  
   for j in range(0, _GLV):
    jj = _IEN[e][j]

    #MSC 2007 pag.84
    K[ii,jj] += 2.0*mini.kxx[i][j] + mini.kyy[i][j] #K11
    K[ii,jj + _NV] += mini.kxy[i][j] #K12
    K[ii + _NV,jj] += mini.kyx[i][j] #K21
    K[ii + _NV,jj + _NV] += mini.kxx[i][j] + 2.0*mini.kyy[i][j] #K22
   
    M[ii,jj] += mini.mass[i][j]
    M[ii + _NV,jj + _NV] += mini.mass[i][j]
    
    MLump[ii,ii] += mini.mass[i][j]
    MLump[ii + _NV,ii + _NV] += mini.mass[i][j]


   for k in range(0, _GLP):
    kk = _IEN[e][k]

    G[ii,kk] += mini.gx[i][k]
    G[ii + _NV,kk] += mini.gy[i][k]

    D[kk,ii] += mini.dx[k][i]
    D[kk,ii + _NV] += mini.dy[k][i]


 return K, M, MLump, G, D



