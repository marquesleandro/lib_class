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
import trigauss
import scipy.sparse as sps
from tqdm import tqdm

def assembly_linear(_npoints, _nelem, _IEN, _x, _y):
 
 Kxx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kxy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 MLump = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gy = sps.lil_matrix((_npoints,_npoints), dtype = float)


 linear = trigauss.Linear(_x, _y, _IEN)

 for e in tqdm(range(0, _nelem)):
  linear.numerical(e)

  for i in range(0,3): 
   ii = _IEN[e][i]
  
   for j in range(0,3):
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

