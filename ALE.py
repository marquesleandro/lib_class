# ==========================================
# Code created by Leandro Marques at 02/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to apply ALE scheme

import sys
import numpy as np


def rotate(_npoints, _x, _y):
 vx_mesh = np.zeros([_npoints,1], dtype = float)
 vy_mesh = np.zeros([_npoints,1], dtype = float)

 for i in range(0, _npoints): 
  vx_mesh[i] = _y[i]*0.1
  vy_mesh[i] = - _x[i]*0.1

 return vx_mesh, vy_mesh
