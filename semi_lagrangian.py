# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used for to use semi-lagrangian scheme

# ------------------------------------------------------------------
# Use:
# scalar_d = semi_lagrangian(mesh.npoints, mesh.IEN, mesh.x, mesh.y, x_d, y_d, mesh.neighbors_elements, bc_dirichlet, bc_neumann, scalar_n)
# ------------------------------------------------------------------

import numpy as np

def Linear_2D(_npoints, _IEN, _xn, _yn, _xd, _yd, _neighbors_elements, _bc_dirichlet, _bc_neumann, _scalar):
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for p in range(0,_npoints):
  x = float(_xd[p])
  y = float(_yd[p])

  node = p
  length = []
  ww = 1
  #print ""
  #print p

  while ww == 1:
   for e in _neighbors_elements[p]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    y1 = float(_yn[v1])
    y2 = float(_yn[v2])
    y3 = float(_yn[v3])
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])
 
    alpha = np.linalg.solve(A,b)
 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     #ee = e + 81
     #print "elemento dominio %s" %ee
     #print "fazer interpolacao triangular" 
     
     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     Li = A1/At
     Lj = A2/At
     Lk = A3/At
     
     Ni = Li
     Nj = Lj
     Nk = Lk

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]

     scalar[node] = Ni*scalar1 + Nj*scalar2 + Nk*scalar3
     #print Ni + Nj + Nk
     #print scalar[node]

     ww = 0
     break

    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
  
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     ww = 1

   # first neighbor is element found 
   if ww == 0:
     break
  
   # coordinate doesn't found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = p
    p = length_min[0]
    #print p

    # outside domain
    if p == node1 and ww == 1:
     #p = p + 1
     #print "elemento contorno proximo ao no %s" %p
     #print "fazer regra da alavanca"

     scalar[node] = _scalar[node]
     
     #print Ni + Nj
     #print scalar[node]

     ww = 0
     break


 return scalar

