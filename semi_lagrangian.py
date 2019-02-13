# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used for to use semi-lagrangian scheme

# ------------------------------------------------------------------------------------
# Use:
# scalar_d = semi_lagrangian.Linear2D(
# mesh.npoints, mesh.IEN, mesh.x, mesh.y, x_d, y_d, mesh.neighbors_elements, scalar_n)
# ------------------------------------------------------------------------------------

import numpy as np


# 1D Semi-Lagrangian using npoints x nelem to find departure node
def Linear1D_v2(_npoints, _nelem, _IEN, _xn, _xd, _scalar):
 
 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(_xd[i])

  breaking = 0
  length = []

  for e in range(0,_nelem):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]

   x1 = float(_xn[v1])
   x2 = float(_xn[v2])

   len1 = abs(x2 - x)
   len2 = abs(x1 - x)
   lent = abs(x1 - x2)

   Li = len1/lent
   Lj = len2/lent

   alpha = [Li,Lj]
   alpha = np.array(alpha)

   if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
    Ni = Li
    Nj = Lj

    scalar1 = _scalar[v1]
    scalar2 = _scalar[v2]

    scalar[i] = Ni*scalar1 + Nj*scalar2
    breaking = 1
    break

   else:
    x_a = x1 - x
    x_b = x2 - x
  
    length1 = np.sqrt(x_a**2)
    length2 = np.sqrt(x_b**2)

    a_1 = [v1,length1]
    a_2 = [v2,length2]
 
    length.append(a_1)
    length.append(a_2)
   
  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
  
 return scalar  


# 1D Semi-Lagrangian using npoints x neighbors_elements to find departure node
def Linear1D(_npoints, _IEN, _xn, _xd, _neighbors_elements, _scalar):
 
 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(_xd[i])

  node = i
  length = []
  breaking = 0
  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])

    len1 = abs(x2 - x)
    len2 = abs(x1 - x)
    lent = abs(x1 - x2)

    Li = len1/lent
    Lj = len2/lent

    alpha = [Li,Lj]
    alpha = np.array(alpha)
 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     Ni = Li
     Nj = Lj
 
     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]

     scalar[i] = Ni*scalar1 + Nj*scalar2
     breaking = 1
     break

    else:
     x_a = x1 - x
     x_b = x2 - x
  
     length1 = np.sqrt(x_a**2)
     length2 = np.sqrt(x_b**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
 
     length.append(a_1)
     length.append(a_2)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate doesn't found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar  



# 2D Semi-Lagrangian using npoints x nelem to find departure node
def Linear2D_v2(_npoints, _nelem, _IEN, _xn, _yn, _xd, _yd, _scalar):
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(_xd[i])
  y = float(_yd[i])
  
  breaking = 0
  length = []

  for e in range(0,_nelem):
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

    scalar[i] = Ni*scalar1 + Nj*scalar2 + Nk*scalar3
    breaking = 1
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
   
  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
     
 return scalar  



# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
def Linear2D(_npoints, _IEN, _xn, _yn, _xd, _yd, _neighbors_elements, _scalar):
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(_xd[i])
  y = float(_yd[i])

  node = i
  length = []
  breaking = 0

  while breaking == 0:
   for e in _neighbors_elements[node]:
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

     scalar[i] = Ni*scalar1 + Nj*scalar2 + Nk*scalar3

     breaking = 1
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
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate doesn't found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar

