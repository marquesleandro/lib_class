# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used for import .msh
# for n governament equations 


# ------------------------------------------------------------------
# Use:
# mesh = import_msh.Linear2D(directory,mesh_name,equation_number)
# mesh.coord()
# mesh.ien()
# ------------------------------------------------------------------


import numpy as np


class Linear1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 15:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,2], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
  
   _self.IEN[e] = [v1,v2]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.length_min = min(length)


  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))



class Quad1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 15:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,3], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))

   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.length_min = min(length)


  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))




class Linear2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.dirichlet_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 1:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = int(_self.gmsh[(jj + 10 + ii)][6])
    a_4 = [a_1,a_2,a_3]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_4)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_4)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.y = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.y[i] = _self.gmsh[_self.nphysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,3], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = []

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v3]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.length_min = min(length)

  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))


class Mini:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.dirichlet_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 1:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = int(_self.gmsh[(jj + 10 + ii)][6])
    a_4 = [a_1,a_2,a_3]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_4)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_4)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1
  _self.NP = _self.npoints
  _self.NV = _self.npoints + _self.nelem

  for i in range(0, _self.NP):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj

    
 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,4], dtype = int)
  _self.GLV = len(_self.IEN[0,:])
  _self.GLP = len(_self.IEN[0,:]) - 1
  
  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
   v4 = _self.NP + e
  
   _self.IEN[e] = [v1,v2,v3,v4]
 
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   _self.neighbors_nodes[v4].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v3]))
   _self.neighbors_nodes[v4] = list(set(_self.neighbors_nodes[v4]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  
   _self.neighbors_elements[v4].append(e)  
 
  for i in range(0, _self.NP):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))
 
 def coord(_self):
  _self.x = np.zeros([_self.NV,1], dtype = float)
  _self.y = np.zeros([_self.NV,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.NP):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.y[i] = _self.gmsh[_self.nphysical + 8 + i][2]
   _self.npts.append(i)

  for e in range(0, _self.nelem):  
   v1 = _self.IEN[e][0]
   v2 = _self.IEN[e][1]
   v3 = _self.IEN[e][2]
   v4 = _self.IEN[e][3]

   _self.x[v4] = (_self.x[v1] + _self.x[v2] + _self.x[v3])/3.0
   _self.y[v4] = (_self.y[v1] + _self.y[v2] + _self.y[v3])/3.0
   _self.npts.append(v4)

class Quad:
 def __init__(_self, _dir, _file):
  _self.name = _dir + '/' + _file
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.dirichlet_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}


 def number_equations(_self, _neq):
  _self.neq = _neq
  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 8:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = int(_self.gmsh[(jj + 10 + ii)][6])
    a_4 = int(_self.gmsh[(jj + 10 + ii)][7])
    a_5 = [a_1,a_2,a_3,a_4]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_5)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_5)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,6], dtype = int)
  
  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
   v4 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][8]) - 1
   v5 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][9]) - 1
   v6 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][10]) - 1
  
   _self.IEN[e] = [v1,v2,v3,v4,v5,v6]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   _self.neighbors_nodes[v4].extend(_self.IEN[e])  
   _self.neighbors_nodes[v5].extend(_self.IEN[e])  
   _self.neighbors_nodes[v6].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v3]))
   _self.neighbors_nodes[v4] = list(set(_self.neighbors_nodes[v4]))
   _self.neighbors_nodes[v5] = list(set(_self.neighbors_nodes[v5]))
   _self.neighbors_nodes[v6] = list(set(_self.neighbors_nodes[v6]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  
   _self.neighbors_elements[v4].append(e)  
   _self.neighbors_elements[v5].append(e)  
   _self.neighbors_elements[v6].append(e)  
 
  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))

 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.y = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.y[i] = _self.gmsh[_self.nphysical + 8 + i][2]
   _self.npts.append(i)


