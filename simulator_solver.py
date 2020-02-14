# =========================================
# Code created by Leandro Marques at 01/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to select solver type


import sys
import numpy as np
import scipy.sparse as sps
from tqdm import tqdm
from time import time
import solver
import export_vtk




def Element1D(_simulator_problem, _scheme_option, _polynomial_option, _x, _IEN, _npoints, _nelem, _scalar, _vx, _dt, _nt, _Re, _Sc, _M, _K, _G, _LHS, _bc_dirichlet, _bc_neumann, _bc_2, _neighbors_nodes, _neighbors_elements, _directory_name):

 simulator_problem = _simulator_problem
 scheme_option = _scheme_option 
 polynomial_option = _polynomial_option
 x = _x
 IEN = _IEN
 npoints = _npoints
 nelem = _nelem
 c = _scalar
 vx = _vx
 dt = _dt
 nt = _nt
 Re = _Re
 Sc = _Sc
 M = _M
 K = _K
 G = _G
 LHS = _LHS
 bc_dirichlet = _bc_dirichlet
 bc_neumann = _bc_neumann
 bc_2 = _bc_2
 neighbors_nodes = _neighbors_nodes
 neighbors_elements = _neighbors_elements
 directory_name = _directory_name


 if simulator_problem == 1:
  # Taylor Galerkin
  if scheme_option == 1:
 
   for t in tqdm(range(0, nt)):
 
    # ------------------------ Export VTK File --------------------------------------
    save = export_vtk.Linear1D(x,IEN,npoints,nelem,c,c,c,vx,vx)
    save.create_dir(directory_name)
    save.saveVTK(directory_name + str(t))
    # -------------------------------------------------------------------------------
 
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation1D(scheme_option)
    scheme.taylor_galerkin(c, vx, dt, M, K, G, LHS, bc_dirichlet, bc_2)
    c = scheme.c
   return c, scheme.scheme_name
    # -------------------------------------------------------------------------------


  # Semi Lagrangian Linear
  elif scheme_option == 2:

   if polynomial_option == 1: #Linear Element
    for t in tqdm(range(0, nt)):
 
     # ------------------------ Export VTK File --------------------------------------
     save = export_vtk.Linear1D(x,IEN,npoints,nelem,c,c,c,vx,vx)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # -------------------------------------------------------------------------------

     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation1D(scheme_option)
     scheme.semi_lagrangian_linear(npoints, neighbors_elements, IEN, x, vx, dt, c, M, LHS, bc_dirichlet, bc_2)
     c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------

   elif polynomial_option == 2: #Quad Element
    for t in tqdm(range(0, nt)):
 
     # ------------------------ Export VTK File --------------------------------------
     save = export_vtk.Linear1D(x,IEN,npoints,nelem,c,c,c,vx,vx)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # -------------------------------------------------------------------------------
 
     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation1D(scheme_option)
     scheme.semi_lagrangian_quad(npoints, nelem, neighbors_elements, IEN, x, vx, dt, c, M, LHS, bc_dirichlet, bc_2)
     c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------




def Element2D(_simulator_problem, _scheme_option, _polynomial_option, _x, _y, _IEN, _npoints, _nelem, _scalar, _vx, _vy, _dt, _nt, _Re, _Sc, _M, _Kxx, _Kyx, _Kxy, _Kyy, _Gx, _Gy, _LHS, _bc_dirichlet, _bc_neumann, _bc_2, _neighbors_nodes, _neighbors_elements, _directory_name):

 simulator_problem = _simulator_problem
 scheme_option = _scheme_option 
 polynomial_option = _polynomial_option
 x = _x
 y = _y
 IEN = _IEN
 npoints = _npoints
 nelem = _nelem
 c = _scalar
 vx = _vx
 vy = _vy
 dt = _dt
 nt = _nt
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
 neighbors_nodes = _neighbors_nodes
 neighbors_elements = _neighbors_elements
 directory_name = _directory_name


 # Poiseuille Problem
 if simulator_problem == 1:

  # Taylor Galerkin Scheme
  if scheme_option == 1:
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
    scheme.taylor_galerkin(c, vx, vy, dt, Re, Sc, M, Kxx, Kyx, Kxy, Kyy, Gx, Gy, LHS, bc_dirichlet, bc_neumann, bc_2)
    c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------


  # Semi Lagrangian Scheme
  elif scheme_option == 2:
 
   if polynomial_option == 1: #Linear Element
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
    scheme.semi_lagrangian_linear(npoints, neighbors_nodes, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
    c = scheme.c
    
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------

   elif polynomial_option == 2: #Mini Element
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
    scheme.semi_lagrangian_mini(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
    c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------


   elif polynomial_option == 3: #Quadratic Element
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
    scheme.semi_lagrangian_quad(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
    c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------

   elif polynomial_option == 4: #Cubic Element
    # -------------------------------- Solver ---------------------------------------
    scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
    scheme.semi_lagrangian_cubic(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
    c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------






 # Convection Problem
 if simulator_problem == 4:
  # Taylor Galerkin Scheme
  if scheme_option == 1:

    for t in tqdm(range(0, nt)):

     # ------------------------ Export VTK File ---------------------------------------
     save = export_vtk.Linear2D(x,y,IEN,npoints,nelem,c,c,c,vx,vy)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # --------------------------------------------------------------------------------

     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
     scheme.taylor_galerkin(c, vx, vy, dt, Re, Sc, M, Kxx, Kyx, Kxy, Kyy, Gx, Gy, LHS, bc_dirichlet, bc_neumann, bc_2)
     c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------


  # Semi Lagrangian Scheme
  elif scheme_option == 2:
 
   if polynomial_option == 1: #Linear Element
    for t in tqdm(range(0, nt)):
     # ------------------------ Export VTK File ---------------------------------------
     save = export_vtk.Linear2D(x,y,IEN,npoints,nelem,c,c,c,vx,vy)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # --------------------------------------------------------------------------------

     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
     scheme.semi_lagrangian_linear(npoints, neighbors_nodes, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
     c = scheme.c
    
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------

   elif polynomial_option == 2: #Mini Element
    for t in tqdm(range(0, nt)):
 
     # ------------------------ Export VTK File ---------------------------------------
     save = export_vtk.Linear2D(x,y,IEN,npoints,nelem,c,c,c,vx,vy)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # --------------------------------------------------------------------------------

     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
     scheme.semi_lagrangian_mini(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
     c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------


   elif polynomial_option == 3: #Quadratic Element
    #for t in tqdm(range(0, nt)):
 
     # ------------------------ Export VTK File ---------------------------------------
     #save = export_vtk.Linear2D(x,y,IEN,npoints,nelem,c,c,c,vx,vy)
     #save.create_dir(directory_name)
     #save.saveVTK(directory_name + str(t))
     # --------------------------------------------------------------------------------
 
     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
     scheme.semi_lagrangian_quad(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
     c = scheme.c
     return c, scheme.scheme_name
    # -------------------------------------------------------------------------------

   elif polynomial_option == 4: #Cubic Element
    for t in tqdm(range(0, nt)):
 
     # ------------------------ Export VTK File ---------------------------------------
     save = export_vtk.Linear2D(x,y,IEN,npoints,nelem,c,c,c,vx,vy)
     save.create_dir(directory_name)
     save.saveVTK(directory_name + str(t))
     # --------------------------------------------------------------------------------

     # -------------------------------- Solver ---------------------------------------
     scheme = solver.SemiImplicit_concentration_equation2D(scheme_option)
     scheme.semi_lagrangian_cubic(npoints, nelem, neighbors_elements, IEN, x, y, vx, vy, dt, Re, Sc, c, M, LHS, bc_dirichlet, bc_neumann, bc_2)
     c = scheme.c
    return c, scheme.scheme_name
    # -------------------------------------------------------------------------------



