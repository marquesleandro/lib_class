# ==========================================
# Code created by Leandro Marques at 01/2020
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to apply boundaries conditions in simulator

import sys
import boundary_condition


def Element1D(_nphysical, _npoints, _x, _neumann_pts, _dirichlet_pts, _neighbors_nodes, _LHS, simulator_option):

 # Convection1D
 if simulator_option == 1:
  # --------- Boundaries conditions --------------------
  condition = boundary_condition.Convection1D(_nphysical, _npoints, _x)
  condition.neumann_condition(_neumann_pts)
  condition.dirichlet_condition(_dirichlet_pts)
  condition.gaussian_elimination(_LHS,_neighbors_nodes)
 
  # --------- Initial condition ------------------------
  condition.initial_condition()
  # ----------------------------------------------------

  return condition.bc_dirichlet, condition.bc_neumann, condition.bc_2, condition.LHS, condition.c, condition.vx


 else:
  print ""
  print " Error: Boundary Condition not found"
  print ""
  sys.exit()





def Element2D(_nphysical, _npoints, _x, _y, _neumann_edges, _dirichlet_pts, _neighbors_nodes, _LHS, simulator_option):

 # Convection2D
 if simulator_option == 4:
  # --------- Boundaries conditions --------------------
  condition = boundary_condition.Convection2D(_nphysical, _npoints, _x, _y)
  condition.neumann_condition(_neumann_edges)
  condition.dirichlet_condition(_dirichlet_pts)
  condition.gaussian_elimination(_LHS,_neighbors_nodes)
 
  # --------- Initial condition ------------------------
  condition.initial_condition()
  # ----------------------------------------------------

  return condition.bc_dirichlet, condition.bc_neumann, condition.bc_2, condition.LHS, condition.c, condition.vx, condition.vy


 else:
  print ""
  print " Error: Boundary Condition not found"
  print ""
  sys.exit()



