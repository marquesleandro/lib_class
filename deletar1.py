class Quad:
 def __init__(_self,_x,_y,_IEN,_npoints,_nelem,_scalar1,_scalar2,_scalar3,_vec1,_vec2):
  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN
  _self.npoints = _npoints
  _self.nelem = _nelem
  _self.scalar1 = _scalar1
  _self.scalar2 = _scalar2
  _self.scalar3 = _scalar3
  _self.vec1 = _vec1
  _self.vec2 = _vec2


 def create_dir(_self,_dir):

  _self.path = '/home/marquesleandro'
  if not 'results' in os.listdir(_self.path):
   os.mkdir('/home/marquesleandro/results')

  _self.path = '/home/marquesleandro/results'
  if not _dir in os.listdir(_self.path):
   #ndir = len(os.listdir(_self.path))
   #_self.dir = _dir + str(ndir)
   _self.dir = _dir
   create_dir = os.mkdir(_self.path + '/' + _self.dir)

  else:
   _self.dir = _dir


 def saveVTK(_self,_file):

  vtkFile = open(_self.path + '/' + _self.dir + '/' + _file + '.vtk', 'w')
  _self.vtkHeader(vtkFile)

  _self.vtkCoords(vtkFile)
  _self.vtkCellArray(vtkFile)
  _self.vtkCellType(vtkFile)
  _self.vtkScalarScalarHeader(vtkFile)
  _self.vtkVector(vtkFile,'vector', _self.vec1, _self.vec2)

  if _self.scalar1 is not None:
   _self.vtkScalar(vtkFile,"scalar1",_self.scalar1);
  
  if _self.scalar2 is not None:
   _self.vtkScalar(vtkFile,"scalar2",_self.scalar2);

  if _self.scalar3 is not None:
   _self.vtkScalar(vtkFile,"scalar3",_self.scalar3);


  vtkFile.close()


 def vtkHeader(_self,_file,_iter=None):
  _file.write( "# vtk DataFile Version 2.0\n" )
  _file.write( "2D Simulation Python\n" )
  _file.write( "ASCII\n" )
  _file.write( "DATASET UNSTRUCTURED_GRID\n" )

  _file.write( "FIELD FieldData 1\n" )
  _file.write( "NODES 1 2 int\n" )
  _file.write( str(_self.npoints) + " " + \
               str(_self.nelem) + "\n" )

  _file.write( "\n" )

 def vtkCoords(_self,_file):
  _file.write( "POINTS " + str(_self.npoints) + " double\n" )
  for i in range(0,_self.npoints):
   _file.write( str(_self.x[i][0]) + " " + \
                str(_self.y[i][0]) + " 0.0\n" )

  _file.write( "\n" )

 def vtkCellArray(_self,_file):
  _file.write( "CELLS " + str(_self.nelem) \
                        + " " + str(4*_self.nelem) + "\n" )
  for i in range(0,_self.nelem):
   _file.write( "6 " + str(_self.IEN[i][0]) + " " + \
                       str(_self.IEN[i][1]) + " " + \
                       str(_self.IEN[i][2]) + " " + \
                       str(_self.IEN[i][3]) + " " + \
                       str(_self.IEN[i][4]) + " " + \
                       str(_self.IEN[i][5]) + "\n" )

  _file.write( "\n" )

 def vtkCellType(_self,_file):
  _file.write( "CELL_TYPES " + str(_self.nelem) + "\n" )
  for i in range(0,_self.nelem):
   _file.write( "22 ")
   _file.write( "\n" )
 
  _file.write( "\n" )
   
 def vtkScalarHeader(_self,_file):
  _file.write( "POINT_DATA " + str(_self.npoints) + "\n" )
  _file.write( "\n" )

 def vtkScalarScalarHeader(_self,_file):
  _file.write( "POINT_DATA " + str(_self.npoints) + "\n" )

 def vtkScalar(_self,_file,_name,_scalar):
  _file.write( "SCALARS " + _name + " double\n" )
  _file.write( "LOOKUP_TABLE default\n" )

  for i in range(0,_self.npoints):
   _file.write( str(_scalar.item(i)) + "\n" )

  _file.write( "\n" )

 def vtkVector(_self,_file,_name,_vec1,_vec2):
  _file.write( "VECTORS " + _name + " double\n" )

  for i in range(0,_self.npoints):
   _file.write( str(_vec1.item(i)) + " " + \
                str(_vec2.item(i)) + " 0.0\n" )

  _file.write( "\n" )


# def printMeshReport(_self):
#  """
#   Print mesh report for lineMesh and Mesh
#  """
#  print ""
#  print ""
#  print "|" + "-"*30 + " Mesh Report " + "-"*30 + "|"
#  print " "*5 + "number of 2D points (npoints):          " + \
#        str(_self.npoints)
#  print " "*5 + "number of triangles (nelem):          " + \
#        str(_self.nelem)
#--------------------------------------------------
#   print ""
#   for nb in range(0,_self.mesh.elemIdRegion.max()+1):
#    print " "*5 + "line (" + str(nb) + ")" 
#    print " "*5 + "  |length (averageLength):              " + \
#         str(_self.mesh.averageEdgeLength)
#-------------------------------------------------- 
#
#  print "|" + "-"*73 + "|"
#  print ""
#  print ""


