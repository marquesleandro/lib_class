# ==========================================
# Code created by Leandro Marques at 02/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code searches a file until 10 directories


# ------------------------------------------------------------------
# Use:
#name_file = 'file'
#path = search_file.Find(name_file)
# ------------------------------------------------------------------




import os

def Find(_file):
  path = '/home/marquesleandro'
  file = _file
  files1 = filter(os.path.isdir,os.listdir(os.curdir))
  breaking = 0

  if file in os.listdir(path):
   print path

  elif files1 == []:
   print "Error"

  else:
   for dir1 in files1:
    os.chdir(dir1)
    path = path + '/' + dir1
    files2 = filter(os.path.isdir,os.listdir(os.curdir))
    if file in os.listdir(path):
     breaking = 1
     break

    elif files2 == []:
     aa = path.split('/')
     directory = ''
   
     for i in range(1,len(aa)-1):
      bb = aa[i]
      directory = directory + '/' + bb
 
     path = directory
     os.chdir(path)
    
    else:
     for dir2 in files2:
      os.chdir(dir2)
      path = path + '/' + dir2
      files3 = filter(os.path.isdir,os.listdir(os.curdir))
      if file in os.listdir(path):
       breaking = 1
       break
 
      elif files3 == []:
       aa = path.split('/')
       directory = ''
    
       for i in range(1,len(aa)-1):
        bb = aa[i]
        directory = directory + '/' + bb
 
       path = directory
       os.chdir(path)
   
      else:
       for dir3 in files3:
        os.chdir(dir3)
        path = path + '/' + dir3
        files4 = filter(os.path.isdir,os.listdir(os.curdir))
        if file in os.listdir(path):
         breaking = 1
         break
 
        elif files4 == []:
         aa = path.split('/')
         directory = ''
    
         for i in range(1,len(aa)-1):
          bb = aa[i]
          directory = directory + '/' + bb
 
         path = directory
         os.chdir(path)
  
        else:
         for dir4 in files4:
          os.chdir(dir4)
          path = path + '/' + dir4
          files5 = filter(os.path.isdir,os.listdir(os.curdir))
          if file in os.listdir(path):
           breaking = 1
           break
 
          elif files5 == []:
           aa = path.split('/')
           directory = ''
    
           for i in range(1,len(aa)-1):
            bb = aa[i]
            directory = directory + '/' + bb
 
           path = directory
           os.chdir(path)
   
          else:
           for dir5 in files5:
            os.chdir(dir5)
            path = path + '/' + dir5
            files6 = filter(os.path.isdir,os.listdir(os.curdir))
            if file in os.listdir(path):
             breaking = 1
             break
 
            elif files6 == []:
             aa = path.split('/')
             directory = ''
    
             for i in range(1,len(aa)-1):
              bb = aa[i]
              directory = directory + '/' + bb
 
             path = directory
             os.chdir(path)
  
            else:
             for dir6 in files6:
              os.chdir(dir6)
              path = path + '/' + dir6
              files7 = filter(os.path.isdir,os.listdir(os.curdir))
              if file in os.listdir(path):
               breaking = 1
               break
 
              elif files7 == []:
               aa = path.split('/')
               directory = ''
    
               for i in range(1,len(aa)-1):
                bb = aa[i]
                directory = directory + '/' + bb
 
               path = directory
               os.chdir(path)
  
              else:
               for dir7 in files7:
                os.chdir(dir7)
                path = path + '/' + dir7
                files8 = filter(os.path.isdir,os.listdir(os.curdir))
                if file in os.listdir(path):
                 breaking = 1
                 break
 
                elif files8 == []:
                 aa = path.split('/')
                 directory = ''
    
                 for i in range(1,len(aa)-1):
                  bb = aa[i]
                  directory = directory + '/' + bb
 
                 path = directory
                 os.chdir(path)
  
                else:
                 for dir8 in files8:
                  os.chdir(dir8)
                  path = path + '/' + dir8
                  files9 = filter(os.path.isdir,os.listdir(os.curdir))
                  if file in os.listdir(path):
                   breaking = 1
                   break
 
                  elif files9 == []:
                   aa = path.split('/')
                   directory = ''
    
                   for i in range(1,len(aa)-1):
                    bb = aa[i]
                    directory = directory + '/' + bb
 
                   path = directory
                   os.chdir(path)
  

                  else:
                   for dir9 in files9:
                    os.chdir(dir9)
                    path = path + '/' + dir9
                    files10 = filter(os.path.isdir,os.listdir(os.curdir))
                    if file in os.listdir(path):
                     breaking = 1
                     break
  
                    elif files10 == []:
                     aa = path.split('/')
                     directory = ''
    
                     for i in range(1,len(aa)-1):
                      bb = aa[i]
                      directory = directory + '/' + bb
  
                     path = directory
                     os.chdir(path)
   
                    else:
                     for dir10 in files10:
                      os.chdir(dir10)
                      path = path + '/' + dir10
                      files11 = filter(os.path.isdir,os.listdir(os.curdir))
                      if file in os.listdir(path):
                       breaking = 1
                       break
  
                      elif files11 == []:
                       aa = path.split('/')
                       directory = ''
     
                       for i in range(1,len(aa)-1):
                        bb = aa[i]
                        directory = directory + '/' + bb
  
                       path = directory
                       os.chdir(path)
   
                      else:
                       for dir11 in files11:
                        os.chdir(dir11)
                        path = path + '/' + dir11
                        files12 = filter(os.path.isdir,os.listdir(os.curdir))
                        if file in os.listdir(path):
                         breaking = 1
                         break
 
                        elif files12 == []:
                         aa = path.split('/')
                         directory = ''
    
                         for i in range(1,len(aa)-1):
                          bb = aa[i]
                          directory = directory + '/' + bb
  
                         path = directory
                         os.chdir(path)
  

                        else:
                         aa = path.split('/')
                         directory = ''
     
                         for i in range(1,len(aa)-1):
                          bb = aa[i]
                          directory = directory + '/' + bb
  
                         path = directory
                         os.chdir(path)
  
  
                       if breaking == 1:
                        break
    
                       else:  
                        aa = path.split('/')
                        directory = ''
    
                        for i in range(1,len(aa)-1):
                         bb = aa[i]
                         directory = directory + '/' + bb
 
                        path = directory
                        os.chdir(path)
 
 
                     if breaking == 1:
                       break
    
                     else:  
                      aa = path.split('/')
                      directory = ''
    
                      for i in range(1,len(aa)-1):
                       bb = aa[i]
                       directory = directory + '/' + bb
 
                      path = directory
                      os.chdir(path)
 
                   if breaking == 1:
                    break
    
                   else:  
                    aa = path.split('/')
                    directory = ''
    
                    for i in range(1,len(aa)-1):
                     bb = aa[i]
                     directory = directory + '/' + bb
 
                    path = directory
                    os.chdir(path)
 
                 if breaking == 1:
                  break
    
                 else:  
                  aa = path.split('/')
                  directory = ''
    
                  for i in range(1,len(aa)-1):
                   bb = aa[i]
                   directory = directory + '/' + bb
 
                  path = directory
                  os.chdir(path)
 
               if breaking == 1:
                break
    
               else:  
                aa = path.split('/')
                directory = ''
    
                for i in range(1,len(aa)-1):
                 bb = aa[i]
                 directory = directory + '/' + bb
 
                path = directory
                os.chdir(path)
 
             if breaking == 1:
              break
   
             else:  
              aa = path.split('/')
              directory = ''
    
              for i in range(1,len(aa)-1):
               bb = aa[i]
               directory = directory + '/' + bb

              path = directory
              os.chdir(path)
 
           if breaking == 1:
            break
   
           else:  
            aa = path.split('/')
            directory = ''
   
            for i in range(1,len(aa)-1):
             bb = aa[i]
             directory = directory + '/' + bb

            path = directory
            os.chdir(path)

         if breaking == 1:
          break
   
         else:  
          aa = path.split('/')
          directory = ''
   
          for i in range(1,len(aa)-1):
           bb = aa[i]
           directory = directory + '/' + bb

          path = directory
          os.chdir(path)

       if breaking == 1:
        break
   
       else:  
        aa = path.split('/')
        directory = ''
   
        for i in range(1,len(aa)-1):
         bb = aa[i]
         directory = directory + '/' + bb

        path = directory
        os.chdir(path)

     if breaking == 1:
      break
   
     else:  
      aa = path.split('/')
      directory = ''
   
      for i in range(1,len(aa)-1):
       bb = aa[i]
       directory = directory + '/' + bb

      path = directory
      os.chdir(path)




  return path 
