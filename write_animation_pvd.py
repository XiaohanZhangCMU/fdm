#!/usr/bin/python

import sys, getopt
import os

def main(argv):
  start_step = 1
  max_step = 1000
  delta_step = 1

  try:
    opts, args = getopt.getopt(argv,"hn:i:s:d:",["total=","step=","start=","directory"])
  except getopt.GetoptError:
    print 'wrie_animation_pvd.py -s <start step> -n <end step> -i <every i steps> -d <dir>'
    sys.exit()

    
  for opt, arg in opts:
    if opt == '-h':
      print 'wrie_animation_pvd.py -s <start step> -n <end step> -i <every i steps>'
      sys.exit()
    elif opt in ("-n", "--total"):
      max_step = int(arg)
    elif opt in ("-i", "--every"):
      delta_step = int(arg)
    elif opt in ("-s", "--start"):
      start_step = int(arg)
    elif opt in ("-d", "--directory"):
      directory = arg

#  filename = "./outputs/ag_dislocation.pvd"
  filename = directory + "/ag_dislocation.pvd"

  if os.path.isfile(filename):
    os.remove(filename)

  f = open( filename, 'w' )

  step = start_step

  f.write( "<?xml version=\"1.0\"?>" + "\n")
  f.write( "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"" + "\n")
  f.write( "compressor=\"vtkZLibDataCompressor\">" + "\n")
  f.write( "<Collection>" + "\n")

  while ( step <= max_step ):
    f.write("  <DataSet timestep=\"" + str(step) + " \" file=\"alpha_solution-" +  str(step).zfill(5) + ".00000.pvtu\"/>" + "\n")
    step = step + delta_step

  f.write( " </Collection>" + "\n")
  f.write( " </VTKFile> " + "\n")

  f.close()


###################################################

#  filename = "./outputs/ag_stress.pvd"
  filename = directory+"/ag_stress.pvd"

  if os.path.isfile(filename):
    os.remove(filename)

  f = open( filename, 'w' )

  step = start_step

  f.write( "<?xml version=\"1.0\"?>" + "\n")
  f.write( "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"" + "\n")
  f.write( "compressor=\"vtkZLibDataCompressor\">" + "\n")
  f.write( "<Collection>" + "\n")

  while ( step <= max_step ):
    f.write("  <DataSet timestep=\"" + str(step) + " \" file=\"fsvr_solution-" +  str(step).zfill(5) + ".00000.pvtu\"/>" + "\n")
    step = step + delta_step

  f.write( " </Collection>" + "\n")
  f.write( " </VTKFile> " + "\n")

  f.close()


##################################################


#  filename = "./outputs/ag_xchi.pvd"
  filename = directory+"/ag_xchi.pvd"

  if os.path.isfile(filename):
    os.remove(filename)

  f = open( filename, 'w' )

  step = start_step

  f.write( "<?xml version=\"1.0\"?>" + "\n")
  f.write( "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"" + "\n")
  f.write( "compressor=\"vtkZLibDataCompressor\">" + "\n")
  f.write( "<Collection>" + "\n")

  while ( step <= max_step ):
    f.write("  <DataSet timestep=\"" + str(step) + " \" file=\"xchi_solution-" +  str(step).zfill(5) + ".00000.pvtu\"/>" + "\n")
    step = step + delta_step

  f.write( " </Collection>" + "\n")
  f.write( " </VTKFile> " + "\n")

  f.close()



##################################################


#  filename = "./outputs/ag_xchi.pvd"
  filename = directory+"/ag_velocity.pvd"

  if os.path.isfile(filename):
    os.remove(filename)

  f = open( filename, 'w' )

  step = start_step

  f.write( "<?xml version=\"1.0\"?>" + "\n")
  f.write( "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"" + "\n")
  f.write( "compressor=\"vtkZLibDataCompressor\">" + "\n")
  f.write( "<Collection>" + "\n")

  while ( step <= max_step ):
    f.write("  <DataSet timestep=\"" + str(step) + " \" file=\"vsvr_solution-" +  str(step).zfill(5) + ".00000.pvtu\"/>" + "\n")
    step = step + delta_step

  f.write( " </Collection>" + "\n")
  f.write( " </VTKFile> " + "\n")

  f.close()




  print "\n... file written ...\n\n"


if __name__ == "__main__":
   main(sys.argv[1:])
