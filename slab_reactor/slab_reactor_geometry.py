# pyexamples/slab_reactor/slab_reactor_geometry.py
#
# This implements the materials for the test cases published 
# in Scott Mosher's Ph.D. thesis, "A Variational Coarse Mesh 
# Transport Method". All data and reference values are from 
# Appendix A of that work.

from detran import *

def get_mesh(geometry_id) :

  # Define the assembly coarse mesh edges.
  cm_assembly = [1.1580, 4.4790, 7.8000, 11.1210, 14.4420, 15.6000]

  # Define the base fine mesh count per coarse mesh.
  #fm_assembly = [1, 2, 2, 2, 2, 1]
  fm_assembly = [2, 4, 4, 4, 4, 2]
   
  # A core is composed of 7 adjacent assemblies.
  cm_core = [0.0]
  fm_core = []
  for i in range(0, 7) :
    for j in range(0, len(cm_assembly)) :
      cm_core.append(cm_assembly[j] + 15.6*float(i))
      fm_core.append(fm_assembly[j])

  # Add the 0.0 edge to assembly coarse mesh
  cm_assembly.insert(0, 0.0)

  # Define the coarse mesh material maps for each assembly type.  
  # There are four of these assemblies.  These correspond, in order,
  # to types A, B, C, and D in the thesis.
  assem = [[ 0, 1, 2, 2, 1, 0 ], \
           [ 0, 1, 1, 1, 1, 0 ], \
           [ 0, 1, 3, 3, 1, 0 ], \
           [ 0, 3, 3, 3, 3, 0 ]]

  # Cores 0, 1 and 2
  core_0 = assem[0]+assem[1]+assem[0]+assem[1]+assem[0]+assem[1]+assem[0]
  core_1 = assem[0]+assem[2]+assem[0]+assem[2]+assem[0]+assem[2]+assem[0]
  core_2 = assem[0]+assem[3]+assem[0]+assem[3]+assem[0]+assem[3]+assem[0]

  # Create the 1D mesh
  if   geometry_id == "assembly0" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[0])
  elif geometry_id == "assembly1" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[1])
  elif geometry_id == "assembly2" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[2])
  elif geometry_id == "assembly3" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[3])
  elif geometry_id == "core0" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_0)
  elif geometry_id == "core1" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_1)
  elif geometry_id == "core2" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_2)
  else :
    print "invalid geometry selected."
    exit()
  
  return mesh
