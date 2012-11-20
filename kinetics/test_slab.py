from detran import *
from twogroup_materials import *
import sys

def get_input() :
  inp = InputDB.Create()
  inp.put_int("number_groups",                2)
  inp.put_int("dimension",                    1)
  inp.put_str("equation",                     "sc")
  inp.put_str("bc_west",                      "reflect")
  inp.put_str("bc_east",                      "reflect")
  inp.put_int("quad_number_polar_octant",     1)
  return inp

def get_mesh(fmm) :
  # define the coarse meshes and coarse mesh material map.  
  cm  = [0, 30, 200, 220, 390, 420]
  mt  = [  2,  0,   4,   0,   2    ]
  # define the base fine mesh count per coarse mesh
  fm = [ 3,  17,  2,  17,  3  ]
  for i in range(0, len(fm)):
    fm[i] = fm[i] * fmm
  # define the meshes
  mesh = Mesh1D.Create(fm, cm, mt)
  return mesh

def get_steady_state(inp, mat, mesh) :
  solver = Eigen1D(inp, mat, mesh)
  solver.solve()
  state = solver.state()
  plot_mesh_function(mesh, state.phi(1))
  kcrit = state.eigenvalue()
  print " eigenvalue is ", kcrit
  return state

if __name__ == "__main__":
    
  Manager.initialize(sys.argv)

  # get the input, mesh, and materials
  inp  = get_input()
  mesh = get_mesh(10) # 1 cm reference

  # get materials
  mat1 = get_rod_in()
  mat2 = get_rod_out()
  mats = vec_material([mat1, mat2, mat2, mat1])
  mat  = LinearMaterial.Create([2.0, 4.0, 10.0, 12.0], mats, "LINEARMAT")

  # get the steady state solution
  mat_ss = downcast(mat1)
  state = get_steady_state(inp, mat_ss, mesh)

  #dt = converged_step(inp, tdmat, mesh, state, 50)
  # run the converged time step on different meshes
  #converged_mesh(inp, tdmat, 1.0/8.0, 50)

  Manager.finalize()
