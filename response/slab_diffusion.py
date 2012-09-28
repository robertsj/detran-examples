import sys
from detran import *
from ResponseDB import *
from Response import *
Manager.initialize(sys.argv)
#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_int("number_groups",              1)
inp.put_str("bc_west",                    "vacuum")
inp.put_str("bc_east",                    "vacuum")
inp.put_int("diffusion_fixed_type",       1)
inp.put_dbl("diffusion_tolerance",        1e-10)
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
mat = Material.Create(1, 1, False)
mat.set_sigma_t(0, 0,     1.0) 
mat.set_sigma_f(0, 0,     0.5) 
mat.set_sigma_s(0, 0, 0,  0.5) 
mat.set_diff_coef(0, 0,   1.0/3.0)
mat.compute_sigma_a()
mat.finalize()
#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
cm = [0.0, 10.0]
fm = vec_int(1, 4)
mesh = Mesh1D.Create(fm, cm, [0])
#-----------------------------------------------------------------------------#
# Run Response Generation
#-----------------------------------------------------------------------------#

# create new response db
db = ResponseDB("test.h5")
# create new response (i.e. a new node)
r = Response(inp, mat, mesh, db, True, [1.0, 1.1])
r.run()




Manager.finalize()

