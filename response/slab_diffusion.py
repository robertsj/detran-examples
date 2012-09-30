import sys
from detran import *
from ResponseDB import *
from Response import *
Manager.initialize(sys.argv)
#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_int("number_groups",              2)
inp.put_int("diffusion_fixed_type",       1)
inp.put_int("response_spatial_order",     1)
inp.put_int("diffusion_monitor_level",    0)
inp.put_dbl("diffusion_tolerance",        1e-10)
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
#mat = Material.Create(1, 1, False)
#mat.set_sigma_t(0, 0,     1.0) 
#mat.set_sigma_f(0, 0,     0.5) 
#mat.set_chi(0, 0,         1.0) 
#mat.set_sigma_s(0, 0, 0,  0.5) 
#mat.set_diff_coef(0, 0,   1.0/3.0)
#mat.compute_sigma_a()
#mat.finalize()

mat = Material.Create(4, 2, False)
# Material 0
mat.set_sigma_t(0, vec_dbl([0.030120, 0.080032]))
mat.set_sigma_s(0, 1, 0, 0.02) 
mat.set_diff_coef(0, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(0, 1, 0.135)
mat.set_chi(0, 0, 1.0)
# Material 1
mat.set_sigma_t(1, vec_dbl([0.030120, 0.085032]))
mat.set_sigma_s(1, 1, 0, 0.02) 
mat.set_diff_coef(1, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(1, 1, 0.135)
mat.set_chi(1, 0, 1.0)
# Material 2
mat.set_sigma_t(2, vec_dbl([0.030120, 0.130032]))
mat.set_sigma_s(2, 1, 0, 0.02) 
mat.set_diff_coef(2, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(2, 1, 0.135)
mat.set_chi(2, 0, 1.0)
## Material 3
mat.set_sigma_t(3, vec_dbl([0.040160, 0.010024]))
mat.set_sigma_s(3, 1, 0, 0.04) 
mat.set_diff_coef(3, vec_dbl([2.0, 0.3]))
mat.compute_sigma_a()
mat.finalize()

#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
cm = [0.0, 20.0]
fm = vec_int(1, 20)
mesh = Mesh1D.Create(fm, cm, [0])
#mesh = Mesh2D.Create(fm, fm, cm, cm, [0])
#mesh = Mesh3D.Create(fm, fm, fm, cm, cm, cm, [0])
#-----------------------------------------------------------------------------#
# Run Response Generation
#-----------------------------------------------------------------------------#

# create new response db
db = ResponseDB("test.h5")
# create new response (i.e. a new node)
r = Response(inp, mat, mesh, db, True, [1.120069900326722])#np.linspace(0.8, 1.2, 100))
r.run()

#poo = db.file["/fuel"]["R"]
#print poo[0,:,:]



Manager.finalize()

