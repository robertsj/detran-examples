"""
1-D Slab Reactor

Slab reactor with five assemblies.  Reflected on 
both sides. 

| R  | A1 | A2 | A3 | A4 | A5 | R  |
0    10   20   30   40   50   60   70 cm

"""

from detran import Mesh1D

def get_core_mesh(cells_per_assembly=10) :
    """ Get the 5-assembly (with reflectors) mesh.  Check the cells per assembly
        is sane for the coarse mesh factor.
    """
    cm = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    fm = [cells_per_assembly]*7
    mt = [0, 1, 2, 3, 4, 5, 0]
    mesh = Mesh1D(fm, cm, mt)
    return mesh

def get_assembly_mesh(cells_per_assembly=10) :
    """ Get mesh for just one assembly. 
    """
    cm = [0.0, 10.0]
    fm = [cells_per_assembly]
    mt = [2]
    mesh = Mesh1D(fm, cm, mt)
    return mesh


if __name__ == "__main__":

    get_core_mesh(3).display()

from detran import *
import pickle

from slab_materials import get_ramp_materials
from input_parameters import get_input

def get_core_mesh(cells_per_assembly=10) :
    """ Return a mesh for the whole core.
    """
    cm = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    fm = [cells_per_assembly]*7
    mt = [0, 1, 2, 3, 4, 5, 0]
    mesh = Mesh1D(fm, cm, mt)
    return mesh

def get_assembly_mesh(cells_per_assembly=10) :
    """ Return the mesh for a single fuel assembly.
    """
    cm = [0.0, 10.0]
    fm = [cells_per_assembly]*7
    mt = [0, 1, 2, 3, 4, 5, 0]
    mesh = Mesh1D(fm, cm, mt)
    return mesh


def compute_worth():
    """ Compute reactivity difference between 
    the rods in and rods out cases.
    """

    inp = get_input()
    mesh = get_mesh(50)
    mat = get_materials()
    
    mat.update(0.0, 0, 1, False)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    k_in = manager.state().eigenvalue()
    
    mat.update(10.0, 0, 1, False)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    k_out = manager.state().eigenvalue()
    return (k_out-k_in)/k_in, mat.beta_total(0)
    
def compute_beta_eff():
    
    mat = get_materials()
    mat.update(0.0, 0, 1, False)
    mesh = get_mesh(50)

    # Compute phi
    inp = get_input()
    inp.put_int("adjoint", 0)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    F = manager.state() 
    
    # Compute phi
    inp.put_int("adjoint", 1)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    A = manager.state()

    num = 0
    den = 0
    mat_map = mesh.mesh_map("MATERIAL")
    for i in range(mesh.number_cells()):
        m = mat_map[i]
        tmp = A.phi(0)[i]*mat.nu_sigma_f(m, 0)*F.phi(0)[i]
        num +=  beta
              
    beta = mat.beta_total(0)
   
    
def compute_power(mesh, mat, state):
    matmap = mesh.mesh_map("MATERIAL");
    F = 0.0;
    for i in range(0, mesh.number_cells()) :
        m = matmap[i]
        F += mesh.volume(0) * \
            (state.phi(0)[i] * mat.sigma_f(m, 0) + \
             state.phi(1)[i] * mat.sigma_f(m, 1))
    F *=  200 * 1.609e-13 # (f/s) * (MeV/f) * (J/MeV)
    return F



    
def run_transient():
    
    # STEADY STATE
    inp = get_input()
    inp.put_int("adjoint", 0)
    mesh = get_mesh(10)
    mat = get_ramp_materials()
    mat.update(0.0, 0, 1, False)
    print(type(mat))
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 

    #return
    ic = manager.state()   
    print("keff =", ic.eigenvalue())

    mat.set_eigenvalue(ic.eigenvalue())
    mat.update(0, 0, 1, False)
    # normalize state.
    F = compute_power(mesh, mat, ic)
    P0 =  100 # W
    ic.scale(P0/F)

    # TRANSIENT
    inp.put_int("ts_max_steps",                   100000)
    inp.put_int("ts_scheme",                      Time1D.BDF2)
    inp.put_int("ts_output",                      0)
    inp.put_dbl("ts_step_size",                   0.1)
    inp.put_dbl("ts_final_time",                  60.0)
    inp.put_int("ts_max_iters",                   1)
    inp.put_dbl("ts_tolerance",                   1.0e-5)
    
    
    data = {}
    data['times'] = []
    data['powers'] = []
    data['flux'] = []
    data['precursors'] = []
    data['sigma_t_thermal'] = []
    data['info'] = \
    """This file includes the times, fluxes (two group),
precursors (8 group), power (core)"""

    x = np.zeros(mesh.number_cells())
    x[0] = mesh.dx(0)/2
    for i in range(1, mesh.number_cells()):
        x[i] = x[i-1] + 0.5*(mesh.dx(i-1)+mesh.dx(i))
    data['x'] = x

    #data['phi]
    times, powers = [], []
    
    def monitor(ts, step, t, dt, it, conv) :
        
        state = ts.state()
        precursor = ts.precursor()
        mesh = ts.mesh()
        mat = ts.material()
        P = compute_power(mesh, mat, state)

        data['times'].append(t)
        data['powers'].append(P)
        data['flux'].append([np.array(state.phi(i)) for i in range(2)])
        prec = [np.array(precursor.C(i)) for i in range(8)]
        data['precursors'].append(prec)

        sig_t = np.zeros(mesh.number_cells())
        mat_map = mesh.mesh_map('MATERIAL')
        for i in range(mesh.number_cells()):
            
            sig_t[i] = mat.sigma_t(mat_map[i], 1)
        data['sigma_t_thermal'].append(sig_t)

        times.append(t)
        powers.append(P)
        print( "{:2} {:2} {:8.4} {:8.4} {:8.4}".format(step, it, t, dt, P))
    

    ts = Time1D(inp, mat, mesh, True)
    ts.set_monitor(monitor)
    ts.solve(ic)
    print("MAX POWER = ", max(powers))
    print("FINAL POWER = ", powers[-1], times[-1])
    pickle.dump(data, open('transient_0_01s.p', 'wb'))
    plt.semilogy(times, powers)
    #plt.show()
    
if __name__ == "__main__":

    solvers.Manager.initialize(sys.argv)
    run_transient()
   #mat = get_materials()
   # mat.display()
   # make_material_table()
    #print(compute_worth())

    
    solvers.Manager.finalize()
