from detran import *
import pickle

from slab_materials import *
from input_parameters import get_input

def get_ramp_materials() :
    """ Create the two group cross sections for the reference slab problem
    """
    #   0   0   1    1    0    
    #   0---2^^^4---10vvv12----\\\
    times = [2.0, 4.0, 10.0, 12.0]
    mat0 = get_rods_in_materials()
    mat1 = get_rods_out_materials()
    mats = vec_material()
    mats.push_back(mat0)
    mats.push_back(mat1)
    mats.push_back(mat1)
    mats.push_back(mat0)
    mat = LinearMaterial.Create(times, mats)
    return mat

def get_mesh(cells_per_assembly=10) :
    """ Get the 5-assembly (with reflectors) mesh.  Check the cells per assembly
        is sane for the coarse mesh factor.
    """
    cm = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    fm = [cells_per_assembly]*7
    mt = [0, 1, 2, 3, 4, 5, 0]
    mesh = Mesh1D.Create(fm, cm, mt)
    return mesh

def compute_power(mesh, mat, state):
    """ Integrate the core fission rate, and multiply by the energy per fission
        to compute the core power.
    """
    matmap = mesh.mesh_map("MATERIAL");
    F = 0.0;
    for i in range(0, mesh.number_cells()) :
        m = matmap[i]
        F += mesh.volume(0) * \
            (state.phi(0)[i] * mat.sigma_f(m, 0) + \
             state.phi(1)[i] * mat.sigma_f(m, 1))
    F *=  200 * 1.609e-13 # (f/s) * (MeV/f) * (J/MeV)
    return F

def run_eigen(case='diffusion'):
    inp = get_input(case)
    mesh = get_mesh(3)
    mat = get_rods_in_materials()
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    return manager.state(), mat

def run_steady_state(inp, mesh, mat):
    mat.update(0.0, 0, 1, False)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 
    ic = manager.state()
    mat.set_eigenvalue(ic.eigenvalue())
    mat.update(0, 0, 1, False)
    F = compute_power(mesh, mat, ic)
    P0 =  100 # W
    ic.scale(P0/F)
    return ic

def run_transient(case='diffusion'):
    
    inp = get_input(case)
    # There is no feedback, so there is no need to iterate.  If the solver 
    # tolerances are looser than the time stepper tolerance, additional
    # iterations *could* be done, so just stop that and adjust the right tolerance!
    inp.put_int("ts_max_iters", 1)
    mesh = get_mesh(10)
    mat = get_ramp_materials()

    ic = run_steady_state(inp, mesh, mat)

    
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
    np.savetxt('powers_'+case+'.txt', np.array([times,powers]).T)
    pickle.dump(data, open('transient_0_01s.p', 'wb'))
    plt.semilogy(times, powers)
    plt.show()

if __name__ == '__main__':
    #ic, mat = run_eigen('transport')
    run_transient('transport_iso')