import numpy as np
import matplotlib.pyplot as plt

from detran import *
import pickle

from slab_materials import *
from slab_geometry import get_core_mesh
from input_parameters import get_input


def get_linear_source(mesh, quad) :
    """ Create the two group cross sections for the reference slab problem
    """
    #   0   0   1    
    #   0---1---1.1
    times = [1.0, 1.1]
    source0 =  ConstantSource(2, mesh, 1e-14, quad)
    source1 =  ConstantSource(2, mesh, 1e8, quad)
    sources = vec_source()
    sources.push_back(source0)
    sources.push_back(source1)
    mat = LinearExternalSource(2, mesh, times, sources)
    return mat

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

def run_eigen(case='diffusion', fractions=[0, 0]):
    inp = get_input(case)
    mesh = get_core_mesh(3)
    mat = get_materials(*fractions)
    manager = Eigen1D(inp, mat, mesh)
    manager.solve() 
    return manager.state(), mat

def run_steady_state(inp, mesh, mat):
    mat.update(0.0, 0, 1, False)
    manager = Eigen1D(inp, mat, mesh)
    manager.solve() 
    ic = manager.state()
    mat.set_eigenvalue(ic.eigenvalue())
    mat.update(0, 0, 1, False)
    F = compute_power(mesh, mat, ic)
    P0 =  100 # W
    ic.scale(P0/F)
    return ic

def run_transient(inp, mesh, mat):
    
    #inp = get_input(case)
    inp.put_int("ts_max_iters", 1)
    # mesh = get_core_mesh(10)
    # mat = get_materials(*fractions)
    tmat = TimeIndependentMaterial(mat)


    
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
        if step % 10 == 0:
            print( "{:8.4} {:8.4}".format(t, P))
    

    ts = Time1D(inp, tmat, mesh, True)
    ic = ts.state()
    ic.clear() # Ensures set to zero (refactor!)
    ts.set_monitor(monitor)
    source = get_linear_source(mesh, ts.quadrature())
    ts.add_source(source)
    ts.solve(ic)
    print("MAX POWER = ", max(powers))
    print("FINAL POWER = ", powers[-1], times[-1])
    return times, powers

if __name__ == '__main__':
    #ic, mat = run_eigen("diffusion")
    run_transient("diffusion", [0.0, 0.0])