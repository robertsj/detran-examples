from model import *


def run_steady_state():
    """ Run steady state problem in forward and adjoint mode.


    Try 100% for a step insertion.  Do same for transient.
    """

    inp = get_input()
    mesh = get_mesh(10)
    mat = get_control_materials(*[1., 1., 1., 1., 1., 1.])

    ## For k-inf with 100% withdrawn fuel, uncomment these lines
    # inp.put_str("bc_west", "reflect")
    # inp.put_str("bc_east", "reflect")
    # mat = get_control_materials(*[1., 1., 1., 1., 1., "fuel"])


    inp.put_int("adjoint", 0)
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve()
    ic_f = manager.state()
    k_f = ic_f.eigenvalue()

    inp.put_int("adjoint", 1)
    manager_a = Eigen1D(inp, downcast(mat), mesh)
    manager_a.solve()
    ic_a = manager_a.state()
    k_a = ic_a.eigenvalue()

    print("Forward keff = ", k_f)
    print("Adjoint keff = ", k_a)


    plt.figure(1)
    plt.plot(ic_f.phi(0))
    plt.plot(ic_f.phi(1))
    plt.legend(['$\phi_{f}$', '$\phi_{t}$'])
    
    #plt.figure(2)
    #plt.plot(ic_a.phi(0))
    #plt.plot(ic_a.phi(1))
    #plt.legend(['$\phi^+_{f}$', '$\phi^+_{t}$'])

    plt.show()
    mat.display()
    return ic_f


def run_transient():
    
    # STEADY STATE
    inp = get_input()
    inp.put_int("adjoint", 0)
    mesh = get_mesh(10)
    mat = get_materials()
    mat.update(0.0, 0, 1, False)
    print(type(mat))
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 

    #return
    ic = manager.state()   
    print("keff =", ic.eigenvalue())
    #return
    mat.set_eigenvalue(ic.eigenvalue())
    mat.update(0, 0, 1, False)
    # normalize state.
    F = compute_power(mesh, mat, ic)
    P0 =  100 # W
    ic.scale(P0/F)

    # TRANSIENT
    inp.put_int("ts_max_steps",                   10000)
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
    plt.plot(times, powers)
    plt.show()
    

if __name__ == "__main__":

    run_steady_state()
    #solvers.Manager.finalize()

    k0 = 1.12196509e+0
    k1 = 1.09362654e+00
    k3 = 1.03935009e+00
    k1and3 = 1.00797253e+00

    rho1 = (k1-1)/k1 - (k0-1)/k0 
    rho3 = (k3-1)/k3 - (k0-1)/k0 

    print(rho1, rho3, rho3/rho1)
