from model import *


def run_steady_state():
    inp = get_input()
    inp.put_int("adjoint", 0)
    mesh = get_mesh(20)
    mat = get_control_materials(*[1., 1., 1., 1., 1.])
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve()
    ic_f = manager.state()

    inp.put_int("adjoint", 1)
    mat = get_control_materials(*[1., 1., 1., 1., 1.])
    manager_a = Eigen1D(inp, downcast(mat), mesh)
    manager_a.solve()
    ic_a = manager_a.state()

    #plot_mesh_function(mesh, ic.phi(0))
    #plot_mesh_function(mesh, ic.phi(1))
    plt.figure(1)
    #plt.plot(ic_f.phi(0))
    #plt.plot(ic_f.phi(1))
    plt.plot(np.array(ic_a.phi(0)))
    plt.plot(np.array(ic_a.phi(1)))
    #plt.legend(['$\phi_{f}$', '$\phi_{t}$',
    #            '$\phi^+_{f}$', '$\phi^+_{t}$'])
   # plt.show()
    #mat.display()
    return ic_f
    

if __name__ == "__main__":

    run_steady_state()
    #solvers.Manager.finalize()
