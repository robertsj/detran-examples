"""
1-D Slab Reactor

Slab reactor with five assemblies.  Reflected on 
both sides. 

| R  | A1 | A2 | A3 | A4 | A5 | R  |
0    10   20   30   40   50   60   70

"""
from detran import *
import pickle

def get_input() :
    inp = utilities.InputDB.Create()
    inp.put_int("number_groups",                  2)
    inp.put_int("dimension",                      2)
    inp.put_str("equation",                       "diffusion")
    inp.put_int("adjoint",                        0)
    inp.put_str("bc_west",                        "vacuum")
    inp.put_str("bc_east",                        "vacuum")
    inp.put_str("eigen_solver",                   "diffusion")
    inp.put_dbl("eigen_tolerance",                1e-12)
    inp.put_int("eigen_max_iters",                1000)
    inp.put_str("outer_solver",                   "GS")
    inp.put_dbl("outer_tolerance",                1e-12)
    inp.put_int("outer_max_iters",                1000)
    inp.put_int("outer_print_level",              0)
    inp.put_int("outer_krylov_group_cutoff",      0)
    inp.put_str("outer_pc_type",                  "mgdsa")
    inp.put_str("inner_solver",                   "SI")
    inp.put_dbl("inner_tolerance",                1e-12)
    inp.put_int("inner_max_iters",                1000)
    inp.put_int("inner_print_level",             1)
    # gmres parameters
    db = utilities.InputDB.Create("callow_db")
   # db.put_dbl("linear_solver_atol",              0.0);
    db.put_dbl("linear_solver_rtol",              1e-13);
    db.put_str("linear_solver_type",              "gmres");
    db.put_int("linear_solver_maxit",             1000);
    db.put_int("linear_solver_gmres_restart",     30);
    db.put_int("linear_solver_monitor_level",     0);
    db.put_str("pc_type",                         "petsc");
    db.put_str("petsc_pc_type",                   "lu");
    db.put_int("petsc_pc_factor_levels",          3);
    db.put_str("eigen_solver_type",               "power");
    db.put_int("eigen_solver_maxit",              1000);
    db.put_int("eigen_solver_monitor_level",      2);
    db.put_dbl("eigen_solver_tol",                1.0e-10)
    inp.put_spdb("inner_solver_db", db)
    inp.put_spdb("inner_pc_db", db)
    inp.put_spdb("outer_solver_db", db) 
    inp.put_spdb("eigen_solver_db", db)
    inp.put_int("ts_max_steps",                   10000)
    inp.put_int("ts_scheme",                      Time1D.IMP)
    inp.put_int("ts_output",                      0)
    inp.put_dbl("ts_step_size",                   0.01)
    inp.put_dbl("ts_final_time",                  3.0)
    #inp.put_int("ts_no_extrapolation",            1)
    inp.put_int("ts_max_iters",                   10)
    inp.put_dbl("ts_tolerance",                   1.0e-5)
    #
    preconditioner_db = utilities.InputDB.Create("preconditioner_db")
    preconditioner_db.put_dbl("linear_solver_atol",              0.0);
    preconditioner_db.put_dbl("linear_solver_rtol",              0.1);
    preconditioner_db.put_str("linear_solver_type",              "gmres");
    preconditioner_db.put_int("linear_solver_maxit",             5000);
    preconditioner_db.put_int("linear_solver_gmres_restart",     30);
    preconditioner_db.put_int("linear_solver_monitor_level",     0);
    preconditioner_db.put_str("pc_type",                         "ilu0");
    preconditioner_db.put_str("petsc_pc_type",                   "ilu");
    preconditioner_db.put_int("petsc_pc_factor_levels",          2);
    #
    return inp

def get_mesh(ffm, rods=False) :
    """ Return the 1-D Mesh
    """
    cm = [10.0*i for i in range(8)]
    fm = [ffm]*7
    mt = [1, 1, 2, 3, 4, 5, 1]
    mesh = Mesh1D.Create(fm, cm, mt)
    return mesh

def fill_water(m, mat):
    mat.set_diff_coef(m, 0,   1.500)
    mat.set_diff_coef(m, 1,   0.500)
    mat.set_sigma_t(m, 0,     0.0002 + 0.032)
    mat.set_sigma_t(m, 1,     0.010)
    mat.set_sigma_s(m, 1, 0,  0.032)
    mat.set_sigma_s(m, 0, 1,  0.000)
    mat.set_sigma_f(m, 0,     0.000)
    mat.set_sigma_f(m, 1,     0.0000)
    
def fill_fuel(m, mat, withdrawn_fraction):
    mat.set_diff_coef(m, 0,   1.300)
    mat.set_diff_coef(m, 1,   0.500)
    mat.set_sigma_t(m, 0,     0.0105 + 0.022)
    mat.set_sigma_t(m, 1,     0.164  - withdrawn_fraction*(.05))
    mat.set_sigma_s(m, 1, 0,  0.022)
    mat.set_sigma_f(m, 0,     0.003)
    mat.set_sigma_f(m, 1,     0.190)
    mat.set_chi(m, 0,         1.0)

def fill_kinetics(mat):
    beta = np.array([.000218, .001023, .000605, .00131, .00220, .00060, .000540, .000152])
    lam  = np.array([0.012467, 0.028292, 0.042524, 0.133042, 0.292467, 0.666488, 1.634781, 3.554601])
    velocity = [2200 * 100 * np.sqrt(0.1e4 / 0.0253),  2200 * 100 * np.sqrt(0.1 / 0.0253)]
    for g in range(2):
        mat.set_velocity(g, velocity[g])
    for i in range(8):
        mat.set_lambda(i, lam[i])
        mat.set_beta(i, beta[i])
        for m in range(mat.number_materials()):
            mat.set_chi_d(m, i, 0, 1.0)
            mat.set_chi_d(m, i, 1, 0.0)
    return mat
 
    
def get_control_materials(a1=1.0, a2=0.25, a3=1.0, a4=0.25, a5=1.0):
    mat = KineticsMaterial.Create(6, 2, 8, "CONTROL")
    fill_water(0, mat)
    fill_fuel(1, mat, withdrawn_fraction=a1)  
    fill_fuel(2, mat, withdrawn_fraction=a2)  
    fill_fuel(3, mat, withdrawn_fraction=a3) 
    fill_fuel(4, mat, withdrawn_fraction=a4)
    fill_fuel(5, mat, withdrawn_fraction=a5)
    fill_kinetics(mat)
    mat.finalize()

    return mat

def get_rods_in_materials():
    mat = KineticsMaterial.Create(6, 2, 8, "RODS_IN")
    fill_water(0, mat)
    fill_fuel(1, mat, withdrawn_fraction=1.00)  
    fill_fuel(2, mat, withdrawn_fraction=0.25)  
    fill_fuel(3, mat, withdrawn_fraction=1.00) 
    fill_fuel(4, mat, withdrawn_fraction=0.25)
    fill_fuel(5, mat, withdrawn_fraction=1.00)
    fill_kinetics(mat)
    mat.finalize()

    return mat

def get_rods_out_materials():
    mat = KineticsMaterial.Create(6, 2, 8, "RODS_OUT")
    fill_water(0, mat)
    fill_fuel(1, mat, withdrawn_fraction=1.00)  
    fill_fuel(2, mat, withdrawn_fraction=0.25)  
    fill_fuel(3, mat, withdrawn_fraction=1.00) 
    fill_fuel(4, mat, withdrawn_fraction=0.30)
    fill_fuel(5, mat, withdrawn_fraction=1.00)
    fill_kinetics(mat)
    mat.finalize()
    return mat

def get_materials() :
    """ Create the two group cross sections for the reference slab problem
    """
    times = [2.0, 4.0, 10.0, 12.0]
    mat0 = get_rods_in_materials()
    mat1 = get_rods_out_materials()
    mat2 = get_rods_in_materials()
    mats = vec_material()
    mats.push_back(mat0)
    mats.push_back(mat1)
    mats.push_back(mat1)
    mats.push_back(mat2)
    mat = LinearMaterial.Create(times, mats)
    return mat

def make_material_table():
    mat = KineticsMaterial.Create(4, 2, 8, "RODS_OUT")
    fill_water(0, mat)
    fill_fuel(1, mat, 1.0)
    fill_fuel(2, mat, 0.0)
    fill_kinetics(mat)
    mat.finalize()
    
    s = r"\hline" + '\n'
    s += r"material & $D_1$  & $D_2$  & $\Sigma_{a1}$  & $\Sigma_{a2}$  & $\Sigma_{s2\gets 1}$ & $\nu\Sigma_{f1}$ &   $\nu\Sigma_{f2}$ \\"
    s += r"\hline" + '\n'
    #                 D1         D2         A1         A2           S21      nuF1     nuF2
    tmpl = r" {:15} & {:8.4f} &  {:8.4f} &  {:8.4f} &   {:8.4f} &  {:8.4f} & {:8.4f} & {:8.4f} \\" 
    names = ['reflector', 'fuel (all out)', 'fuel (all in)']
    for m in range(3):
        D0 = mat.diff_coef(m, 0)
        D1 = mat.diff_coef(m, 1)
        A1 = mat.sigma_t(m, 0) - mat.sigma_s(m, 1, 0)
        A2 = mat.sigma_t(m, 1)
        S21 = mat.sigma_s(m, 1, 0)
        nF1 = mat.nu_sigma_f(m, 0)
        nF2 = mat.nu_sigma_f(m, 1)
        s += tmpl.format(names[m], D0, D1, A1, A2, S21, nF1, nF2) + '\n'
        s += r"\hline" + '\n'
    print(s)
    print("")
    s = ""
    s = r"\hline" + '\n'
    s += r"group & $\beta_i$ & $\lambda_i$ \\" + '\n'
    s += r"\hline" + '\n'
    for i in range(mat.number_precursor_groups()):
        beta = mat.beta(0, i)
        lam = mat._lambda(i)
        #mat.
        s += r" {:8} & {:12.2e} & {:12.6e} \\".format(i, beta, lam) + '\n'
        s += r"\hline" + '\n'
    print(s)
    

    
   # mat.set_eigenvalue(ic.eigenvalue())
   # mat.update(0, 0, 1, False) 

def compute_worth():
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
    mat = get_materials()
    mat.update(0.0, 0, 1, False)
    print(type(mat))
    manager = Eigen1D(inp, downcast(mat), mesh)
    manager.solve() 

    #return
    ic = manager.state()   
    print("keff =", ic.eigenvalue())
    return
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
    solvers.Manager.initialize(sys.argv)
    #run_transient()
   #mat = get_materials()
   # mat.display()
   # make_material_table()
    print(compute_worth())
    solvers.Manager.finalize()
