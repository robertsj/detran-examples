# pyexamples/run_assembly_fixed.py
#
# Test the efficacy of preconditioners for multigroup fixed 
# source multiplying problems.  Moreover, we also look to see
# what sort of convergence criteria yields 1%, 0.1%, and 0.01% 
# on pin powers.  Parameters to investigate:
#   SI-GS, converging inners or not
#   Krylov, number of restarts [10, 20, 30, 40, 50]
#   quadrature  1x1, 2x2, 3x3, ..., 10x10
#   DSA solver tolerance LU vs ILU(2) @ 0.001, 0.01, 0.1

import numpy as np
import time
import sys
from detran import *
sys.path.append('../c5g7/')
from assemblies_c5g7 import get_assemblies
from pins_c5g7 import get_pins
from material_c5g7 import get_materials
import time
import cPickle as pickle

q = 8 # number of polar/azimuth per octant

# Iterations for 1e-3 on fission rates


def get_input() :
    inp = utilities.InputDB.Create()
    inp.put_str("equation", "dd")
    inp.put_int("store_angular_flux", 0)
    inp.put_int("store_current", 0)
    inp.put_dbl("cmfd_relaxation", 1.0)
    inp.put_str("problem_type", "fixed")
    inp.put_int("number_groups", 7)
    inp.put_str("inner_solver", "SI")
    inp.put_int("inner_max_iters", 1000000)
    inp.put_dbl("inner_tolerance", 1.0e-16)
    inp.put_int("inner_print_level", 0)
    inp.put_str("outer_solver", "GS")
    inp.put_int("outer_max_iters", 100000)
    inp.put_int("outer_print_level", 2)
    inp.put_int("outer_print_interval", 1)
    inp.put_int("mgdsa_disable_fission", 0)
    #inp.put_str("bc_west", "fixed")
    inp.put_str("quad_type", "u-dgl")
    inp.put_int("quad_number_polar_octant",   q)
    inp.put_int("quad_number_azimuth_octant", q)
    inp.put_int("inner_pc_side", 2)
    inp.put_int("outer_pc_side", 2)
    inp.put_int("outer_krylov_group_cutoff", 0)
    solver_db = utilities.InputDB.Create("solver_db")
    solver_db.put_dbl("linear_solver_atol",              0.0);
    solver_db.put_str("linear_solver_type",              "gmres");
    solver_db.put_int("linear_solver_maxit",             50000);
    solver_db.put_int("linear_solver_gmres_restart",     30);
    solver_db.put_int("linear_solver_monitor_level",     2);
    preconditioner_db = utilities.InputDB.Create("preconditioner_db")
    preconditioner_db.put_dbl("linear_solver_atol",              0.0);
    preconditioner_db.put_dbl("linear_solver_rtol",              0.01);
    preconditioner_db.put_str("linear_solver_type",              "gmres");
    preconditioner_db.put_int("linear_solver_maxit",             5000);
    preconditioner_db.put_int("linear_solver_gmres_restart",     30);
    preconditioner_db.put_int("linear_solver_monitor_level",     0);
    preconditioner_db.put_str("pc_type",                         "ilut");
    preconditioner_db.put_str("petsc_pc_type",                   "ilut");
    preconditioner_db.put_int("petsc_pc_factor_levels",          2);
    inp.put_spdb("inner_solver_db", solver_db)
    inp.put_spdb("inner_pc_db",     preconditioner_db)
    inp.put_spdb("outer_solver_db", solver_db)
    inp.put_spdb("outer_pc_db",     preconditioner_db)
    return inp

def get_mesh() :
  assemblies = get_assemblies(7, True)
  mesh = assemblies[0].mesh()
  return mesh

def get_data() :
  inp = get_input()
#   inp.get_spdb("inner_pc_db").put_dbl("linear_solver_rtol",  pctol);
#   inp.get_spdb("outer_pc_db").put_dbl("linear_solver_rtol",  pctol);
#   inp.get_spdb("inner_solver_db").put_int("linear_solver_gmres_restart", r);
#   inp.get_spdb("outer_solver_db").put_int("linear_solver_gmres_restart", r);
  return inp, get_materials(), get_mesh()

def get_solver(type, tol) :
    """ Gets the solver
    """ 
    inp, mat, mesh = get_data()
    set_solver(type, inp, tol)
    solver = Fixed2D(inp, mat, mesh, True)
    solver.setup()
    q = ConstantSource.Create(7, solver.mesh(), 1.0)
    solver.set_source(q)
    solver.set_solver()
    #set_boundary(solver)
    return solver, mat, mesh

def set_boundary(solver) :
    #b = cast_sn(solver.boundary())
    #bc = cast_fixed(b.bc(0)) 
    #quad = solver.quadrature()
    #g = 0 
    #for o in range(0, 2) :
    #    for a in range(0, quad.number_angles_octant()) :
    #        for j in range(0, solver.mesh().number_cells_y()) :
    #            bc[o, a, g, j] = 1.0  
    pass

def set_solver(type, inp, tol) :
    """ Set database parameters for a specific solver/preconditioner 
        scheme.  Tolerance defaults to ~0 so that the number of iterations
        is selected.
    """
    assert(len(type) == 4)
    print 'TYPE=', type
    # SETUP INNER
    if type[0] == 0 :
        inp.put_str("inner_solver", "SI")
        inp.put_int("inner_max_iters", 20);
        inp.put_dbl("inner_tolerance", tol);
        inp.get_spdb("inner_solver_db").put_dbl("linear_solver_rtol", tol);
    elif type[0] == 1 :
        inp.put_str("inner_solver", "SI")
        inp.put_int("inner_max_iters", 1); # flattened inners for SI
    elif type[0] == 2 :
        inp.put_str("inner_solver", "GMRES")
        inp.put_int("inner_max_iters", 20);
        inp.get_spdb("inner_solver_db").put_int("linear_solver_maxit", 20);
        inp.put_dbl("inner_tolerance", tol);
        inp.get_spdb("inner_solver_db").put_dbl("linear_solver_rtol", tol);
    else :
        error("invalid inner")  

    # SETUP OUTER
    if type[1] == 0 :
        inp.put_str("outer_solver", "GS")
    elif type[1] == 1 :
        inp.put_str("outer_solver", "GMRES")
    elif type[1] == 2 :
        inp.put_str("outer_solver", "CMFD")
        inp.put_int("store_current", 1)
    else :
        error("invalid outer")  
    inp.put_dbl("outer_tolerance", tol);
    inp.get_spdb("outer_solver_db").put_dbl("linear_solver_rtol", tol);

    # SETUP INNER PC
    if type[2] == 0 :
        inp.put_str("inner_pc_type", "none")
    elif type[2] == 1 :
        inp.put_str("inner_pc_type", "DSA") 
        inp.get_spdb("inner_pc_db").put_str("petsc_pc_type", "lu");
    elif type[2] == 2 :
        inp.put_str("inner_pc_type", "DSA") 
        inp.get_spdb("inner_pc_db").put_str("petsc_pc_type", "ilu");
    else :
        error("invalid inner pc")
        
    # SETUP OUTER PC
    if type[3] == 0 :
        inp.put_str("outer_pc_type", "none")       
    elif type[3] == 1 :
        inp.put_str("outer_pc_type", "mgdsa")
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "lu");
    elif type[3] == 2 :
        inp.put_str("outer_pc_type", "mgdsa")
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "ilu");
    elif type[3] == 3 :
        inp.put_str("outer_pc_type", "mgcmdsa")
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "ilu");
        inp.put_int('mgpc_coarse_mesh_level',     7)
        inp.put_int('mgpc_condensation_option',   1)
        inp.put_int('mgpc_cmdsa_use_smoothing',   0)      
    elif type[3] == 4 :
        inp.put_str("outer_pc_type", "mgcmdsa")
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "ilu");
        inp.put_int('mgpc_coarse_mesh_level',     7)
        inp.put_int('mgpc_condensation_option',   1)
        inp.put_int('mgpc_cmdsa_use_smoothing',   1)
        inp.put_int('mgpc_cmdsa_smoothing_iters', 3)
        inp.put_dbl('mgpc_cmdsa_smoothing_relax', 0.7)  
    elif type[3] == 5 :
        inp.put_str("outer_pc_type", "mgtcdsa")
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "ilu")
        inp.put_int('mgpc_coarse_mesh_level',               7)
        inp.put_int('mgpc_condensation_option',             1)
        inp.put_int('mgpc_cmdsa_use_smoothing',             1)
        inp.put_int('mgpc_cmdsa_smoothing_iters',           3)
        inp.put_dbl('mgpc_cmdsa_smoothing_relax',           0.7)
        inp.put_int('mgpc_tcdsa_number_coarse_corrections', 3)
        inp.put_int('mgpc_tcdsa_use_fine_correction',       1)
        inp.put_int('mgpc_number_azimuth_octant',           2)
        inp.put_int('mgpc_number_polar_octant',             2)      
    elif type[3] == 6 :
        inp.get_spdb("outer_pc_db").put_str("petsc_pc_type", "ilu");
        inp.get_spdb("outer_pc_db").put_dbl("linear_solver_rtol",  1.0e-12);
        inp.put_int('mgpc_coarse_mesh_level',     7) 
    else :
        error("invalid outer pc")
    
def run_reference() :
    """ Compute reference solution for determining convergence criteria 
        of all solver schemes.  Uses unpreconditioned MG-GMRES
    """
    solver, mat, mesh = get_solver([0, 1, 0, 0], 1.0e-12)
    t = time.time()
    solver.solve()
    print " reference elapsed = ", time.time()-t
    state = solver.state()
    rates = ReactionRates(mat, mesh, state)
    reference_data = {}
    reference_data['fission'] = np.asarray(rates.region_power(("PINS"), 1.0))
    for g in range(0, 7) :
        reference_data['phi'+str(g)] = np.asarray(state.phi(g))
    pickle.dump(reference_data, open('reference_data.p', 'wb'))


def flux_norm(state, mat, mesh, data) :
    """ Compute ||phi-ref||_2 / ||ref||_2
    """
    nrm = 0.0
    nrm_phi_r = 0.0
    for g in range(0, state.number_groups()) :
        phi        = vec_asarray(state.phi(g))
        phi_r      = data['phi'+str(g)]
        assert(np.linalg.norm(phi) > 0.0)
        nrm       += np.linalg.norm(phi - phi_r) ** 2
        nrm_phi_r += np.linalg.norm(phi_r) ** 2
    nrm_phi = np.sqrt(nrm) / np.sqrt(nrm_phi_r)
    
    mask = data['fission'] > 0.0
    ref_fiss = data['fission'][mask]
    rates = ReactionRates(mat, mesh, state)
    fiss = np.asarray(rates.region_power(("PINS"), 1.0))[mask]
    nrm_fiss = np.max(np.abs(fiss-ref_fiss)/ref_fiss)
    return nrm_phi, nrm_fiss

def gs_compare_residuals() :
    """
    """
    pass

def flux_norm_vs_pin_power() :
    """ Compute a comparison of flux norm vs pin power.
    """
    data = pickle.load(open('reference_data.p', 'rb'))
    n = 10
    nrm_phi = np.zeros(n)
    nrm_fis = np.zeros(n)
    tols    = np.zeros(n)
    tol = 1.0
    for i in range(0, n) :
        tol *= 0.5
        solver, mat, mesh = get_solver([0,0,0,0], tol)
        solver.solve()
        a, b = flux_norm(solver.state(), mat, mesh, data)
        nrm_phi[i] = a 
        nrm_fis[i] = b
        tols[i] = tol
    x = np.array(range(1, n+1))
    tol_data = np.array([tols, nrm_phi, nrm_fis])
    np.savetxt('tol_data.txt', tol_data.T)

def plot_flux_norm_vs_pin_power():

    tols, nrm_phi, nrm_fis = np.loadtxt('tol_data.txt', unpack=True)

    print(tols)
    plt.figure(1)
    plt.loglog(tols, nrm_phi, 'k-o', tols, nrm_fis, 'b-*')
    #plt.legend(['$  \\frac{ || \phi - \phi_{\mathrm{ref}} ||_2 }{ || \phi_{\mathrm{ref}} ||_2 }$', 
    #            '$  \mathrm{max}   \left | \\frac{d - d_{\mathrm{ref}} } { d_{\mathrm{ref}} } \\right | $'],2)
    plt.legend(['a','b'])
    plt.xlabel('$ \mathrm{max}| \phi_{n+1}-\phi_{n} | $ ')
    plt.ylabel('residual norm')
    plt.grid(True)
    plt.savefig('flux_norm_vs_pin_power.pdf')

    plt.figure(2)
    plt.semilogx(tols, nrm_fis / nrm_phi, 'k-o')
    plt.grid(True)
    plt.xlabel('$ \mathrm{max}| \phi_{n+1}-\phi_{n} | $ ')
    plt.ylabel('ratio of norms')
    #plt.show()
    plt.savefig('flux_norm_vs_pin_power_ratio.pdf')
          
def run_convergence(re) :
    """ Decide unpreconditioned solver iterations yielding desired 
        accuracy.  For GS+SI, we use the same iteration count for 
        both GS and SI.
    """ 
  
    # Load reference data
    data = pickle.load(open('reference_data.p', 'rb'))
    mask = data['fission'] > 0.0
    ref_fiss = data['fission'][mask]
  
    # Define cases
    cases = [[0,0,0,0], [1,0,0,0], [2,0,0,0], [0,1,0,0], [1,2,0,6]]
    tols = {}
    sweeps = []
    for n in range(3, 5) :
        # first guess 
        tol0 = 0.01
        print cases[n]
        solver, mat, mesh = get_solver(cases[n], tol0) 
        solver.input().display()
        
        solver.solve()
        re0, tmp = flux_norm(solver.state(), mat, mesh, data)
        # second guess
        tol1 = 0.001
        solver, mat, mesh = get_solver(cases[n], tol1) 
        solver.solve()
        re1, tmp = flux_norm(solver.state(), mat, mesh, data)
        assert(re0 != re1)
        iter = 0
        print re0, re1
        done = False
        while iter < 50 :
            print "iter = ", iter
            a  = (tol1-tol0) / (re1-re0)
            tol = a * (re-re1) + tol1
            if re1 == re0 :
                tol = 1.0e-14
            tol = max(tol, 0.5*tol1)
            print "--->" ,tol1, tol, a
            print tol
            re0  = re1
            tol0 = tol1
            tol1 = tol
            solver, mat, mesh = get_solver(cases[n], tol1)
            solver.solve() 
            re1, tmp = flux_norm(solver.state(), mat, mesh, data)
            print "re =  %12.8e, its=  %4i" % (re1, tol1)
            if (tol0 == tol1 or re1 < re) : 
                print "yeah!", re-re1, re
                tols[cases[n][0], cases[n][1]] = tol1
                done = True
                break
            iter += 1
        assert(done)
        sweeps.append(solver.number_sweeps())
    print sweeps
    
    data = {}
    data['tols'] = tols
    data['sweeps'] = sweeps
    pickle.dump(data, open('convergence_data.p', 'wb'))
    return data

def run(iters, tol) :
    print iters
    #exit('lala')
    cases = []
    # add GS+SI cases
    #for i in range(0, 2) : 
    #    cases.append([i,0,0,0])
    # add GS+GMRES cases
    #for i in range(0, 3) : 
    #    cases.append([2,0,i,0])
    # add MG-GMRES cases
    for i in range(0, 6) : 
        cases.append([0,1,0,i])
    # add CMFD
    cases.append([1,2,0,6])
    num = len(cases)
    
    data = {}
    data["times"]  = np.zeros(num)
    data["sweeps"] = np.zeros(num)
    
    for n in range(0, num) :
        case = cases[n]
        print "case:",n , case
        solver, mat, mesh = get_solver(case, iters[(case[0],case[1])]) 
        t = time.time()
        solver.solve()
        data["times"][n]  = time.time() - t
        data["sweeps"][n] = solver.number_sweeps()
        print data["times"][n]
        print data["sweeps"][n]
    f = "Xdata_" + str(tol) + "_" + str(q) + ".p"
    pickle.dump(data, open(f, 'wb'))

def print_table(data, tol) :

    solvers = ['GMRES', 'GMRES+DSA(LU)', 'GMRES+DSA(ILU)', 
               'GMRES+CMDSA', 'GMRES+CMDSA+S', 'GMRES+TC-CMDSA+S', 'CMFD']

    out = \
    "\\begin{table}[ht]\n" + \
    "\\begin{center}\n" + \
    "\\begin{tabular*}{0.7\\textwidth}{@{\\extracolsep{\\fill}} rccc}\n" + \
    "\\toprule\n" + \
    "                        & sweeps & sweeps/groups & time (s) \\\\ \n" + \
    "      \midrule \n" 

    for i in range(0, len(solvers)) :
        out += solvers[i] + " & " + "%4i" % (data['sweeps'][i]) + " & " + "%4i" % (data['sweeps'][i]/7) + \
               ' & ' + "%8.1f" % (data['times'][i]) + ' \\\\ \n'

    out += \
    "\\bottomrule \n" + \
    "\\end{tabular*} \n" + \
    "\\end{center} \n" + \
    "\\caption{Sweeps and timings for a C5G7 assembly with $\\tau = $" + "%8.4e" % (tol) + " } \n" + \
    "\\end{table} \n" 

    print out


if __name__ == "__main__":            
    tol = 1.0e-6
    #run_reference() 
    #flux_norm_vs_pin_power()
    #plot_flux_norm_vs_pin_power()
    #data = run_convergence(tol)
    #print iters
    data = pickle.load(open('convergence_data.p', 'rb'))
    #print data
    #run(data["tols"], tol)
    data = pickle.load(open('Xdata_1e-06_8.p', 'rb'))
    #print data
    print_table(data, tol)
