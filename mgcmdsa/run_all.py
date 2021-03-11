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

from detran import *
import sys
sys.path.append('../c5g7')
from assemblies_c5g7 import get_assemblies
from pins_c5g7 import get_pins
from material_c5g7 import get_materials
import time
import cPickle as pickle
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams.update({'figure.autolayout': True})
q = 8 # number of polar/azimuth per octant

# Iterations for 1e-3 on fission rates


def get_input() :
    inp = utilities.InputDB.Create()
    inp.put_str("equation", "dd")
    inp.put_str("problem_type", "eigenvalue")
    inp.put_int("number_groups", 7)
    inp.put_str("inner_solver", "SI")
    inp.put_int("inner_max_iters", 1000000)
    inp.put_dbl("inner_tolerance", 1.0e-16)
    inp.put_int("inner_print_level", 0)
    inp.put_str("outer_solver", "GS")
    inp.put_int("outer_max_iters", 1000000)
    inp.put_int("outer_print_level", 2)
    inp.put_int("outer_print_interval", 1)
    inp.put_int("mgdsa_disable_fission", 0)
    inp.put_str("bc_west", "fixed")
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
    preconditioner_db.put_str("linear_solver_type",              "petsc");
    preconditioner_db.put_int("linear_solver_maxit",             5000);
    preconditioner_db.put_int("linear_solver_gmres_restart",     30);
    preconditioner_db.put_int("linear_solver_monitor_level",     0);
    preconditioner_db.put_str("pc_type",                         "ilu0");
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
    solver.set_solver()
    set_boundary(solver)
    return solver, mat, mesh

def set_boundary(solver) :
    b = cast_sn(solver.boundary())
    bc = cast_fixed(b.bc(0))
    quad = solver.quadrature()
    g = 0
    for o in range(0, 2) :
        for a in range(0, quad.number_angles_octant()) :
            for j in range(0, solver.mesh().number_cells_y()) :
                bc[o, a, g, j] = 1.0

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
    else :
        error("invalid outer")
    inp.put_dbl("outer_tolerance", tol);
    inp.get_spdb("outer_solver_db").put_dbl("linear_solver_rtol", tol);
    inp.get_spdb("outer_pc_db").put_dbl("linear_solver_rtol",     tol);
    inp.get_spdb("inner_pc_db").put_dbl("linear_solver_rtol",     tol);

    # SETUP INNER PC
    # Note, in the paper, we use PETSc and the ILU preconditioner for the
    # diffusion solvers.  Here, we've got only ILU0.
    if type[2] == 0 :
        inp.put_str("inner_pc_type", "none")
    elif type[2] == 1 :
        inp.put_str("inner_pc_type", "DSA")
    else :
        error("invalid inner pc")

    # SETUP OUTER PC
    # Same, only ILU0 with built-in GMRES.
    if type[3] == 0 :
        inp.put_str("outer_pc_type", "none")
    elif type[3] == 1 :
        inp.put_str("outer_pc_type", "mgdsa")
    elif type[3] == 2 :
        inp.put_str("outer_pc_type", "mgcmdsa")
        inp.put_int('mgpc_coarse_mesh_level',     7)
        inp.put_int('mgpc_condensation_option',   1)
        inp.put_int('mgpc_cmdsa_use_smoothing',   0)
    elif type[3] == 3 :
        inp.put_str("outer_pc_type", "mgcmdsa")
        inp.put_int('mgpc_coarse_mesh_level',     7)
        inp.put_int('mgpc_condensation_option',   1)
        inp.put_int('mgpc_cmdsa_use_smoothing',   1)
        inp.put_int('mgpc_cmdsa_smoothing_iters', 3)
        inp.put_dbl('mgpc_cmdsa_smoothing_relax', 0.7)
    elif type[3] == 5 :
        inp.put_str("outer_pc_type", "mgtcdsa")
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
        inp.put_str("outer_pc_type", "mgtcdsa")
        inp.put_int('mgpc_coarse_mesh_level',               7)
        inp.put_int('mgpc_condensation_option',             1)
        inp.put_int('mgpc_cmdsa_use_smoothing',             1)
        inp.put_int('mgpc_cmdsa_smoothing_iters',           3)
        inp.put_dbl('mgpc_cmdsa_smoothing_relax',           0.7)
        inp.put_int('mgpc_tcdsa_number_coarse_corrections', 1)
        inp.put_int('mgpc_tcdsa_use_fine_correction',       0)
        inp.put_int('mgpc_number_azimuth_octant',           2)
        inp.put_int('mgpc_number_polar_octant',             2)
    else :
        error("invalid outer pc")

def run_reference(keff) :
    """ Compute reference solution for determining convergence criteria
        of all solver schemes.  Uses unpreconditioned MG-GMRES
    """
    solver, mat, mesh = get_solver([0, 1, 0, 0], 1.0e-12)
    #solver, mat, mesh = get_solver([0, 0, 0, 0], 0.001)

    t = time.time()
    solver.solve(keff)
    print " reference elapsed = ", time.time()-t
    state = solver.state()
    rates = ReactionRates(mat, mesh, state)
    reference_data = {}
    reference_data['fission'] = np.asarray(rates.region_power(("PINS"), 1.0))
    for g in range(0, 7) :
        reference_data['phi'+str(g)] = np.asarray(state.phi(g))
    pickle.dump(reference_data, open('reference_data_'+str(keff)+'_keff.p', 'wb'))


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

def flux_norm_vs_pin_power(keff) :
    """ Compute a comparison of flux norm vs pin power.
    """
    data = pickle.load(open('reference_data_'+str(keff)+'.p', 'rb'))
    n = 12
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

    data2 = {}
    data2['tols']    = tols
    data2['nrm_phi'] = nrm_phi
    data2['nrm_fis'] = nrm_fis
    pickle.dump(data2, open('flux_norm_vs_pin_power_'+str(keff)+'.p', 'wb'))

def plot_flux_norm_vs_pin_power(keff) :

    data2 = pickle.load(open('flux_norm_vs_pin_power_'+str(keff)+'.p', 'rb'))

    tols = data2['tols']
    nrm_phi = data2['nrm_phi']
    nrm_fis = data2['nrm_fis']

    plt.figure(1, figsize=(8., 6.))
    plt.loglog(tols, nrm_phi, 'k-o', tols, nrm_fis, 'b-s',markersize=10)
    plt.legend(['$  \\frac{ || \phi - \phi_{\mathrm{ref}} ||_2 }{ || \phi_{\mathrm{ref}} ||_2 }$',
                '$  \mathrm{max}   \left | \\frac{d - d_{\mathrm{ref}} } { d_{\mathrm{ref}} } \\right | $'],loc=2,prop={'size':22},numpoints=1)
    plt.xlabel('$ \mathrm{max}| \phi_{n+1}-\phi_{n} | $ ', fontsize=22)
    plt.ylabel('residual norm', fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.savefig('flux_norm_vs_pin_power.pdf')

    plt.figure(2, figsize=(8., 6.))
    plt.semilogx(tols, nrm_fis / nrm_phi, 'k-o',markersize=10)
    plt.grid(True)
    plt.xlabel('$ \mathrm{max}| \phi_{n+1}-\phi_{n} | $ ',fontsize=22)
    plt.ylabel('ratio of norms', fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.show()
    plt.savefig('flux_norm_vs_pin_power_ratio.pdf')

def run_convergence2(re, keff) :
    # Load reference data
    data = pickle.load(open('reference_data_'+str(keff)+'_keff.p', 'rb'))
    mask = data['fission'] > 0.0
    ref_fiss = data['fission'][mask]
    sweeps = []
    tols = {}
    # Define cases  7.72044953288e-06 4.75762866421e-07
    cases = [[0,0,0,0], [1,0,0,0], [2,0,0,0], [0,1,0,0]]
    tols0 = [2.18466246478e-05, 2.18466246478e-05, 0.0005, 0.001]
    for n in range(0, 4) :
        tol = tols0[n]
        iter = 0
        while iter < 50 :
            print "iter = ", iter, tol
            solver, mat, mesh = get_solver(cases[n], tol)
            solver.solve(keff)
            re1, tmp = flux_norm(solver.state(), mat, mesh, data)
            print "re =  %12.8e, its=  %4i" % (re1, tol)
            if (re1 < re) :
                print "yeah!", re-re1, re, tol
                tols[cases[n][0], cases[n][1]] = tol
                done = True
                break
            tol = tol * 0.707
            iter += 1
        sweeps.append(solver.number_sweeps())
    print sweeps
    # 7.51820723957e-05
    data = {}
    data['tols']   = tols
    data['sweeps'] = sweeps
    pickle.dump(data, open('convergence_data_'+str(re)+'_'+str(keff)+'_keff.p', 'wb'))
    return data

def run_convergence(re, keff) :
    """ Decide unpreconditioned solver iterations yielding desired
        accuracy.  For GS+SI, we use the same iteration count for
        both GS and SI.
    """

    # Load reference data
    data = pickle.load(open('reference_data_'+str(keff)+'_keff.p', 'rb'))
    mask = data['fission'] > 0.0
    ref_fiss = data['fission'][mask]

    # Define cases
    cases = [[0,0,0,0], [1,0,0,0], [2,0,0,0], [0,1,0,0]]
    tols0 = [9.50309158e-04,  0.01, 0.01, 0.01]
    tols1 = [0.001,  0.001, 0.001, 0.001]
    tols = {}
    sweeps = []
    toler = []
    rerr = []
    for n in range(0, 4) :
        print "CONV n =",n
        # first guess
        tol0 = tols0[n]
        print cases[n]
        solver, mat, mesh = get_solver(cases[n], tol0)
        solver.solve(keff)
        re0, tmp = flux_norm(solver.state(), mat, mesh, data)
        print "RE0 = ", re0
        # second guess
        tol1 = tols1[n]
        solver, mat, mesh = get_solver(cases[n], tol1)
        solver.solve()
        re1, tmp = flux_norm(solver.state(), mat, mesh, data)
        print "RE1 = ", re1
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
            solver.solve(keff)
            re1, tmp = flux_norm(solver.state(), mat, mesh, data)
            print "re =  %12.8e, its=  %4i" % (re1, tol1)



            if (tol0 == tol1 or re1 < re) :
                print "yeah!", re-re1, re, tol1
                tols[cases[n][0], cases[n][1]] = tol1
                done = True
                break
            iter += 1
        assert(done)
        sweeps.append(solver.number_sweeps())
        toler.append(tol1)
        rerr.append(re1)
    print sweeps
    # 7.51820723957e-05
    data = {}
    data['tols']   = tols
    data['sweeps'] = sweeps
    data['toler'] = toler
    data['rerr']  = rerr
    print data
    pickle.dump(data, open('convergence_data_'+str(re)+'_'+str(keff)+'_keff.p', 'wb'))
    return data

def run(iters, tol, keff) :

    refdata = pickle.load(open('reference_data_'+str(keff)+'.p', 'rb'))

    cases = []
    # add GS+SI cases
    for i in range(0, 2) :
        cases.append([i,0,0,0])
    # add GS+GMRES cases
    for i in range(0, 3) :
        cases.append([2,0,i,0])
    # add MG-GMRES cases
    for i in range(0, 7) :
        cases.append([0,1,0,i])
    num = len(cases)

    data = {}
    data["times"]  = np.zeros(num)
    data["sweeps"] = np.zeros(num)
    data["tol"]    = np.zeros(num)
    data["error"]  = np.zeros(num)

    for n in range(0, num) :
        case = cases[n]
        print "case:",n , case
        solver, mat, mesh = get_solver(case, iters[(case[0],case[1])])
        t = time.time()
        solver.solve(keff)
        data["times"][n]  = time.time() - t
        data["sweeps"][n] = solver.number_sweeps()
        data["tol"][n]    = iters[(case[0],case[1])]
        data["error"][n], tmp = flux_norm(solver.state(), mat, mesh, refdata)
        print data["times"][n]
        print data["sweeps"][n]
        print data["tol"][n]
        print data["error"][n]

    f = "data_" + str(tol) + "_" + str(keff) + "_desktop.p"
    pickle.dump(data, open(f, 'wb'))

def print_table(data, tol) :

    solvers = ['GS, SI(20)', 'GS, SI(1)',
               'GS, GMRES', 'GS, GMRES+DSA(LU)', 'GS, GMRES+DSA(ILU)',
               'GMRES', 'GMRES+DSA(LU)', 'GMRES+DSA(ILU)',
               'GMRES+CMDSA', 'GMRES+CMDSA+S', 'GMRES+TC-CMDSA+S-f', 'GMRES+TC-CMDSA+S']

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
    "\\caption{Sweeps and timings for a C5G7 assembly with $\\tau =  " + "%8.4e" % (tol) + " $ } \n" + \
    "\\end{table} \n"

    print out

def plot_stuff() :
    """ Plot the number of iterations and the time as a function of k """


    kv = 0.534164595
    keffs = [kv+1e-2, kv+1e-1, 0.7, 1.0, 1.3]
    tol = 1.0e-3



if __name__ == "__main__":

    # parameters of the study
    kv = 0.534164595
    keffs = [kv+1e-8,kv+1e-7,kv+1e-6,kv+1e-5,kv+1e-4,kv+1e-3,kv+1e-2, kv+1e-1, 0.7, 1.0, 1.3]



    tols = [1.0e-4, 1.0e-6, 1.0e-8]

    #plot_flux_norm_vs_pin_power(1.0)

    for keff in [0.7] :

        # compute the reference for this keff
        #run_reference(keff)

        for tol in [1.0e-4] : #tols :
            #pass
            # compute and plot norm comparison
            #flux_norm_vs_pin_power(keff)
            #plot_flux_norm_vs_pin_power(keff)

            # compute the tolerances needed to achieve the error tolerance
            data = run_convergence2(tol, keff)

            #data = pickle.load(open('convergence_data_'+str(tol)+'_'+str(keff)+'c.p', 'rb'))

            #print data
            #run(data["tols"], tol, keff)
            #data = pickle.load(open("data_" + str(tol) + "_" + str(keff) + ".p", 'rb'))
            #print data
            #print_table(data, tol)
