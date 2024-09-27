"""
Materials for 1-D Transient Slab Problem
"""


import numpy as np 
from detran import KineticsMaterial

def compute_sigma(D, A, S_out, F):
    T = 1.0/(3.0*D)
    S_in = T - A - S_out 
    return D, T, A, S_in, S_out, F 

def fill_water(m, mat):

    D, T, A, S_in, S_out, F = compute_sigma(D=1.5, A=0.0002, S_out=0.032, F=0.0)
    mat.set_diff_coef(m, 0,   D)
    mat.set_sigma_t(m, 0,     T) 
    mat.set_sigma_a(m, 0,     A)
    mat.set_sigma_s(m, 0, 0,  S_in)
    mat.set_sigma_s(m, 1, 0,  S_out)
    mat.set_sigma_f(m, 0,     F)

    D, T, A, S_in, S_out, F = compute_sigma(D=0.5, A=0.01, S_out=0.0, F=0.0)
    mat.set_diff_coef(m, 1,   D)
    mat.set_sigma_t(m, 1,     T) 
    mat.set_sigma_a(m, 1,     A) 
    mat.set_sigma_s(m, 0, 1,  S_out)
    mat.set_sigma_s(m, 1, 1,  S_in)
    mat.set_sigma_f(m, 1,     F)
    
def fill_fuel(m, mat, withdrawn_fraction):

    D, T, A, S_in, S_out, F = compute_sigma(D=1.3, A=0.0105, S_out=0.022, F=0.003)
    mat.set_diff_coef(m, 0,   D)
    mat.set_sigma_t(m, 0,     T) 
    mat.set_sigma_a(m, 0,     A)
    mat.set_sigma_s(m, 0, 0,  S_in)
    mat.set_sigma_s(m, 1, 0,  S_out)
    mat.set_sigma_f(m, 0,     F)
    mat.set_chi(m, 0,         1.0)
    mat.set_nu(m,  0,         1.0)

    D, T, A, S_in, S_out, F = compute_sigma(D=0.5, A=0.164-0.05*withdrawn_fraction, S_out=0.0, F=0.19)
    mat.set_diff_coef(m, 1,   D)
    mat.set_sigma_t(m, 1,     T) 
    mat.set_sigma_a(m, 1,     A) 
    mat.set_sigma_s(m, 0, 1,  S_out)
    mat.set_sigma_s(m, 1, 1,  S_in)
    mat.set_sigma_f(m, 1,     F)
    mat.set_chi(m, 1,         0.0)
    mat.set_nu(m,  1,         1.0)

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

def get_materials(fraction_A2, fraction_A4):
    assert 0.0 <= fraction_A2 <= 1.0, "fraction out of range!"
    assert 0.0 <= fraction_A4 <= 1.0, "fraction out of range!"
    mat = KineticsMaterial(6, 2, 8, "SLAB_MATERIAL")
    fill_water(0, mat)
    fill_fuel(1, mat, withdrawn_fraction=1.00)  
    fill_fuel(2, mat, withdrawn_fraction=fraction_A2)  
    fill_fuel(3, mat, withdrawn_fraction=1.00) 
    fill_fuel(4, mat, withdrawn_fraction=fraction_A4)
    fill_fuel(5, mat, withdrawn_fraction=1.00)
    fill_kinetics(mat)
    mat.finalize()
    return mat


def get_rods_in_materials():
    mat = KineticsMaterial(6, 2, 8, "RODS_IN")
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
    mat = KineticsMaterial(6, 2, 8, "RODS_OUT")
    fill_water(0, mat)
    fill_fuel(1, mat, withdrawn_fraction=1.00)  
    fill_fuel(2, mat, withdrawn_fraction=0.25)  
    fill_fuel(3, mat, withdrawn_fraction=1.00) 
    fill_fuel(4, mat, withdrawn_fraction=0.30)
    fill_fuel(5, mat, withdrawn_fraction=1.00)
    fill_kinetics(mat)
    mat.finalize()
    return mat



def make_material_table():
    mat = KineticsMaterial(3, 2, 8, "Material Tables")
    fill_water(0, mat)
    fill_fuel(1, mat, 1.0) # all out
    fill_fuel(2, mat, 0.0) # all in
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
        A1 = mat.sigma_a(m, 0) 
        A2 = mat.sigma_a(m, 1)
        S21 = mat.sigma_s(m, 1, 0)
        nF1 = mat.nu_sigma_f(m, 0)
        nF2 = mat.nu_sigma_f(m, 1)
        s += tmpl.format(names[m], D0, D1, A1, A2, S21, nF1, nF2) + '\n'
        s += r"\hline" + '\n'
    f = open('cross_section_table.txt', 'w')
    f.write(s)
    f.close()

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
    f = open('kinetics_table.txt', 'w')
    f.write(s)
    f.close()

    mat.display()
    

if __name__ == '__main__':

    make_material_table()