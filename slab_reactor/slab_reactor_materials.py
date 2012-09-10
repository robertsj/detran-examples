# pyexamples/slab_reactor/slab_reactor_materials.py
#
# This implements the materials for the test cases published 
# in Scott Mosher's Ph.D. thesis, "A Variational Coarse Mesh 
# Transport Method". All data and reference values are from 
# Appendix A of that work.

from detran import *

def get_materials() :
    # Two-group data from 1-d coarse mesh benchmarks (Mosher, Ilas etc.)
    
    # Create the Materials object.
    mat = Material.Create(4, 2, False);

    # ---------------------------
    # Material 0: Water           
    # ---------------------------

    # Total
    mat.set_sigma_t(0, 0, 0.1890);       # (obj, matid, g, value);
    mat.set_sigma_t(0, 1, 1.4633);       

    # Fission 
    mat.set_sigma_f(0, 0, 0.0);         # Note, default is zero
    mat.set_sigma_f(0, 1, 0.0);   
    mat.set_chi(0, 0, 0.0); 
    mat.set_chi(0, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(0, 0, 0, 0.1507);    # 1 <- 1
    mat.set_sigma_s(0, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(0, 1, 0, 0.0380);    # 2 <- 1
    mat.set_sigma_s(0, 1, 1, 1.4536);    # 2 <- 2

    # ---------------------------
    # Material 1: Fuel I           
    # ---------------------------

    # Total
    mat.set_sigma_t(1, 0, 0.2263);       # (obj, matid, g, value);
    mat.set_sigma_t(1, 1, 1.0119);       

    # Fission 
    mat.set_sigma_f(1, 0, 0.0067);
    mat.set_sigma_f(1, 1, 0.1241);   
    mat.set_chi(1, 0, 1.0); 
    mat.set_chi(1, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(1, 0, 0, 0.2006);    # 1 <- 1
    mat.set_sigma_s(1, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(1, 1, 0, 0.0161);    # 2 <- 1
    mat.set_sigma_s(1, 1, 1, 0.9355);    # 2 <- 2

    # ---------------------------
    # Material 3: Fuel II          
    # ---------------------------

    # Total
    mat.set_sigma_t(2, 0, 0.2252);       # (obj, matid, g, value);
    mat.set_sigma_t(2, 1, 0.9915);       

    # Fission 
    mat.set_sigma_f(2, 0, 0.0078);
    mat.set_sigma_f(2, 1, 0.1542);   
    mat.set_chi(2, 0, 1.0); 
    mat.set_chi(2, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(2, 0, 0, 0.1995);    # 1 <- 1
    mat.set_sigma_s(2, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(2, 1, 0, 0.0156);    # 2 <- 1
    mat.set_sigma_s(2, 1, 1, 0.9014);    # 2 <- 2

    # ---------------------------
    # Material 4: Fuel II + Gd          
    # ---------------------------
    
    # Total
    mat.set_sigma_t(3, 0, 0.2173);       # (obj, matid, g, value);
    mat.set_sigma_t(3, 1, 1.0606);       

    # Fission 
    mat.set_sigma_f(3, 0, 0.0056);
    mat.set_sigma_f(3, 1, 0.0187);   
    mat.set_chi(3, 0, 1.0); 
    mat.set_chi(3, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(3, 0, 0, 0.1902);	   # 1 <- 1
    mat.set_sigma_s(3, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(3, 1, 0, 0.0136);    # 2 <- 1
    mat.set_sigma_s(3, 1, 1, 0.5733);    # 2 <- 2

    # ---------------------------
    # FINALIZE     
    # ---------------------------

    mat.finalize();

    return mat

