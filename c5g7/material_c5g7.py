# pyexamples/c5g7/material_c5g7.py
#
# C5G7 materials

from detran import Material

def get_materials() :

  # Create the new database.
  mat = Material.Create(7, 7, False)

  # --------------------------------------------
  # Material 0: UO2 fuel-clad
  # --------------------------------------------
  m = 0
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.77949E-01)
  mat.set_sigma_t(m, 1, 3.29805E-01)
  mat.set_sigma_t(m, 2, 4.80388E-01)
  mat.set_sigma_t(m, 3, 5.54367E-01)
  mat.set_sigma_t(m, 4, 3.11801E-01)
  mat.set_sigma_t(m, 5, 3.95168E-01)
  mat.set_sigma_t(m, 6, 5.64406E-01)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 8.02480E-03)
  # mat.set_sigma_a(m, 1, 3.71740E-03)
  # mat.set_sigma_a(m, 2, 2.67690E-02)
  # mat.set_sigma_a(m, 3, 9.62360E-02)
  # mat.set_sigma_a(m, 4, 3.00200E-02)
  # mat.set_sigma_a(m, 5, 1.11260E-01)
  # mat.set_sigma_a(m, 6, 2.82780E-01)
  # Fission times nu
  mat.set_sigma_f(m, 0, 7.21206E-03)
  mat.set_sigma_f(m, 1, 8.19301E-04)
  mat.set_sigma_f(m, 2, 6.45320E-03)
  mat.set_sigma_f(m, 3, 1.85648E-02)
  mat.set_sigma_f(m, 4, 1.78084E-02)
  mat.set_sigma_f(m, 5, 8.30348E-02)
  mat.set_sigma_f(m, 6, 2.16004E-01)
  mat.set_nu(m, 0, 2.78145E+00)
  mat.set_nu(m, 1, 2.47443E+00)
  mat.set_nu(m, 2, 2.43383E+00)
  mat.set_nu(m, 3, 2.43380E+00)
  mat.set_nu(m, 4, 2.43380E+00)
  mat.set_nu(m, 5, 2.43380E+00)
  mat.set_nu(m, 6, 2.43380E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 1.27537E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 4.23780E-02)
  mat.set_sigma_s(m, 1, 1, 3.24456E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 9.43740E-06)
  mat.set_sigma_s(m, 2, 1, 1.63140E-03)
  mat.set_sigma_s(m, 2, 2, 4.50940E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 5.51630E-09)
  mat.set_sigma_s(m, 3, 1, 3.14270E-09)
  mat.set_sigma_s(m, 3, 2, 2.67920E-03)
  mat.set_sigma_s(m, 3, 3, 4.52565E-01)
  mat.set_sigma_s(m, 3, 4, 1.25250E-04)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 3, 5.56640E-03)
  mat.set_sigma_s(m, 4, 4, 2.71401E-01)
  mat.set_sigma_s(m, 4, 5, 1.29680E-03)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 4, 1.02550E-02)
  mat.set_sigma_s(m, 5, 5, 2.65802E-01)
  mat.set_sigma_s(m, 5, 6, 8.54580E-03)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 4, 1.00210E-08)
  mat.set_sigma_s(m, 6, 5, 1.68090E-02)
  mat.set_sigma_s(m, 6, 6, 2.73080E-01)

  # --------------------------------------------
  # Material 1: 4.3w/o MOX fuel-clad
  # --------------------------------------------
  m = 1
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.78731E-01)
  mat.set_sigma_t(m, 1, 3.30849E-01)
  mat.set_sigma_t(m, 2, 4.83772E-01)
  mat.set_sigma_t(m, 3, 5.66922E-01)
  mat.set_sigma_t(m, 4, 4.26227E-01)
  mat.set_sigma_t(m, 5, 6.78997E-01)
  mat.set_sigma_t(m, 6, 6.82852E-01)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 8.43390E-03)
  # mat.set_sigma_a(m, 1, 3.75770E-03)
  # mat.set_sigma_a(m, 2, 2.79700E-02)
  # mat.set_sigma_a(m, 3, 1.04210E-01)
  # mat.set_sigma_a(m, 4, 1.39940E-01)
  # mat.set_sigma_a(m, 5, 4.09180E-01)
  # mat.set_sigma_a(m, 6, 4.09350E-01)
  # Fission times nu
  mat.set_sigma_f(m, 0, 7.62704E-03)
  mat.set_sigma_f(m, 1, 8.76898E-04)
  mat.set_sigma_f(m, 2, 5.69835E-03)
  mat.set_sigma_f(m, 3, 2.28872E-02)
  mat.set_sigma_f(m, 4, 1.07635E-02)
  mat.set_sigma_f(m, 5, 2.32757E-01)
  mat.set_sigma_f(m, 6, 2.48968E-01)
  mat.set_nu(m, 0, 2.85209E+00)
  mat.set_nu(m, 1, 2.89099E+00)
  mat.set_nu(m, 2, 2.85486E+00)
  mat.set_nu(m, 3, 2.86073E+00)
  mat.set_nu(m, 4, 2.85447E+00)
  mat.set_nu(m, 5, 2.86415E+00)
  mat.set_nu(m, 6, 2.86780E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 1.28876E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 4.14130E-02)
  mat.set_sigma_s(m, 1, 1, 3.25452E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 8.22900E-06)
  mat.set_sigma_s(m, 2, 1, 1.63950E-03)
  mat.set_sigma_s(m, 2, 2, 4.53188E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 5.04050E-09)
  mat.set_sigma_s(m, 3, 1, 1.59820E-09)
  mat.set_sigma_s(m, 3, 2, 2.61420E-03)
  mat.set_sigma_s(m, 3, 3, 4.57173E-01)
  mat.set_sigma_s(m, 3, 4, 1.60460E-04)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 3, 5.53940E-03)
  mat.set_sigma_s(m, 4, 4, 2.76814E-01)
  mat.set_sigma_s(m, 4, 5, 2.00510E-03)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 4, 9.31270E-03)
  mat.set_sigma_s(m, 5, 5, 2.52962E-01)
  mat.set_sigma_s(m, 5, 6, 8.49480E-03)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 4, 9.16560E-09)
  mat.set_sigma_s(m, 6, 5, 1.48500E-02)
  mat.set_sigma_s(m, 6, 6, 2.65007E-01)

  # --------------------------------------------
  # Material 2: 7.0 w/o MOX fuel-clad
  # --------------------------------------------
  m = 2
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.81323E-01)
  mat.set_sigma_t(m, 1, 3.34368E-01)
  mat.set_sigma_t(m, 2, 4.93785E-01)
  mat.set_sigma_t(m, 3, 5.91216E-01)
  mat.set_sigma_t(m, 4, 4.74198E-01)
  mat.set_sigma_t(m, 5, 8.33601E-01)
  mat.set_sigma_t(m, 6, 8.53603E-01)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 9.06570E-03)
  # mat.set_sigma_a(m, 1, 4.29670E-03)
  # mat.set_sigma_a(m, 2, 3.28810E-02)
  # mat.set_sigma_a(m, 3, 1.22030E-01)
  # mat.set_sigma_a(m, 4, 1.82980E-01)
  # mat.set_sigma_a(m, 5, 5.68460E-01)
  # mat.set_sigma_a(m, 6, 5.85210E-01)
  # Fission times nu
  mat.set_sigma_f(m, 0, 8.25446E-03)
  mat.set_sigma_f(m, 1, 1.32565E-03)
  mat.set_sigma_f(m, 2, 8.42156E-03)
  mat.set_sigma_f(m, 3, 3.28730E-02)
  mat.set_sigma_f(m, 4, 1.59636E-02)
  mat.set_sigma_f(m, 5, 3.23794E-01)
  mat.set_sigma_f(m, 6, 3.62803E-01)
  mat.set_nu(m, 0, 2.88498E+00)
  mat.set_nu(m, 1, 2.91079E+00)
  mat.set_nu(m, 2, 2.86574E+00)
  mat.set_nu(m, 3, 2.87063E+00)
  mat.set_nu(m, 4, 2.86714E+00)
  mat.set_nu(m, 5, 2.86658E+00)
  mat.set_nu(m, 6, 2.87539E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 1.30457E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 4.17920E-02)
  mat.set_sigma_s(m, 1, 1, 3.28428E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 8.51050E-06)
  mat.set_sigma_s(m, 2, 1, 1.64360E-03)
  mat.set_sigma_s(m, 2, 2, 4.58371E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 5.13290E-09)
  mat.set_sigma_s(m, 3, 1, 2.20170E-09)
  mat.set_sigma_s(m, 3, 2, 2.53310E-03)
  mat.set_sigma_s(m, 3, 3, 4.63709E-01)
  mat.set_sigma_s(m, 3, 4, 1.76190E-04)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 3, 5.47660E-03)
  mat.set_sigma_s(m, 4, 4, 2.82313E-01)
  mat.set_sigma_s(m, 4, 5, 2.27600E-03)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 4, 8.72890E-03)
  mat.set_sigma_s(m, 5, 5, 2.49751E-01)
  mat.set_sigma_s(m, 5, 6, 8.86450E-03)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 4, 9.00160E-09)
  mat.set_sigma_s(m, 6, 5, 1.31140E-02)
  mat.set_sigma_s(m, 6, 6, 2.59529E-01)

  # --------------------------------------------
  # Material 3: 8.7 w/o MOX fuel-clad
  # --------------------------------------------
  m = 3
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.83045E-01)
  mat.set_sigma_t(m, 1, 3.36705E-01)
  mat.set_sigma_t(m, 2, 5.00507E-01)
  mat.set_sigma_t(m, 3, 6.06174E-01)
  mat.set_sigma_t(m, 4, 5.02754E-01)
  mat.set_sigma_t(m, 5, 9.21028E-01)
  mat.set_sigma_t(m, 6, 9.55231E-01)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 9.48620E-03)
  # mat.set_sigma_a(m, 1, 4.65560E-03)
  # mat.set_sigma_a(m, 2, 3.62400E-02)
  # mat.set_sigma_a(m, 3, 1.32720E-01)
  # mat.set_sigma_a(m, 4, 2.08400E-01)
  # mat.set_sigma_a(m, 5, 6.58700E-01)
  # mat.set_sigma_a(m, 6, 6.90170E-01)
  # Fission times nu
  mat.set_sigma_f(m, 0, 8.67209E-03)
  mat.set_sigma_f(m, 1, 1.62426E-03)
  mat.set_sigma_f(m, 2, 1.02716E-02)
  mat.set_sigma_f(m, 3, 3.90447E-02)
  mat.set_sigma_f(m, 4, 1.92576E-02)
  mat.set_sigma_f(m, 5, 3.74888E-01)
  mat.set_sigma_f(m, 6, 4.30599E-01)
  mat.set_nu(m, 0, 2.90426E+00)
  mat.set_nu(m, 1, 2.91795E+00)
  mat.set_nu(m, 2, 2.86986E+00)
  mat.set_nu(m, 3, 2.87491E+00)
  mat.set_nu(m, 4, 2.87175E+00)
  mat.set_nu(m, 5, 2.86752E+00)
  mat.set_nu(m, 6, 2.87808E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 1.31504E-01)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 4.20460E-02)
  mat.set_sigma_s(m, 1, 1, 3.30403E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 8.69720E-06)
  mat.set_sigma_s(m, 2, 1, 1.64630E-03)
  mat.set_sigma_s(m, 2, 2, 4.61792E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 5.19380E-09)
  mat.set_sigma_s(m, 3, 1, 2.60060E-09)
  mat.set_sigma_s(m, 3, 2, 2.47490E-03)
  mat.set_sigma_s(m, 3, 3, 4.68021E-01)
  mat.set_sigma_s(m, 3, 4, 1.85970E-04)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 3, 5.43300E-03)
  mat.set_sigma_s(m, 4, 4, 2.85771E-01)
  mat.set_sigma_s(m, 4, 5, 2.39160E-03)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 4, 8.39730E-03)
  mat.set_sigma_s(m, 5, 5, 2.47614E-01)
  mat.set_sigma_s(m, 5, 6, 8.96810E-03)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 4, 8.92800E-09)
  mat.set_sigma_s(m, 6, 5, 1.23220E-02)
  mat.set_sigma_s(m, 6, 6, 2.56093E-01)

  # --------------------------------------------
  # Material 4: fission chamber
  # --------------------------------------------
  m = 4
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.26032E-01)
  mat.set_sigma_t(m, 1, 2.93160E-01)
  mat.set_sigma_t(m, 2, 2.84250E-01)
  mat.set_sigma_t(m, 3, 2.81020E-01)
  mat.set_sigma_t(m, 4, 3.34460E-01)
  mat.set_sigma_t(m, 5, 5.65640E-01)
  mat.set_sigma_t(m, 6, 1.17214E+00)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 5.11320E-04)
  # mat.set_sigma_a(m, 1, 7.58130E-05)
  # mat.set_sigma_a(m, 2, 3.16430E-04)
  # mat.set_sigma_a(m, 3, 1.16750E-03)
  # mat.set_sigma_a(m, 4, 3.39770E-03)
  # mat.set_sigma_a(m, 5, 9.18860E-03)
  # mat.set_sigma_a(m, 6, 2.32440E-02)
  # Fission times nu
  mat.set_sigma_f(m, 0, 4.79002E-09)
  mat.set_sigma_f(m, 1, 5.82564E-09)
  mat.set_sigma_f(m, 2, 4.63719E-07)
  mat.set_sigma_f(m, 3, 5.24406E-06)
  mat.set_sigma_f(m, 4, 1.45390E-07)
  mat.set_sigma_f(m, 5, 7.14972E-07)
  mat.set_sigma_f(m, 6, 2.08041E-06)
  mat.set_nu(m, 0, 2.76283E+00)
  mat.set_nu(m, 1, 2.46239E+00)
  mat.set_nu(m, 2, 2.43380E+00)
  mat.set_nu(m, 3, 2.43380E+00)
  mat.set_nu(m, 4, 2.43380E+00)
  mat.set_nu(m, 5, 2.43380E+00)
  mat.set_nu(m, 6, 2.43380E+00)
  # Fission spectrum
  mat.set_chi(m, 0, 5.87819E-01)
  mat.set_chi(m, 1, 4.11760E-01)
  mat.set_chi(m, 2, 3.39060E-04)
  mat.set_chi(m, 3, 1.17610E-07)
  mat.set_chi(m, 4, 0.00000E+00)
  mat.set_chi(m, 5, 0.00000E+00)
  mat.set_chi(m, 6, 0.00000E+00)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 6.61659E-02)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 5.90700E-02)
  mat.set_sigma_s(m, 1, 1, 2.40377E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 2.83340E-04)
  mat.set_sigma_s(m, 2, 1, 5.24350E-02)
  mat.set_sigma_s(m, 2, 2, 1.83425E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 1.46220E-06)
  mat.set_sigma_s(m, 3, 1, 2.49900E-04)
  mat.set_sigma_s(m, 3, 2, 9.22880E-02)
  mat.set_sigma_s(m, 3, 3, 7.90769E-02)
  mat.set_sigma_s(m, 3, 4, 3.73400E-05)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 0, 2.06420E-08)
  mat.set_sigma_s(m, 4, 1, 1.92390E-05)
  mat.set_sigma_s(m, 4, 2, 6.93650E-03)
  mat.set_sigma_s(m, 4, 3, 1.69990E-01)
  mat.set_sigma_s(m, 4, 4, 9.97570E-02)
  mat.set_sigma_s(m, 4, 5, 9.17420E-04)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 1, 2.98750E-06)
  mat.set_sigma_s(m, 5, 2, 1.07900E-03)
  mat.set_sigma_s(m, 5, 3, 2.58600E-02)
  mat.set_sigma_s(m, 5, 4, 2.06790E-01)
  mat.set_sigma_s(m, 5, 5, 3.16774E-01)
  mat.set_sigma_s(m, 5, 6, 4.97930E-02)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 1, 4.21400E-07)
  mat.set_sigma_s(m, 6, 2, 2.05430E-04)
  mat.set_sigma_s(m, 6, 3, 4.92560E-03)
  mat.set_sigma_s(m, 6, 4, 2.44780E-02)
  mat.set_sigma_s(m, 6, 5, 2.38760E-01)
  mat.set_sigma_s(m, 6, 6, 1.09910E+00)

  # --------------------------------------------
  # Material 5: guide tube
  # --------------------------------------------
  m = 5
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.26032E-01)
  mat.set_sigma_t(m, 1, 2.93160E-01)
  mat.set_sigma_t(m, 2, 2.84240E-01)
  mat.set_sigma_t(m, 3, 2.80960E-01)
  mat.set_sigma_t(m, 4, 3.34440E-01)
  mat.set_sigma_t(m, 5, 5.65640E-01)
  mat.set_sigma_t(m, 6, 1.17215E+00)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 5.11320E-04)
  # mat.set_sigma_a(m, 1, 7.58010E-05)
  # mat.set_sigma_a(m, 2, 3.15720E-04)
  # mat.set_sigma_a(m, 3, 1.15820E-03)
  # mat.set_sigma_a(m, 4, 3.39750E-03)
  # mat.set_sigma_a(m, 5, 9.18780E-03)
  # mat.set_sigma_a(m, 6, 2.32420E-02)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 6.61659E-02)
  # 2 <- g'
  mat.set_sigma_s(m, 1, 0, 5.90700E-02)
  mat.set_sigma_s(m, 1, 1, 2.40377E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 2.83340E-04)
  mat.set_sigma_s(m, 2, 1, 5.24350E-02)
  mat.set_sigma_s(m, 2, 2, 1.83297E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 1.46220E-06)
  mat.set_sigma_s(m, 3, 1, 2.49900E-04)
  mat.set_sigma_s(m, 3, 2, 9.23970E-02)
  mat.set_sigma_s(m, 3, 3, 7.88511E-02)
  mat.set_sigma_s(m, 3, 4, 3.73330E-05)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 0, 2.06420E-08)
  mat.set_sigma_s(m, 4, 1, 1.92390E-05)
  mat.set_sigma_s(m, 4, 2, 6.94460E-03)
  mat.set_sigma_s(m, 4, 3, 1.70140E-01)
  mat.set_sigma_s(m, 4, 4, 9.97372E-02)
  mat.set_sigma_s(m, 4, 5, 9.17260E-04)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 1, 2.98750E-06)
  mat.set_sigma_s(m, 5, 2, 1.07900E-03)
  mat.set_sigma_s(m, 5, 3, 2.58600E-02)
  mat.set_sigma_s(m, 5, 4, 2.06790E-01)
  mat.set_sigma_s(m, 5, 5, 3.16774E-01)
  mat.set_sigma_s(m, 5, 6, 4.97930E-02)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 1, 4.21400E-07)
  mat.set_sigma_s(m, 6, 2, 2.05430E-04)
  mat.set_sigma_s(m, 6, 3, 4.92560E-03)
  mat.set_sigma_s(m, 6, 4, 2.44780E-02)
  mat.set_sigma_s(m, 6, 5, 2.38760E-01)
  mat.set_sigma_s(m, 6, 6, 1.09910E+00)

  # --------------------------------------------
  # Material 6: moderator
  # --------------------------------------------
  m = 6
  # Transport cross section
  mat.set_sigma_t(m, 0, 1.59206E-01)
  mat.set_sigma_t(m, 1, 4.12970E-01)
  mat.set_sigma_t(m, 2, 5.90310E-01)
  mat.set_sigma_t(m, 3, 5.84350E-01)
  mat.set_sigma_t(m, 4, 7.18000E-01)
  mat.set_sigma_t(m, 5, 1.25445E+00)
  mat.set_sigma_t(m, 6, 2.65038E+00)
  # Absorption cross section
  # mat.set_sigma_a(m, 0, 6.01050E-04)
  # mat.set_sigma_a(m, 1, 1.57930E-05)
  # mat.set_sigma_a(m, 2, 3.37160E-04)
  # mat.set_sigma_a(m, 3, 1.94060E-03)
  # mat.set_sigma_a(m, 4, 5.74160E-03)
  # mat.set_sigma_a(m, 5, 1.50010E-02)
  # mat.set_sigma_a(m, 6, 3.72390E-02)
  # Scattering
  # 1 <- g'
  mat.set_sigma_s(m, 0, 0, 4.44777E-02)
  # 2 <- g
  mat.set_sigma_s(m, 1, 0, 1.13400E-01)
  mat.set_sigma_s(m, 1, 1, 2.82334E-01)
  # 3 <- g'
  mat.set_sigma_s(m, 2, 0, 7.23470E-04)
  mat.set_sigma_s(m, 2, 1, 1.29940E-01)
  mat.set_sigma_s(m, 2, 2, 3.45256E-01)
  # 4 <- g'
  mat.set_sigma_s(m, 3, 0, 3.74990E-06)
  mat.set_sigma_s(m, 3, 1, 6.23400E-04)
  mat.set_sigma_s(m, 3, 2, 2.24570E-01)
  mat.set_sigma_s(m, 3, 3, 9.10284E-02)
  mat.set_sigma_s(m, 3, 4, 7.14370E-05)
  # 5 <- g'
  mat.set_sigma_s(m, 4, 0, 5.31840E-08)
  mat.set_sigma_s(m, 4, 1, 4.80020E-05)
  mat.set_sigma_s(m, 4, 2, 1.69990E-02)
  mat.set_sigma_s(m, 4, 3, 4.15510E-01)
  mat.set_sigma_s(m, 4, 4, 1.39138E-01)
  mat.set_sigma_s(m, 4, 5, 2.21570E-03)
  # 6 <- g'
  mat.set_sigma_s(m, 5, 1, 7.44860E-06)
  mat.set_sigma_s(m, 5, 2, 2.64430E-03)
  mat.set_sigma_s(m, 5, 3, 6.37320E-02)
  mat.set_sigma_s(m, 5, 4, 5.11820E-01)
  mat.set_sigma_s(m, 5, 5, 6.99913E-01)
  mat.set_sigma_s(m, 5, 6, 1.32440E-01)
  # 7 <- g'
  mat.set_sigma_s(m, 6, 1, 1.04550E-06)
  mat.set_sigma_s(m, 6, 2, 5.03440E-04)
  mat.set_sigma_s(m, 6, 3, 1.21390E-02)
  mat.set_sigma_s(m, 6, 4, 6.12290E-02)
  mat.set_sigma_s(m, 6, 5, 5.37320E-01)
  mat.set_sigma_s(m, 6, 6, 2.48070E+00)

  # Finalize
  mat.finalize()

  return mat
