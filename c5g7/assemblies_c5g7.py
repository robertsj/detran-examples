# pyexamples/assemblies_c5g7.py
#
# 2-D assemblies definitions for the C5G7 benchmark

from pins_c5g7 import *
from detran import *

def get_assemblies(number, flag) :
  """ Return the assemblies for the C5G7 benchmark.

  See get_pincells for parameter definition.
  """

  # Shared things
  G = 4 # guide tube
  F = 5 # fission chamber
  pin0, pin1, pin2, pin3, pin4, pin5, pin6 = get_pins(number, flag)

  # Assembly 1 -- UO2
  assem1 = Assembly.Create(17)
  assem1.add_pincell(pin0)
  assem1.add_pincell(pin1)
  assem1.add_pincell(pin2)
  assem1.add_pincell(pin3)
  assem1.add_pincell(pin4)
  assem1.add_pincell(pin5)
  assem1.add_pincell(pin6)
  pin_map1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
              0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,F,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
              0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  assem1.finalize(pin_map1)

  # Assembly 2 - MOX
  assem2 = Assembly.Create(17)
  assem2.add_pincell(pin0)
  assem2.add_pincell(pin1)
  assem2.add_pincell(pin2)
  assem2.add_pincell(pin3)
  assem2.add_pincell(pin4)
  assem2.add_pincell(pin5)
  assem2.add_pincell(pin6)
  pin_map2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
              1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
              1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
              1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
              1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
              1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,G,3,3,G,3,3,F,3,3,G,3,3,G,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
              1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
              1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
              1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
              1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
  assem2.finalize(pin_map2)

  # Assembly 3 - Moderator
  assem3 = Assembly.Create(17)
  assem3.add_pincell(pin0)
  assem3.add_pincell(pin1)
  assem3.add_pincell(pin2)
  assem3.add_pincell(pin3)
  assem3.add_pincell(pin4)
  assem3.add_pincell(pin5)
  assem3.add_pincell(pin6)
  pin_map3 = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];
  assem3.finalize(pin_map3)

  return assem1, assem2, assem3

