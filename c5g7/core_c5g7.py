# pyexamples/core_c5g7.py
#
# 2-D core definition for the C5G7 benchmark

from assemblies_c5g7 import get_assemblies
from detran import *

def get_core(number, flag) :
  """ Return the core for the C5G7 benchmark.

  See get_pincells for parameter definition.
  """

  assemblies = get_assemblies(number, flag)
  core       = Core.Create(3)
  core.add_assembly(assemblies[0])
  core.add_assembly(assemblies[1])
  core.add_assembly(assemblies[2])
  core_map = [0,1,2,
              1,0,2,
              2,2,2]
  core.finalize(core_map)
  return core


