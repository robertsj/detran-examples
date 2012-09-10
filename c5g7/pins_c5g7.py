# pyexamples/pins_c5g7.py
#
# 2-D pin cell definitions for the C5G7 benchmark

from detran import *

def get_pins(number, flag) :
  """ Return the pins for the C5G7 benchmark.

  @param number   The number of meshes per dimension
  @param flag     The mesh type (false for uniform mesh with 
                  cell-center material, true for staircase)
  """

  # Shared
  pitch = 1.26
  radii = [0.54]
  # Pin 0 - UO2 
  pin0 = PinCell.Create(pitch, radii, [0,6])
  pin0.meshify(number, flag)
  # Pin 1 - 4.3% MOX
  pin1 = PinCell.Create(pitch, radii, [1,6])
  pin1.meshify(number, flag)
  # Pin 2 - 7.0% MOX
  pin2 = PinCell.Create(pitch, radii, [2,6])
  pin2.meshify(number, flag)
  # Pin 3 - 8.7% MOX
  pin3 = PinCell.Create(pitch, radii, [3,6])
  pin3.meshify(number, flag)
  # Pin 4 - Guide Tube
  pin4 = PinCell.Create(pitch, radii, [4,6])
  pin4.meshify(number, flag)
  # Pin 5 - Fission Chamber
  pin5 = PinCell.Create(pitch, radii, [5,6])
  pin5.meshify(number, flag)
  # Pin 6 - Moderator
  pin6 = PinCell.Create(pitch, radii, [6,6])
  pin6.meshify(number, flag)

  return pin0, pin1, pin2, pin3, pin4, pin5, pin6

