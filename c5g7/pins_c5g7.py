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
  pitch = Point(1.26, 1.26)
  radii = [0.54]
  # Pin 0 - UO2 
  pin0 = PinCell.Create(pitch, [0,6], radii)
  pin0.meshify(number, flag)
  # Pin 1 - 4.3% MOX
  pin1 = PinCell.Create(pitch, [1,6], radii)
  pin1.meshify(number, flag)
  # Pin 2 - 7.0% MOX
  pin2 = PinCell.Create(pitch, [2,6], radii)
  pin2.meshify(number, flag)
  # Pin 3 - 8.7% MOX
  pin3 = PinCell.Create(pitch,[3,6], radii)
  pin3.meshify(number, flag)
  # Pin 4 - Guide Tube
  pin4 = PinCell.Create(pitch,[4,6], radii)
  pin4.meshify(number, flag)
  # Pin 5 - Fission Chamber
  pin5 = PinCell.Create(pitch, [5,6], radii)
  pin5.meshify(number, flag)
  # Pin 6 - Moderator
  pin6 = PinCell.Create(pitch, [6,6], radii)
  pin6.meshify(number, flag)


  return pin0, pin1, pin2, pin3, pin4, pin5, pin6

