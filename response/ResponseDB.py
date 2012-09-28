##
#   @file   ResponseDB.py
#   @author Jeremy Roberts
#   @brief  ResponseDB class definition
#
import h5py as h5

class ResponseDB(object) :
  """ Supports HDF5 I/O of response functions
  """
  
  def __init__(self, filename) :
    # name of hdf5 file
    self.filename = filename
    # hdf5 file
    self.file = h5.File(filename, 'w')
    
  def initialize_node(self, response, responsename) :
    self.nodegroup = self.file.create_group(responsename)
    # if we're interpolating, set the keffs
    if response._flag :
      self.nodegroup['keffs'] = response._keffs

  def put(self, response, keffname) :
    # size of response block and number of sides
    n = response._number_responses
    s = response._number_sides
    # create /nodename/ki
    print "keffname = ", keffname
    keffgroup = self.nodegroup.create_group(keffname)
    # create datasets
    keffgroup['R'] = response._R
    keffgroup['F'] = response._F
    keffgroup['A'] = response._A
    keffgroup['L'] = response._L

