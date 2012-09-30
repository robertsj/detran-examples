##
#   @file   ResponseDB.py
#   @author Jeremy Roberts
#   @brief  ResponseDB class definition
#
import h5py as h5
import matplotlib.pyplot as plt

class ResponseDB(object) :
  """ Supports HDF5 I/O of response functions
  """
  
  def __init__(self, filename) :
    # name of hdf5 file
    self.filename = filename
    # hdf5 file
    self.file = h5.File(filename, 'w')
    
  def initialize_node(self, response, responsename) :
    """ Add a group to the database for this node 
    """
    # set the group for this node
    self.nodegroup = self.file.create_group(responsename)
    # get the node size
    size = response._number_responses
    # get number sides
    sides = response._number_sides
    # if we're interpolating, set the keffs and presize
    if response._flag :
      self.nodegroup['keffs'] = response._keffs
      n = len(response._keffs)
      dset = self.nodegroup.create_dataset("R", (n, size, size),  '=f8')
      dset = self.nodegroup.create_dataset("L", (n, sides, size), '=f8')
      dset = self.nodegroup.create_dataset("F", (n, size), '=f8')
      dset = self.nodegroup.create_dataset("A", (n, size), '=f8')

  def put(self, response, kindex) :

    # append datasets
    self.nodegroup['R'][kindex, :, :] = response._R
    self.nodegroup['F'][kindex, :   ] = response._F
    self.nodegroup['A'][kindex, :   ] = response._A
    self.nodegroup['L'][kindex, :, :] = response._L

