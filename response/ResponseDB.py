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

        Each node can be based on a completely different
        approximation; the end user must ensure nodes are
        coupled appropriately.  Moreover, each node can 
        be expanded or interpolated.
    """
    # set the group for this node
    self.nodegroup = self.file.create_group(responsename)
    # get the node size
    size = response.number_responses
    self.nodegroup.attrs['response_size'] = size
    # get number sides
    sides = response.number_sides
    self.nodegroup.attrs['number_surfaces'] = sides
    if response.flag :
      # we're interpolating, so set the keffs
      self.nodegroup['keffs'] = response.keffs
      n = len(response.keffs)
      self.nodegroup.attrs['scheme'] = 1
    else :
      # we're expanding, so set the number of terms
      n = response.keffs
      self.nodegroup.attrs['scheme'] = 0
    # set the number of terms (applies to both schemes)
    self.nodegroup.attrs['number_keffs'] = n
    dset = self.nodegroup.create_dataset("R", (n, size, size),  '=f8')
    dset = self.nodegroup.create_dataset("L", (n, sides, size), '=f8')
    dset = self.nodegroup.create_dataset("F", (n, size), '=f8')
    dset = self.nodegroup.create_dataset("A", (n, size), '=f8')

  def put(self, response, kindex) :

    # append datasets
    self.nodegroup['R'][kindex, :, :] = response.R
    self.nodegroup['F'][kindex, :   ] = response.F
    self.nodegroup['A'][kindex, :   ] = response.A
    self.nodegroup['L'][kindex, :, :] = response.L
    print response.F

  def close(self) :
   print "closing hdf5"
   self.file.close()
