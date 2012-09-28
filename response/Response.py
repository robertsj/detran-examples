##
#   @file   Response.py
#   @author Jeremy Roberts
#   @brief  Response class definition
#
from detran import *
from h5py import *

class Response(object) :
  """ Computes response functions

      This is a scoping study for computing multidimensional responses
      using diffusion.
  """
  
  def __init__(self, inp, material, mesh, db, flag=False, keffs=[]) :
    """ Initialize the response generator.

    Args:
           inp (InputDB):           Input database
           material (Material):     Material database
           mesh (Mesg):             Geometry definition

    Kwargs:
           state (bool): Current state to be in.
    """
    self._inp       = inp                                
    self._material  = material                         
    self._mesh      = mesh                                
    self._number_groups = self._material.number_groups()  
    self._dimension = self._mesh.dimension()
    self._number_sides = self._dimension * 2
    self._flag = flag
    if self._flag :
      self._keffs = keffs
    self._db = db # response database

    # obtain the spatial order.
    self._order = 0
    if self._inp.check("response_spatial_order") :
      self._order = self._inp.get_int("response_spatial_order")
    if self._dimension == 1 :
      self._order = 0 # 1D has no spatial dependence

    # create the state vector
    self._state = State.Create(self._inp, self._mesh)

    # compute the number of responses
    self._number_responses = self._number_sides             * \
                             (1 + self._order)              * \
                             self._material.number_groups() 
    print "number responses = ", self._number_responses 

    # allocate the response data
    self._R = np.zeros((self._number_responses, self._number_responses), 'd')
    self._F = np.zeros(self._number_responses, 'd')
    self._A = np.zeros(self._number_responses, 'd')
    self._L = np.zeros((2*self._dimension, self._number_responses), 'd')

    # define a response index
    

  def run(self) :
    """ Run the response generator for a given keff

        Note, once solver/loss is updated to recompute for keff, we can 
        make the solver once.  Also, need a way to pull out fission
        iteration coefficients
    """ 

    # initialize the inp for this node
    self._db.initialize_node(self, 'fuel')

    print "RESPONSE GENERATOR: "
    k = 0
    for keff in self._keffs :
        print "keff = ", keff
        self._inp.put_dbl("fission_scaling", 11000)
        for side in range(0, self._number_sides) :
            print "  side = ", side
            for group in range(0, self._number_groups) :
                print "    groups = ", group
                for order in range(0, self._order + 1) :
                    print "    order = ", order
                    # clear the state
                    self._state.clear()
                    # create solver
                    solver = DiffusionFixed1D(self._inp, self._material, self._mesh, self._state)
                    solver.Q().set(0)
                    # set boundary
                    b = solver.boundary()
                    b[side, group, BoundBASE1D.IN, 0, 0] = 1.0
                    solver.build_source()
                    # solve
                    solver.solve()
                    print np.asarray(self._state.phi(0))
                    # compute responses
                    self.compute(side, group, order, b)
        # put the response for this keff and increment k
        self._db.put(self, 'k'+str(k))
        k = k + 1
        # clear responses
        self._R[:, :] = 0.0
        self._F[:]    = 0.0
        self._A[:]    = 0.0
        self._L[:, :] = 0.0

  def compute(self, side, group, order, b) :

    # compute cardinal incident response index
    index = order + \
            group * (self._order + 1) * self._number_sides + \
            side * (self._number_groups * (self._order + 1))
    print "      index = ", index

    # get mesh map
    mt = self._mesh.mesh_map("MATERIAL")

    # loop over all outgoing energy groups
    for g in range(0, self._material.number_groups()) :

      # loop over all cells to compute integrated response contributions
      for i in range(0, self._mesh.number_cells()) :
        phi_times_volume = self._mesh.volume(i) * self._state.phi(g)[i]
        self._F[index] = self._F[index] + phi_times_volume * self._material.sigma_f(mt[i], g)
        self._A[index] = self._A[index] + phi_times_volume * self._material.sigma_a(mt[i], g) 

      # loop over all sides and compute current response and leakage
      for s in range(0, self._number_sides) :
        index_o = 0 + \
                  g * (self._order + 1) * self._number_sides + \
                  s * (self._number_groups * (self._order + 1))
        self._L[s][index] = self._L[s][index] + b[s, g, BoundBASE1D.OUT] - b[s, g, BoundBASE1D.IN]
        self._R[index_o][index] = b[s, g, BoundBASE1D.OUT, 0, 0]
    print self._F


