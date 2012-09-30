##
#   @file   Response.py
#   @author Jeremy Roberts
#   @brief  Response class definition
#
from detran import *
from h5py import *

class Response(object) :
  """ Computes response functions based on diffusion theory.

      This is a scoping study for computing multidimensional responses
      using diffusion.  It generates HDF5 response databases to help
      in testing serment.

      Responses are defined in one of two ways.  The first way is 
      to define responses at discrete values of the k-eigenvalue.  For
      arbitrary k values, the resulting response can be found by 
      interpolation.  A second approach is to expand the responses
      in the form 
        R(k) ~ R_0 + R_1 / k + R_2 / k^2 + ...
      which be analytically summed by noting
        R_{n+1}/R_{n} -> k_{vac} = node eigenvalue in vacuum
      
      For the first approach, the user gives the values of k for which
      to evaluate responses.  For the second approach, the user specifies
      the number of terms to use.  The eigenvalue in vacuum is estimated
      by taking the ratio of the final two terms.  There is probably a
      sensible tolerance that can be specified to estimate the number of
      terms needed a priori.

      The second approach is expected to be more computationally and
      efficient if only a few terms are required. 
  """
  
  def __init__(self, inp, material, mesh, db, flag=False, keffs=[]) :
    """ Initialize the response generator.

    Args:
           inp (InputDB):           Input database
           material (Material):     Material database
           mesh (Mesh):             Geometry definition
           db (ResponseDB):         HDF5 database for responses
           flag (bool):             Switches between interpolation and expansion
           keffs (...):             Either the number of terms or the discrete k values

    """
    self.inp       = inp                                
    self.material  = material                         
    self.mesh      = mesh                                
    self.number_groups = self.material.number_groups()  
    self.dimension     = self.mesh.dimension()
    self.number_sides  = self.dimension * 2
    self.flag = flag
    if self.flag :
      self.keffs = keffs
    self.db = db

    # obtain the spatial order.  assumed constant for both dimensions in 3D.
    self.order = 0
    if self.inp.check("response_spatial_order") :
      so = self.inp.get_int("response_spatial_order")
      self.so0 = so
      self.so1 = so
    if self.dimension == 1 :
      self.so0 = 0 # 1D has no spatial dependence
      self.so1 = 0
    if self.dimension == 2 :
      self.so1 = 0 # 2D has only one spatial dof

    # create the state vector
    self.state = State.Create(self.inp, self.mesh)

    # define a response index.  this gives the energy and spatial orders.
    self.indices = []
    index = {}
    for s in range(0, self.number_sides) :
      index['s'] = s
      for g in range(0, self.number_groups) :
        index['g'] = g
        for o0 in range(0, self.so0 + 1) :
          index['s0'] = o0
          for o1 in range(0, self.so1 + 1) :
            index['s1'] = o1
            index['eo'] = (o0 + o1) % 2
            self.indices.append(index)
    print self.indices

    # compute the number of responses
    self.number_responses = len(self.indices)
    print "number responses = ", self.number_responses 

    # allocate the response data
    self.R = np.zeros((self.number_responses, self.number_responses), 'd')
    self.F = np.zeros(self.number_responses, 'd')
    self.A = np.zeros(self.number_responses, 'd')
    self.L = np.zeros((2*self.dimension, self.number_responses), 'd')
    

  def run(self) :
    """ Run the response generator for a given keff

        Note, once solver/loss is updated to recompute for keff, we can 
        make the solver once.  Also, need a way to pull out fission
        iteration coefficients
    """ 

    # initialize the inp for this node
    self.db.initialize_node(self, 'fuel')

    dims = [0, 0, 1, 1, 2, 2]

    print "RESPONSE GENERATOR: "
    k = 0
    for keff in self.keffs :
        print "keff = ", keff
        self.inp.put_dbl("fission_scaling", keff)
        for index in self.indices :
          print " indices = ", index
          self.state.clear()
          # create solver
          solver = DiffusionFixed1D(self.inp, self.material, self.mesh, self.state)
          solver.Q().set(0)
          b = solver.boundary()
          # set boundary and build source
          set_boundary(b)
          solver.build_source()
          # solve
          solver.solve()
          print np.asarray(self.state.phi(0))
          # compute responses
          self.compute(side, group, order, b)
        # put the response for this keff and increment k
        self.db.put(self, k)
        k = k + 1
        # clear responses
        print self.R 
        e = np.linalg.eigvals(self.R)
        print " e = ", np.max(abs(e))
        self.R[:, :] = 0.0
        self.F[:]    = 0.0
        self.A[:]    = 0.0
        self.L[:, :] = 0.0
  
  def compute(self, side, group, order, b) :
    """ Compute the responses for incident side, group, and order
    """

    # compute cardinal incident response index
    index = order + \
            group * (self.order + 1) * self.number_sides + \
            side * (self.number_groups * (self.order + 1))
    print "      index = ", index

    # get mesh map
    mt = self.mesh.mesh_map("MATERIAL")

    # loop over all outgoing energy groups
    for g in range(0, self.material.number_groups()) :

      # loop over all cells to compute integrated response contributions
      for i in range(0, self.mesh.number_cells()) :
        phi_times_volume = self.mesh.volume(i) * self.state.phi(g)[i]
        self.F[index] = self.F[index] + phi_times_volume * self.material.sigma_f(mt[i], g)
        self.A[index] = self.A[index] + phi_times_volume * self.material.sigma_a(mt[i], g) 

      # loop over all sides and compute current response and leakage
      for s in range(0, self.number_sides) :
        index_o = 0 + \
                  g * (self.order + 1) * self.number_sides + \
                  s * (self.number_groups * (self.order + 1))
        # compute leakage (this will be best done with a specialized function)
        self.L[s][index] = self.L[s][index] + b[s, g, BoundBASE1D.OUT] - b[s, g, BoundBASE1D.IN]
        # current response (specialized)
        self.R[index_o][index] = b[s, g, BoundBASE1D.OUT, 0, 0]
    
  def set_boundary(self, side, group, order, b) :
    
    b[side, group, BoundBASE1D.IN, 0, 0] = 1.0

