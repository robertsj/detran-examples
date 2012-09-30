##
#   @file   Response.py
#   @author Jeremy Roberts
#   @brief  Response class definition
#
from detran import *
from h5py import *
from DLP import *

class Response(object) :
  """ Computes response functions based on diffusion theory.

      This is a scoping study for computing multidimensional responses
      using diffusion.  It generates HDF5 response databases to help
      in testing serment and to support our paper on convergence of
      the RMM.

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
    """ Constructor.

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
    else :
      self.number_terms = keffs
    self.db = db

    # obtain the expansion orders.
    self.ao  = 0
    self.po  = 0
    self.so0 = 0
    self.so1 = 0
    if self.inp.check("response_spatial_order") :
      so = self.inp.get_int("response_spatial_order")
      self.so0 = so
      self.so1 = so # using constant expansion in all directions
    if self.inp.check("response_polar_order") :
      self.po = self.inp.get_int("response_polar_order")
    if self.inp.check("response_azimuth_order") :
      self.ao = self.inp.get_int("response_azimuth_order")
    if self.dimension == 1 :
      self.so0 = 0 # 1D has no spatial or azimuthal dependence
      self.so1 = 0 
      self.ao  = 0
    if self.dimension == 2 :
      self.so1 = 0 # 2D has only one spatial dof

    # create the state vector
    self.state = State.Create(self.inp, self.mesh)

    # define a response index.  this gives the energy and spatial orders.
    self.indices = []
    index = {}
    cardinal = 0
    for s in range(0, self.number_sides) :
      index['s'] = s
      for e in range(0, self.number_groups) :
        index['e'] = e
        for p in range(0, self.po + 1) :
          index['p'] = p
          for a in range(0, self.ao + 1) :
            index['a'] = a
            for o0 in range(0, self.so0 + 1) :
              index['s0'] = o0
              for o1 in range(0, self.so1 + 1) :
                index['s1'] = o1
                # Compute the polarity. 1D has no odd polarity combinations.  
                # 2D can have odd polarity via space and azimuth. In 3D, the
                # polarity can switch due to space, azimuth, and polar at top
                # and bottom boundaries and due to space and azimuth for other
                # boundaries.  Here, we assume bottom and top are sides 4 and 5.
                if self.dimension == 1 :
                  index['eo'] = 0                        
                elif self.dimension == 2 :
                  index['eo'] = (o0 + a) % 2            
                elif self.dimension == 3 :
                  if s > 3 :
                    index['eo'] = (o0 + o1 + a + p) % 2 
                  else :
                    index['eo'] = (o0 + o1 + a) % 2  
                index['index'] = cardinal
                self.indices.append(index.copy())
                cardinal = cardinal + 1

    # compute the number of responses
    self.number_responses = len(self.indices)
    print "number responses = ", self.number_responses 

    # allocate the response data
    self.R = np.zeros((self.number_responses, self.number_responses), 'd')
    self.F = np.zeros(self.number_responses, 'd')
    self.A = np.zeros(self.number_responses, 'd')
    self.L = np.zeros((2*self.dimension, self.number_responses), 'd')
    
    # initialize basis functions for each phase space variable
    self.basis_s = DLP(self.so0, self.mesh.number_cells_x())
    self.basis_s.build()

  def run(self) :
    """ Run the response generator for a given keff

        Note, once solver/loss is updated to recompute for keff, we can 
        make the solver once.  Also, need a way to pull out fission
        iteration coefficients
    """ 

    # initialize the inp for this node
    self.db.initialize_node(self, 'fuel')

    print "RESPONSE GENERATOR: "
    k = 0
    for keff in self.keffs :
        print "keff = ", keff
        self.inp.put_dbl("fission_scaling", keff)
        for index in self.indices :
          print " indices = ", index
          self.state.clear()
          # create solver
          if self.dimension == 1 :
            solver = DiffusionFixed1D(self.inp, self.material, self.mesh, self.state)
          elif self.dimension == 2 :
            solver = DiffusionFixed2D(self.inp, self.material, self.mesh, self.state)
          else :
            solver = DiffusionFixed3D(self.inp, self.material, self.mesh, self.state)
          solver.Q().set(0)
          b = solver.boundary()
          # set boundary and build source
          self.set_boundary(b, index)
          solver.build_source()
          # solve
          solver.solve()
          # compute responses
          self.compute(b, index)
        # put the response for this keff and increment k
        self.db.put(self, k)
        k = k + 1
        # sanity check: for diffusion, k=kinf should yield e=lambda=1.0
        e = np.linalg.eigvals(self.R)
        print " e = ", np.max(abs(e))
        # clear responses
        self.R[:, :] = 0.0
        self.F[:]    = 0.0
        self.A[:]    = 0.0
        self.L[:, :] = 0.0
  
  def set_boundary(self, b, index) :
    """ Set the boundary condition for a given incident order.  Add angular
        dependence later.
    """

    if self.dimension == 1 :

      # For 1D diffusion, we use an incident unit source on one side in one group
      b[index['s'], index['e'], BoundBASE1D.IN, 0, 0] = 1.0

    elif self.dimension == 2 :
  
      for i in range(0, len(self.basis_s.P[index['s0'], :])) :
        b[index['s'], index['e'], BoundBASE2D.IN, i, 0] = self.basis_s.P[index['s0'], i]

    else :
      for i in range(0, len(self.basis_s.P[index['s0'], :])) :
        for j in range(0, len(self.basis_s.P[index['s1'], :])) :
          b[index['s'], index['e'], BoundBASE3D.IN, i, j] = \
            self.basis_s.P[index['s0'], i] * self.basis_s.P[index['s1'], j]

  def compute(self, b, index) :
    """ Compute the responses produced from an incident order.
    """

    # incident cardinal response index
    r_index = index['index']

    # get mesh map
    mt = self.mesh.mesh_map("MATERIAL")

    # loop over all outgoing energy groups
    for g in range(0, self.material.number_groups()) :

      # loop over all cells to compute integrated response contributions
      for i in range(0, self.mesh.number_cells()) :
        phi_times_volume = self.mesh.volume(i) * self.state.phi(g)[i]
        self.F[r_index] = self.F[r_index] + \
                          phi_times_volume * self.material.sigma_f(mt[i], g)
        self.A[r_index] = self.A[r_index] + \
                          phi_times_volume * self.material.sigma_a(mt[i], g) 

    # compute the current response
    self.compute_leakage(b, index)

    # compute the boundary response
    self.compute_boundary(b, index)
    
  def compute_leakage(self, b, index) :
    """ Compute the leakage response
    """

    # incident cardinal response index
    r_index = index['index']

    # loop over all outgoing energy groups
    for g in range(0, self.material.number_groups()) :
      # loop over outgoing sides
      for s in range(0, self.number_sides) :

        # note, the leakage response is always a *net* leakage.  that's why the 
        # incident boundary is subtracted.  In many cases, the incident boundary
        # is zero, but this is a simple way to make sure it's accounted for.

        if self.dimension == 1 :
          self.L[s][r_index] = self.L[s][r_index] + \
                               b[s, g, BoundBASE1D.OUT] - b[s, g, BoundBASE1D.IN]

        elif self.dimension == 2 :

          dim = 1     # surfaces 0 and 1 have y-dependent boundaries
          if s > 1 :
            dim = 0   # surfaces 2 and 3 have x-dependent boundaries
          for i in range(0, self.mesh.number_cells(dim)) :
            self.L[s][r_index] = self.L[s][r_index] + \
                                 b[s, g, BoundBASE2D.OUT, i] - b[s, g, BoundBASE2D.IN, i]

        else :
          if s < 2 :
            dim = [1, 2]          # surfaces 0 and 1 have yz-dependent boundaries
          elif s == 2 or s == 3 :
            dim = [0, 2]          # surfaces 2 and 3 have xz-dependent boundaries
          else :
            dim = [0, 1]          # surfaces 4 and 5 have xy-dependent boundaries

          for i in range(0, self.mesh.number_cells(dim[0])) :
            for j in range(0, self.mesh.number_cells(dim[1])) :
              self.L[s][r_index] = self.L[s][r_index] + \
                                   b[s, g, BoundBASE3D.OUT, i, j] - b[s, g, BoundBASE3D.IN, i, j]

  def compute_boundary(self, b, index_i) :
    """ Compute the boundary response
    """

    dims = [ [1,2], [0,2], [0,1] ]

    # incident cardinal response index
    r_index_i = index_i['index']

    # loop over all outgoing moment indices
    for index_o in self.indices :
      r_index_o = index_o['index']
      
      s = index_o['s']
      g = index_o['e']
      s0 = index_o['s0']
      s1 = index_o['s1']

      if self.dimension == 1 :

        # 1D has no expansion
        self.R[r_index_o][r_index_i] = b[s, g, BoundBASE1D.OUT]

      else :

        # 2D/3D require we extract the boundary function and expand
        bf = np.asarray(b(s, g, BoundBASE1D.OUT))
        if self.dimension == 2 :
          # in 2D, compute the coefficient directly by dotting with the basic function
          self.R[r_index_o][r_index_i] = np.dot(self.basis_s.P[s0, :], bf)
        else :
          # in 3D, first collapse the boundary matrix from bf(i, j) --> bf(i), a vector
          # of coefficients retaining dependence in one dimension.  Then dot this with its
          # basis to recover the coefficient.
          bf0 = np.zeros(len(bf[0]))
          for i in range(0, len(bf0)) :
            bf0[i] = np.dot(self.basis_s.P[s0, :], bf[i, :])
          self.R[r_index_o][r_index_i] = np.dot(self.basis_s.P[s1, :], bf0)

        

  
