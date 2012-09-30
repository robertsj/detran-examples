import numpy as np
from OrthogonalPoly import OrthogonalPoly

class DLP(OrthogonalPoly) :
  """ Discrete Legendre Polynomials
  """

  def __init__(self, order, length) :
    """ Constructor.

    Args:
           order (int):   order of expansion (maximum)
           length (int):  length of basis vectors
    """
    OrthogonalPoly.__init__(self, order, length)

  def build(self) :
    """ Build the basis set.
    """     
    # define P0
    self.P[0, :] = 1.0
    # define P1
    if self.order > 0 :
      self.P[1, :] = 1.0 - 2.0*np.asarray(range(0, self.length))/(self.length-1);
    # recursively define the rest
    N = self.length
    for o in range(2, self.order + 1) :
      for j in range(0, self.length) :
        c0 = (o-1)*(N-1+o)
        c1 = (2*(o-1)+1)*(N-1-2*j)
        c2 = o * (N - 1 - (o-1))
        self.P[o, j] = (c1 * self.P[o-1, j] - c0 * self.P[o-2, j])/c2
    # normalize the basis
    for o in range(0, self.order + 1) :
      self.P[o, :] = self.P[o, :] / np.linalg.norm(self.P[o, :])
    print self.P
