import numpy as np
from numpy.polynomial.legendre import *
from OrthogonalPoly import OrthogonalPoly

class CLP(OrthogonalPoly) :
  """ Continuous Legendre Polynomials

      Unlike the discrete variant, to apply CLP's in expansions, we need
      to integrate over a discretization of a continuous field.  Hence,
      the function evaluation points are needed rather than just the
      number of points.  We build the integration into the basic directly.
  """

  def __init__(self, order, abscissa, width) :
    """ Constructor.

    Args:
           order (int):       order of expansion (maximum)
           abscissa (array):  points at which function is evaluated
           width (float):     width of domain (starting at zero)
    """
    OrthogonalPoly.__init__(self, order, len(abscissa))
    # shift the abscissa to [-1, 1]
    self.abscissa = 2.0 * abscissa / width - 1.0

  def build(self) :
    """ Build the basis set.
    """     
    self.P = (legvander(self.abscissa, self.order)).transpose

