import numpy as np

class OrthogonalPoly(object) :
  """ Base class for orthogonal polynomials.

      An orthogonal polynomial basis can be used to 
      approximate functions by finite expansions.
  """

  def __init__(self, order, length) :
    """ Constructor.

    Args:
           order (int):   order of expansion (maximum)
           length (int):  length of basis vectors
    """
    self.order  = order
    self.length = length
    self.P = np.zeros((order + 1, length), 'd')

  def build(self) :
    """ Build the basis set.
    """     
    raise NotImplementedError   

  def expand(self, v) :
    """ Expand one or more vectors in the basis, returning the coefficients.

    In general, we project a set of # vectors into the subspace by
      Y = PV     sizes: (order+1, length) * (length, #vecs) -> (order+1, #vecs)
    from which we can get the approximation
      V ~ P'Y   sizes: (length, order+1) * (order+1, #vecs) --> (length, #vecs)
    We recover an exact (excluding roundoff) solution if order+1 = length.

    """
    try :
      nvecs = np.size(v[:, 0])
      V = v
    except :
      nvecs = 1
      V = np.zeros((1, len(v)))
      V[0, :] = v[:]
    coefs = np.zeros((self.order + 1, nvecs), 'd') 
    for o in range(0, self.order + 1) :
      for n in range(0, nvecs) :
        coefs[o, n] = np.dot(self.P[o, :], V[n, :])
    return coefs

  def recover(self, coefs) :
    """ Recover the approximation.
    """
    try :
      nvecs = np.size(coef[0, :])
    except :
      nvecs = 1
    V = np.zeros((nvecs, self.length))
    for n in range(0, nvecs) :
      for o in range(0, self.order + 1) :
        V[n, :] = V[n, :] + self.P[o, :]*coefs[o, n]
    if nvecs == 1 :
      v = np.zeros(self.length)
      v[:] = V[0, :]
      return v
    else :
      return V
