import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def plot_mesh_map(obj, map_key, edges=False) :
  """ Plot the mesh map.
  """
  # Get the map.
  m_map = np.asarray(obj.mesh_map(map_key))
  if len(m_map) == 0 :
    print "Oops! Map has size zero!"
    return
        
  # Reshape it.
  m_map2 = m_map.reshape(obj.number_cells_x(), obj.number_cells_y())
   
  # Define random color map.
  colormap = matplotlib.colors.ListedColormap ( np.random.rand ( 10000,3))
    
  # Get the mesh axes and then make a grid of them for plotting.
  x, y = mesh_axes(obj)
  X, Y = np.meshgrid(x, y)
    
  if edges == False :
    edges = None
  else :
    edges = "black"
 
  # Create the plot a
  plt.pcolor(X, Y, m_map2, cmap=colormap)
  plt.show()


def plot_flux(obj, f) :
  """ Plot a flux or flux-shaped vector on the mesh."
  """
  # Reshape flux.
  f = f.reshape(obj.number_cells_x(), obj.number_cells_y())
  #print f
  # Plot me. 
  x, y = mesh_axes(obj)
  X, Y = np.meshgrid(x, y)
  plt.pcolor(X, Y, f, cmap="hot")
  plt.colorbar()
  plt.show()

#  plt.imshow(f, interpolation='nearest', cmap=plt.cm.hot)
#  plt.title('flux map')
#  plt.colorbar()
#  plt.show()

def mesh_axes(obj) :
  """ Get the fine mesh esges defining cells for plotting.
  """
  x = np.zeros(obj.number_cells_x()+1)
  y = np.zeros(obj.number_cells_y()+1)
  for i in range(0, obj.number_cells_x()) :
    x[i + 1] = x[i] + obj.dx(i)
  for j in range(0, obj.number_cells_y()) :
    y[j + 1] = y[j] + obj.dy(j)
  return (x, y)
