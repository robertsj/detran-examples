# test construction of linear material

from detran import *
from twogroup_materials import *

mat1 = get_rod_in()
mat2 = get_rod_out()
#mat1.display()
#mat2.display()

mats = vec_material([mat1, mat2, mat2, mat1])

mat = LinearMaterial.Create([2.0, 4.0, 10.0, 12.0], mats, "LINEARMAT")

n  = 40
sa = np.zeros(n)
times = np.zeros(n)
t = 0.0
for i in range(0, n) :
  mat.update(t, 1000000000.0, 1)
  times[i] = t
  sa[i] = mat.sigma_t(4, 1)
  t += 0.5

plt.plot(times, sa, 'k-o')
plt.show()

mat_ss = downcast(mat1)
print mat_ss
