import sys, os
import numpy as np
import igl

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../tutorial/shared/beetle.off", V, F)

# Convert the mesh to numpy matrices (without copying it)
Vn = np.array(V, copy=False)
Fn = np.array(F, copy=False)

# Plot using matplotlib
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_trisurf(Vn[:,0], Vn[:,1], Vn[:,2], triangles=Fn, cmap=plt.cm.Spectral)
plt.show()
