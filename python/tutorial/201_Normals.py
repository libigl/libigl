# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

N_vertices = igl.eigen.MatrixXd()
N_faces = igl.eigen.MatrixXd()
N_corners = igl.eigen.MatrixXd()

# This function is called every time a keyboard button is pressed
def key_pressed(viewer, key, modifier):
    if key == ord('1'):
      viewer.data.set_normals(N_faces)
      return True
    elif key == ord('2'):
      viewer.data.set_normals(N_vertices)
      return True
    elif key == ord('3'):
      viewer.data.set_normals(N_corners)
      return True
    return False

# Load a mesh in OFF format
igl.readOFF("../../tutorial/shared/fandisk.off", V, F);

# Compute per-face normals
N_faces = igl.eigen.MatrixXd()
igl.per_face_normals(V,F,N_faces)

# Compute per-vertex normals
N_vertices = igl.eigen.MatrixXd()
igl.per_vertex_normals(V,F,igl.PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,N_vertices)

# Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
N_corners = igl.eigen.MatrixXd()
igl.per_corner_normals(V,F,20,N_corners)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.callback_key_pressed = key_pressed
viewer.core.show_lines = False
viewer.data.set_mesh(V, F)
viewer.data.set_normals(N_faces)
print("Press '1' for per-face normals.")
print("Press '2' for per-vertex normals.")
print("Press '3' for per-corner normals.")
viewer.launch()
