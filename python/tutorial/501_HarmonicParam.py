import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
V_uv = igl.eigen.MatrixXd()


def key_down(viewer, key, modifier):
    if key == ord('1'):
        # Plot the 3D mesh
        viewer.data.set_mesh(V, F)
        viewer.core.align_camera_center(V, F)
    elif key == ord('2'):
        # Plot the mesh in 2D using the UV coordinates as vertex coordinates
        viewer.data.set_mesh(V_uv, F)
        viewer.core.align_camera_center(V_uv, F)
    viewer.data.compute_normals()
    return False


# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "camelhead.off", V, F)

# Find the open boundary
bnd = igl.eigen.MatrixXi()
igl.boundary_loop(F, bnd)

# Map the boundary to a circle, preserving edge proportions
bnd_uv = igl.eigen.MatrixXd()
igl.map_vertices_to_circle(V, bnd, bnd_uv)

# Harmonic parametrization for the internal vertices
igl.harmonic(V, F, bnd, bnd_uv, 1, V_uv)

# Scale UV to make the texture more clear
V_uv *= 5

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.data.set_uv(V_uv)
viewer.callback_key_down = key_down

# Disable wireframe
viewer.core.show_lines = False

# Draw checkerboard texture
viewer.core.show_texture = True

# Launch the viewer
viewer.launch()
