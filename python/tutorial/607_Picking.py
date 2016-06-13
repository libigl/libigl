import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl


from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


# Mesh with per-face color
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
C = igl.eigen.MatrixXd()

viewer = igl.viewer.Viewer()

def mouse_down(viewer, a, b):
    bc = igl.eigen.MatrixXd()

    # Cast a ray in the view direction starting from the mouse position
    fid = igl.eigen.MatrixXi([-1])
    coord = igl.eigen.MatrixXd([viewer.current_mouse_x, viewer.core.viewport[3] - viewer.current_mouse_y])
    hit = igl.unproject_onto_mesh(coord, viewer.core.view * viewer.core.model,
      viewer.core.proj, viewer.core.viewport, V, F, fid, bc)
    if hit:
        C.setRow(fid[0, 0], igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data.set_colors(C)
        return True

    return False


print("Usage: [LeftMouseClick] to select a face")

# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "fertility.off", V, F)

# Initialize white
C.setConstant(F.rows(), 3, 1.0)

viewer.data.set_mesh(V, F)
viewer.data.set_colors(C)
viewer.core.show_lines = False
viewer.callback_mouse_down = mouse_down
viewer.launch()
