import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl


from shared import TUTORIAL_SHARED_PATH, check_dependencies, print_usage

dependencies = ["viewer"]
check_dependencies(dependencies)


def mouse_down(viewer, a, b):
    bc = igl.eigen.MatrixXd()

    # Cast a ray in the view direction starting from the mouse position
    fid = igl.eigen.MatrixXi([-1])
    coord = igl.eigen.MatrixXd([viewer.current_mouse_x, viewer.core.viewport[3] - viewer.current_mouse_y])
    hit = igl.unproject_onto_mesh(coord, viewer.core.view * viewer.core.model,
      viewer.core.proj, viewer.core.viewport, V, F, fid, bc)
    if hit:
        # paint hit red
        C.setRow(fid[0, 0], igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data.set_colors(C)
        return True

    return False


if __name__ == "__main__":
    keys = {"click": "Pick face on shape"}
    print_usage(keys)

    # Mesh with per-face color
    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()
    C = igl.eigen.MatrixXd()

    # Load a mesh in OFF format
    igl.readOFF(TUTORIAL_SHARED_PATH + "fertility.off", V, F)

    # Initialize white
    C.setConstant(F.rows(), 3, 1.0)

    # Show mesh
    viewer = igl.viewer.Viewer()
    viewer.data.set_mesh(V, F)
    viewer.data.set_colors(C)
    viewer.core.show_lines = False
    viewer.callback_mouse_down = mouse_down
    viewer.launch()
