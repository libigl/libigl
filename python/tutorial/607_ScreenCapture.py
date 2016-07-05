import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["png", "viewer"]
check_dependencies(dependencies)

temp_png = os.path.join(os.getcwd(),"out.png")

def key_down(viewer, key, modifier):
    if key == ord('1'):
        # Allocate temporary buffers
        R = igl.eigen.MatrixXuc(1280, 800)
        G = igl.eigen.MatrixXuc(1280, 800)
        B = igl.eigen.MatrixXuc(1280, 800)
        A = igl.eigen.MatrixXuc(1280, 800)

        # Draw the scene in the buffers
        viewer.core.draw_buffer(viewer.data, viewer.opengl, False, R, G, B, A)

        # Save it to a PNG
        igl.png.writePNG(R, G, B, A, temp_png)
    elif key == ord('2'):
        # Allocate temporary buffers
        R = igl.eigen.MatrixXuc()
        G = igl.eigen.MatrixXuc()
        B = igl.eigen.MatrixXuc()
        A = igl.eigen.MatrixXuc()

        # Read the PNG
        igl.png.readPNG(temp_png, R, G, B, A)

        # Replace the mesh with a triangulated square
        V = igl.eigen.MatrixXd([[-0.5, -0.5, 0],
                                [0.5, -0.5, 0],
                                [0.5, 0.5, 0],
                                [-0.5, 0.5, 0]])

        F = igl.eigen.MatrixXi([[0, 1, 2], [2, 3, 0]])

        UV = igl.eigen.MatrixXd([[0, 0], [1, 0], [1, 1], [0, 1]])

        viewer.data.clear()
        viewer.data.set_mesh(V, F)
        viewer.data.set_uv(UV)
        viewer.core.align_camera_center(V)
        viewer.core.show_texture = True

        # Use the image as a texture
        viewer.data.set_texture(R, G, B)

    else:
        return False

    return True


if __name__ == "__main__":
    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()

    # Load meshes in OFF format
    igl.readOFF(TUTORIAL_SHARED_PATH + "bunny.off", V, F)

    viewer = igl.viewer.Viewer()

    print(
        "Usage: Press 1 to render the scene and save it in a png. \nPress 2 to load the saved png and use it as a texture.")

    viewer.callback_key_down = key_down
    viewer.data.set_mesh(V, F)
    viewer.launch()

    os.remove(temp_png)
