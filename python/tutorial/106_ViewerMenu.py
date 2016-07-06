# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import pyigl as igl
import nanogui

V1 = igl.eigen.MatrixXd()
F1 = igl.eigen.MatrixXi()

V2 = igl.eigen.MatrixXd()
F2 = igl.eigen.MatrixXi()

float_variable = 0.1
bool_variable = True
dir = 0


def make_accessors(name):
    def setter(value):
        globals()[name] = value

    def getter():
        return globals()[name]
    return setter, getter


def viewer_init(viewer):
    # add new group
    viewer.ngui.addGroup("New Group")

    # Expose the using general callback
    viewer.ngui.addDoubleVariable("double", *make_accessors("float_variable"))

    def setter(val):
        global bool_variable
        bool_variable = val

    def getter():
        global bool_variable
        return bool_variable

    # ... or using a custom callback
    viewer.ngui.addBoolVariable("bool", setter, getter)

    viewer.ngui.addEnumVariable("Direction", *make_accessors("dir")) \
        .setItems(["Up", "Down", "Left", "Right"])

    # Add a button
    def cb():
        print("Hello")
    viewer.ngui.addButton("Print Hello", cb)

    #Add an additional menu window
    viewer.ngui.addWindow(nanogui.Vector2i(220, 10), "New Window")

    # add accessor
    viewer.ngui.addDoubleVariable("double", *make_accessors("float_variable"))

    #Generate menu
    viewer.screen.performLayout()

    return False



def main():
    # Load a mesh in OFF format
    igl.readOFF("../../tutorial/shared/bunny.off", V1, F1)

    # Init the viewer
    viewer = igl.viewer.Viewer()

    # Extend viewer menu
    viewer.callback_init = viewer_init

    # Plot the mesh
    viewer.data.set_mesh(V1, F1)
    viewer.launch()


if __name__ == "__main__":
    main()
