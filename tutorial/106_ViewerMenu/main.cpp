#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

bool boolVariable = true;
float floatVariable = 0.1f;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  // Init the viewer
  igl::viewer::Viewer viewer;

  // Extend viewer menu
  viewer.callback_init = [&](igl::viewer::Viewer& viewer)
  {
    // Add new group
    viewer.ngui->addGroup("New Group");

    // Expose variable directly ...
    viewer.ngui->addVariable("float",floatVariable);

    // ... or using a custom callback
    viewer.ngui->addVariable<bool>("bool",[&](bool val) {
      boolVariable = val; // set
    },[&]() {
      return boolVariable; // get
    });

    // Add a button
    viewer.ngui->addButton("Print Hello",[](){ std::cout << "Hello\n"; });

    // Add an additional bar
    viewer.ngui->addWindow(Eigen::Vector2i(220,10),"New Window");

    // Expose the same variable directly ...
    viewer.ngui->addVariable("float",floatVariable);

    // Generate menu
    viewer.screen->performLayout();

    return false;
  };

  // Plot the mesh
  viewer.data.set_mesh(V, F);
  viewer.launch();
}
