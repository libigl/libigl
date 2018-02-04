#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  bool boolVariable = true;
  float floatVariable = 0.1f;
  enum Orientation { Up=0,Down,Left,Right } dir = Up;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

  // Extend viewer menu
  viewer.callback_init = [&](igl::opengl::glfw::Viewer& viewer)
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

    // Expose an enumaration type
    viewer.ngui->addVariable<Orientation>("Direction",dir)->setItems({"Up","Down","Left","Right"});

    // Add a button
    viewer.ngui->addButton("Print Hello",[](){ std::cout << "Hello\n"; });

    // Add an additional menu window
    viewer.ngui->addWindow(Eigen::Vector2i(220,10),"New Window");

    // Expose the same variable directly ...
    viewer.ngui->addVariable("float",floatVariable);

    // Generate menu
    viewer.screen->performLayout();

    return false;
  };

  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.launch();
}
