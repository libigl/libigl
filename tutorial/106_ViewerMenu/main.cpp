#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

bool boolVariable = true;
float floatVariable = 0.1f;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/bunny.off", V, F);

  // Init the viewer
  igl::viewer::Viewer viewer;

  // Extend viewer menu
  viewer.callback_init = [&](igl::viewer::Viewer& viewer)
  {
    // Add new group
    viewer.ngui->addNewGroup("New Group");

    // Expose variable directly ...
    viewer.ngui->addVariable(floatVariable,"float");

    // ... or using a custom callback
    viewer.ngui->addVariable([&](bool val) {
      boolVariable = val; // set
    },[&]() {
      return boolVariable; // get
    },"bool",true);

    // Add a button
    viewer.ngui->addButton("Print Hello",[](){ cout << "Hello\n"; });

    // Add an additional bar
    viewer.ngui->addNewWindow(Eigen::Vector2i(220,10),"New Window");

    // Expose the same variable directly ...
    viewer.ngui->addVariable(floatVariable,"float");

    // Generate menu
    viewer.ngui->layout();

    return false;
  };

  // Plot the mesh
  viewer.data.set_mesh(V, F);
  viewer.launch();
}
