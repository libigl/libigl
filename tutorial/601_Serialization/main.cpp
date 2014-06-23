#include <igl/readOFF.h>
#include <iostream>

#include <igl/xml/XMLSerializer.h>
#include <igl/xml/XMLSerialization.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

class State : public ::igl::XMLSerialization
{
public:
  State() : XMLSerialization("dummy") {}

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::vector<int> ids;

  // You have to define this function to
  // register the fields you want to serialize
  void InitSerialization()
  {
    xmlSerializer->Add(V  , "V");
    xmlSerializer->Add(F  , "F");
    xmlSerializer->Add(ids, "ids");
  }

};

int main(int argc, char *argv[])
{
  State state;

  // Load a mesh in OFF format
  igl::readOFF("../shared/cube.off", state.V, state.F);

  // Save some integers in a vector
  state.ids.push_back(6);
  state.ids.push_back(7);

  // Serialize to XML the state of the application
  ::igl::XMLSerializer serializer_save("601_Serialization");
  serializer_save.Add(state,"State");
  serializer_save.Save("temp.xml",true);

  // Load the state from the same XML file
  State loaded_state;
  ::igl::XMLSerializer serializer_load("601_Serialization");
  serializer_load.Add(loaded_state,"State");
  serializer_load.Load("temp.xml");

  // Plot the state
  std::cout << "Vertices: " << std::endl << loaded_state.V << std::endl;
  std::cout << "Faces:    " << std::endl << loaded_state.F << std::endl;
  std::cout << "ids:      " << std::endl
            << loaded_state.ids[0] << " " << loaded_state.ids[1] << std::endl;

}
