#include <igl/readOFF.h>
#include <igl/serialize.h>
#include <igl/xml/serialize_xml.h>
#include <iostream>


Eigen::MatrixXd V;
Eigen::MatrixXi F;

// derive from igl::Serializable to serialize your own type
struct State : public igl::Serializable
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::vector<int> ids;

  // You have to define this function to
  // register the fields you want to serialize
  virtual void InitSerialization()
  {
    this->Add(V  , "V");
    this->Add(F  , "F");
    this->Add(ids, "ids");
  }
};

//// alternatively you can do it like the following to have
//// a non-intrusive serialization:
////
//struct State
//{
//  Eigen::MatrixXd V;
//  Eigen::MatrixXi F;
//  std::vector<int> ids;
//};
//
//
//namespace igl
//{
//  namespace serialization
//  {
//    // the `template <>` is essential
//    template <> inline void serialize(const State& obj,std::vector<char>& buffer){
//      ::igl::serialize(obj.V,std::string("V"),buffer);
//      ::igl::serialize(obj.F,std::string("F"),buffer);
//      ::igl::serialize(obj.ids,std::string("ids"),buffer);
//    }
//    template <> inline void deserialize(State& obj,const std::vector<char>& buffer){
//      ::igl::deserialize(obj.V,std::string("V"),buffer);
//      ::igl::deserialize(obj.F,std::string("F"),buffer);
//      ::igl::deserialize(obj.ids,std::string("ids"),buffer);
//    }
//  }
//}
//
////OR:
//
//SERIALIZE_TYPE(State,
// SERIALIZE_MEMBER(V)
//  SERIALIZE_MEMBER(F)
//  SERIALIZE_MEMBER_NAME(ids,"ids")
//)

int main(int argc, char *argv[])
{
  std::string binaryFile = "binData";
  std::string xmlFile = "data.xml";

  bool b = true;
  unsigned int num = 10;
  std::vector<float> vec = {0.1f,0.002f,5.3f};

  // use overwrite = true for the first serialization to create or overwrite an existing file
  igl::serialize(b,"B",binaryFile,true);
  // append following serialization to existing file
  igl::serialize(num,"Number",binaryFile);
  igl::serialize(vec,"VectorName",binaryFile);

  // deserialize back to variables
  igl::deserialize(b,"B",binaryFile);
  igl::deserialize(num,"Number",binaryFile);
  igl::deserialize(vec,"VectorName",binaryFile);

  State stateIn, stateOut;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/2triangles.off",stateIn.V,stateIn.F);

  // Save some integers in a vector
  stateIn.ids.push_back(6);
  stateIn.ids.push_back(7);

  // Serialize the state of the application
  igl::serialize(stateIn,"State",binaryFile,true);

  // Load the state from the same file
  igl::deserialize(stateOut,"State",binaryFile);

  // Plot the state
  std::cout << "Vertices: " << std::endl << stateOut.V << std::endl;
  std::cout << "Faces:    " << std::endl << stateOut.F << std::endl;
  std::cout << "ids:      " << std::endl
            << stateOut.ids[0] << " " << stateOut.ids[1] << std::endl;

  // XML serialization

  // binary = false, overwrite = true
  igl::xml::serialize_xml(vec,"VectorXML",xmlFile,false,true);
  // binary = true, overwrite = false
  igl::xml::serialize_xml(vec,"VectorBin",xmlFile,true,false);
  igl::xml::deserialize_xml(vec,"VectorXML",xmlFile);
  igl::xml::deserialize_xml(vec,"VectorBin",xmlFile);
}
