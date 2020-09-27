#include <igl/opengl/glfw/Viewer.h>
#include <igl/barycenter.h>
#include <igl/colormap.h>

#include <igl/readMSH.h>
#include <igl/readMESH.h>


Eigen::MatrixXd X,B;
Eigen::MatrixXi Tri;
Eigen::MatrixXi Tet;
Eigen::VectorXi TriTag;
Eigen::VectorXi TetTag;

Eigen::VectorXd D;

std::vector<std::string> XFields;
std::vector<std::string> EFields;

std::vector<Eigen::MatrixXd> XF;
std::vector<Eigen::MatrixXd> TriF;
std::vector<Eigen::MatrixXd> TetF;


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key >= '1' && key <= '9')
  {
    double t = double((key - '1')+1) / 9.0;

    VectorXd v = B.col(2).array() - B.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    vector<int> s;

    for (unsigned i=0; i<v.size();++i)
      if (v(i) < t && v(i)>(t-0.1)) // select a thick slab
        s.push_back(i);

    MatrixXd V_temp(s.size()*4,3);
    MatrixXi F_temp(s.size()*4,3);
    VectorXd D_temp(s.size()*4);

    for (unsigned i=0; i<s.size();++i)
    {
      V_temp.row(i*4+0) = X.row(Tet(s[i],0));
      V_temp.row(i*4+1) = X.row(Tet(s[i],1));
      V_temp.row(i*4+2) = X.row(Tet(s[i],2));
      V_temp.row(i*4+3) = X.row(Tet(s[i],3));

      F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
      F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
      F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
      F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;

      D_temp(i*4+0) = D(s[i]);
      D_temp(i*4+1) = D(s[i]);
      D_temp(i*4+2) = D(s[i]);
      D_temp(i*4+3) = D(s[i]);
    }

    viewer.data().clear();
    viewer.data().set_mesh(V_temp, F_temp);

    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, D_temp, true, C);
    viewer.data().set_face_based(true);
    viewer.data().set_colors(C);

  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  igl::readMSH(argc > 1 ? argv[1] : TUTORIAL_SHARED_PATH "/hand.msh", X, Tri, Tet, TriTag, TetTag, XFields, XF, EFields, TriF, TetF);

  for(auto i:EFields)
    std::cout<<i<<"\t";
  std::cout<<std::endl;

  // search for a predefined field name "E"
  for(int i=0;i<EFields.size();++i)
  {
    if(EFields[i]=="E")
      D = TetF[i].rowwise().norm(); // take a row-wise norm
  }
  std::cout<<"D:"<<D.rows()<<"x"<<D.cols()<<std::endl;

  // generate fake data
  if(D.rows()==0)
    D = TetTag.cast<double>();

  // Compute barycenters
  igl::barycenter(X, Tet, B);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;

  viewer.callback_key_down = &key_down;
  key_down(viewer,'5',0);
  viewer.launch();
}
