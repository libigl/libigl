#define IGL_HEADER_ONLY
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/comiso/miq.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <sstream>


void line_texture(Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
  {
    unsigned size = 128;
    unsigned size2 = size/2;
    unsigned lineWidth = 3;
    texture_R.setConstant(size, size, 255);
    for (unsigned i=0; i<size; ++i)
      for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
        texture_R(i,j) = 0;
    for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
      for (unsigned j=0; j<size; ++j)
        texture_R(i,j) = 0;

    texture_G = texture_R;
    texture_B = texture_R;
  }


bool readPolyVf(const char *fname,
                Eigen::VectorXi &isConstrained,
                std::vector<Eigen::MatrixXd> &polyVF)
{
  FILE *fp = fopen(fname,"r");
  if (!fp)
    return false;
  int degree, numF;
  if (fscanf(fp,"%d %d", &degree, &numF) !=2)
    return false;
  polyVF.resize(degree, Eigen::MatrixXd::Zero(numF, 3));
  isConstrained.setZero(numF,1);
  int vali; float u0,u1,u2;
  for (int i = 0; i<numF; ++i)
  {
    if (fscanf(fp,"%d", &vali)!=1)
      return false;
    isConstrained[i] = vali;
    for (int j = 0; j<degree; ++j)
    {
      if (fscanf(fp,"%g %g %g", &u0, &u1, &u2) !=3)
        return false;
      polyVF[j](i,0) = u0;
      polyVF[j](i,1) = u1;
      polyVF[j](i,2) = u2;
    }
  }
  fclose(fp);
  return true;
}

void writePolyVf(const char *fname,
                 const Eigen::VectorXi &isConstrained,
                 const std::vector<Eigen::MatrixXd> &polyVF)
{
  int numF = polyVF[0].rows();
  int degree = polyVF.size();
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"%d %d\n", degree,numF);
  for (int i = 0; i<numF; ++i)
  {
    fprintf(fp,"%d ", isConstrained[i]);
    for (int j = 0; j<degree; ++j)
      fprintf(fp,"%.15g %.15g %.15g ", polyVF[j](i,0), polyVF[j](i,1), polyVF[j](i,2));
    fprintf(fp, "\n");
  }
  fclose(fp);

}


int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Load a mesh in OFF format
  igl::readOBJ("../shared/lilium.obj", V, F);


  // Load a frame field
  Eigen::VectorXi isConstrained;
  std::vector<Eigen::MatrixXd> polyVF;
  readPolyVf("../shared/lilium.crossfield", isConstrained, polyVF);

  Eigen::MatrixXd UV;
  Eigen::MatrixXi FUV;

  double gradientSize = 50;
  double quadIter = 0;
  double stiffness = 5.0;
  bool directRound = 1;
  igl::miq(V,
           F,
           polyVF[0],
           polyVF[1],
           UV,
           FUV,
           gradientSize,
           stiffness,
           directRound,
           quadIter);


  // Face barycenters
  Eigen::MatrixXd MF;
  igl::barycenter(V, F, MF);

  double scale =  .5*igl::avg_edge_length(V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Plot the field
  viewer.add_edges (MF, MF+scale*polyVF[0],Eigen::RowVector3d(1,0,1));
  viewer.add_edges (MF, MF+scale*polyVF[1],Eigen::RowVector3d(1,0,1));
  viewer.set_uv(UV,FUV);
  viewer.options.show_texture = true;

  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;
  line_texture(texture_R, texture_G, texture_B);
  viewer.set_texture(texture_R, texture_B, texture_G);
  // Increase the thickness of the lines
  viewer.options.line_width = 2.0f;

  // Launch the viewer
  viewer.launch();
}
