#include <igl/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>


int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MatrixXi F;
  MatrixXd V;
  // Read in inputs as double precision floating point meshes
  read_triangle_mesh(
      TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  cout<<"Creating grid..."<<endl;
  // number of vertices on the largest side
  const int s = 100;
  // create grid
  MatrixXd GV;
  Eigen::RowVector3i res;
  igl::voxel_grid(V,0,s,1,GV,res);
 
  // compute values
  cout<<"Computing distances..."<<endl;

  // Batch based function for evaluating implicit
  const auto batch_implicit = [&V,&F](const MatrixXd & Q)->VectorXd
  {
    VectorXd S;
    {
      VectorXi I;
      MatrixXd C,N;
      signed_distance(Q,V,F,SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,S,I,C,N);
      // Extremely flatten out near zero
      S = S.array().sign() * S.array().abs().exp();
    }
    return S;
  };

  VectorXd S = batch_implicit(GV);
  cout<<"Marching cubes..."<<endl;
  MatrixXd SV;
  MatrixXi SF;
  igl::marching_cubes(S,GV,res(0),res(1),res(2),0,SV,SF);

  std::unordered_map<std::int64_t,int> E2V;
  igl::marching_cubes(S,GV,res(0),res(1),res(2),0,SV,SF,E2V);

  // Initialize min and max for root finding bisection
  assert(E2V.size() == SV.rows());
  Eigen::MatrixXd T(SV.rows(),2);
  Eigen::VectorXi I(SV.rows());

  // Precompute slices for end points
  Eigen::MatrixXd GVi(SV.rows(),3);
  Eigen::MatrixXd GVj(SV.rows(),3);

  // This is only used for the assertion below
  const auto ij2key = [](std::int32_t i,std::int32_t j)
  {
    if(i>j){ std::swap(i,j); }
    std::int64_t ret = 0;
    ret |= i;
    ret |= static_cast<std::int64_t>(j) << 32;
    return ret;
  };
  const auto key2ij = [](const std::int64_t & key, std::int32_t & i, std::int32_t & j)
  {
    i = key & 0xFFFFFFFF;
    j = key >> 32;
  };
  for (const auto& e2v: E2V)
  {
    const std::int64_t key = e2v.first;
    const std::int32_t v = e2v.second;
    std::int32_t i,j;
    key2ij(key,i,j);
    const std::int64_t key0 = ij2key(i,j);
    assert(key0 != key);
    T.row(v) << 0,1;
    // (i,j) is ordered so that i<j, but let's order so that S(i)<S(j)
    if(S(i)>S(j)) { std::swap(i,j); }
    GVi.row(v) = GV.row(i);
    GVj.row(v) = GV.row(j);
  }

  const auto root_find_iteration = [&SV,&I,&T,&GVi,&GVj,&batch_implicit]()
  {
    Eigen::VectorXd T_mid = (T.col(0)+T.col(1))/2;
    // Use this midpoint guess for the current visualization
    SV = GVi.array().colwise() * (1.0 - T_mid.array()) + GVj.array().colwise() * T_mid.array();
    // Compute values at midpoints
    VectorXd S_mid = batch_implicit(SV);
    // Update bounds
    T.col(1) = (S_mid.array() >  0).select(T_mid, T.col(1));
    T.col(0) = (S_mid.array() <= 0).select(T_mid, T.col(0));
  };

  cout<<R"(Usage:
' '  Conduct a bisection iteration
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(SV,SF);
  viewer.data().show_lines = false;
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case ' ':
          root_find_iteration();
          viewer.data().set_vertices(SV);
          viewer.data().compute_normals();
          break;
      }
      return true;
    };
  viewer.launch();
}
