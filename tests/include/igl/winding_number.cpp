#include <test_common.h>
#include <igl/winding_number.h>
#include <igl/PI.h>

TEST_CASE("winding_number: bunny", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("bunny.off"), V, F);
  Eigen::RowVectorXd V_mean = V.colwise().mean();
  Eigen::MatrixXd P = (V.rowwise() - V_mean) * 0.8;
  P.rowwise() += V_mean;
  Eigen::VectorXd W;
  igl::winding_number(V, F, P, W);
}

TEST_CASE("winding_number: near_planar_node_boundary", "[igl]")
{
  // Closed UV sphere scaled to 4:2:1 (#2579). Its rings put vertices exactly
  // on the symmetry planes, so the hierarchy used by the batch winding_number
  // has nodes whose boundary loops lie in those planes. Query points a hair
  // off these planes lie just outside such a node's box and almost in the
  // plane of its cap. The hierarchical result must still match the direct sum
  // over all faces.
  const int nlat = 16, nlon = 32;
  Eigen::MatrixXd V((nlat-1)*nlon+2,3);
  V.row(0) << 0,0,1;
  for(int i = 1;i<nlat;i++)
  {
    for(int j = 0;j<nlon;j++)
    {
      const double th = igl::PI*i/nlat, ph = 2.*igl::PI*j/nlon;
      Eigen::RowVector3d x(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
      for(int k = 0;k<3;k++){ if(std::abs(x(k))<1e-12){ x(k) = 0; } }
      V.row(1+(i-1)*nlon+j) = x;
    }
  }
  V.row(V.rows()-1) << 0,0,-1;
  V = V*Eigen::Vector3d(2,1,0.5).asDiagonal();
  const auto idx = [&](const int i,const int j){ return 1+i*nlon+(j%nlon); };
  Eigen::MatrixXi F(2*nlon*(nlat-1),3);
  int f = 0;
  for(int j = 0;j<nlon;j++){ F.row(f++) << 0,idx(0,j),idx(0,j+1); }
  for(int i = 0;i<nlat-2;i++)
  {
    for(int j = 0;j<nlon;j++)
    {
      F.row(f++) << idx(i,j),idx(i+1,j),idx(i,j+1);
      F.row(f++) << idx(i,j+1),idx(i+1,j),idx(i+1,j+1);
    }
  }
  const int s = V.rows()-1;
  for(int j = 0;j<nlon;j++){ F.row(f++) << s,idx(nlat-2,j+1),idx(nlat-2,j); }
  REQUIRE(f == F.rows());

  // 5x5 grid in each symmetry plane (including the symmetry lines), moved
  // off the plane by +-1e-16 and +-1e-12.
  const double half[3] = {2,1,0.5};
  const double offsets[4] = {-1e-12,-1e-16,1e-16,1e-12};
  const int g = 5;
  Eigen::MatrixXd P(3*g*g*4,3);
  int r = 0;
  for(int d = 0;d<3;d++)
  {
    for(int a = 0;a<g;a++)
    {
      for(int b = 0;b<g;b++)
      {
        for(const double e : offsets)
        {
          P(r,d) = e;
          P(r,(d+1)%3) = half[(d+1)%3]*(-0.8+1.6*a/(g-1));
          P(r,(d+2)%3) = half[(d+2)%3]*(-0.8+1.6*b/(g-1));
          r++;
        }
      }
    }
  }
  Eigen::VectorXd W;
  igl::winding_number(V,F,P,W);
  double max_err = 0;
  for(int i = 0;i<P.rows();i++)
  {
    const double w = igl::winding_number(V,F,Eigen::RowVector3d(P.row(i)));
    max_err = std::max(max_err,std::abs(W(i)-w));
  }
  REQUIRE(max_err < 1e-10);
}
