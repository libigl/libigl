#include "update.h"
#include <igl/heat_geodesics.h>
#include <igl/unproject_onto_mesh.h>

bool update(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double t,
  const double x,
  const double y,
  const Eigen::Matrix4f& view,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::HeatGeodesicsData<double>& data,
  Eigen::VectorXd& D)
{
  int fid;
  Eigen::Vector3f bc;
  // Cast a ray in the view direction starting from the mouse position
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), view,
    proj, viewport, V, F, fid, bc))
  {
    // if big mesh, just use closest vertex. Otherwise, blend distances to
    // vertices of face using barycentric coordinates.
    if(F.rows()>100000)
    {
      // 3d position of hit
      const Eigen::RowVector3d m3 =
        V.row(F(fid,0))*bc(0) + V.row(F(fid,1))*bc(1) + V.row(F(fid,2))*bc(2);
      int cid = 0;
      Eigen::Vector3d(
          (V.row(F(fid,0))-m3).squaredNorm(),
          (V.row(F(fid,1))-m3).squaredNorm(),
          (V.row(F(fid,2))-m3).squaredNorm()).minCoeff(&cid);
      const int vid = F(fid,cid);
      igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<vid).finished(),D);
    }else
    {
      D = Eigen::VectorXd::Zero(V.rows());
      for(int cid = 0;cid<3;cid++)
      {
        const int vid = F(fid,cid);
        Eigen::VectorXd Dc;
        igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<vid).finished(),Dc);
        D += Dc*bc(cid);
      }
    }
    return true;
  }
  return false;
}
