#include "in_element.h"

IGL_INLINE void igl::in_element(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & Ele,
  const Eigen::MatrixXd & Q,
  const InElementAABB & aabb,
  Eigen::VectorXi & I)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  const int Qr = Q.rows();
  I.setConstant(Qr,1,-1);
#pragma omp parallel for
  for(int e = 0;e<Qr;e++)
  {
    // find all
    const auto R = aabb.find(V,Ele,Q.row(e),true);
    if(!R.empty())
    {
      I(e) = R[0];
    }
  }
}

IGL_INLINE void igl::in_element(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & Ele,
  const Eigen::MatrixXd & Q,
  const InElementAABB & aabb,
  Eigen::SparseMatrix<double> & I)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  const int Qr = Q.rows();
  std::vector<Triplet<double> > IJV;
  IJV.reserve(Qr);
#pragma omp parallel for
  for(int e = 0;e<Qr;e++)
  {
    // find all
    const auto R = aabb.find(V,Ele,Q.row(e),false);
    for(const auto r : R)
    {
#pragma omp critical
      IJV.push_back(Triplet<double>(e,r,1));
    }
  }
  I.resize(Qr,Ele.rows());
  I.setFromTriplets(IJV.begin(),IJV.end());
}
