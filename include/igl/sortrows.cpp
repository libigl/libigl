#include "sortrows.h"

#include "SortableRow.h"
#include "sort.h"

#include <vector>

template <typename DerivedX, typename DerivedIX>
IGL_INLINE void igl::sortrows(
  const Eigen::PlainObjectBase<DerivedX>& X,
  const bool ascending,
  Eigen::PlainObjectBase<DerivedX>& Y,
  Eigen::PlainObjectBase<DerivedIX>& IX)
{
  using namespace std;
  using namespace Eigen;
  typedef Eigen::Matrix<typename DerivedX::Scalar, Eigen::Dynamic, 1> RowVector;
  vector<SortableRow<RowVector> > rows;
  rows.resize(X.rows());
  // Loop over rows
  for(int i = 0;i<X.rows();i++)
  {
    RowVector ri = X.row(i);
    rows[i] = SortableRow<RowVector>(ri);
  }
  vector<SortableRow<RowVector> > sorted;
  std::vector<size_t> index_map;
  // Perform sort on rows
  igl::sort(rows,ascending,sorted,index_map);
  // Resize output
  Y.resize(X.rows(),X.cols());
  IX.resize(X.rows(),1);
  // Convert to eigen
  for(int i = 0;i<X.rows();i++)
  {
    Y.row(i) = sorted[i].data;
    IX(i) = index_map[i];
  }
}

#ifndef IGL_HEADER_ONLY
template void igl::sortrows<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sortrows<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
