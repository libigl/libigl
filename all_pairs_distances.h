#ifndef IGL_ALL_PAIRS_DISTANCES_H
#define IGL_ALL_PAIRS_DISTANCES_H

namespace igl
{
  // ALL_PAIRS_DISTANCES compute distances between each point i in V and point j
  // in U
  // 
  // D = all_pairs_distances(V,U)
  // 
  // Templates:
  //   Mat  matrix class like MatrixXd
  // Inputs:
  //   V  #V by dim list of points
  //   U  #U by dim list of points
  //   squared  whether to return squared distances
  // Outputs:
  //   D  #V by #U matrix of distances, where D(i,j) gives the distance or
  //     squareed distance between V(i,:) and U(j,:)
  // 
  template <typename Mat>
  void all_pairs_distances(
    const Mat & V,
    const Mat & U,
    const bool squared, 
    Mat & D);
}

// Implementation

template <typename Mat>
void igl::all_pairs_distances(
  const Mat & V,
  const Mat & U,
  const bool squared,
  Mat & D)
{
  // dimension should be the same
  assert(V.cols() == U.cols());
  // resize output
  D.resize(V.rows(),U.rows());
  for(int i = 0;i<V.rows();i++)
  {
    for(int j=0;j<U.rows();j++)
    {
      D(i,j) = (V.row(i)-U.row(j)).array().pow(2).sum();
      if(!squared)
      {
        D(i,j) = sqrt(D(i,j));
      }
    }
  }
}

#endif
