#include "boundary_faces.h"

// IGL includes
#include "sort.h"

// STL includes
#include <map>

template <typename IntegerT, typename IntegerF>
IGL_INLINE void igl::boundary_faces(
  const std::vector<std::vector<IntegerT> > & T,
  std::vector<std::vector<IntegerF> > & F)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;

  // Get a list of all faces
  vector<vector<IntegerF> > allF(T.size()*4,vector<IntegerF>(3));
  // Gather faces, loop over tets
  for(int i = 0; i< (int)T.size();i++)
  {
    assert(T[i].size() == 4);
    // get face in correct order
    allF[i*4+0][0] = T[i][0];
    allF[i*4+0][1] = T[i][1];
    allF[i*4+0][2] = T[i][2];
    // get face in correct order
    allF[i*4+1][0] = T[i][0];
    allF[i*4+1][1] = T[i][2];
    allF[i*4+1][2] = T[i][3];
    // get face in correct order
    allF[i*4+2][0] = T[i][0];
    allF[i*4+2][1] = T[i][3];
    allF[i*4+2][2] = T[i][1];
    // get face in correct order
    allF[i*4+3][0] = T[i][1];
    allF[i*4+3][1] = T[i][3];
    allF[i*4+3][2] = T[i][2];
  }
  // Get a list of sorted faces
  vector<vector<IntegerF> > sortedF = allF;
  for(int i = 0; i < (int)allF.size();i++)
  {
    sort(sortedF[i].begin(),sortedF[i].end());
  }
  // Count how many times each sorted face occurs
  map<vector<IntegerF>,int> counts;
  int twos = 0;
  for(int i = 0; i < (int)sortedF.size();i++)
  {
    if(counts.find(sortedF[i]) == counts.end())
    {
      // initialize to count of 1
      counts[sortedF[i]] = 1;
    }else
    {
      // increment count
      counts[sortedF[i]]++;
      assert(counts[sortedF[i]] == 2);
      // increment number of twos
      twos++;
    }
  }

  // Resize output to fit number of ones
  F.resize(allF.size() - twos*2);
  int j = 0;
  for(int i = 0;i< (int)allF.size();i++)
  {
    // sorted face should definitely be in counts map
    assert(counts.find(sortedF[i]) != counts.end());
    if(counts[sortedF[i]] == 1)
    {
      assert(j<(int)F.size());
      F[j] = allF[i];
      j++;
    }
  }
}

#ifndef IGL_NO_EIGEN
#include "list_to_matrix.h"
#include "matrix_to_list.h"

template <typename DerivedT, typename DerivedF>
IGL_INLINE void igl::boundary_faces(
  const Eigen::PlainObjectBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedF>& F)
{
  assert(T.cols() == 0 || T.cols() == 4);
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  // Cop out: use vector of vectors version
  vector<vector<typename Eigen::PlainObjectBase<DerivedT>::Scalar> > vT;
  matrix_to_list(T,vT);
  vector<vector<typename Eigen::PlainObjectBase<DerivedF>::Scalar> > vF;
  boundary_faces(vT,vF);
  list_to_matrix(vF,F);
}
#endif


#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::boundary_faces<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif

