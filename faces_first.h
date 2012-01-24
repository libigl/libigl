#ifndef IGL_FACES_FIRST_H
#define IGL_FACES_FIRST_H
namespace igl
{
  // FACES_FIRST Reorder vertices so that vertices in face list come before
  // vertices that don't appear in the face list. This is especially useful if
  // the face list contains only surface faces and you want surface vertices
  // listed before internal vertices
  //
  // [RV,RT,RF,IM] = faces_first(V,T,F);
  //
  // Templates:
  //   MatV  matrix for vertex positions, e.g. MatrixXd
  //   MatF  matrix for vertex positions, e.g. MatrixXi
  //   VecI  matrix for vertex positions, e.g. VectorXi
  // Input:
  //  V  # vertices by 3 vertex positions
  //  F  # faces by 3 list of face indices
  // Output: 
  //  RV  # vertices by 3 vertex positions, order such that if the jth vertex is
  //    some face in F, and the kth vertex is not then j comes before k
  //  RF  # faces by 3 list of face indices, reindexed to use RV
  //  IM  # faces by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  //    and RV(IM,:) = V
  //
  template <typename MatV, typename MatF, typename VecI>
  inline void faces_first(
    const MatV & V, 
    const MatF & F, 
    MatV & RV, 
    MatF & RF, 
    VecI & IM);
}

// Implementation
#include <vector>
#include <Eigen/Dense>

template <typename MatV, typename MatF, typename VecI>
inline void igl::faces_first(
  const MatV & V, 
  const MatF & F, 
  MatV & RV, 
  MatF & RF, 
  VecI & IM)
{
  using namespace std;
  using namespace Eigen;
  vector<bool> in_face(V.rows());
  for(int i = 0; i<F.rows(); i++)
  {
    for(int j = 0; j<F.cols(); j++)
    {
      in_face[F(i,j)] = true;
    }
  }
  // count number of vertices not in faces
  int num_in_F = 0;
  for(int i = 0;i<V.rows();i++)
  {
    num_in_F += (in_face[i]?1:0);
  }
  // list of unique vertices that occur in F
  VectorXi U(num_in_F);
  // list of unique vertices that do not occur in F
  VectorXi NU(V.rows()-num_in_F);
  int Ui = 0;
  int NUi = 0;
  // loop over vertices
  for(int i = 0;i<V.rows();i++)
  {
    if(in_face[i])
    {
      U(Ui) = i;
      Ui++;
    }else
    {
      NU(NUi) = i;
      NUi++;
    }
  }
  IM.resize(V.rows());
  // reindex vertices that occur in faces to be first
  for(int i = 0;i<U.size();i++)
  {
    IM(U(i)) = i;
  }
  // reindex vertices that do not occur in faces to come after those that do
  for(int i = 0;i<NU.size();i++)
  {
    IM(NU(i)) = i+U.size();
  }
  RF.resize(F.rows(),F.cols());
  // Reindex faces
  for(int i = 0; i<F.rows(); i++)
  {
    for(int j = 0; j<F.cols(); j++)
    {
      RF(i,j) = IM(F(i,j));
    }
  }
  RV.resize(V.rows(),V.cols());
  // Reorder vertices
  for(int i = 0;i<V.rows();i++)
  {
    RV.row(IM(i)) = V.row(i);
  }
}

#endif
