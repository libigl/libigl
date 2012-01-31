#ifndef IGL_LIMIT_FACES_H
#define IGL_LIMIT_FACES_H
namespace igl
{
  // LIMIT_FACES limit given faces F to those which contain (only) indices found
  // in L.
  //
  // [LF] = limit_faces(F,L,exclusive);
  // [LF,in] = limit_faces(F,L,exclusive);
  //
  // Templates:
  //   MatF matrix type of faces, matrixXi
  //   VecL  matrix type of vertex indices, VectorXi
  // Inputs:
  //   F  #F by 3 list of face indices
  //   L  #L by 1 list of allowed indices
  //   exclusive  flag specifying whether a face is included only if all its
  //     indices are in L, default is false
  // Outputs:
  //   LF  #LF by 3 list of remaining faces after limiting
  //   in  #F list of whether given face was included
  //
  template <typename MatF, typename VecL>
  inline void limit_faces(
    const MatF & F, 
    const VecL & L, 
    const bool exclusive,
    MatF & LF);
}

// Implementation
#include <vector>

template <typename MatF, typename VecL>
inline void igl::limit_faces(
  const MatF & F, 
  const VecL & L, 
  const bool exclusive,
  MatF & LF)
{
  using namespace std;
  using namespace Eigen;
  vector<bool> in(F.rows(),false);
  int num_in = 0;
  // loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    bool all = true;
    bool any = false;
    for(int j = 0;j<F.cols();j++)
    {
      bool found = false;
      // loop over L
      for(int l = 0;l<L.size();l++)
      {
        if(F(i,j) == L(l))
        {
          found = true;
          break;
        }
      }
      any |= found;
      all &= found;
    }
    in[i] = (exclusive?all:any);
    num_in += (in[i]?1:0);
  }

  LF.resize(num_in,F.cols());
  // loop over faces
  int lfi = 0;
  for(int i = 0;i<F.rows();i++)
  {
    if(in[i])
    {
      LF.row(lfi) = F.row(i);
      lfi++;
    }
  }
}

#endif
