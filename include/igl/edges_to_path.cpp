#include "edges_to_path.h"
#include "ismember.h"
#include "adjacency_matrix.h"

template <
  typename DerivedE,
  typename DerivedI,
  typename DerivedJ,
  typename DerivedK>
IGL_INLINE void igl::edges_to_path(
  const Eigen::MatrixBase<DerivedE> & OE,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  assert(OE.rows()>=1);
  if(OE.rows() == 1)
  {
    I.resize(2);
    I(0) = OE(0);
    I(1) = OE(1);
    J.resize(1);
    J(0) = 0;
    K.resize(1);
    K(0) = 0;
  }

  // Compute on reduced graph
  DerivedI U;
  DerivedE E;
  {
    Eigen::VectorXi IA;
    unique(OE,U,IA,E);
  }

  Eigen::VectorXi V = Eigen::VectorXi::Zero(E.maxCoeff()+1);
  for(int e = 0;e<E.size();e++)
  {
    V(E(e))++;
    assert(V(E(e))<=2);
  }
  // Try to find a vertex with valence = 1
  int c = 2;
  int s = E(0);
  for(int v = 0;v<V.size();v++)
  {
    if(V(v) == 1)
    {
      c = V(v);
      s = v;
      break;
    }
  }
  assert(V(s) == c);
  assert(c == 2 || c == 1);

  // reshape E to be #E by 2
  E = Eigen::Map<DerivedE>(E.data(),OE.rows(),OE.cols()).eval();
  {
    std::vector<std::vector<int> > A;
    igl::adjacency_list(E,A);
    Eigen::VectorXi P,C;
    dfs(A,s,I,P,C);
  }
  if(c == 2)
  {
    I.conservativeResize(I.size()+1);
    I(I.size()-1) = I(I.size()-2);
  }

  DerivedE sE;
  Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,2> sEI;
  {
    Eigen::MatrixXi _;
    sort(E,2,true,sE,_);
    Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,2> EI(I.size()-1,2);
    EI.col(0) = I.head(I.size()-1);
    EI.col(1) = I.tail(I.size()-1);
    sort(EI,2,true,sEI,_);
  }
  {
    Eigen::Array<bool,Eigen::Dynamic,1> F;
    ismember_rows(sEI,sE,F,J);
  }
  K.resize(I.size()-1);
  for(int k = 0;k<K.size();k++)
  {
    K(k) = 1 + (E(J(k),0) != I(k) ? 1 : 0);
  }

  // Map vertex indices onto original graph
  slice(U,DerivedI(I),1,I);
}
