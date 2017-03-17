#include "crouzeix_raviart_cotmatrix.h"
#include "unique_simplices.h"
#include "all_edges.h"
#include "is_edge_manifold.h"
#include "cotmatrix_entries.h"

template <typename DerivedV, typename DerivedF, typename LT, typename DerivedE, typename DerivedEMAP>
void igl::crouzeix_raviart_cotmatrix(
  const Eigen::MatrixBase<DerivedV> & V, 
  const Eigen::MatrixBase<DerivedF> & F, 
  Eigen::SparseMatrix<LT> & L,
  Eigen::PlainObjectBase<DerivedE> & E,
  Eigen::PlainObjectBase<DerivedEMAP> & EMAP)
{
  // All occurances of directed edges
  Eigen::MatrixXi allE;
  all_edges(F,allE);
  Eigen::VectorXi _1;
  unique_simplices(allE,E,_1,EMAP);
  return crouzeix_raviart_cotmatrix(V,F,E,EMAP,L);
}

template <typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP, typename LT>
void igl::crouzeix_raviart_cotmatrix(
  const Eigen::MatrixBase<DerivedV> & V, 
  const Eigen::MatrixBase<DerivedF> & F, 
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedEMAP> & EMAP,
  Eigen::SparseMatrix<LT> & L)
{
  assert(F.cols() == 3);
  // number of rows
  const int m = F.rows();
  // Mesh should be edge-manifold
  assert(is_edge_manifold(F));
  typedef Eigen::Matrix<LT,Eigen::Dynamic,Eigen::Dynamic> MatrixXS;
  MatrixXS C;
  cotmatrix_entries(V,F,C);
  Eigen::MatrixXi F2E(m,3);
  {
    int k =0;
    for(int c = 0;c<3;c++)
    {
      for(int f = 0;f<m;f++)
      {
        F2E(f,c) = k++;
      }
    }
  }
  std::vector<Eigen::Triplet<LT> > LIJV(12*m);
  Eigen::VectorXi LI(12),LJ(12),LV(12);
  LI<<0,1,2,1,2,0,0,1,2,1,2,0;
  LJ<<1,2,0,0,1,2,0,1,2,1,2,0;
  LV<<2,0,1,2,0,1,2,0,1,2,0,1;

  for(int f=0;f<m;f++)
  {
    for(int c = 0;c<12;c++)
    {
      LIJV.emplace_back(
        EMAP(F2E(f,LI(c))),
        EMAP(F2E(f,LJ(c))),
        (c<6?-1.:1.) * 4. *C(f,LV(c)));
    }
  }
  L.resize(E.rows(),E.rows());
  L.setFromTriplets(LIJV.begin(),LIJV.end());
}
