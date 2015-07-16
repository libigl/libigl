#include <igl/colon.h>
#include <algorithm>
#include <deque>
#include <igl/polyvector_field_comb_from_matchings_and_cuts.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>

template <typename DerivedV, typename DerivedF, typename DerivedTT, typename DerivedS, typename DerivedM, typename DerivedC>
IGL_INLINE void igl::polyvector_field_comb_from_matchings_and_cuts(
  const Eigen::PlainObjectBase<DerivedV> &V,
  const Eigen::PlainObjectBase<DerivedF> &F,
  const Eigen::PlainObjectBase<DerivedTT> &TT,
  const Eigen::MatrixXi &E2F,
  const Eigen::MatrixXi &F2E,
  const Eigen::PlainObjectBase<DerivedS> &sol3D,
  const Eigen::PlainObjectBase<DerivedM> &match_ab,
  const Eigen::PlainObjectBase<DerivedM> &match_ba,
  const Eigen::PlainObjectBase<DerivedC> &cuts,
  Eigen::PlainObjectBase<DerivedS> &sol3D_combed)
{

  int half_degree = sol3D.cols()/3;
  int full_degree = 2*half_degree;

  Eigen::MatrixXi used; used.setConstant(F.rows(),half_degree,-1);

  //  Eigen::VectorXi mark;
  std::deque<int> d;
  sol3D_combed.setZero(sol3D.rows(), sol3D.cols());

  int start = 0;
  d.push_back(start);
  used.row(start) = igl::colon<int>(0, half_degree-1);
  sol3D_combed.row(start) = sol3D.row(start);
  while (!d.empty())
  {
    int f0 = d.at(0);
    d.pop_front();
    for (int k=0; k<3; k++)
    {
      int f1 = TT(f0,k);
      if (f1==-1) continue;
      //do not comb across cuts
      if ((used.row(f1).array()!=-1).any()||cuts(f0,k))
        continue;


      // look at the edge between the two faces
      const int &current_edge = F2E(f0,k);

      // check its orientation to determine whether match_ab or match_ba should be used
      if ((E2F(current_edge,0) == f0) &&
          (E2F(current_edge,1) == f1) )
      {
        //look at match_ab
        for(int i=0; i<half_degree; ++i)
        {
          int ii = used(f0,i);
          used(f1,i) = (match_ab(current_edge,ii%half_degree) + (ii>=half_degree)*half_degree)%full_degree;
        }
      }
      else
      {
        assert((E2F(current_edge,1) == f0) &&
               (E2F(current_edge,0) == f1));
        //look at match_ba
        for(int i=0; i<half_degree; ++i)
        {
          int ii = used(f0,i);
          used(f1,i) = (match_ba(current_edge,ii%half_degree)+ (ii>=half_degree)*half_degree)%full_degree;
        }
      }

      for (int i = 0; i<half_degree; ++i)
      {
        int sign = (used(f1,i)<half_degree)?1:-1;
        sol3D_combed.block(f1,i*3,1,3) = sign*sol3D.block(f1,(used(f1,i)%half_degree)*3,1,3);
      }
      d.push_back(f1);
    }
  }

  // everything should be marked
  assert((used.rowwise().minCoeff().array()>=0).all());
}


template <typename DerivedV, typename DerivedF, typename DerivedS, typename DerivedM, typename DerivedC>
IGL_INLINE void igl::polyvector_field_comb_from_matchings_and_cuts(
  const Eigen::PlainObjectBase<DerivedV> &V,
  const Eigen::PlainObjectBase<DerivedF> &F,
  const Eigen::PlainObjectBase<DerivedS> &sol3D,
  const Eigen::PlainObjectBase<DerivedM> &match_ab,
  const Eigen::PlainObjectBase<DerivedM> &match_ba,
  const Eigen::PlainObjectBase<DerivedC> &cuts,
  Eigen::PlainObjectBase<DerivedS> &sol3D_combed)
  {
    Eigen::MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(V,F,TT,TTi);

    Eigen::MatrixXi E, E2F, F2E;
    igl::edge_topology(V,F,E,F2E,E2F);

    igl::polyvector_field_comb_from_matchings_and_cuts(V, F, TT, E2F, F2E, sol3D, match_ab, match_ba, cuts, sol3D_combed);
  }


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::polyvector_field_comb_from_matchings_and_cuts<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::polyvector_field_comb_from_matchings_and_cuts<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
