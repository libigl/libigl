#include "outer_hull.h"
#include "outer_facet.h"
#include "facet_components.h"
#include "winding_number.h"
#include "triangle_triangle_adjacency.h"
#include "unique_edge_map.h"
#include "barycenter.h"
#include "per_face_normals.h"
#include "all_edges.h"
#include "colon.h"
#include "get_seconds.h"

#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <queue>
#include <iostream>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedG,
  typename DerivedJ,
  typename Derivedflip>
IGL_INLINE void igl::outer_hull(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedG> & G,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<Derivedflip> & flip)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  typedef typename DerivedF::Index Index;
  Matrix<Index,DerivedF::RowsAtCompileTime,1> C;
  typedef Matrix<typename DerivedV::Scalar,Dynamic,DerivedV::ColsAtCompileTime> MatrixXV;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,DerivedF::ColsAtCompileTime> MatrixXF;
  typedef Matrix<typename DerivedG::Scalar,Dynamic,DerivedG::ColsAtCompileTime> MatrixXG;
  typedef Matrix<typename DerivedJ::Scalar,Dynamic,DerivedJ::ColsAtCompileTime> MatrixXJ;
  const Index m = F.rows();

  typedef Matrix<typename DerivedF::Scalar,Dynamic,2> MatrixX2I;
  typedef Matrix<typename DerivedF::Index,Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI EMAP;
  vector<vector<typename DerivedF::Index> > uE2E;
  unique_edge_map(F,E,uE,EMAP,uE2E);

  vector<vector<vector<Index > > > TT,_1;
  triangle_triangle_adjacency(E,EMAP,uE2E,false,TT,_1);
  VectorXI counts;
  facet_components(TT,C,counts);
  assert(C.maxCoeff()+1 == counts.rows());
  const size_t ncc = counts.rows();
  G.resize(0,F.cols());
  J.resize(0,1);
  flip.setConstant(m,1,false);
  // precompute face normals
  typename Eigen::Matrix<
    typename DerivedV::Scalar,
    DerivedF::RowsAtCompileTime,
    3> N;
  per_face_normals(V,F,N);
  // H contains list of faces on outer hull;
  vector<bool> FH(m,false);
  vector<bool> EH(3*m,false);
  vector<MatrixXG> vG(ncc);
  vector<MatrixXJ> vJ(ncc);
  vector<MatrixXJ> vIM(ncc);
  for(size_t id = 0;id<ncc;id++)
  {
    vIM[id].resize(counts[id],1);
  }
  // current index into each IM
  vector<size_t> g(ncc,0);
  // place order of each face in its respective component
  for(Index f = 0;f<m;f++)
  {
    vIM[C(f)](g[C(f)]++) = f;
  }

  // assumes that "resolve" has handled any coplanar cases correctly and nearly
  // coplanar cases can be sorted based on barycenter.
  MatrixXV BC;
  barycenter(V,F,BC);

  for(Index id = 0;id<(Index)ncc;id++)
  {
    auto & IM = vIM[id];
    // starting face that's guaranteed to be on the outer hull and in this
    // component
    int f;
    bool f_flip;
    outer_facet(V,F,N,IM,f,f_flip);
    int FHcount = 0;
    // Q contains list of face edges to continue traversing upong
    queue<int> Q;
    Q.push(f+0*m);
    Q.push(f+1*m);
    Q.push(f+2*m);
    flip[f] = f_flip;
    while(!Q.empty())
    {
      // face-edge
      const int e = Q.front();
      Q.pop();
      // face
      const int f = e%m;
      // corner
      const int c = e/m;
      // Should never see edge again...
      if(EH[e] == true)
      {
        continue;
      }
      EH[e] = true;
      // first time seeing face
      if(!FH[f])
      {
        FH[f] = true;
        FHcount++;
      }
      // find overlapping face-edges
      const auto & neighbors = uE2E[EMAP(e)];
      const auto & fN = (flip[f]?-1.:1.)*N.row(f);
      // source of edge according to f
      const int fs = flip[f]?F(f,(c+2)%3):F(f,(c+1)%3);
      // destination of edge according to f
      const int fd = flip[f]?F(f,(c+1)%3):F(f,(c+2)%3);
      const auto & eV = (V.row(fd)-V.row(fs)).normalized();
      // Loop over and find max dihedral angle
      typename DerivedV::Scalar max_di = -1;
      int max_ne = -1;
      typename Eigen::Matrix< typename DerivedV::Scalar, 1, 3> max_nN;
      for(const auto & ne : neighbors)
      {
        const int nf = ne%m;
        if(nf == f)
        {
          continue;
        }
        const int nc = ne/m;
        // are faces consistently oriented
        //const int ns = F(nf,(nc+1)%3);
        const int nd = F(nf,(nc+2)%3);
        const bool cons = (flip[f]?fd:fs) == nd;
        const auto & nN = (cons? (flip[f]?-1:1.) : (flip[f]?1.:-1.) )*N.row(nf);
        const auto & ndi = M_PI - atan2( fN.cross(nN).dot(eV), fN.dot(nN));
        if(ndi>=max_di)
        {
          max_ne = ne;
          max_di = ndi;
          max_nN = nN;
        }
      }
      if(max_ne>=0)
      {
        const int nf = max_ne%m;
        const int nc = max_ne/m;
        const int nd = F(nf,(nc+2)%3);
        const bool cons = (flip[f]?fd:fs) == nd;
        flip[nf] = (cons ? flip[f] : !flip[f]);
        const int ne1 = nf+((nc+1)%3)*m;
        const int ne2 = nf+((nc+2)%3)*m;
        if(!EH[ne1])
        {
          Q.push(ne1);
        }
        if(!EH[ne2])
        {
          Q.push(ne2);
        }
      }
    }
    
    {
      vG[id].resize(FHcount,3);
      vJ[id].resize(FHcount,1);
      //nG += FHcount;
      size_t h = 0;
      assert(counts(id) == IM.rows());
      for(int i = 0;i<counts(id);i++)
      {
        const size_t f = IM(i);
        //if(f_flip)
        //{
        //  flip[f] = !flip[f];
        //}
        if(FH[f])
        {
          vG[id].row(h) = (flip[f]?F.row(f).reverse().eval():F.row(f));
          vJ[id](h,0) = f;
          h++;
        }
      }
      assert(h == FHcount);
    }
  }

  // Is A inside B? Assuming A and B are consistently oriented but closed and
  // non-intersecting.
  const auto & is_component_inside_other = [](
    const Eigen::PlainObjectBase<DerivedV> & V,
    const MatrixXV & BC,
    const MatrixXG & A,
    const MatrixXJ & AJ,
    const MatrixXG & B)->bool
  {
    const auto & bounding_box = [](
      const Eigen::PlainObjectBase<DerivedV> & V,
      const MatrixXG & F)->
      MatrixXV
    {
      MatrixXV BB(2,3);
      BB<<
         1e26,1e26,1e26,
        -1e26,-1e26,-1e26;
      const size_t m = F.rows();
      for(size_t f = 0;f<m;f++)
      {
        for(size_t c = 0;c<3;c++)
        {
          const auto & vfc = V.row(F(f,c));
          BB.row(0) = BB.row(0).array().min(vfc.array()).eval();
          BB.row(1) = BB.row(1).array().max(vfc.array()).eval();
        }
      }
      return BB;
    };
    // A lot of the time we're dealing with unrelated, distant components: cull
    // them.
    MatrixXV ABB = bounding_box(V,A);
    MatrixXV BBB = bounding_box(V,B);
    if( (BBB.row(0)-ABB.row(1)).maxCoeff()>0  ||
        (ABB.row(0)-BBB.row(1)).maxCoeff()>0 )
    {
      // bounding boxes do not overlap
      return false;
    }
    ////////////////////////////////////////////////////////////////////////
    // POTENTIAL ROBUSTNESS WEAK AREA
    ////////////////////////////////////////////////////////////////////////
    //
    // q could be so close (<~1e-16) to B that the winding number is not a robust way to
    // determine inside/outsideness. We could try to find a _better_ q which is
    // farther away, but couldn't they all be bad?
    MatrixXV q = BC.row(AJ(1));
    // In a perfect world, it's enough to test a single point.
    double w;
    winding_number_3(
      V.data(),V.rows(),
      B.data(),B.rows(),
      q.data(),1,&w);
    return fabs(w)>0.5;
  };

  // Reject components which are completely inside other components
  vector<bool> keep(ncc,true);
  size_t nG = 0;
  // This is O( ncc * ncc * m)
  for(size_t id = 0;id<ncc;id++)
  {
    for(size_t oid = 0;oid<ncc;oid++)
    {
      if(id == oid)
      {
        continue;
      }
      keep[id] = keep[id] && 
        !is_component_inside_other(V,BC,vG[id],vJ[id],vG[oid]);
    }
    if(keep[id])
    {
      nG += vJ[id].rows();
    }
  }

  // collect G and J across components
  G.resize(nG,3);
  J.resize(nG,1);
  {
    size_t off = 0;
    for(Index id = 0;id<(Index)ncc;id++)
    {
      if(keep[id])
      {
        assert(vG[id].rows() == vJ[id].rows());
        G.block(off,0,vG[id].rows(),vG[id].cols()) = vG[id];
        J.block(off,0,vJ[id].rows(),vJ[id].cols()) = vJ[id];
        off += vG[id].rows();
      }
    }
  }
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::outer_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#endif
