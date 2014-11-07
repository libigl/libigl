#include "outer_hull.h"
#include "outer_facet.h"
#include "facet_components.h"

#include "triangle_triangle_adjacency.h"
#include "unique_edge_map.h"
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
  typedef Matrix<typename DerivedF::Scalar,Dynamic,DerivedF::ColsAtCompileTime> MatrixXF;
  typedef Matrix<typename DerivedG::Scalar,Dynamic,DerivedG::ColsAtCompileTime> MatrixXG;
  typedef Matrix<typename DerivedJ::Scalar,Dynamic,DerivedJ::ColsAtCompileTime> MatrixXJ;
  typedef Matrix<typename Derivedflip::Scalar,Dynamic,Derivedflip::ColsAtCompileTime> MatrixXflip;
  const Index m = F.rows();

  typedef Matrix<typename DerivedF::Scalar,Dynamic,2> MatrixX2I;
  typedef Matrix<typename DerivedF::Index,Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI EMAP;
  vector<vector<typename DerivedF::Index> > uE2E;
  unique_edge_map(F,E,uE,EMAP,uE2E);

  vector<vector<vector<Index > > > TT,_1;
  triangle_triangle_adjacency(E,EMAP,uE2E,false,TT,_1);
  facet_components(TT,C);
  const Index ncc = C.maxCoeff()+1;
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

  // Total size of G (and J)
  size_t nG = 0;
  // This is O( (n+m) * nnc) !
  for(Index id = 0;id<ncc;id++)
  {
    // Determine number of facets in component id
    Index num_id = 0;
    for(Index f = 0;f<m;f++)
    {
      if(C(f) == id)
      {
        num_id++;
      }
    }
    //MatrixXF Fc(num_id,F.cols());
    MatrixXJ IM(num_id,1);
    if(C.maxCoeff() == 0)
    {
      assert(num_id == F.rows());
      //Fc = F;
      IM = MatrixXJ::LinSpaced(F.rows(),0,F.rows()-1);
    }else
    {
      int g = 0;
      for(Index f = 0;f<m;f++)
      {
        if(C(f) == id)
        {
          //Fc.row(g) = F.row(f);
          IM(g) = f;
          g++;
        }
      }
    }
    {
      // starting face that's guaranteed to be on the outer hull and in this
      // component
      int f;
      bool f_flip;
      outer_facet(V,F,N,IM,f,f_flip);
      // Q contains list of face edges to continue traversing upong
      queue<int> Q;
      Q.push(f+0*m);
      Q.push(f+1*m);
      Q.push(f+2*m);
      flip[f] = f_flip;
      int FHcount = 0;
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
        nG += FHcount;
        size_t h = 0;
        for(int i = 0;i<num_id;i++)
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
  }
  // collect G and J across components
  G.resize(nG,3);
  J.resize(nG,1);
  {
    size_t off = 0;
    for(Index id = 0;id<ncc;id++)
    {
      assert(vG[id].rows() == vJ[id].rows());
      G.block(off,0,vG[id].rows(),vG[id].cols()) = vG[id];
      J.block(off,0,vJ[id].rows(),vJ[id].cols()) = vJ[id];
      off += vG[id].rows();
    }
  }
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::outer_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#endif
