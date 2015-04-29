#include "outer_hull.h"
#include "outer_facet.h"
#include "sortrows.h"
#include "facet_components.h"
#include "winding_number.h"
#include "triangle_triangle_adjacency.h"
#include "unique_edge_map.h"
#include "barycenter.h"
#include "per_face_normals.h"

#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <queue>
#include <iostream>
//#define IGL_OUTER_HULL_DEBUG

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedG,
  typename DerivedJ,
  typename Derivedflip>
IGL_INLINE void igl::outer_hull(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedN> & N,
  Eigen::PlainObjectBase<DerivedG> & G,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<Derivedflip> & flip)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedF::Index Index;
  Matrix<Index,DerivedF::RowsAtCompileTime,1> C;
  typedef Matrix<typename DerivedV::Scalar,Dynamic,DerivedV::ColsAtCompileTime> MatrixXV;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,DerivedF::ColsAtCompileTime> MatrixXF;
  typedef Matrix<typename DerivedG::Scalar,Dynamic,DerivedG::ColsAtCompileTime> MatrixXG;
  typedef Matrix<typename DerivedJ::Scalar,Dynamic,DerivedJ::ColsAtCompileTime> MatrixXJ;
  typedef Matrix<typename DerivedN::Scalar,1,3> RowVector3N;
  const Index m = F.rows();

  const auto & duplicate_simplex = [&F](const int f, const int g)->bool
  {
    return 
      (F(f,0) == F(g,0) && F(f,1) == F(g,1) && F(f,2) == F(g,2)) ||
      (F(f,1) == F(g,0) && F(f,2) == F(g,1) && F(f,0) == F(g,2)) ||
      (F(f,2) == F(g,0) && F(f,0) == F(g,1) && F(f,1) == F(g,2)) ||
      (F(f,0) == F(g,2) && F(f,1) == F(g,1) && F(f,2) == F(g,0)) ||
      (F(f,1) == F(g,2) && F(f,2) == F(g,1) && F(f,0) == F(g,0)) ||
      (F(f,2) == F(g,2) && F(f,0) == F(g,1) && F(f,1) == F(g,0));
  };

#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"outer hull..."<<endl;
#endif

#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"edge map..."<<endl;
#endif
  typedef Matrix<typename DerivedF::Scalar,Dynamic,2> MatrixX2I;
  typedef Matrix<typename DerivedF::Index,Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI EMAP;
  vector<vector<typename DerivedF::Index> > uE2E;
  unique_edge_map(F,E,uE,EMAP,uE2E);

  // TODO:
  // uE --> face-edge index, sorted CCW around edge according to normal
  // uE --> sorted order index 
  // uE --> bool, whether needed to flip face to make "consistent" with unique
  //   edge
  // Place order of each half-edge in its corresponding sorted list around edge
  VectorXI diIM(3*m);
  // Whether face's edge used for sorting is consistent with unique edge
  VectorXI dicons(3*m);
  // dihedral angles of faces around edge with face of edge in dicons
  vector<vector<typename DerivedV::Scalar> > di(uE2E.size());
  // For each list of face-edges incide on a unique edge
  for(size_t ui = 0;ui<(size_t)uE.rows();ui++)
  {
    // Base normal vector to orient against
    const auto fe0 = uE2E[ui][0];
    const RowVector3N & eVp = N.row(fe0%m);
    MatrixXd di_I(uE2E[ui].size(),2);

    const typename DerivedF::Scalar d = F(fe0%m,((fe0/m)+2)%3);
    const typename DerivedF::Scalar s = F(fe0%m,((fe0/m)+1)%3);
    // Edge vector
    const auto & eV = (V.row(d)-V.row(s)).normalized();

    vector<bool> cons(uE2E[ui].size());
    // Loop over incident face edges
    for(size_t fei = 0;fei<uE2E[ui].size();fei++)
    {
      const auto & fe = uE2E[ui][fei];
      const auto f = fe % m;
      const auto c = fe / m;
      // source should match destination to be consistent
      cons[fei] = (d == F(f,(c+1)%3));
      assert( cons[fei] ||  (d == F(f,(c+2)%3)));
      assert(!cons[fei] || (s == F(f,(c+2)%3)));
      assert(!cons[fei] || (d == F(f,(c+1)%3)));
      // Angle between n and f
      const RowVector3N & n = N.row(f);
      di_I(fei,0) = M_PI - atan2( eVp.cross(n).dot(eV), eVp.dot(n));
#ifdef IGL_OUTER_HULL_DEBUG
      if(di_I(fei,0) != di_I(fei,0) )
      {
        cout<<"NaN from face: "<<(f+1)<<endl;
        cout<<"  n: "<<n<<endl;
        cout<<"  eVp: "<<eVp<<endl;
        cout<<"  eV: "<<eV<<endl;
        cout<<"  eVp x n . eV: "<<(eVp.cross(n).dot(eV))<<endl;
        cout<<"  eVp . n: "<<(eVp.dot(n))<<endl;
      }
#endif
      assert(di_I(fei,0) == di_I(fei,0) && "NaN Alert!");
      if(!cons[fei])
      {
        di_I(fei,0) = di_I(fei,0) + M_PI;
        if(di_I(fei,0)>=2.*M_PI)
        {
          di_I(fei,0) = di_I(fei,0) - 2.*M_PI;
        }
      }
      // This signing is very important to make sure different edges sort
      // duplicate faces the same way, regardless of their orientations
      di_I(fei,1) = (cons[fei]?1.:-1.)*f;
    }

    // Despite the effort to get stable normals the atan2 up doesn't
    // compute (exactly) -θ for -n if it computes θ for n. So just
    // explicitly check if there's a duplicate face
    // Shitty O(val^2) implementation
    for(size_t fei = 0;fei<uE2E[ui].size();fei++)
    {
      const auto & fe = uE2E[ui][fei];
      const auto f = fe % m;
      for(size_t gei = fei+1;gei<uE2E[ui].size();gei++)
      {
        const auto & ge = uE2E[ui][gei];
        const auto g = ge % m;
        if(duplicate_simplex(f,g))
        {
#ifdef IGL_OUTER_HULL_DEBUG
          cout<<"Forcing duplicate: "<<(f+1)<<","<<(g+1)<<endl;
#endif
          di_I(gei,0) = di_I(fei,0);
        }
      }
    }
    VectorXi IM;
    //igl::sort(di[ui],true,di[ui],IM);
    // Sort, but break ties using "signed index" to ensure that duplicates
    // always show up in same order.
    MatrixXd s_di_I;
    igl::sortrows(di_I,true,s_di_I,IM);
    di[ui].resize(uE2E[ui].size());
    for(size_t i = 0;i<di[ui].size();i++)
    {
      di[ui][i] = s_di_I(i,0);
    }

    // copy old list
    vector<typename DerivedF::Index> temp = uE2E[ui];
    for(size_t fei = 0;fei<uE2E[ui].size();fei++)
    {
      uE2E[ui][fei] = temp[IM(fei)];
      const auto & fe = uE2E[ui][fei];
      diIM(fe) = fei;
      dicons(fe) = cons[IM(fei)];
    }

  }

  vector<vector<vector<Index > > > TT,_1;
  triangle_triangle_adjacency(E,EMAP,uE2E,false,TT,_1);
  VectorXI counts;
#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"facet components..."<<endl;
#endif
  facet_components(TT,C,counts);
  assert(C.maxCoeff()+1 == counts.rows());
  const size_t ncc = counts.rows();
  G.resize(0,F.cols());
  J.resize(0,1);
  flip.setConstant(m,1,false);

#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"reindex..."<<endl;
#endif
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

#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"barycenters..."<<endl;
#endif
  // assumes that "resolve" has handled any coplanar cases correctly and nearly
  // coplanar cases can be sorted based on barycenter.
  MatrixXV BC;
  barycenter(V,F,BC);

#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"loop over CCs (="<<ncc<<")..."<<endl;
#endif
  for(Index id = 0;id<(Index)ncc;id++)
  {
    auto & IM = vIM[id];
    // starting face that's guaranteed to be on the outer hull and in this
    // component
    int f;
    bool f_flip;
#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"outer facet..."<<endl;
#endif
    outer_facet(V,F,N,IM,f,f_flip);
#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"outer facet: "<<f<<endl;
#endif
    int FHcount = 1;
    FH[f] = true;
    // Q contains list of face edges to continue traversing upong
    queue<int> Q;
    Q.push(f+0*m);
    Q.push(f+1*m);
    Q.push(f+2*m);
    flip(f) = f_flip;
    //cout<<"flip("<<f<<") = "<<(flip(f)?"true":"false")<<endl;
#ifdef IGL_OUTER_HULL_DEBUG
  cout<<"BFS..."<<endl;
#endif
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
      // source of edge according to f
      const int fs = flip(f)?F(f,(c+2)%3):F(f,(c+1)%3);
      // destination of edge according to f
      const int fd = flip(f)?F(f,(c+1)%3):F(f,(c+2)%3);
      // edge valence
      const size_t val = uE2E[EMAP(e)].size();
      //// find overlapping face-edges
      //const auto & neighbors = uE2E[EMAP(e)];
      //// normal after possible flipping 
      //const auto & fN = (flip(f)?-1.:1.)*N.row(f);
      //// Edge vector according to f's (flipped) orientation.
      ////const auto & eV = (V.row(fd)-V.row(fs)).normalized();

//#warning "EXPERIMENTAL, DO NOT USE"
      //// THIS IS WRONG! The first face is---after sorting---no longer the face
      //// used for orienting the sort.
      //const auto ui = EMAP(e);
      //const auto fe0 = uE2E[ui][0];
      //const auto es = F(fe0%m,((fe0/m)+1)%3);

      // is edge consistent with edge of face used for sorting
      const int e_cons = (dicons(e) ? 1: -1);
      int nfei = -1;
      // Loop once around trying to find suitable next face
      for(size_t step = 1; step<val+2;step++)
      {
        const int nfei_new = (diIM(e) + 2*val + e_cons*step*(flip(f)?-1:1))%val;
        const int nf = uE2E[EMAP(e)][nfei_new] % m;
        // Don't consider faces with identical dihedral angles
        if((di[EMAP(e)][diIM(e)] != di[EMAP(e)][nfei_new]))
//#warning "THIS IS HACK, FIX ME"
//        if( abs(di[EMAP(e)][diIM(e)] - di[EMAP(e)][nfei_new]) < 1e-16 )
        {
#ifdef IGL_OUTER_HULL_DEBUG
        cout<<"Next facet: "<<(f+1)<<" --> "<<(nf+1)<<", |"<<
          di[EMAP(e)][diIM(e)]<<" - "<<di[EMAP(e)][nfei_new]<<"| = "<<
            abs(di[EMAP(e)][diIM(e)] - di[EMAP(e)][nfei_new])
            <<endl;
#endif
        


          // Only use this face if not already seen
          if(!FH[nf])
          {
            nfei = nfei_new;
          }
          break;
        }
//#ifdef IGL_OUTER_HULL_DEBUG
//        cout<<"Skipping co-planar facet: "<<(f+1)<<" --> "<<(nf+1)<<endl;
//#endif
      }

      int max_ne = -1;
      //// Loop over and find max dihedral angle
      //typename DerivedV::Scalar max_di = -1;
      //for(const auto & ne : neighbors)
      //{
      //  const int nf = ne%m;
      //  if(nf == f)
      //  {
      //    continue;
      //  }
      //  // Corner of neighbor
      //  const int nc = ne/m;
      //  // Is neighbor oriented consistently with (flipped) f?
      //  //const int ns = F(nf,(nc+1)%3);
      //  const int nd = F(nf,(nc+2)%3);
      //  const bool cons = (flip(f)?fd:fs) == nd;
      //  // Normal after possibly flipping to match flip or orientation of f
      //  const auto & nN = (cons? (flip(f)?-1:1.) : (flip(f)?1.:-1.) )*N.row(nf);
      //  // Angle between n and f
      //  const auto & ndi = M_PI - atan2( fN.cross(nN).dot(eV), fN.dot(nN));
      //  if(ndi>=max_di)
      //  {
      //    max_ne = ne;
      //    max_di = ndi;
      //  }
      //}

      ////cout<<(max_ne != max_ne_2)<<" =?= "<<e_cons<<endl;
      //if(max_ne != max_ne_2)
      //{
      //  cout<<(f+1)<<" ---> "<<(max_ne%m)+1<<" != "<<(max_ne_2%m)+1<<" ... "<<e_cons<<" "<<flip(f)<<endl;
      //  typename DerivedV::Scalar max_di = -1;
      //  for(size_t nei = 0;nei<neighbors.size();nei++)
      //  {
      //    const auto & ne = neighbors[nei];
      //    const int nf = ne%m;
      //    if(nf == f)
      //    {
      //      cout<<"  "<<(ne%m)+1<<":\t"<<0<<"\t"<<di[EMAP[e]][nei]<<" "<<diIM(ne)<<endl;
      //      continue;
      //    }
      //    // Corner of neighbor
      //    const int nc = ne/m;
      //    // Is neighbor oriented consistently with (flipped) f?
      //    //const int ns = F(nf,(nc+1)%3);
      //    const int nd = F(nf,(nc+2)%3);
      //    const bool cons = (flip(f)?fd:fs) == nd;
      //    // Normal after possibly flipping to match flip or orientation of f
      //    const auto & nN = (cons? (flip(f)?-1:1.) : (flip(f)?1.:-1.) )*N.row(nf);
      //    // Angle between n and f
      //    const auto & ndi = M_PI - atan2( fN.cross(nN).dot(eV), fN.dot(nN));
      //    cout<<"  "<<(ne%m)+1<<":\t"<<ndi<<"\t"<<di[EMAP[e]][nei]<<" "<<diIM(ne)<<endl;
      //    if(ndi>=max_di)
      //    {
      //      max_ne = ne;
      //      max_di = ndi;
      //    }
      //  }
      //}
      if(nfei >= 0)
      {
        max_ne = uE2E[EMAP(e)][nfei];
      }

      if(max_ne>=0)
      {
        // face of neighbor
        const int nf = max_ne%m;
#ifdef IGL_OUTER_HULL_DEBUG
        if(!FH[nf])
        {
          // first time seeing face
          cout<<(f+1)<<" --> "<<(nf+1)<<endl;
        }
#endif
        FH[nf] = true;
        FHcount++;
        // corner of neighbor
        const int nc = max_ne/m;
        const int nd = F(nf,(nc+2)%3);
        const bool cons = (flip(f)?fd:fs) == nd;
        flip(nf) = (cons ? flip(f) : !flip(f));
        //cout<<"flip("<<nf<<") = "<<(flip(nf)?"true":"false")<<endl;
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
        //  flip(f) = !flip(f);
        //}
        if(FH[f])
        {
          vG[id].row(h) = (flip(f)?F.row(f).reverse().eval():F.row(f));
          vJ[id](h,0) = f;
          h++;
        }
      }
      assert((int)h == FHcount);
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
    MatrixXV q = BC.row(AJ(0));
    // In a perfect world, it's enough to test a single point.
    double w;

    // winding_number_3 expects colmajor
    const typename DerivedV::Scalar * Vdata;
    Vdata = V.data();
    Matrix<
      typename DerivedV::Scalar,
      DerivedV::RowsAtCompileTime,
      DerivedV::ColsAtCompileTime,
      ColMajor> Vcol;
    if(DerivedV::IsRowMajor)
    {
      // copy to convert to colmajor
      Vcol = V;
      Vdata = Vcol.data();
    }
    winding_number_3(
      Vdata,V.rows(),
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
      const bool inside = is_component_inside_other(V,BC,vG[id],vJ[id],vG[oid]);
#ifdef IGL_OUTER_HULL_DEBUG
      cout<<id<<" is inside "<<oid<<" ? "<<inside<<endl;
#endif
      keep[id] = keep[id] && !inside;
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
  Eigen::Matrix<typename DerivedV::Scalar,DerivedF::RowsAtCompileTime,3> N;
  per_face_normals_stable(V,F,N);
  return outer_hull(V,F,N,G,J,flip);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::outer_hull<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
template void igl::outer_hull<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::outer_hull<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#endif
