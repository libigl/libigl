#include "remesh_at_points.h"
#include "STR.h"
#include "write_triangle_mesh.h"
#include "unique_rows.h"
#include "unique_edge_map.h"
#include "triangle/triangulate.h"
#include "PlainMatrix.h"
#include "PlainVector.h"
#include "find.h"
#include "sort.h"
#include <iostream>

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedB,
  typename DerivedFI,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedJ,
  typename DerivedK
  >
IGL_INLINE void igl::remesh_at_points(
  const Eigen::MatrixBase<DerivedV> & V, 
  const Eigen::MatrixBase<DerivedF> & F, 
  const Eigen::MatrixBase<DerivedB> & B,
  const Eigen::MatrixBase<DerivedFI> & FI,
  Eigen::PlainObjectBase<DerivedVV> & VV,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  using BScalar = typename DerivedB::Scalar;
  // static assert that B has 3 columns or Eigen::Dynamic and if Dynamic then
  // runtime assert that it's 3
  static_assert(DerivedB::ColsAtCompileTime == 3 || DerivedB::ColsAtCompileTime == Eigen::Dynamic, "B must have 3 columns");
  assert(B.cols() == 3);

  // 1. Remove any duplicate points in B and FI, track the mapping from old to
  // new indices
  assert(B.minCoeff() >= 0 && B.maxCoeff() <= 1);

  // 2. Compute whether on vertex, edge or face based on counting the number of
  // nonzero barycentric coordinates
  const auto nnz = (B.array()>0).rowwise().count();

  Eigen::VectorXi EMAP, uEC, uEE;
  Eigen::MatrixXi E, uE;
  igl::unique_edge_map(F,E,uE,EMAP,uEC,uEE);

  VV.resize(V.rows() + B.rows(), V.cols());
  VV.topRows(V.rows()) = V;

  std::vector<int> vBI;
  std::vector<int> vFI;
  std::vector<Eigen::Matrix<BScalar,1,3>> vB;
  K.resize(B.rows());
  Eigen::Array<bool,Eigen::Dynamic,1> empty
      = Eigen::Array<bool,Eigen::Dynamic,1>::Constant(F.rows(),true);
  {
    int insert_count = 0;
    for(int i = 0;i<B.rows();i++)
    {
      if(nnz(i) == 1)
      {
        int j;
        B.row(i).maxCoeff(&j);
        K(i) = F(FI(i),j);
        continue; //don't increment insert_count
      }
      else if(nnz(i) == 2)
      {
        // Vertex that has minimum barycentric coordinate is the vertex that is
        // opposite the edge in question.
        int c;
        B.row(i).minCoeff(&c);
        const int f = FI(i);
        const int ue = EMAP(F.rows()*c + f);

        const int a = F(f,(c+1)%3);
        const BScalar Ba = B(i,(c+1)%3);
        const int b = F(f,(c+2)%3);
        const BScalar Bb = B(i,(c+2)%3);
       for(int k = uEC(ue);k<uEC(ue+1);k++)
       {
         const int e = uEE(k); // i = j-uEC(u);
         const int f2 = e % F.rows();
         const int c2 = e / F.rows();
         const int a2 = F(f2,(c2+1)%3);
         const int b2 = F(f2,(c2+2)%3);
         
         BScalar Ba2 = Ba;
         BScalar Bb2 = Bb;
         if((a2 == b && b2 == a))
         {
           // swap Ba2 and Bb2
           std::swap(Ba2,Bb2);
         }else
         {
           assert(a2 == a && b2 == b);
         }
         Eigen::Matrix<BScalar,1,3> B2(0,0,0);
         B2((c2+1)%3) = Ba2;
         B2((c2+2)%3) = Bb2;
         vBI.push_back(i);
         vFI.push_back(f2);
         vB.push_back(B2);
         empty(f2) = false;
       }
      }else if(nnz(i) == 3)
      {
        vBI.push_back(i);
        vFI.push_back(FI(i));
        vB.push_back(B.row(i));
        empty(FI(i)) = false;
      }
      K(i) = V.rows() + insert_count;
      VV.row(K(i)) = 
          B(i,0) * V.row(F(FI(i),0)) +
          B(i,1) * V.row(F(FI(i),1)) +
          B(i,2) * V.row(F(FI(i),2));
      insert_count++;
    }
    VV.conservativeResize(V.rows() + insert_count, V.cols());
  }

  Eigen::VectorXi BI_rep(vBI.size());
  igl::PlainMatrix<DerivedB,Eigen::Dynamic> B_rep(vBI.size(),3);
  igl::PlainVector<DerivedFI,Eigen::Dynamic> FI_rep(vFI.size());
  for(int i = 0;i<vBI.size();i++)
  {
    BI_rep(i) = vBI[i];
    B_rep.row(i) = vB[i];
    FI_rep(i) = vFI[i];
  }

  //Eigen::VectorXi BI = Eigen::VectorXi::LinSpaced(B.rows(),0,B.rows()-1);
  //Eigen::VectorXi BI_rep = BI;
  //igl::PlainMatrix<DerivedB,Eigen::Dynamic> B_rep = B;
  //igl::PlainVector<DerivedFI,Eigen::Dynamic> FI_rep = FI;

  // 3. Handle on vertex cases: just track mapping into original vertices based
  // on F(FI(i),j) where j is the index of the nonzero barycentric coordinate

  // 4. Handle on edge cases into replicated face cases. Use igl::unique_edge_map (uEE,uEC) varient. And
  // replicate every edge case for all incident faces (mapping given barycentric
  // coordinates to each face). Track indices of replication for de-replication
  // later.

  igl::PlainVector<DerivedFI> sFI;
  Eigen::VectorXi sI;
  igl::sort(FI_rep,1,true,sFI,sI);
  igl::PlainMatrix<DerivedB> sB = B_rep(sI, Eigen::placeholders::all);
  Eigen::VectorXi sBI = BI_rep(sI);
  

  // 5. Gather points on each face and retriangulate them independently using
  // igl::delaunay_triangulation or 
  //
  // Assert that on-edge points ended up on edges.
  int start = 0;
  int iter = 0;


  FF.resize(F.rows() + 2*sB.rows(), 3);
  J.resize(FF.rows());
  //printf("guessing that |FF| = %d = %d + 2*%d\n", FF.rows(), F.rows(), B.rows());
  const int empty_count = empty.count();
  const std::vector<int> J_empty = igl::find(empty);
  FF.topRows( empty_count ) = F(J_empty, Eigen::placeholders::all);
  for(int i = 0;i<empty_count;i++) { J(i) = J_empty[i]; }
  int nf = empty_count;
  while(start < sFI.rows())
  {
    int end = start+1;
    while(end < sFI.rows() && sFI(end) == sFI(start) ) { end++; }
    //printf("%d: ",iter);
    //std::cout<<sFI.segment(start, end-start).transpose()<<std::endl;

    igl::PlainMatrix<DerivedB,Eigen::Dynamic> BBi_in(3 + (end-start),3);
    BBi_in.row(0) << 1,0,0;
    BBi_in.row(1) << 0,1,0;
    BBi_in.row(2) << 0,0,1;
    BBi_in.block(3,0,end-start,3) = sB.block(start,0,end-start,3);

    std::vector<Eigen::RowVector2i> E_in;
    for(int z = 0;z<3;z++)
    {
      std::vector<std::pair<BScalar,int>> x_with_index;
      for(int j = 0; j < BBi_in.rows(); j++)
      {
        if(BBi_in(j,z) == 0)
        {
          x_with_index.emplace_back(BBi_in(j,(z+1)%3),j);
        }
      }
      // sort
      std::sort(x_with_index.begin(), x_with_index.end());
      for(int j = 0; j < x_with_index.size()-1; j++)
      {
        Eigen::RowVector2i e(x_with_index[j].second, x_with_index[j+1].second);
        E_in.push_back(e);
      }
    }
    Eigen::MatrixXi E(E_in.size(),2);
    for(int i = 0; i < E_in.size(); i++) { E.row(i) = E_in[i]; }


    igl::PlainMatrix<DerivedB,Eigen::Dynamic> BBi;
    igl::PlainMatrix<DerivedF,Eigen::Dynamic> FFi;
    //igl::write_triangle_mesh(STR("remesh_at_points_input-" << sBI(start) << ".ply"),BBi_in,E);
    igl::triangle::triangulate(BBi_in.leftCols(2).eval(),E,Eigen::MatrixXi(),"pQY",BBi,FFi);
    //Eigen::MatrixXd BBi3(BBi.rows(),3);
    //BBi3<<BBi, (1.0 - BBi.array().rowwise().sum()).eval();
    //igl::write_triangle_mesh(STR("remesh_at_points_output-" << sBI(start) << ".ply"),BBi3,FFi);
    assert(BBi_in.rows() == BBi.rows());

    //printf("  %d\n",FFi.rows());
    if(nf + FFi.rows() > FF.rows())
    {
      FF.conservativeResize(2*FF.rows() + FFi.rows(),3);
      J.conservativeResize(FF.rows());
    }
    for(int i = 0; i < FFi.rows(); i++)
    {
      //FF.row(nf + i) << sI(start + FFi(i,0) - 3), sI(start + FFi(i,1) - 3), sI(start + FFi(i,2) - 3);
      for(int j = 0; j < 3; j++)
      {
        if(FFi(i,j) < 3)
        {
          FF(nf + i,j) = F(sFI(start),FFi(i,j));
        }else
        {
          FF(nf + i,j) = K(sBI(start + FFi(i,j) - 3));
        }
      }
    }
    J.segment(nf, FFi.rows()).setConstant(sFI(start));
    //printf("%d (%d | %d) %d new faces\n", iter,sFI(start),empty[sFI(start)], FFi.rows());
    nf += FFi.rows();

    start = end;
    iter++;
  }
  FF.conservativeResize(nf,3);
  J.conservativeResize(FF.rows());
  //printf("FF.rows() = %d\n",FF.rows());
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::remesh_at_points<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>>&);
#endif
