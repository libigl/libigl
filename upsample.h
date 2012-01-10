#ifndef IGL_UPSAMPLE_H
#define IGL_UPSAMPLE_H
namespace igl
{
  // Subdivide a mesh without moving vertices: loop subdivision but odd
  // vertices stay put and even vertices are just edge midpoints
  // 
  // Templates:
  //   MatV  matrix for vertex positions, e.g. MatrixXd
  //   MatF  matrix for vertex positions, e.g. MatrixXi
  // Inputs:
  //   V  #V by dim  mesh vertices
  //   F  #F by 3  mesh triangles
  // Outputs:
  //   NV new vertex positions, V is guaranteed to be at top
  //   NF new list of face indices
  //
  // NOTE: V should not be the same as NV,
  // NOTE: F should not be the same as NF, use other proto
  template <typename MatV, typename MatF>
  void upsample( const MatV & V, const MatF & F, MatV & NV, MatF & NF);
  // Virtually in place wrapper
  template <typename MatV, typename MatF>
  void upsample( MatV & V,MatF & F);
}

// Implementation
#include <tt.h>
#include <adjacency_list.h>
#include <Eigen/Dense>

template <typename MatV, typename MatF>
void igl::upsample( const MatV & V, const MatF & F, MatV & NV, MatF & NF)
{
  // Use "in place" wrapper instead
  assert(&V != &NV);
  assert(&F != &NF);
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  MatF FF, FFi;
  tt<double>(V,F,FF,FFi);

  // TODO: Cache optimization missing from here, it is a mess
  
  // Compute the number and positions of the vertices to insert (on edges)
  MatF NI = MatF::Constant(FF.rows(),FF.cols(),-1);
  int counter = 0;
  
  for(int i=0;i<FF.rows();++i)
  {
    for(int j=0;j<3;++j)
    {
      if(NI(i,j) == -1)
      {
        NI(i,j) = counter;
        if (FF(i,j) != -1) // If it is not a border
          NI(FF(i,j),FFi(i,j)) = counter;
        ++counter;
      }
    }
  }

  int n_odd = V.rows();
  int n_even = counter;

  Eigen::DynamicSparseMatrix<double> SUBD(V.rows()+n_even,V.rows());
  SUBD.reserve(15 * (V.rows()+n_even));
  
  // Preallocate NV and NF
  NV = MatV(V.rows()+n_even,V.cols());
  NF = MatF(F.rows()*4,3);
  
  // Fill the odd vertices position
  NV.block(0,0,V.rows(),V.cols()) = V;

  // Fill the even vertices position
  for(int i=0;i<FF.rows();++i)
  {
    for(int j=0;j<3;++j)
    {
      NV.row(NI(i,j) + n_odd) = 0.5 * V.row(F(i,j)) + 0.5 * V.row(F(i,(j+1)%3));
    }
  }

  // Build the new topology (Every face is replaced by four)
  for(int i=0; i<F.rows();++i)
  {
    VectorXi VI(6);
    VI << F(i,0), F(i,1), F(i,2), NI(i,0) + n_odd, NI(i,1) + n_odd, NI(i,2) + n_odd;
    
    VectorXi f0(3), f1(3), f2(3), f3(3);
    f0 << VI(0), VI(3), VI(5);
    f1 << VI(1), VI(4), VI(3);
    f2 << VI(3), VI(4), VI(5);
    f3 << VI(4), VI(2), VI(5);
    
    NF.row((i*4)+0) = f0;
    NF.row((i*4)+1) = f1;
    NF.row((i*4)+2) = f2;
    NF.row((i*4)+3) = f3;
  }
  
}

template <typename MatV, typename MatF>
void igl::upsample( MatV & V,MatF & F)
{
  const MatV V_copy = V;
  const MatF F_copy = F;
  return upsample(V_copy,F_copy,V,F);
}

#endif
