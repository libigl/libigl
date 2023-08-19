#include "precomputation.h"

#include <igl/arap.h>
#include <igl/arap_dof.h>
#include <igl/lbs_matrix.h>
#include <igl/columnize.h>
#include <igl/partition.h>
#include <igl/max.h>


void precomputation(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & W,
  Eigen::MatrixXd & M,
  Eigen::VectorXi & b,
  Eigen::MatrixXd & L,
  igl::ARAPData & arap_data,
  igl::ARAPData & arap_grouped_data,
  igl::ArapDOFData<Eigen::MatrixXd,double> & arap_dof_data)
{
  using namespace Eigen;
  using namespace std;
  igl::lbs_matrix_column(V,W,M);

  // Cluster according to weights
  VectorXi G;
  {
    VectorXi S;
    VectorXd D;
    igl::partition(W,50,G,S,D);
  }

  // vertices corresponding to handles (those with maximum weight)
  {
    VectorXd maxW;
    igl::max(W,1,maxW,b);
  }

  // Precomputation for FAST
  cout<<"Initializing Fast Automatic Skinning Transformations..."<<endl;
  // number of weights
  const int m = W.cols();
  Eigen::SparseMatrix<double> Aeq;
  Aeq.resize(m*3,m*3*(3+1));
  vector<Triplet<double> > ijv;
  for(int i = 0;i<m;i++)
  {
    RowVector4d homo;
    homo << V.row(b(i)),1.;
    for(int d = 0;d<3;d++)
    {
      for(int c = 0;c<(3+1);c++)
      {
        ijv.push_back(Triplet<double>(3*i + d,i + c*m*3 + d*m, homo(c)));
      }
    }
  }
  Aeq.setFromTriplets(ijv.begin(),ijv.end());
  igl::arap_dof_precomputation(V,F,M,G,arap_dof_data);
  igl::arap_dof_recomputation(VectorXi(),Aeq,arap_dof_data);
  // Initialize
  MatrixXd Istack = MatrixXd::Identity(3,3+1).replicate(1,m);
  igl::columnize(Istack,m,2,L);

  // Precomputation for ARAP
  cout<<"Initializing ARAP..."<<endl;
  arap_data.max_iter = 1;
  igl::arap_precomputation(V,F,V.cols(),b,arap_data);
  // Grouped arap
  cout<<"Initializing ARAP with grouped edge-sets..."<<endl;
  arap_grouped_data.max_iter = 2;
  arap_grouped_data.G = G;
  igl::arap_precomputation(V,F,V.cols(),b,arap_grouped_data);

}
