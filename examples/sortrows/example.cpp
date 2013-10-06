#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <igl/sort.h>
#include <igl/unique.h>
#include <igl/get_seconds.h>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/sortrows.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  //MatrixXd X;
  //if(!readDMAT("X.dmat",X))
  //{
  //  X = MatrixXd::Random(766443,2);
  //  for(int x = 0;x<X.size();x++)
  //  {
  //    X(x) = floor(X(x)*680);
  //  }
  //}
  MatrixXd X(4,2);
  //X << 1,2,1,2,1,1,2,2;
  X << 1,3,2,4,6,1,7,2;
  cout<<X<<endl<<endl;
  typedef MatrixXi::Scalar Scalar;
  bool ascending = false;
  const int dim = 2;
  double t;

//#define NUM_RUNS 100
//  MatrixXd Y;
//  MatrixXi IX;
//  t = get_seconds();
//  for(int r = 0;r<NUM_RUNS;r++)
//  sort(X,dim,ascending,Y,IX);
//  cout<<(get_seconds()-t)/(double)NUM_RUNS<<endl;
//  //cout<<Y<<endl<<endl;
//  //cout<<IX<<endl<<endl;
//  MatrixXd Y2;
//  MatrixXi IX2;
//  t = get_seconds();
//  for(int r = 0;r<NUM_RUNS;r++)
//  sort2(X,dim,ascending,Y2,IX2);
//  cout<<(get_seconds()-t)/(double)NUM_RUNS<<endl;
//  //cout<<Y2<<endl<<endl;
//  //cout<<IX2<<endl<<endl;

#define NUM_RUNS 1
  MatrixXd Y;
  VectorXi IA,IC;

  t = get_seconds();
  for(int r = 0;r<NUM_RUNS;r++)
  unique_rows(X,Y,IA,IC);
  cout<<(get_seconds()-t)/(double)NUM_RUNS<<endl;

  t = get_seconds();
  for(int r = 0;r<NUM_RUNS;r++)
  sortrows(X,true,Y,IC);
  cout<<(get_seconds()-t)/(double)NUM_RUNS<<endl;
  writeDMAT("X.dmat",X,true);
  writeDMAT("Y.dmat",Y,true);
  writeDMAT("IC.dmat",IC,true);

  //cout<<Y<<endl<<endl;
  //cout<<IA<<endl<<endl;
  //cout<<IC<<endl<<endl;
  cout<<Y.rows()<<endl;

}
