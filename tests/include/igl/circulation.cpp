#include <test_common.h>
#include <igl/circulation.h>
#include <igl/edge_flaps.h>
#include <igl/unique_edge_map.h>
#include <igl/matlab_format.h>

TEST(circulation,single_edge)
{
  //       7    
  //     /₆|₇\
  //   4 - 5 - 6
  //   |₂/₃|₄\₅|
  //   1 - 2 - 3
  //     \₀|₁/  
  //       0     
  const Eigen::MatrixXi F = (Eigen::MatrixXi(8,3)<<
    0,2,1,
    0,3,2,
    1,5,4,
    1,2,5,
    2,3,5,
    3,6,5,
    4,5,7,
    5,6,7).finished();
  Eigen::MatrixXi E,uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<int> > uE2E;
  igl::unique_edge_map(F, E, uE, EMAP, uE2E);
  Eigen::MatrixXi EI,EF;
  {
    const auto & cuE = uE;
    const auto & cEMAP = EMAP;
    igl::edge_flaps(F,cuE,cEMAP,EF,EI);
  }
  // Find (2,5) in uE
  int ei = 0;
  for(;ei<E.rows();ei++)
  {
    if(uE(ei,0) == 2 && uE(ei,1) == 5) break;
    if(uE(ei,1) == 2 && uE(ei,0) == 5) break;
  }
  Eigen::VectorXi Nccw;
  igl::circulation(ei,true,EMAP,EF,EI,Nccw);
  const Eigen::VectorXi Nccwgt = 
    (Eigen::VectorXi(6)<<
     4,
     5,
     7,
     6,
     2,
     3).finished();
  test_common::assert_eq(Nccw,Nccwgt);
  Eigen::VectorXi Ncw;
  igl::circulation(ei,false,EMAP,EF,EI,Ncw);
  const Eigen::VectorXi Ncwgt = 
    (Eigen::VectorXi(4)<<
     4,
     1,
     0,
     3).finished();
  test_common::assert_eq(Ncw,Ncwgt);
}
