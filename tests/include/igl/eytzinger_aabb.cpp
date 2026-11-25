#include <test_common.h>
#include <igl/eytzinger_aabb.h>
#include <random>

TEST_CASE("eytzinger_aabb: tied_centers_preserve_every_primitive", "[igl]")
{
  std::mt19937 generator(836);
  for(int count:{10,128,8192})
  {
    Eigen::MatrixXd minimum(count,3),maximum(count,3),lower,upper;
    for(int i=0;i<count;++i) for(int axis=0;axis<3;++axis)
    {
      minimum(i,axis)=double(generator()%4);
      maximum(i,axis)=minimum(i,axis)+0.25;
    }
    if(count==10)
    {
      const double tied[10]={2,3,1,0,0,0,2,0,0,0};
      minimum.setZero(); maximum.setZero();
      for(int i=0;i<count;++i) minimum(i,0)=maximum(i,0)=tied[i];
    }
    Eigen::VectorXi leaves;
    igl::eytzinger_aabb(minimum,maximum,lower,upper,leaves);
    std::vector<int> visited(count,0);
    for(int node=0;node<leaves.size();++node)
    {
      const int primitive=leaves(node);
      if(primitive>=0)
      {
        REQUIRE(primitive<count);
        ++visited[primitive];
        REQUIRE(lower.row(node)==minimum.row(primitive));
        REQUIRE(upper.row(node)==maximum.row(primitive));
      }else if(primitive==-1)
      {
        REQUIRE(2*node+2<leaves.size());
        REQUIRE(leaves(2*node+1)!=-2);
        REQUIRE(leaves(2*node+2)!=-2);
      }
    }
    for(int occurrences:visited) REQUIRE(occurrences==1);
    std::vector<int> subtree_size(leaves.size(),0);
    for(int node=leaves.size()-1;node>=0;--node)
    {
      if(leaves(node)>=0) subtree_size[node]=1;
      else if(leaves(node)==-1)
      {
        const int left=subtree_size[2*node+1],right=subtree_size[2*node+2];
        REQUIRE(std::abs(left-right)<=1);
        subtree_size[node]=left+right;
      }
    }
  }
}
