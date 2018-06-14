#include "build_octree.h"
#include <vector>
#include <queue>


namespace igl {
  template <typename DerivedP, typename IndexType, typename CentersType,
  	typename WidthsType>
  IGL_INLINE void build_octree(const Eigen::MatrixBase<DerivedP>& P,
    std::vector<std::vector<IndexType> > & point_indices,
    std::vector<Eigen::Matrix<IndexType,8,1>,
    	Eigen::aligned_allocator<Eigen::Matrix<IndexType,8,1> > > & children,
    std::vector<Eigen::Matrix<CentersType,1,3>,
    	Eigen::aligned_allocator<Eigen::Matrix<CentersType,1,3> > > & centers,
    std::vector<WidthsType> & widths)
  {
    const int MAX_DEPTH = 30000;

    typedef Eigen::Matrix<IndexType,8,1> Vector8i;
    typedef Eigen::Matrix<typename DerivedP::Scalar, 1, 3> RowVector3PType;
    typedef Eigen::Matrix<CentersType, 1, 3> RowVector3CentersType;
    
    auto get_octant = [](RowVector3PType location,
                         RowVector3CentersType center){
      // We use a binary numbering of children. Treating the parent cell's
      // center as the origin, we number the octants in the following manner:
      // The first bit is 1 iff the octant's x coordinate is positive
      // The second bit is 1 iff the octant's y coordinate is positive
      // The third bit is 1 iff the octant's z coordinate is positive
      //
      // For example, the octant with negative x, positive y, positive z is:
      // 110 binary = 6 decimal
      IndexType index = 0;
      if( location(0) >= center(0)){
        index = index + 1;
      }
      if( location(1) >= center(1)){
        index = index + 2;
      }
      if( location(2) >= center(2)){
        index = index + 4;
      }
      return index;
    };

    
    auto translate_center = [](const RowVector3CentersType & parent_center,
                               const CentersType h,
                               const IndexType child_index){
      RowVector3CentersType change_vector;
      change_vector << -h,-h,-h;
      //positive x chilren are 1,3,4,7
      if(child_index % 2){
        change_vector(0) = h;
      }
      //positive y children are 2,3,6,7
      if(child_index == 2 || child_index == 3 ||
        child_index == 6 || child_index == 7){
        change_vector(1) = h;
      }
      //positive z children are 4,5,6,7
      if(child_index > 3){
        change_vector(2) = h;
      }
      return parent_center + change_vector;
    };
  
    // How many cells do we have so far?
    IndexType m = 0;
  
    // Useful list of number 0..7
    const Vector8i zero_to_seven = (Vector8i()<<0,1,2,3,4,5,6,7).finished();
    const Vector8i neg_ones = (Vector8i()<<-1,-1,-1,-1,-1,-1,-1,-1).finished();
  
    std::function< void(const int, const int) > helper;
    helper = [&helper,&translate_center,&get_octant,&m,
              &zero_to_seven,&neg_ones,&P,
              &point_indices,&children,&centers,&widths]
    (const IndexType index, const int depth)-> void
    {
      if(point_indices.at(index).size() > 1 && depth < MAX_DEPTH){
        //give the parent access to the children
        children.at(index) = zero_to_seven.array() + m;
        //make the children's data in our arrays
      
        //Add the children to the lists, as default children
        WidthsType h = widths.at(index)/2;
        RowVector3CentersType curr_center = centers.at(index);
        for(int i = 0; i < 8; i++){
          children.emplace_back(neg_ones);
          point_indices.emplace_back(std::vector<IndexType>());
          centers.emplace_back(translate_center(curr_center,h,i));
          widths.emplace_back(h);
        }
      
        //Split up the points into the corresponding children
        for(int j = 0; j < point_indices.at(index).size(); j++){
          IndexType curr_point_index = point_indices.at(index).at(j);
          IndexType cell_of_curr_point =
          	get_octant(P.row(curr_point_index),curr_center)+m;
          point_indices.at(cell_of_curr_point).emplace_back(curr_point_index);
        }
      
        //Now increase m
        m += 8;
      
        // Look ma, I'm calling myself.
        for(int i = 0; i < 8; i++){
          helper(children.at(index)(i),depth+1);
        }
      }
    };
  
    {
      std::vector<IndexType> all(P.rows());
      for(IndexType i = 0;i<all.size();i++) all[i]=i;
      point_indices.emplace_back(all);
    }
    children.emplace_back(neg_ones);
  
    //Get the minimum AABB for the points
    RowVector3PType backleftbottom(P.col(0).minCoeff(),
                                   P.col(1).minCoeff(),
                                   P.col(2).minCoeff());
    RowVector3PType frontrighttop(P.col(0).maxCoeff(),
                                  P.col(1).maxCoeff(),
                                  P.col(2).maxCoeff());
    RowVector3CentersType aabb_center = (backleftbottom+frontrighttop)/2.0;
    WidthsType aabb_width = std::max(std::max(
                                          frontrighttop(0) - backleftbottom(0),
                                          frontrighttop(1) - backleftbottom(1)),
                                 					frontrighttop(2) - backleftbottom(2));
    centers.emplace_back( aabb_center );
  
    //Widths are the side length of the cube, (not half the side length):
    widths.emplace_back( aabb_width );
    m++;
    // then you have to actually call the function
    helper(0,0);
  }
}

