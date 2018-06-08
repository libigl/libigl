#include "build_octree.h"
#include <vector>
#include <queue>

const int MAX_DEPTH = 1000;

int get_octant(const Eigen::RowVector3d & location, const Eigen::RowVector3d & center){
    //Binary numbering of children:
    //if we treat location as the origin,
    //first bit is 1 if x is positive, 0 if not
    //second bit is 1 if y is positive, 0 if not
    //third bit is 1 if z is positive, 0 if not
    //for example, the quadrant with negative x, positive y, positive z is:
    //110 binary = 6 decimal
    int index = 0;
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
}


Eigen::RowVector3d translate_center(const Eigen::RowVector3d & parent_center, double h, int child_index){
    Eigen::RowVector3d change_vector;
    change_vector << h,h,h;
    if(child_index % 2){ //positive x chilren are 1,3,4,7
        change_vector(0) = -h;
    }
    if(child_index == 2 || child_index == 3 ||
       child_index == 6 || child_index == 7){ //positive y children are 2,3,6,7
        change_vector(1) = -h;
    }
    if(child_index > 3){ //positive z children are 4,5,6,7
        change_vector(2) = -h;
    }
    return parent_center - change_vector;
}

void build_octree(const Eigen::MatrixXd & P,
                  std::vector<std::vector<int> > & point_indices,
                  std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                  std::vector<Eigen::RowVector3d, Eigen::aligned_allocator<Eigen::RowVector3d>> & centers,
                  std::vector<double> & widths
                  ){
    typedef Eigen::Matrix<int,8,1> Vector8i;
    typedef Eigen::Matrix<double,8,1> Vector8d;

    // How many cells do we have so far?
    int m = 0;
    
    // Useful list of number 0..7
    const Vector8i zero_to_seven = (Vector8i()<<0,1,2,3,4,5,6,7).finished();
    const Vector8i neg_ones = (Vector8i()<<-1,-1,-1,-1,-1,-1,-1,-1).finished();
    
    // This function variable needs to be declared before it is defined so that
    // you can capture it and call it recursively. Annoyingly this means you need
    // to fill in the function prototype here, too.
    std::function< void(const int, const int) > helper;
    helper = [
              // These variables from the parent scope are "captured" you can read and
              // write to them
              &helper,
              &zero_to_seven,&neg_ones,&P,&point_indices,&children,&centers,&widths,&m]
    // The "-> bool" means that this function will return bool (I don't think
    // you need this, but in case you do.)
    (const int index, const int depth)-> void
    {
        if(point_indices.at(index).size() > 1 && depth < MAX_DEPTH){
            //give the parent access to the children
            children.at(index) = zero_to_seven.array() + m;
            //make the children's data in our arrays
            
            //Add the children to the lists, as default children
            double h = widths.at(index)/2;
            Eigen::RowVector3d curr_center = centers.at(index);
            for(int i = 0; i < 8; i++){
                children.emplace_back(-1 * Eigen::MatrixXi::Ones(8,1));
                point_indices.emplace_back(std::vector<int>());
                centers.emplace_back(translate_center(curr_center,h,i));
                widths.emplace_back(h);
            }
            
            //Split up the points into the corresponding children
            for(int j = 0; j < point_indices.at(index).size(); j++){
                int curr_point_index = point_indices.at(index).at(j);
                int cell_of_curr_point = get_octant(P.row(curr_point_index),curr_center)+m;
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
        std::vector<int> all(P.rows());
        for(int i = 0;i<all.size();i++) all[i]=i;
        point_indices.emplace_back(all);
    }
    children.emplace_back(-1 * Eigen::MatrixXi::Ones(8,1));
    
    //Get the minimum AABB for the points
    Eigen::RowVector3d backleftbottom(P.col(0).minCoeff(),P.col(1).minCoeff(),P.col(2).minCoeff());
    Eigen::RowVector3d frontrighttop(P.col(0).maxCoeff(),P.col(1).maxCoeff(),P.col(2).maxCoeff());
    Eigen::RowVector3d aabb_center = (backleftbottom+frontrighttop)/2.0;
    double aabb_width = std::max(std::max(frontrighttop(0) - backleftbottom(0),
                                          frontrighttop(1) - backleftbottom(1)),
                                 frontrighttop(2) - backleftbottom(2));
    centers.emplace_back( aabb_center );
    widths.emplace_back( aabb_width ); //Widths are the side length of the cube, (not half the side length)
    m++;
    // then you have to actually call the function
    helper(0,0);
    
    assert(point_indices.size() == m);
    assert(children.size() == m);
    assert(centers.size() == m);
    assert(widths.size() == m);
}
//}
//#ifndef IGL_STATIC_LIBRARY
//#  include "point_areas_and_normals.cpp"
//#endif
//#endif

