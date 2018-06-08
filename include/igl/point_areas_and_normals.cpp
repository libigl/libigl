#include "point_areas_and_normals.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <igl/copyleft/cgal/delaunay_triangulation.h>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>
#include <igl/parallel_for.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;

template <typename DerivedP, typename DerivedI, typename DerivedO>
IGL_INLINE void point_areas_and_normals(
  const Eigen::MatrixBase<DerivedP>& P,
  const Eigen::MatrixBase<DerivedI>& I,
  const Eigen::MatrixBase<DerivedO>& O,
  Eigen::PlainObjectBase<DerivedA> & A,
  Eigen::PlainObjectBase<DerivedA> & N);
{
//namespace igl {
//    void point_areas_and_normals(const Eigen::MatrixXd & P,
//                                 const Eigen::MatrixXi & I,
//                                 const Eigen::MatrixXd & O,
//                                 Eigen::VectorXd & A,
//                                 Eigen::MatrixXd & N
//                                 )
//    {
        const int n = P.rows();
        A.setZero(n,1);
        N.setZero(n,3);
        igl::parallel_for(P.rows(),[&](int i)
      {
          Eigen::MatrixXi neighbour_index = I.row(i);
          Eigen::MatrixXd neighbours;
          igl::slice(P,neighbour_index,1,neighbours);
          if(O.rows() && neighbours.rows() > 1){
              Eigen::MatrixXd neighbour_normals;
              igl::slice(O,neighbour_index,1,neighbour_normals);
              Eigen::Vector3d poi_normal = neighbour_normals.row(0);
              Eigen::VectorXd dotprod = poi_normal(0)*neighbour_normals.col(0)
              + poi_normal(1)*neighbour_normals.col(1)
              + poi_normal(2)*neighbour_normals.col(2);
              Eigen::Array<bool,Eigen::Dynamic,1> keep = dotprod.array() > 0;
              igl::slice_mask(Eigen::MatrixXd(neighbours),keep,1,neighbours);
          }
          if(neighbours.rows() <= 2){
              A(i) = 0;
              N.row(i) = Eigen::RowVector3d::Zero();
          } else {
              //subtract the mean from neighbours, then take svd, the scores will be U*S;
              //This is our pca plane fitting
              Eigen::RowVector3d mean = neighbours.colwise().mean();
              Eigen::MatrixXd mean_centered = neighbours.rowwise() - mean;
              Eigen::JacobiSVD<Eigen::MatrixXd> svd(mean_centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
              Eigen::MatrixXd scores = svd.matrixU() * svd.singularValues().asDiagonal();
              N.row(i) = svd.matrixV().col(2).transpose();
              if(N.row(i).dot(O.row(i)) < 0){
                  N.row(i) *= -1;
              }
              Eigen::MatrixXd plane;
              igl::slice(scores,igl::colon<int>(0,scores.rows()-1),igl::colon<int>(0,1),plane);
              std::vector< std::pair<Point,unsigned> > points;
              
              //This is where we obtain a delaunay triangulation of the points
              for(unsigned iter = 0; iter < plane.rows(); iter++){
                  points.push_back( std::make_pair( Point(plane(iter,0),plane(iter,1)), iter ) );
              }
              Delaunay triangulation;
              triangulation.insert(points.begin(),points.end());
              Eigen::MatrixXi F(triangulation.number_of_faces(),3);
              int f_row = 0;
              for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
                  fit != triangulation.finite_faces_end(); ++fit) {
                  Delaunay::Face_handle face = fit;
                  F.row(f_row) = Eigen::RowVector3i((int)face->vertex(0)->info(),
                                                    (int)face->vertex(1)->info(),
                                                    (int)face->vertex(2)->info());
                  f_row++;
              }
              
              //Here we calculate the voronoi area of the point
              double area_accumulator = 0;
              for(int f = 0; f < F.rows(); f++){
                  int X = -1;
                  for(int face_iter = 0; face_iter < 3; face_iter++){
                      if(F(f,face_iter) == 0){
                          X = face_iter;
                      }
                  }
                  if(X >= 0){
                      //Triangle XYZ with X being the point we want the area of
                      int Y = (X+1)%3;
                      int Z = (X+2)%3;
                      double x = (plane.row(F(f,Y))-plane.row(F(f,Z))).norm();
                      double y = (plane.row(F(f,X))-plane.row(F(f,Z))).norm();
                      double z = (plane.row(F(f,Y))-plane.row(F(f,X))).norm();
                      double cosX = (z*z + y*y - x*x)/(2*y*z);
                      double cosY = (z*z + x*x - y*y)/(2*x*z);
                      double cosZ = (x*x + y*y - z*z)/(2*y*x);
                      Eigen::Vector3d barycentric;
                      barycentric << x*cosX, y*cosY, z*cosZ;
                      barycentric /= (barycentric(0) + barycentric(1) + barycentric(2));
                      //TODO: to make numerically stable, reorder so that x≥y≥z:
                      double full_area = 0.25*std::sqrt((x+(y+z))*(z-(x-y))*(z+(x-y))*(x+(y-z)));
                      Eigen::Vector3d partial_area = barycentric * full_area;
                      if(cosX < 0){
                          area_accumulator += 0.5*full_area;
                      } else if (cosY < 0 || cosZ < 0){
                          area_accumulator += 0.25*full_area;
                      } else {
                          area_accumulator += (partial_area(1) + partial_area(2))/2;
                      }
                  }
              }
              if(std::isfinite(area_accumulator)){
                  A(i) = area_accumulator;
              } else {
                  A(i) = 0;
              }
          }
      },1000);
    }
}




