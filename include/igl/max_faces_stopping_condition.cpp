#include "max_faces_stopping_condition.h"

IGL_INLINE void igl::max_faces_stopping_condition(
  int & m,
  const int max_m,
  std::function<bool(
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const std::set<std::pair<double,int> > &,
    const std::vector<std::set<std::pair<double,int> >::iterator > &,
    const Eigen::MatrixXd &,
    const int,
    const int,
    const int,
    const int,
    const int)> & stopping_condition)
{
  stopping_condition = 
    [&max_m,&m](
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const std::set<std::pair<double,int> > &,
    const std::vector<std::set<std::pair<double,int> >::iterator > &,
    const Eigen::MatrixXd &,
    const int,
    const int,
    const int,
    const int,
    const int)->bool
    {
      m-=2;
      return m<=(int)max_m;
    };
}

IGL_INLINE 
  std::function<bool(
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const std::set<std::pair<double,int> > &,
    const std::vector<std::set<std::pair<double,int> >::iterator > &,
    const Eigen::MatrixXd &,
    const int,
    const int,
    const int,
    const int,
    const int)> 
  igl::max_faces_stopping_condition(
    int & m,
    const int max_m)
{
  std::function<bool(
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const std::set<std::pair<double,int> > &,
    const std::vector<std::set<std::pair<double,int> >::iterator > &,
    const Eigen::MatrixXd &,
    const int,
    const int,
    const int,
    const int,
    const int)> stopping_condition;
  max_faces_stopping_condition(
      m,max_m,stopping_condition);
  return stopping_condition;
}
