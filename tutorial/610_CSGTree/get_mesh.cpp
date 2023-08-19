#include "get_mesh.h"
#include <igl/copyleft/cgal/CSGTree.h>

void get_mesh(
    const Eigen::MatrixXd &VA,
    const Eigen::MatrixXi &FA,
    const Eigen::MatrixXd &VB,
    const Eigen::MatrixXi &FB,
    const Eigen::MatrixXd &VC,
    const Eigen::MatrixXi &FC,
    const Eigen::MatrixXd &VD,
    const Eigen::MatrixXi &FD,
    const Eigen::MatrixXd &VE,
    const Eigen::MatrixXi &FE,
    const int view_id,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    Eigen::VectorXd &I)
{
  using namespace Eigen;
  const auto & set_mesh = [&V,&F,&I](const Eigen::MatrixXd &V_, const Eigen::MatrixXi &F_, const int i)
  {
    V = V_;
    F = F_;
    I = VectorXd::Constant(F.rows(),1,i);
  };
  switch(view_id)
  {
    case 0:
      set_mesh(VA,FA,5);
      break;
    case 1:
      set_mesh(VB,FB,4);
      break;
    case 2:
      set_mesh(VC,FC,3);
      break;
    case 3:
      set_mesh(VD,FD,2);
      break;
    case 4:
      set_mesh(VE,FE,1);
      break;
    default:
    {
      igl::copyleft::cgal::CSGTree M;
      VectorXi J;
      switch(view_id)
      {
        case 5:
          // Compute result of (A ∩ B)
          M = {{VA,FA},{VB,FB},"i"};
          J = M.J().array()+0;
          break;
        case 6:
          // Compute result of (C ∪ D)
          M = {{VC,FC},{VD,FD},"u"};
          J = M.J().array()+FA.rows()+FB.rows();
          break;
        case 7:
          // Compute result of (C ∪ D) ∪ E
          M = {{{VC,FC},{VD,FD},"u"},{VE,FE},"u"};
          J = M.J().array()+FA.rows()+FB.rows();
          break;
        case 8:
          // Compute result of (A ∩ B) \ ((C ∪ D) ∪ E)
          M = {{{VA,FA},{VB,FB},"i"},{{{VC,FC},{VD,FD},"u"},{VE,FE},"u"},"m"};
          J = M.J().array()+0;
          break;
        default:
          assert(false && "unknown view id");
      }
      V = M.cast_V<MatrixXd>();
      F = M.F();
      I.resize(M.F().rows(),1);
      // Compute colors based on original facets
      for(int f = 0;f<M.F().rows();f++)
      {
        const int j = J(f);
        I(f) = 
          (int)(j<FA.rows())+
          (int)(j<FA.rows()+FB.rows())+
          (int)(j<FA.rows()+FB.rows()+FC.rows())+
          (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows())+
          (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows()+FE.rows());
      }
    }
  }
}
