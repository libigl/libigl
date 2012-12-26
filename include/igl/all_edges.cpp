#include "all_edges.h"

IGL_INLINE void igl::all_edges(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & E)
{
  E.resize(F.rows()*F.cols(),F.cols()-1);
  switch(F.cols())
  {
    case 4:
      E.block(0*F.rows(),0,F.rows(),1) = F.col(1);
      E.block(0*F.rows(),1,F.rows(),1) = F.col(3);
      E.block(0*F.rows(),2,F.rows(),1) = F.col(2);

      E.block(1*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(1*F.rows(),1,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),2,F.rows(),1) = F.col(3);

      E.block(2*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),1,F.rows(),1) = F.col(3);
      E.block(2*F.rows(),2,F.rows(),1) = F.col(1);

      E.block(3*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(3*F.rows(),1,F.rows(),1) = F.col(1);
      E.block(3*F.rows(),2,F.rows(),1) = F.col(2);
      return;
    case 3:
      E.block(0*F.rows(),0,F.rows(),1) = F.col(1);
      E.block(0*F.rows(),1,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),0,F.rows(),1) = F.col(2);
      E.block(1*F.rows(),1,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),0,F.rows(),1) = F.col(0);
      E.block(2*F.rows(),1,F.rows(),1) = F.col(1);
      return;
  }
}
