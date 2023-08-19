#include "get_soft_constraint_for_circle.h"
#include <igl/boundary_loop.h>

void get_soft_constraint_for_circle(Eigen::MatrixXd& V_o, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc) {

    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);
    const int B_STEPS = 22; // constraint every B_STEPS vertices of the boundary

    b.resize(bnd.rows()/B_STEPS);
    bc.resize(b.rows(),2);

    int c_idx = 0;
    for (int i = B_STEPS; i < bnd.rows(); i+=B_STEPS) {
        b(c_idx) = bnd(i);
        c_idx++;
    }

    bc.row(0) = V_o.row(b(0)); // keep it there for now
    bc.row(1) = V_o.row(b(2));
    bc.row(2) = V_o.row(b(3));
    bc.row(3) = V_o.row(b(4));
    bc.row(4) = V_o.row(b(5));


    bc.row(0) << V_o(b(0),0), 0;
    bc.row(4) << V_o(b(4),0), 0;
    bc.row(2) << V_o(b(2),0), 0.1;
    bc.row(3) << V_o(b(3),0), 0.05;
    bc.row(1) << V_o(b(1),0), -0.15;
    bc.row(5) << V_o(b(5),0), +0.15;
}
