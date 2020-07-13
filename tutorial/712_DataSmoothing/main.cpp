#include <igl/read_triangle_mesh.h>
#include <igl/hessian_energy.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/isolines_map.h>
#include <igl/parula.h>
#include <igl/vertex_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <iostream>
#include <set>
#include <limits>
#include <stdlib.h>

#include "tutorial_shared_path.h"



int main(int argc, char * argv[])
{
    typedef Eigen::SparseMatrix<double> SparseMat;

    //Read our mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if(!igl::read_triangle_mesh(
        argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off",V,F)) {
        std::cout << "Failed to load mesh." << std::endl;
    }

    //Constructing an exact function to smooth
    Eigen::VectorXd zexact = V.block(0,2,V.rows(),1).array()
        + 0.5*V.block(0,1,V.rows(),1).array()
        + V.block(0,1,V.rows(),1).array().pow(2)
        + V.block(0,2,V.rows(),1).array().pow(3);

    //Make the exact function noisy
    srand(5);
    const double s = 0.2*(zexact.maxCoeff() - zexact.minCoeff());
    Eigen::VectorXd znoisy = zexact + s*Eigen::VectorXd::Random(zexact.size());

    //Constructing the squared Laplacian and squared Hessian energy
    SparseMat L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    Eigen::SimplicialLDLT<SparseMat> solver(M);
    SparseMat MinvL = solver.solve(L);
    SparseMat QL = L.transpose()*MinvL;
    SparseMat QH;
    igl::hessian_energy(V, F, QH);

    //Solve to find Laplacian-smoothed and Hessian-smoothed solutions
    const double al = 8e-4;
    Eigen::SimplicialLDLT<SparseMat> lapSolver(al*QL + (1.-al)*M);
    Eigen::VectorXd zl = lapSolver.solve(al*M*znoisy);
    const double ah = 5e-6;
    Eigen::SimplicialLDLT<SparseMat> hessSolver(ah*QH + (1.-ah)*M);
    Eigen::VectorXd zh = hessSolver.solve(ah*M*znoisy);

    //Viewer that shows all functions: zexact, znoisy, zl, zh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
    viewer.data().show_lines = false;
    viewer.callback_key_down =
      [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
      {
        //Graduate result to show isolines, then compute color matrix
        const Eigen::VectorXd* z;
        switch(key) {
            case '1':
                z = &zexact;
                break;
            case '2':
                z = &znoisy;
                break;
            case '3':
                z = &zl;
                break;
            case '4':
                z = &zh;
                break;
            default:
                return false;
        }
        viewer.data().set_data(*z);
        return true;
    };
    std::cout << R"(Usage:
1  Show original function
2  Show noisy function
3  Biharmonic smoothing (zero Neumann boundary)
4  Biharmonic smoothing (natural Hessian boundary)

)";
    Eigen::MatrixXd CM;
    igl::parula(Eigen::VectorXd::LinSpaced(21,0,1).eval(),false,CM);
    igl::isolines_map(Eigen::MatrixXd(CM),CM);
    viewer.data().set_colormap(CM);
    viewer.data().set_data(znoisy);
    viewer.launch();

    return 0;
}
