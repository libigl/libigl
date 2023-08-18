#include <igl/readOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/edges.h>
#include <igl/find.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polygons_to_triangles.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OBJ format
  Eigen::MatrixXd OV,V;
  Eigen::VectorXi I,C;
  igl::readOBJ(
    argc<=1?TUTORIAL_SHARED_PATH "/cylinder.obj" :argv[1],V,I,C);
  // center
  V.rowwise() -= V.colwise().mean();
  OV = V;
  // Convert polygon representation to triangles
  Eigen::MatrixXi F;
  Eigen::VectorXi J;
  igl::polygons_to_triangles(I,C,F,J);
  Eigen::SparseMatrix<double> pL,pM,pP;
  igl::cotmatrix(V,I,C,pL,pM,pP);
  Eigen::SparseMatrix<double> tL,tM;
  igl::cotmatrix(V,F,tL);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,tM);
  const double bbd = (V.colwise().maxCoeff()- V.colwise().minCoeff()).norm();
  igl::opengl::glfw::Viewer vr;
  vr.data_list[0].set_mesh(V,F);
  vr.append_mesh();
  vr.selected_data_index = 0;
  vr.data_list[0].set_face_based(true);
  Eigen::MatrixXi E;
  igl::edges(I,C,E);
  bool show_edges = true;
  bool use_poly = true;
  Eigen::MatrixXd pHN;
  Eigen::MatrixXd tHN;
  const auto update = [&]()
  {
    // why does windows need this Eigen::VectorXd(...) ?
    pHN = (pL*V).array().colwise() / Eigen::VectorXd(pM.diagonal()).array();
    tHN = (tL*V).array().colwise() / Eigen::VectorXd(tM.diagonal()).array();
    pHN *= 1.0/pHN.rowwise().norm().maxCoeff();
    tHN *= 1.0/tHN.rowwise().norm().maxCoeff();
    const auto was_face_based  = vr.data_list[0].face_based;
    Eigen::MatrixXd QV(V.rows()*2,3);
    QV.topRows(V.rows()) = V;
    if(use_poly)
    {
      printf("using polygon Laplacian\n");
      QV.bottomRows(V.rows()) = V-pHN;
    }else
    {
      printf("using triangle Laplacian\n");
      QV.bottomRows(V.rows()) = V-tHN;
    }
    Eigen::MatrixXi QE(V.rows(),2);
    for(int i = 0;i<V.rows();i++){ QE(i,0)=i;QE(i,1)=i+V.rows();}
    vr.data_list[1].set_edges(QV,QE,Eigen::RowVector3d(1,1,1));

    if(use_poly)
    {
      vr.data_list[0].show_lines = false;
      if(show_edges)
      {
        vr.data_list[0].set_edges(V,E,Eigen::RowVector3d(0,0,0));
      }else
      {
        vr.data_list[0].clear_edges();
      }
    }else
    {
      vr.data_list[0].clear_edges();
      vr.data_list[0].show_lines = show_edges;
    }
    vr.data_list[0].set_face_based(was_face_based);
  };
  const double original_area = pM.diagonal().sum();
  const auto recompute_M = [&]()
  {
    Eigen::SparseMatrix<double> _1,_2;
    igl::cotmatrix(V,I,C,_1,pM,_2);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,tM);
    V *= sqrt(original_area / pM.diagonal().sum());
  };
  const auto cmcf_step = [&]()
  {
    const Eigen::SparseMatrix<double> S =
      use_poly? ((pM) - 0.05*(pL)): ((tM) - 0.05*(tL));
    const Eigen::MatrixXd rhs = use_poly? pM*V : tM*V;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
    assert(solver.info() == Eigen::Success);
    V = solver.solve(rhs).eval();
    // recompute just mass matrices
    recompute_M();
    // center
    V.rowwise() -= V.colwise().mean();
    vr.data_list[0].set_vertices(V);
    vr.data_list[0].compute_normals();
    update();
  };
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ' ': cmcf_step(); return true;
      case 'R': case 'r': V=OV;recompute_M();vr.data_list[0].set_vertices(V);vr.data_list[0].compute_normals(); update();return true;
      case 'P': case 'p': use_poly=!use_poly; update();return true;
      case 'L': case 'l': show_edges=!show_edges; update();return true;
    }
    return false;
  };
  update();
  vr.launch();
}
