#include <igl/fast_winding_number.h>
#include <igl/read_triangle_mesh.h>
#include <igl/slice_mask.h>
#include <Eigen/Geometry>
#include <igl/octree.h>
#include <igl/barycenter.h>
#include <igl/knn.h>
#include <igl/random_points_on_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/per_face_normals.h>
#include <igl/copyleft/cgal/point_areas.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  const auto time = [](std::function<void(void)> func)->double
  {
    const double t_before = igl::get_seconds();
    func();
    const double t_after = igl::get_seconds();
    return t_after-t_before;
  };

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc>1?argv[1]:TUTORIAL_SHARED_PATH "/bunny.off",V,F);
  // Sample mesh for point cloud
  Eigen::MatrixXd P,N;
  {
    Eigen::VectorXi I;
    Eigen::MatrixXd B;
    igl::random_points_on_mesh(10000,V,F,B,I,P);
    Eigen::MatrixXd FN;
    igl::per_face_normals(V,F,FN);
    N.resize(P.rows(),3);
    for(int p = 0;p<I.rows();p++)
    {
      N.row(p) = FN.row(I(p));
    }
  }
  // Build octree
  std::vector<std::vector<int > > O_PI;
  Eigen::MatrixXi O_CH;
  Eigen::MatrixXd O_CN;
  Eigen::VectorXd O_W;
  igl::octree(P,O_PI,O_CH,O_CN,O_W);
  Eigen::VectorXd A;
  {
    Eigen::MatrixXi I;
    igl::knn(P,20,O_PI,O_CH,O_CN,O_W,I);
    // CGAL is only used to help get point areas
    igl::copyleft::cgal::point_areas(P,I,N,A);
  }

  if(argc<=1)
  {
    // corrupt mesh
    Eigen::MatrixXd BC;
    igl::barycenter(V,F,BC);
    Eigen::MatrixXd OV = V;
    V.resize(F.rows()*3,3);
    for(int f = 0;f<F.rows();f++)
    {
      for(int c = 0;c<3;c++)
      {
        int v = f+c*F.rows();
        // random rotation about barycenter
        Eigen::AngleAxisd R(
          0.5*static_cast <double> (rand()) / static_cast <double> (RAND_MAX),
          Eigen::Vector3d::Random(3,1));
        V.row(v) = (OV.row(F(f,c))-BC.row(f))*R.matrix()+BC.row(f);
        F(f,c) = v;
      }
    }
  }

  // Generate a list of random query points in the bounding box
  Eigen::MatrixXd Q = Eigen::MatrixXd::Random(1000000,3);
  const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
  const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
  const Eigen::RowVector3d Vdiag = Vmax-Vmin;
  for(int q = 0;q<Q.rows();q++)
  {
    Q.row(q) = (Q.row(q).array()*0.5+0.5)*Vdiag.array() + Vmin.array();
  }

  // Positions of points inside of point cloud P
  Eigen::MatrixXd QiP;
  {
    Eigen::MatrixXd O_CM;
    Eigen::VectorXd O_R;
    Eigen::MatrixXd O_EC;
    printf("     point cloud precomputation (% 8ld points):    %g secs\n",
      P.rows(),
      time([&](){igl::fast_winding_number(P,N,A,O_PI,O_CH,2,O_CM,O_R,O_EC);}));
    Eigen::VectorXd WiP;
    printf("        point cloud evaluation  (% 8ld queries):   %g secs\n",
      Q.rows(),
      time([&](){igl::fast_winding_number(P,N,A,O_PI,O_CH,O_CM,O_R,O_EC,Q,2,WiP);}));
    igl::slice_mask(Q,(WiP.array()>0.5).eval(),1,QiP);
  }

  // Positions of points inside of triangle soup (V,F)
  Eigen::MatrixXd QiV;
  {
    igl::FastWindingNumberBVH fwn_bvh;
    printf("triangle soup precomputation    (% 8ld triangles): %g secs\n",
      F.rows(),
      time([&](){igl::fast_winding_number(V.cast<float>().eval(),F,2,fwn_bvh);}));
    Eigen::VectorXf WiV;
    printf("      triangle soup evaluation  (% 8ld queries):   %g secs\n",
      Q.rows(),
      time([&](){igl::fast_winding_number(fwn_bvh,2,Q.cast<float>().eval(),WiV);}));
    igl::slice_mask(Q,WiV.array()>0.5,1,QiV);
  }


  // Visualization
  igl::opengl::glfw::Viewer viewer;
  // For dislpaying normals as little line segments
  Eigen::MatrixXd PN(2*P.rows(),3);
  Eigen::MatrixXi E(P.rows(),2);
  const double bbd = igl::bounding_box_diagonal(V);
  for(int p = 0;p<P.rows();p++)
  {
    E(p,0) = 2*p;
    E(p,1) = 2*p+1;
    PN.row(E(p,0)) = P.row(p);
    PN.row(E(p,1)) = P.row(p)+bbd*0.01*N.row(p);
  }

  bool show_P = false;
  int show_Q = 0;

  int query_data = 0;
  viewer.data_list[query_data].set_mesh(V,F);
  viewer.data_list[query_data].clear();
  viewer.data_list[query_data].point_size = 2;
  viewer.append_mesh();
  int object_data = 1;
  viewer.data_list[object_data].set_mesh(V,F);
  viewer.data_list[object_data].point_size = 5;

  const auto update = [&]()
  {
    viewer.data_list[query_data].clear();
    switch(show_Q)
    {
      case 1:
        // show all Q
        viewer.data_list[query_data].set_points(Q,Eigen::RowVector3d(0.996078,0.760784,0.760784));
        break;
      case 2:
        // show all Q inside
        if(show_P)
        {
          viewer.data_list[query_data].set_points(QiP,Eigen::RowVector3d(0.564706,0.847059,0.768627));
        }else
        {
          viewer.data_list[query_data].set_points(QiV,Eigen::RowVector3d(0.564706,0.847059,0.768627));
        }
        break;
    }
    
    viewer.data_list[object_data].clear();
    if(show_P)
    {
      viewer.data_list[object_data].set_points(P,Eigen::RowVector3d(1,1,1));
      viewer.data_list[object_data].set_edges(PN,E,Eigen::RowVector3d(0.8,0.8,0.8));
    }else
    {
      viewer.data_list[object_data].set_mesh(V,F);
    }
  };



  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
      default: 
        return false;
      case '1':
        show_P = !show_P;
        break;
      case '2':
        show_Q = (show_Q+1) % 3;
        break;
    }
    update();
    return true;
  };

  std::cout<<R"(
FastWindingNumber
  1  Toggle point cloud and triangle soup
  2  Toggle hiding query points, showing query points, showing inside queries
)";

  update();
  viewer.launch();

}
