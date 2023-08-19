#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <igl/stb/read_image.h>
#include <igl/get_seconds.h>
#include <igl/PI.h>
#include <Eigen/Core>

#include <algorithm>

void floor(const Eigen::MatrixXd & V, 
    Eigen::MatrixXd & fV,
    Eigen::MatrixXd & fU,
    Eigen::MatrixXi & fF)
{
  igl::triangulated_grid(2,2,fU,fF);
  fV = fU;
  fV.array() -= 0.5;
  fV.array() *= 2 * 2 * 
    (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  fV = (fV * (Eigen::Matrix<double,2,3>()<<1,0,0,0,0,-1).finished() ).eval();
  fV.col(0).array() += 0.5*(V.col(0).minCoeff()+ V.col(0).maxCoeff());
  fV.col(1).array() += V.col(1).minCoeff();
  fV.col(2).array() += 0.5*(V.col(2).minCoeff()+ V.col(2).maxCoeff());
}

void checker_texture(const int s, const int f,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & X,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & A)
{
  X.resize(s*f,s*f);
  A.resize(s*f,s*f);
  for(int i = 0;i<s*f;i++)
  {
    const double x = double(i)/double(s*f-1)*2-1;
    for(int j = 0;j<s*f;j++)
    {
      const int u = i/f;
      const int v = j/f;
      const double y = double(j)/double(s*f-1)*2-1;
      const double r1 = std::min(std::max( (1.0 - sqrt(x*x+y*y))*1.0 ,0.0),1.0);
      const double r3 = std::min(std::max( (1.0 - sqrt(x*x+y*y))*3.0 ,0.0),1.0);
      //const double a = 3*r*r - 2*r*r*r;
      const auto smooth_step = [](const double w)
      {
        return ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) ;
      };
      double a3 = smooth_step(r1);
      double a1 = smooth_step(r1);
      X(i,j) = (0.75+0.25*a1) * (u%2 == v%2 ? 245 : 235);
      A(i,j) = a3 * 255;
    }
  }
}

int main(int argc, char *argv[])
{
  igl::opengl::glfw::Viewer vr;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);

  // Create a floor
  Eigen::MatrixXd fV, fU;
  Eigen::MatrixXi fF;
  floor(V,fV,fU,fF);
  const int s = 16;
  const int f = 100;
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> X;
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A;
  checker_texture(s,f,X,A);
  vr.data().set_mesh(fV,fF);
  vr.data().set_uv(fU);
  vr.data().uniform_colors(Eigen::Vector3d(0.3,0.3,0.3),Eigen::Vector3d(0.8,0.8,0.8),Eigen::Vector3d(0,0,0));
  vr.data().set_texture(X,X,X,A);
  vr.data().show_texture = true;
  vr.data().show_lines = false;

  // Move the light a bit off center to cast a more visible shadow.
  vr.core().light_position << 1.0f, 2.0f, 0.0f;
  // For now, the default is a positional light with no shadows. Meanwhile,
  // shadows only support a directional light. To best match the appearance of
  // current lighting use this conversion when turning on shadows. In the
  // future, hopefully this will reduce to just 
  //     core().is_shadow_mapping = true
  vr.core().is_directional_light = true;
  vr.core().light_position = vr.core().light_position + vr.core().camera_eye;
  vr.core().is_shadow_mapping = true;

  // Send the main object to the viewer
  vr.append_mesh();
  vr.data().set_mesh(V,F);
  vr.data().show_lines = false;
  vr.data().set_face_based(true);

  // If a second argument is present read it as a matcap
  if(argc>2)
  {
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
    igl::stb::read_image(argv[2],R,G,B,A);
    // If more args, read them as light direction
    if(argc>2+3)
    {
      Eigen::Vector3f D;
      D << std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]);
      D.normalize();
      vr.core().light_position = D;
      Eigen::Vector3d Ka(0.14,0.14,0.14);
      if(argc>2+3+1 && std::atoi(argv[6]))
      {
        // Assume that color opposite D along umbra boundary is ambient color
        const double s = -D(2)*1.0/sqrt((D(0)*D(0)+D(1)*D(1))*(D(0)*D(0)+D(1)*D(1)+D(2)*D(2))); 
        const int i = ((D(0)*s)*0.5+0.5)*R.cols();
        const int j = ((D(1)*s)*0.5+0.5)*R.rows();
        Ka << double(R(i,j))/255.0, double(G(i,j))/255.0, double(B(i,j))/255.0;
      }
      std::cout<<"Ka : "<<Ka<<std::endl;
      // viewer only exposes ambient color through per-face and per-vertex
      // materials
      vr.data().V_material_ambient.col(0).setConstant( Ka(0) );
      vr.data().V_material_ambient.col(1).setConstant( Ka(1) );
      vr.data().V_material_ambient.col(2).setConstant( Ka(2) );
      vr.data().F_material_ambient.col(0).setConstant( Ka(0) );
      vr.data().F_material_ambient.col(1).setConstant( Ka(1) );
      vr.data().F_material_ambient.col(2).setConstant( Ka(2) );
    }
    vr.data().set_texture(R,G,B,A);
    vr.data().use_matcap = true;
  }
  vr.core().is_animating = true;
  vr.core().camera_zoom *= 1.5;
  vr.callback_pre_draw = [&](decltype(vr)&)
  {
    if(vr.core().is_animating)
    {
      vr.core().trackball_angle = Eigen::AngleAxisf( 
          sin(igl::get_seconds())*igl::PI*0.5,
          Eigen::Vector3f(0,1,0));
    }
    return false;
  };
  vr.launch();
}
