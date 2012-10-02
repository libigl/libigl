#include "draw_mesh.h"

IGL_INLINE void igl::draw_mesh(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & N)
{
  glBegin(GL_TRIANGLES);
  // loop over faces
  for(int i = 0; i<F.rows();i++)
  {
    // loop over corners of triangle
    for(int j = 0;j<3;j++)
    {
      glNormal3d(N(F(i,j),0),N(F(i,j),1),N(F(i,j),2));
      glVertex3d(V(F(i,j),0),V(F(i,j),1),V(F(i,j),2));
    }
  }
  glEnd();
}

IGL_INLINE void igl::draw_mesh(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & N,
  const Eigen::MatrixXd & C)
{
  glBegin(GL_TRIANGLES);
  // loop over faces
  for(int i = 0; i<F.rows();i++)
  {
    // loop over corners of triangle
    for(int j = 0;j<3;j++)
    {
      glColor3d(C(F(i,j),0),C(F(i,j),1),C(F(i,j),2));
      glNormal3d(N(F(i,j),0),N(F(i,j),1),N(F(i,j),2));
      glVertex3d(V(F(i,j),0),V(F(i,j),1),V(F(i,j),2));
    }
  }
  glEnd();
}

IGL_INLINE void igl::draw_mesh(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & N,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXd & TC,
  const Eigen::MatrixXd & W,
  const GLuint W_index,
  const Eigen::MatrixXi & WI,
  const GLuint WI_index)
{
  using namespace std;
  if(F.size() > 0)
  {
    assert(F.maxCoeff() < V.rows());
    assert(V.cols() == 3);
    assert(C.rows() == V.rows() || C.rows() == F.rows()*3 || C.size() == 0);
    assert(TC.rows() == V.rows() || TC.rows() == F.rows()*3 || TC.size() == 0);
    assert(C.cols() == 3 || C.size() == 0);
    assert(
      N.rows() == V.rows() || N.rows() == F.rows()*3 || N.rows() ==F.rows());
    assert(N.cols() == 3);
  }
  if(W.size()>0)
  {
    assert(W.rows() == V.rows());
    assert(WI.rows() == V.rows());
    assert(W.cols() == WI.cols());
  }

  glBegin(GL_TRIANGLES);
  // loop over faces
  for(int i = 0; i<F.rows();i++)
  {
    // loop over corners of triangle
    for(int j = 0;j<3;j++)
    {
      if(W.size()>0 && W_index !=0 && WI_index != 0)
      {
        int weights_left = W.cols();
        while(weights_left != 0)
        {
          int pass_size = std::min(4,weights_left);
          int weights_already_passed = W.cols()-weights_left;
          // Get attribute location of next 4 weights
          int pass_W_index = W_index + weights_already_passed/4;
          int pass_WI_index = WI_index + weights_already_passed/4;
          switch(pass_size)
          {
            case 1:
              glVertexAttrib1d(
                pass_W_index,
                W(F(i,j),0+weights_already_passed));
              glVertexAttrib1d(
                pass_WI_index,
                WI(F(i,j),0+weights_already_passed));
              break;
            case 2:
              glVertexAttrib2d(
                pass_W_index,
                W(F(i,j),0+weights_already_passed),
                W(F(i,j),1+weights_already_passed));
              glVertexAttrib2d(
                pass_WI_index,
                WI(F(i,j),0+weights_already_passed),
                WI(F(i,j),1+weights_already_passed));
              break;
            case 3:
              glVertexAttrib3d(
                pass_W_index,
                W(F(i,j),0+weights_already_passed),
                W(F(i,j),1+weights_already_passed),
                W(F(i,j),2+weights_already_passed));
              glVertexAttrib3d(
                pass_WI_index,
                WI(F(i,j),0+weights_already_passed),
                WI(F(i,j),1+weights_already_passed),
                WI(F(i,j),2+weights_already_passed));
              break;
            default:
              glVertexAttrib4d(
                pass_W_index,
                W(F(i,j),0+weights_already_passed),
                W(F(i,j),1+weights_already_passed),
                W(F(i,j),2+weights_already_passed),
                W(F(i,j),3+weights_already_passed));
              glVertexAttrib4d(
                pass_WI_index,
                WI(F(i,j),0+weights_already_passed),
                WI(F(i,j),1+weights_already_passed),
                WI(F(i,j),2+weights_already_passed),
                WI(F(i,j),3+weights_already_passed));
              break;
          }
          weights_left -= pass_size;
        }
      }
      if(TC.rows() == V.rows())
      {
        glTexCoord2d(TC(F(i,j),0),TC(F(i,j),1));
      }else if(TC.rows() == F.rows()*3)
      {
        glTexCoord2d(TC(F(i,j),0),TC(F(i,j),1));
      }
      if(C.rows() == V.rows())
      {
        glColor3d(C(F(i,j),0),C(F(i,j),1),C(F(i,j),2));
      }else if(C.rows() == F.rows()*3)
      {
        glColor3d(C(i*3+j,0), C(i*3+j,1), C(i*3+j,2));
      }
      if(N.rows() == V.rows())
      {
        glNormal3d(N(F(i,j),0),N(F(i,j),1),N(F(i,j),2));
      }else if(N.rows() == F.rows()*3)
      {
        glNormal3d(N(i*3+j,0),N(i*3+j,1),N(i*3+j,2));
      }else if(N.rows() == F.rows())
      {
        glNormal3d(N(i,0),N(i,1),N(i,2));
      }
      glVertex3d(V(F(i,j),0),V(F(i,j),1),V(F(i,j),2));
    }
  }
  glEnd();
}

