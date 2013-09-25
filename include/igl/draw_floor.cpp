#include "draw_floor.h"
#ifndef IGL_NO_OPENGL

#include "OpenGL_convenience.h"

IGL_INLINE void igl::draw_floor(const float * colorA, const float * colorB)
{
  glDisable(GL_LIGHTING);
  glColorMaterial( GL_FRONT, GL_EMISSION);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  // Set material
  const float black[] = {0.,0.,0.,1.};
  glMaterialfv(GL_FRONT, GL_AMBIENT, black);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, black);
  glMaterialfv(GL_FRONT, GL_SPECULAR, black);
  glMaterialfv(GL_FRONT, GL_EMISSION, black);
  glMaterialf(GL_FRONT, GL_SHININESS,0);
  const bool use_lighting = false;
  if(use_lighting)
  {
    glEnable(GL_LIGHTING);
  }else
  {
    glDisable(GL_LIGHTING);
  }

  int GridSizeX = 100;
  int GridSizeY = 100;
  float SizeX = 0.5f;
  float SizeY = 0.5f;

  glBegin(GL_QUADS);
  glNormal3f(0,1,0);
  for (int x =-GridSizeX/2;x<GridSizeX/2;++x)
  {
    for (int y =-GridSizeY/2;y<GridSizeY/2;++y)
    {
      if ((x+y)&0x00000001) //modulo 2
      {
        glColor4fv(colorA);
      }else
      {
        glColor4fv(colorB);
      }
      glVertex3f(    x*SizeX,0,    y*SizeY);
      glVertex3f((x+1)*SizeX,0,    y*SizeY);
      glVertex3f((x+1)*SizeX,0,(y+1)*SizeY);
      glVertex3f(    x*SizeX,0,(y+1)*SizeY);
    }
  }
  glEnd();
  glDisable(GL_COLOR_MATERIAL);
} 

IGL_INLINE void igl::draw_floor()
{
  const float grey[] = {0.80,0.80,0.80,1.};
  const float white[] = {0.95,0.95,0.95,1.};
  igl::draw_floor(grey,white);
}
#endif
