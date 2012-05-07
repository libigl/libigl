#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <cstdio>

int main(int argc,char * argv[])
{
  // Make GLUT window just to get a valid OpenGL context
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(640, 480);
  glutCreateWindow("Dummy");
  glutCreateMenu(NULL);

  printf("GL_VERSION: %s\n", 
    glGetString(GL_VERSION));
  printf("GL_SHADING_LANGUAGE_VERSION: %s\n", 
    glGetString(GL_SHADING_LANGUAGE_VERSION));
  return 0;
}
