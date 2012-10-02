#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdarg.h>
#include "medit.h"
#include "extern.h"
#include "sproto.h"


static GLfloat IdMatrix[16] = {
   1.0, 0.0, 0.0, 0.0,
   0.0, 1.0, 0.0, 0.0,
   0.0, 0.0, 1.0, 0.0,
   0.0, 0.0, 0.0, 1.0
};

/* set font style and size */
void setFont(char* name,int size) {
  GLvoid *font_style = GLUT_BITMAP_HELVETICA_10;

  if ( !strcmp(name,"helvetica") ) {
    if (size == 12)
      font_style = GLUT_BITMAP_HELVETICA_12;
    else if (size == 18)
      font_style = GLUT_BITMAP_HELVETICA_18;
  }
  else if (strcmp(name, "times roman") == 0) {
    font_style = GLUT_BITMAP_TIMES_ROMAN_10;
    if (size == 24)
      font_style = GLUT_BITMAP_TIMES_ROMAN_24;
  }
  /*
  else if (strcmp(name, "8x13") == 0)
    font_style = GLUT_BITMAP_8_BY_13;
  */
  else if (strcmp(name, "9x15") == 0)
    font_style = GLUT_BITMAP_9_BY_15;
}

/* display string format at pos(x,y) */
void drwstr(GLuint x,GLuint y,char* format, ...) {
  va_list  args;
  char    *s,buffer[255];

  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);

  glRasterPos2i(x,y);
  for (s=buffer; *s; s++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*s);
}

void output2(GLfloat x,GLfloat y,char *format,...) {
  va_list  args;
  char    *s,buffer[255];

  /*strcpy(myerror.procname,"output2");*/
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);

  glRasterPos2f(x,y);
  for (s=buffer; *s; s++) {
    /*glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*s);*/
    /*glutBitmapCharacter(font_style,*s);*/
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,*s);
  }
}

void output3(GLfloat x,GLfloat y,GLfloat z,char *format,...) {
  va_list args;
  char buffer[255], *s;

  /*strcpy(myerror.procname,"output3");*/
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  glRasterPos3f(x,y,z);
  for (s=buffer; *s; s++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,*s);
}

/* color converter */
void hsvrgb(double *hsv,double *rgb) {
  double f,p,q,t;
  int    i;

  hsv[0] = ((int)hsv[0] % 360) / 60.;
  i = (int)floor((double)hsv[0]);    /* largest int <= h     */
  f = hsv[0] - i;                    /* fractional part of h */
  p = hsv[2] * (1.0 - hsv[1]);
  q = hsv[2] * (1.0 - (hsv[1] * f));
  t = hsv[2] * (1.0 - (hsv[1] * (1.0 - f)));

  switch(i) {
  case 0: rgb[0] = hsv[2]; rgb[1] = t;      rgb[2] = p; break;
  case 1: rgb[0] = q;      rgb[1] = hsv[2]; rgb[2] = p; break;
  case 2: rgb[0] = p;      rgb[1] = hsv[2]; rgb[2] = t; break;
  case 3: rgb[0] = p;      rgb[1] = q;      rgb[2] = hsv[2]; break;
  case 4: rgb[0] = t;      rgb[1] = p;      rgb[2] = hsv[2]; break;
  case 5: rgb[0] = hsv[2]; rgb[1] = p;      rgb[2] = q; break;
  }
}

/* transform: u = MxV */
void transformPoint(double u[4],float v[4],float m[16]) {
  u[0] = v[0] * m[0]  + v[1] * m[1]  + v[2] * m[2]  + v[3] * m[3];
  u[1] = v[0] * m[4]  + v[1] * m[5]  + v[2] * m[6]  + v[3] * m[7];
  u[2] = v[0] * m[8]  + v[1] * m[9]  + v[2] * m[10] + v[3] * m[11];
  u[3] = v[0] * m[12] + v[1] * m[13] + v[2] * m[14] + v[3] * m[15];
}

void transformPoint2(double u[4],float v[4],float m[16]) {
  u[0] = v[0] * m[0] + v[1] * m[4] + v[2] * m[8]  + v[3] * m[12];
  u[1] = v[0] * m[1] + v[1] * m[5] + v[2] * m[9]  + v[3] * m[13];
  u[2] = v[0] * m[2] + v[1] * m[6] + v[2] * m[10] + v[3] * m[14];
  u[3] = v[0] * m[3] + v[1] * m[7] + v[2] * m[11] + v[3] * m[15];
}

void transformPointd(double u[3],double v[3],double m[16]) {
  u[0] = v[0] * m[0]  + v[1] * m[1]  + v[2] * m[2];
  u[1] = v[0] * m[4]  + v[1] * m[5]  + v[2] * m[6];
  u[2] = v[0] * m[8]  + v[1] * m[9]  + v[2] * m[10];
}

void transformVector(float u[4],float v[4],float m[16]) {
  u[0] = v[0] * m[0] + v[1] * m[4] + v[2] * m[8];
  u[1] = v[0] * m[1] + v[1] * m[5] + v[2] * m[9];
  u[2] = v[0] * m[2] + v[1] * m[6] + v[2] * m[10];
  u[3] = v[0] * m[3] + v[1] * m[7] + v[2] * m[11];
}


/* p = axb */
void multMatrix(GLfloat *p,GLfloat *a,GLfloat *b) {
  GLint i,row;

  for (i=0; i<4; i++) {
    row = i*4;
    p[row+0] = a[row] * b[0] + a[row+1] * b[4] + a[row+2] * b[8]  + a[row+3] * b[12];
    p[row+1] = a[row] * b[1] + a[row+1] * b[5] + a[row+2] * b[9]  + a[row+3] * b[13];
    p[row+2] = a[row] * b[2] + a[row+1] * b[6] + a[row+2] * b[10] + a[row+3] * b[14];
    p[row+3] = a[row] * b[3] + a[row+1] * b[7] + a[row+2] * b[11] + a[row+3] * b[15];
  }
}

void rotateMatrix(GLfloat angle,GLfloat x,GLfloat y,GLfloat z,GLfloat rm[16]) {
   GLfloat mag,s,c;
   GLfloat xx,yy,zz,xy,yz,zx,xs,ys,zs,one_c;

   if ( angle == 0.0f ) {
     memcpy(rm,IdMatrix,16*sizeof(GLfloat));
     return;
   }
   mag = x*x + y*y + z*z;

   if ( mag == 0.0f ) {
     memcpy(rm,IdMatrix,16*sizeof(GLfloat));
     return;
   }
   mag = 1.0f / sqrt(mag);
   x *= mag;    
   y *= mag;    
   z *= mag;
   s  = sin(angle * DTOR);
   c  = cos(angle * DTOR);
   xx = x*x;  yy = y*y;  zz = z*z;
   xy = x*y;  yz = y*z;  zx = z*x;
   xs = x*s;  ys = y*s;  zs = z*s;
   one_c = 1.0f - c;

   rm[0] = (one_c * xx) + c;
   rm[1] = (one_c * xy) - zs;
   rm[2] = (one_c * zx) + ys;
   rm[3]  = 0.0f;

   rm[4] = (one_c * xy) + zs;
   rm[5] = (one_c * yy) + c;
   rm[6] = (one_c * yz) - xs;
   rm[7] = 0.0f;

   rm[8] = (one_c * zx) - ys;
   rm[9] = (one_c * yz) + xs;
   rm[10] = (one_c * zz) + c;
   rm[11] = 0.0f;

   rm[12] = rm[13] = rm[14] = 0.0f;
   rm[15] = 1.0f;
}

int invertMatrix(float src[16],float inverse[16]) {
  double  t;
  int     i, j, k, swap;
  double  tmp[4][4];

  memcpy(inverse,IdMatrix,16*sizeof(GLfloat));
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      tmp[i][j] = src[i*4+j];

  for (i=0; i<4; i++) {
    /* look for largest element in column. */
    swap = i;
    for (j=i+1; j<4; j++) {
      if ( fabs(tmp[j][i]) > fabs(tmp[i][i]) )
	swap = j;
    }
    if ( swap != i ) {
      /* swap rows. */
      for (k=0; k<4; k++) {
	t            = tmp[i][k];
	tmp[i][k]    = tmp[swap][k];
	tmp[swap][k] = t;
	t                 = inverse[i*4+k];
	inverse[i*4+k]    = inverse[swap*4+k];
	inverse[swap*4+k] = t;
      }
    }
    /* The matrix is singular. */
    if ( tmp[i][i] == 0 )
      return(0);
    
    t = tmp[i][i];
    for (k=0; k<4; k++) {
      tmp[i][k]      /= t;
      inverse[i*4+k] /= t;
    }
    for (j=0; j<4; j++) {
      if ( j != i ) {
	t = tmp[j][i];
	for (k=0; k<4; k++) {
	  tmp[j][k]      -= tmp[i][k]*t;
	  inverse[j*4+k] -= inverse[i*4+k]*t;
	}
      }
    }
  }

  return(1);
}

void print_matrix(const GLfloat m[16],const char *ligne) {
  int i;

  printf("---- %s ----\n",ligne);
  for (i=0; i<4; i++)
    printf("%f %f %f %f\n",m[i],m[4+i],m[8+i],m[12+i]);
  printf("---------------------------------\n");
}

void print_matrixd(const GLdouble m[16],const char *ligne) {
  int i;

  printf("---- %s ----\n",ligne);
  for (i=0; i<4; i++)
    printf("%f %f %f %f\n",m[i],m[4+i],m[8+i],m[12+i]);
  printf("---------------------------------\n");
}
int filnum(char *data,int numdep,char *ext) {
  FILE  *in;
  char   tmpstr[256];

  do {
    sprintf(tmpstr,"%s.%.3d.%s",data,numdep,ext);
    in = fopen(tmpstr,"r");
    if ( !in ) return(numdep);
    fclose(in);
  }
  while ( ++numdep < 999 );
  return(-1);
}

#ifdef __cplusplus
}
#endif
