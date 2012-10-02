#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "extern.h"
#include "sproto.h"


/* OpenGL's GL_3D_COLOR feedback vertex format */
typedef struct Feedback3Dcolor {
  GLfloat   x,y,z;
  GLfloat   r,g,b,a;
} Feedback3Dcolor;

typedef struct Token {
  GLfloat   depth;
  GLfloat  *ptr;
} Token;


static int compare(const void *a, const void *b) {
  Token   *p1 = (Token *)a;
  Token   *p2 = (Token *)b;
  GLfloat  diff = p2->depth - p1->depth;

  if ( diff > 0.0f )
    return 1;
  else if ( diff < 0.0f )
    return -1;
  else
    return 0;
}


void headps(FILE *file) {
  GLfloat  color[4],viewport[4];
  GLfloat  linewidth = 0.1;
  GLfloat  pointsize;
 

  glGetFloatv(GL_VIEWPORT,viewport);
  glGetFloatv(GL_COLOR_CLEAR_VALUE,color);
  glGetFloatv(GL_LINE_WIDTH,&linewidth);
  glGetFloatv(GL_POINT_SIZE,&pointsize);

  /* write EPS header */
  fprintf(file,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(file,"%%LanguageLevel: 1\n");
  fprintf(file,"%%%%Creator: Medit (via OpenGL feedback)\n");
  fprintf(file,"%%%%Pages: (atend)\n");
  fprintf(file,"%%%%BoundingBox: %g %g %g %g\n",
	  viewport[0],viewport[1],viewport[2],viewport[3]);
  fprintf(file,"%%EndComments\n\n");

  /* define macros */
  /* fprintf(file,"100 dict begin\n"); */
  fprintf(file,"%%%%BeginProlog\n");
  fprintf(file,"/sgm  {moveto lineto  stroke } def\n");  
  fprintf(file,"/stg {setgray} def\n");
  fprintf(file,"/stc {setrgbcolor} def\n");
  fprintf(file,"/mt  {moveto} def\n");
  fprintf(file,"/lt  {lineto} def\n");
  fprintf(file,"/dc {currentrgbcolor mark} def\n");
  fprintf(file,"/fc {stc newpath moveto\n");
  fprintf(file,"     counttomark 2 idiv {lineto} repeat\n");
  fprintf(file,"     closepath fill cleartomark  stc} def\n");
  fprintf(file,"%%%%EndProlog\n\n");

  fprintf(file,"gsave\n\n");

  /* set background color */
  fprintf(file,"%g %g %g setrgbcolor\n",color[0],color[1],color[2]);
  fprintf(file,"%g %g %g %g rectfill\n\n",
	  viewport[0],viewport[1],viewport[2],viewport[3]);
  fprintf(file, "%g setlinewidth\n\n",0.01*linewidth);
}

void tailps(FILE *file) {
  /* write EPS tail */
  fprintf(file,"\n");
  fprintf(file,"%%%%Trailer\n");
  fprintf(file,"%%%%Pages:        1\n");
  fprintf(file,"grestore\n\n");
  fprintf(file,"showpage\n\n");
  fclose(file);
}

int coreps(FILE *file,GLsizei size,GLfloat *buffer) {
  Token           *ttoken;
  pScene           sc;
  Feedback3Dcolor *v;
  GLfloat         *ptr;
  int              i,k,idw,token,nit,nbv,sorting=TRUE;

  /* default */
  if ( ddebug ) printf("dump EPS file\n");
  idw = currentScene();
  sc = cv.scene[idw];

  /* Count tokens */
  nit = 0;
  ptr = buffer;
  while ( *ptr && ptr < buffer+size ) {
    token = (int)*ptr;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      nit++;
      ptr += 14;
      break;

    case GL_POLYGON_TOKEN:
      nit++;
      ++ptr;
      nbv = *ptr;
      ptr += (7 * nbv);
      break;

    case GL_POINT_TOKEN:
      ptr += 7;
      break;

    case GL_BITMAP_TOKEN:
      puts("BITMAP");
      nit++;
      ++ptr;
      nbv = *ptr;
      printf("nbv = %d\n",nbv);
      break;

    default:
      fprintf(stdout,"  Unrecognized token %d\n",token);
      break;
    }
    ptr++;
  }

  if ( ddebug )
    printf("size = %d  ptr = %p  buffer = %p  -> size = %d\n",
	   size,ptr,buffer,ptr-buffer);

  /* allocate mem to store tokens */
  if ( ddebug ) printf("%d tokens found\n",nit);
  ttoken = (Token *)malloc((1+nit)*sizeof(struct Token));
  assert(ttoken);

  /* store tokens */
  nit = 0;
  ptr = buffer;
  while ( *ptr && ptr < buffer+size ) {
    ttoken[++nit].ptr = ptr;
    token = (int)*ptr;
    ptr++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      v = (Feedback3Dcolor *)ptr;
      ttoken[nit].depth = 0.5*(v[0].z+v[1].z);
      ptr += 14;
      break;

    case GL_POLYGON_TOKEN:
      nbv = *ptr;
      ptr++;
      v = (Feedback3Dcolor *)ptr;
      ttoken[nit].depth = v[0].z;
      for (k=1; k<nbv; k++)
	ttoken[nit].depth += v[k].z;
      ttoken[nit].depth /= nbv;
      ptr += (7 * nbv);
      break;

    case GL_POINT_TOKEN:
      ptr += 7;
      nit--;
      break;

    default:
      ptr += 4;
      break;
    }
  }

  /* Sort the primitives according to depth */
  if ( sc->mode == WIRE || sc->mode == WIRE+S_MATERIAL )
    sorting = FALSE;
  if ( ddebug ) printf("prim = %d  size = %d, %d\n",nit,ptr-buffer-1,size);

  if ( sorting == TRUE ) {
    if ( ddebug ) printf("start sorting %d tokens...\n",nit);
    qsort(ttoken+1,nit,sizeof(struct Token),compare);
    if ( ddebug ) printf("end sorting\n");
  }

  /* write tokens in EPS file */
  for (k=1; k<=nit; k++) {
    ptr = ttoken[k].ptr;
    if ( *ptr == 0 )  continue;
    token = *ptr;
    ptr++;
    switch(token) {
    case GL_LINE_RESET_TOKEN:
    case GL_LINE_TOKEN:
      v = (Feedback3Dcolor *)ptr;
      fprintf(file,"%g %g %g stc\n",v[0].r,v[0].g,v[0].b);
      fprintf(file,"%g %g %g %g sgm\n",v[0].x,v[0].y,v[1].x,v[1].y);
      ptr += 14;
      break;

    case GL_POLYGON_TOKEN:
      nbv = *ptr;
      ptr++;
      v = (Feedback3Dcolor *)ptr;

      /* draw filled polygon */
      if ( sorting == TRUE ) {
	/* fprintf(file,"1. stg\n"); */
	fprintf(file, "dc ");
	for (i=0; i<nbv; i++)
	  fprintf(file," %g %g",v[i].x,v[i].y);
	fprintf(file,"  %g %g %g fc\n",v[0].r,v[0].g,v[0].b);
      }

      /* draw polygon border */
      if ( sc->mode & S_BDRY ) {
	fprintf(file,"0. stg\n");
	for (i=0; i<nbv-1; i++)
	  fprintf(file,"%g %g %g %g sgm\n",v[i].x,v[i].y,v[i+1].x,v[i+1].y);
	fprintf(file,"%g %g %g %g sgm\n\n",v[nbv-1].x,v[nbv-1].y,v[0].x,v[0].y);
      }

      ptr += (7*nbv);
      break;

    case GL_POINT_TOKEN:
      v = (Feedback3Dcolor *)ptr;
      ptr += 7;
      break; 
    default:
      printf("problem: token %d\n",token);
    }
  }

  return(1);
}


int sftcpy(pScene sc,pMesh mesh) {
  FILE     *file;
  GLfloat  *fbbuffer;
  GLint     nvalues;
  GLsizei   size;
  char     *ptr,data[128];
  static int nfree=0;

  /* default */
  if ( ddebug ) printf("soft copy\n");

  /* get file name */
  strcpy(data,mesh->name);
  ptr = (char*)strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  nfree = filnum(data,nfree,"ps");
  if ( nfree == -1 )  return(0);

  /* open PS file */
  sprintf(data,"%s.%.3d.ps",data,nfree);
  file = fopen(data,"w");
  if ( !file ) {
    fprintf(stdout,"  Unable to open %s\n",data);
    return(0);
  }

  /* size for feedback buffer */
  size    =  0;
  nvalues = -1;
  do {
    size += 1024*1024;
    fbbuffer = (GLfloat *)calloc(1+size,sizeof(GLfloat));
    if ( !fbbuffer ) {
      return(0);
    }
    if ( ddebug ) printf("feedback pointer = %p\n",fbbuffer);

    /* draw scene in back buffer */
    glFeedbackBuffer(size,GL_3D_COLOR,fbbuffer);
    (void)glRenderMode(GL_FEEDBACK);
/*
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.,0.,-sc->persp->depth, 0.,0.,0., 0.0,1.0,0.0);

    setupView(sc);
    glMultMatrixf(sc->view->matrix);
    glTranslatef(sc->cx,sc->cy,sc->cz);
*/
    drawModel(sc);
    if ( sc->type & S_DECO )  redrawStatusBar(sc);

    /*drawModel(sc);*/
    nvalues = (GLint)glRenderMode(GL_RENDER);
    if ( nvalues < 0 )
      free(fbbuffer);
  }
  while ( nvalues < 0 );

  if ( nvalues < 1 ) {
    return(0);
  }
  else if ( ddebug )  printf("nvalues = %d  size = %d\n",nvalues,size);
  
  /* write EPS file */
  glutSetCursor(GLUT_CURSOR_WAIT);

  headps(file);
  coreps(file,size,fbbuffer);
  tailps(file);

  if ( ddebug ) fprintf(stdout,"%s written\n",data);
  glutSetCursor(GLUT_CURSOR_INHERIT);

  return(1);
}


#ifdef __cplusplus
}
#endif
