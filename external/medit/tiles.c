#include "medit.h"
#include "sproto.h"
#include "extern.h"

extern ubyte ddebug;
ubyte  tiling;

#define IN2CM        2.54
#define CM2IN        0.3937


void dumpTile(char *data,int width,int height,GLubyte *buffer) {
  FILE *out2; 

  out2 = fopen(data,"w");
  fprintf(out2,"P6\n");
  fprintf(out2,"# Created using medit %s %s, (c) INRIA\n",ME_VER,ME_REL);
  fprintf(out2,"%d %d\n",width,height);
  fprintf(out2,"255\n");
  fwrite(buffer,sizeof(GLubyte),width*height*3,out2);
  fclose(out2);
}


/* dump big image */
int imgTiling(pScene sc,char *data,char key) {
  FILE     *out;
  pCamera   c = sc->camera;
  GLubyte  *tile,*buffer,*rowPtr;
  GLint     matmode,viewport[4];
  double    xmin,xmax,ymin,ymax,left,right,top,bottom,ratio;
  float     debhaut,finhaut,deblarg,finlarg,look[3];
  int       i,tw,th,col,row,nbcol,nbrow;
  int       imgWidth,imgHeight,tileWidth,tileHeight,tileWidthNB,tileHeightNB;
  int       bitsTileRow,bitsImgOffset,bitsCurTileRow,border;
  int       bitsTileOffset,bitsPixel,imgRowSize;
  char     *ptr,name[256];
  ubyte     bckbyte;
  static GLfloat up[3] = { 0.0, 1.0, 0.0};

  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  bckbyte = (ubyte)(255*(0.30*sc->par.back[0] + 
                         0.59*sc->par.back[1] + 
                         0.11*sc->par.back[2])+ 0.5);

  /* compute image size */
  ratio = (double)sc->par.xs/(double)sc->par.ys;
  ymax  = 0.1f*sc->dmax * tan(sc->persp->fovy * M_PI/360.0);
  ymin  = -ymax;
  xmin  = ymin * ratio;
  xmax  = ymax * ratio;

  bitsPixel  = 3 * sizeof(GLubyte);
  imgWidth   = (int)(sc->par.dpi*CM2IN*sc->par.cm+0.5);
  imgHeight  = (int)((double)imgWidth / ratio +0.5);
  imgRowSize = imgWidth * bitsPixel;

  tileWidth    = sc->par.xs;
  tileHeight   = sc->par.ys;
  border       = sc->mode & S_BDRY ? 1 : 0;
  tileWidthNB  = tileWidth  - 2*border;
  tileHeightNB = tileHeight - 2*border;

  if ( ddebug ) {
    fprintf(stdout,"   Generating %d by %d image\n",imgWidth,imgHeight);
    fprintf(stdout,"   tile %d x %d\n",tileWidth,tileHeight);
    fprintf(stdout,"   image size %f x %f cm\n",
            sc->par.cm,(double)imgHeight/sc->par.dpi/CM2IN);
  }

  /* buffer to store one tile */
  tile = (GLubyte*)calloc(tileWidthNB*tileHeightNB,bitsPixel);
  if ( !tile ) {
    fprintf(stderr,"  ## Unable to store buffer!\n");
    return(0);
  }
  
  /* buffer for a row of tiles */
  buffer = (GLubyte*)calloc(imgWidth*tileHeightNB,bitsPixel);
  if ( !buffer ) {
    free(tile);
    fprintf(stderr,"  ## Unable to store a row of tiles!\n");
    return(0);
  }
  
  /* open EPS file */
  strcpy(name,data);
  ptr  = (char*)strstr(name,".ps");
  if ( !ptr ) strcat(name,".ps");
  out = fopen(name,"w");
  if ( !out ) {
    fprintf(stderr,"  ## Unable to open file %s.\n",name);
    free(tile);
    return(0);
  }
  writeEPSheader(out,name,key,imgWidth,imgHeight,sc->par.cm,sc->par.dpi);
  
  /* save current viewport */
  glGetIntegerv(GL_VIEWPORT,viewport);
  glDrawBuffer(GL_BACK_LEFT);
  glReadBuffer(GL_BACK_LEFT);

  /* dump tiles */
  nbcol  = (int)((float)imgWidth / tileWidthNB)  + 1;
  nbrow  = (int)((float)imgHeight/ tileHeightNB) + 1;
  tiling = 1;  
  th = tileHeightNB;
  finhaut = sc->par.ys;

  for (row=nbrow-1; row>=0; row--) {
    if ( row < nbrow-1 )
      th = tileHeightNB;
    else
      th = imgHeight - row*tileHeightNB;
    debhaut = finhaut - sc->par.ys * th / imgHeight;
    deblarg = 0;
    
    for (col=0; col<nbcol; col++) {
      if ( col < nbcol-1 )
        tw = tileWidthNB;
      else
        tw = imgWidth - col*tileWidthNB;
      if ( th > tileHeightNB || tw > tileWidthNB ) {
        fprintf(stderr,"  %%%% Wrong tile size (%d,%d).\n",th,tw);
        free(buffer);
        free(tile);
        return(0); 
      }
      finlarg = deblarg + sc->par.xs * tw / imgWidth;

      /* set viewport to tilesize (with border) */
      glViewport(0,0,tw+2*border,th+2*border);

      /* current matrix */
      glGetIntegerv(GL_MATRIX_MODE,&matmode);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glMultMatrixf(sc->persp->matrix);

      /* compute projection parameters */
      if ( sc->persp->pmode == PERSPECTIVE ) {
        left   = xmin   + (xmax-xmin) * (col*tileWidthNB) / imgWidth;
        right  = left   + (xmax-xmin) * tw / imgWidth;
        bottom = ymin   + (ymax-ymin) * (row*tileHeightNB) / imgHeight;
        top    = bottom + (ymax-ymin) * th / imgHeight;

        glFrustum(left,right,bottom,top,0.1*sc->dmax,10.0f*sc->dmax);
        glTranslatef(0.0f,0.0,sc->persp->depth);
      }
      else if ( sc->persp->pmode == CAMERA ) {
        left   = xmin   + (xmax-xmin) * (col*tileWidthNB) / imgWidth;
        right  = left   + (xmax-xmin) * tw / imgWidth;
        bottom = ymin   + (ymax-ymin) * (row*tileHeightNB) / imgHeight;
        top    = bottom + (ymax-ymin) * th / imgHeight;

        glFrustum(left,right,bottom,top,0.1f*sc->dmax,10.0f*sc->dmax);

        look[0] = c->eye[0] + 0.001*sc->dmax*c->speed[0];
        look[1] = c->eye[1] + 0.001*sc->dmax*c->speed[1];
        look[2] = c->eye[2] + 0.001*sc->dmax*c->speed[2];
        gluLookAt(c->eye[0],c->eye[1],c->eye[2], 
                  look[0],look[1],look[2],
                  up[0],up[1],up[2]);
        glTranslatef(0.0f,0.0f,0.5*sc->persp->depth);
      }
      else if ( sc->persp->pmode == ORTHO ) {
        glOrtho(-1.,1.,-1.,0.1,0.01,0.01);
        glTranslatef(0.0,0.0,sc->persp->depth);
      }

      /* redraw scene */
      glDisable(GL_LIGHTING);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      /*gluLookAt(0.,0.,-sc->persp->depth, 0.,0.,0., 0.0,1.0,0.0);*/

      glMultMatrixf(sc->view->matrix);
      glTranslatef(sc->cx,sc->cy,sc->cz);
      drawModel(sc);

      /* read buffer */
      glFlush();
      glReadPixels(border,border,tileWidthNB,tileHeightNB,
                   GL_RGB,GL_UNSIGNED_BYTE,tile);

      /* store to row buffer */
      bitsImgOffset  = col * tileWidthNB * bitsPixel;
      bitsTileRow    = tileWidthNB * bitsPixel;
      bitsCurTileRow = tw * bitsPixel;
      bitsTileOffset = border*bitsPixel;
      bitsTileOffset = 0;
      for (i=0; i<th; i++) {
        memcpy(buffer + i*imgRowSize  + bitsImgOffset,
               tile   + i*bitsTileRow + bitsTileOffset,
               bitsCurTileRow);
      }
      deblarg = finlarg+1;
    }
    finhaut = debhaut-1;

    /* modify color */
    if ( sc->par.coeff > 0.0f ) {
      for (i=0; i<imgWidth*tileHeightNB*bitsPixel; i++)
        if ( buffer[i] > 10 )
          buffer[i] += (255-buffer[i]) * sc->par.coeff;
        else
          buffer[i] += buffer[i] * sc->par.coeff;
    }

    /* write row of tiles */
    for (i=0; i<th; i++) {
      /* reverse image */
      rowPtr = buffer + (th-1-i) * imgRowSize;
      writeEPSRow(out,key,rowPtr,imgWidth,bckbyte);
    }
  }

  writeEPStrailer(out);
  fclose(out);
  free(tile);
  free(buffer);
  tiling = 0;
  
  /* restore viewport */
  glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
  glDrawBuffer(GL_FRONT | GL_BACK);

  farclip(GL_TRUE);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if ( ddebug ) fprintf(stdout,"   Tiling completed.\n");
  return(1);
}
