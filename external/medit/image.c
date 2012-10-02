#include "medit.h"
#include "sproto.h"
#include "extern.h"


PPMimage *loadPPM(const char *imgname,int *type) {
  pPPMimage  result;
  FILE      *fp;
  int        i,k,typimg,ret,r,g,b,s,maxval,bitsize;
  char      *ptr,c,buff[1024],data[256];

  /* search for image */
  ptr = strstr(imgname,".ppm");
  if ( !ptr ) {
    strcpy(data,imgname);
    strcat(data,".ppm");
    fp = fopen(data,"rb");
  } 
  else
    fp = fopen(imgname,"rb");
  if ( !fp ) {
    fprintf(stderr,"  ## Unable to open file %s.\n",imgname);
    return(0);
  }

  if ( !fgets(buff,sizeof(buff),fp) ) {
    fprintf(stderr,"  ## Invalid file header.\n");
    return(0);
  }

  /* check header file */
  if ( buff[0] != 'P' ) {
    fprintf(stderr,"  ## Invalid image format.\n");
    return(0);
  }
  
  switch(buff[1]) {
  case '2': typimg = P2;  break;
  case '3': typimg = P3;  break;
  case '5': typimg = P5;  break;
  case '6': typimg = P6;  break;
  default:
    fprintf(stderr,"  ## Invalid image format.\n");
    return(0);
  }
  
  /* allocate memory to store imagee */
  result = malloc(sizeof(PPMimage));
  if ( !result ) {
    fprintf(stderr,"  ## Unable to load image.\n");
    return(0);
  }

  do {
    ret = fscanf(fp,"%s",buff);
    if ( ret == EOF ) break;
    /* check and strip comments */
    if ( buff[0] == '#' )
      do
        c = getc(fp);
      while ( c != '\n' );
    else break;
  }
  while (1);

  /* read columns + lines */
  ret  = sscanf(buff,"%d",&s);
  result->sizeX = (short)s;
  ret += fscanf(fp,"%d",&s);
  result->sizeY = (short)s;
  if ( ret != 2 ) {
    fprintf(stderr,"  ## Error loading image.\n");
    free(result);
    return(0);
  }
  if ( !quiet )
    fprintf(stdout,"  image size:   %d x %d\n",result->sizeX,result->sizeY);

  if ( fscanf(fp,"%d",&maxval) != 1 ) {
    fprintf(stderr,"  ## Invalid image size.\n");
    free(result);
    return(0);
  }

  /* strip line */
  while ( fgetc(fp) != '\n' ) ;

  /* size based on type */
  if ( typimg == P2 || typimg == P5 )
    bitsize = result->sizeX*result->sizeY;
  else
    bitsize = 3*result->sizeX*result->sizeY;
  result->data = (ubyte*)malloc(bitsize*sizeof(ubyte));
  if ( !result ) {
    fprintf(stderr,"  ## Unable to load image.\n");
    free(result);
    return(0);
  }

  /* read data file */
  switch( typimg ) {
  case P2:  /* ascii file (grey)  */
  case P3:  /* ascii file (color) */
    for (i=0; i<bitsize; i++) {
      fscanf(fp,"%d",&r);
      result->data[i] = (ubyte)r;
    }
    break;
    
  case P5:  /* binary file (grey) */
  case P6:  /* binary file (color) */
    ret = fread(result->data,sizeof(ubyte),bitsize,fp);
    if ( ret != bitsize ) {
      fprintf(stderr,"  ## Error loading image.\n");
      free(result->data);
      free(result);
      return(0);
    }
    break;
  }
  fclose(fp);

  if ( *type == DEFAULT )
    switch( typimg ) {
     case P2:
     case P5:
       *type = GREY;  break;
     case P3:
     case P6:
       *type = RGB;  break;
    }

  /* convert to grey levels */
  else if ( *type == GREY && (typimg == P3 || typimg == P6) ) {
    fprintf(stdout,"  converting to grey levels\n");
    for (i=0,k=0; i<bitsize; i+=3,k++) {
      r = (int)result->data[i];
      g = (int)result->data[i+1];
      b = (int)result->data[i+2];
      result->data[k] = (ubyte)(0.3*r+0.59*g+0.11*b);
    }
    result->data = (ubyte*)realloc(result->data,sizeof(ubyte)*bitsize/3);
  }
  
  return(result);
}


int savePPM(const char *imgname,pPPMimage img,int typimg) {
  FILE      *out;
  int        i,c,bitsize;
  char      *ptr,data[512];

  strcpy(data,imgname);
  ptr  = (char*)strstr(data,".ppm");
  if ( !ptr ) strcat(data,".ppm");
  out = fopen(data,"w");
  if ( !out ) {
    fprintf(stderr,"  ## Unable to open file %s.\n",data);
    exit(1);
  }

  bitsize = img->sizeX*img->sizeY;
  switch(typimg) {
  case P2:
    fprintf(out,"P2\n");
    fprintf(out,"# Created using medit %s %s, (c) INRIA\n",ME_VER,ME_REL);
    fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
    fprintf(out,"255\n");
    c = 0;
    for (i=0; i<img->sizeX*img->sizeY; i++) {
      fprintf(out,"%3d ",(int)img->data[i]);
      if ( ++c == 17 ) { 
        c = 0; 
        fprintf(out,"\n");
      }
    }
    fprintf(out,"\n");
    break;
  case P5:
    fprintf(out,"P5\n");
    fprintf(out,"# Created using medit %s %s, (c) INRIA\n",ME_VER,ME_REL);
    fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
    fprintf(out,"255\n");
    fwrite(img->data,sizeof(ubyte),bitsize,out);
    break;
  case P6:
    fprintf(out,"P6\n");
    fprintf(out,"# Created using medit %s %s, (c) INRIA\n",ME_VER,ME_REL);
    fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
    fprintf(out,"255\n");
    bitsize = (img->sizeX*24+7)/8*img->sizeY;
    if ( fwrite(img->data,sizeof(ubyte),bitsize,out) < bitsize )
      fprintf(stderr,"  ## Data file corrupted.\n");
    break;
  }
  fclose(out);

  return(1);
}


int saveTGA(const char *imgname,GLubyte *img,int w,int h) {
  TGAheader  tga;
  int        i;
  char      *ptr,data[256];
  FILE      *out;
  
  strcpy(data,imgname);
  ptr  = (char*)strstr(data,".tga");
  if ( !ptr ) strcat(data,".tga");
  out = fopen(data,"wb");
  if ( !out ) {
    fprintf(stderr,"  ## UNABLE TO OPEN FILE %s.\n",data);
    exit(1);
  }

  tga.idfield_len = 0;
  tga.cmap_type   = 0;
  tga.image_type  = 2;
  for (i=0; i<5; i++)
    tga.cmap_spec[i] = 0;
  for (i=0; i<2; i++) {
    tga.x_orig[i] = 0;
    tga.y_orig[i] = 0;
  }
  /* Lo bits */
  tga.width[0] = w & 0xFF;
  /* Hi bits */
  tga.width[1] = (w >> 8) & 0xFF;
  tga.height[0] = h & 0xFF;
  tga.height[1] = (h >> 8) & 0xFF;
  tga.pixel_size = 24;
  tga.image_desc = 0;
  /* Output header */
  fwrite(&tga,sizeof(TGAheader),1,out);
  
  /* Output image */
  fwrite(img, sizeof(unsigned char),w*h*3,out);
  return(1);
}

void swapPixels(PPMimage *pixels) {
  GLubyte   *row;
  int        i,k,ck,bits;

  bits = 3*pixels->sizeX;
  row  = (GLubyte*)malloc(bits*sizeof(GLubyte));
  if ( !row ) return;

  /* exchange rows */
  for (i=0; i<pixels->sizeY/2; i++) {
    k  = 3*i*pixels->sizeX;
    ck = 3*(pixels->sizeY-i-1)*pixels->sizeX;
    memcpy(row,&pixels->data[k],bits);
    memcpy(&pixels->data[k],&pixels->data[ck],bits);
    memcpy(&pixels->data[ck],row,bits);
  }
  free(row);
}


int imgHard(pScene sc,char *data,char key) {
  PPMimage  *pixels;
  GLint      viewport[4];
  pPersp     p;
  int        xx0,yy0,ww,hh;

  pixels = (PPMimage*)M_malloc(sizeof(PPMimage),"imgHard");
  if ( !pixels ) {
    fprintf(stderr,"  ## UNABLE TO ALLOCATE MEMORY FOR IMAGE.\n");
    return(0);
  }

  p = sc->persp;
  if ( abs(p->rubfx-p->rubix) > 0 ) {
    ww  = abs(p->rubfx-p->rubix)-2;
    hh  = abs(p->rubfy-p->rubiy)-2;
    xx0 = min(p->rubix,p->rubfx)+1;
    yy0 = min(p->rubiy,p->rubfy)+1;
  }
  else {
    glGetIntegerv(GL_VIEWPORT,viewport);
    ww  = viewport[2];
    hh  = viewport[3];
    xx0 = 0;
    yy0 = 0;
  }
  
  /* align to 8 bytes */
  ww = ww & ~7;
  hh = hh & ~7;

  pixels->sizeX = (short)ww;
  pixels->sizeY = (short)hh;
  pixels->data = (ubyte*)M_calloc(3*ww*hh,sizeof(ubyte),"imgHard.data");
  if ( !pixels->data ) {
    fprintf(stderr,"  ## Not enough memory to save image.\n");
    M_free(pixels);
    return(0);
  }

  if ( ddebug ) fprintf(stdout,"size %d x %d\n",ww,hh);
  
  glFinish();
  if ( saveimg )
    glReadBuffer(GL_BACK_LEFT);
  else
    glReadBuffer(GL_FRONT);
  glPixelStorei(GL_PACK_ALIGNMENT,4);
  glPixelStorei(GL_PACK_ROW_LENGTH,0);
  glPixelStorei(GL_PACK_SKIP_ROWS,0);
  glPixelStorei(GL_PACK_SKIP_PIXELS,0);

  glReadPixels(xx0,yy0,ww,hh,GL_RGB,GL_UNSIGNED_BYTE,pixels->data);
  if ( glGetError() != GL_NO_ERROR ) {
    fprintf(stderr,"  ## Unable to save image\n");
    M_free(pixels->data);
    M_free(pixels);
    return(0);
  }

  if ( key == 'H' ) {
    swapPixels(pixels);
    savePPM(data,pixels,P6);
  }
  else if ( key == 'T' )
    saveTGA(data,pixels->data,ww,hh);

  M_free(pixels->data);
  M_free(pixels);
  return(1);
}
