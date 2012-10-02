#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>

#ifndef ubyte
typedef unsigned char  ubyte;
#endif

#define CM2IN     0.3937


void writeEPSheader(FILE *out,char *data,char key,int ww,int hh,float cm,float dpi) {
  int wpad;
  
  fprintf(out,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(out,"%%LanguageLevel: 1\n");
  fprintf(out,"%%%%Title: %s\n",data);
  fprintf(out,"%%%%Creator: medit -  (C) INRIA-Rocquencourt, 1999-2001\n");
  fprintf(out,"%%%%BoundingBox: 50 50 %d %d\n",
               (int)(ww*72./dpi+50.5),(int)(hh*72.0/dpi+50.5));
  fprintf(out,"%%%%Pages: (atend)\n");
  fprintf(out,"%%DocumentFonts:\n");
  fprintf(out,"%%%%EndComments\n");
  fprintf(out,"%%%%EndProlog\n");
  fprintf(out,"\n%%%%Page: 1 1\n");
  fprintf(out,"\n%% remember original state\n");
  fprintf(out,"/origstate save def\n");
  fprintf(out,"\n%% build a temporary dictionary\n");
  fprintf(out,"20 dict begin\n");

  fprintf(out,"\n%% define space for color conversions\n");
  fprintf(out,"/grays %d string def  %% space for gray scale line\n",ww);
  fprintf(out,"/npixls 0 def\n");
  fprintf(out,"/rgbindx 0 def\n");
  fprintf(out,"\n%% lower left corner\n");
  fprintf(out,"50 50 translate\n");
  
  fprintf(out,"\n%% size of image (on paper, in 1/72inch coords)\n");
  fprintf(out,"%g %g scale\n",72.0*CM2IN*cm,72.0*CM2IN*cm*(float)hh/ww);
  fprintf(out,"\n%% define string to hold a scanline's worth of data\n");

  /* write BW image */
  if ( key == 'B') {
    wpad = (ww+7) & ~7;
    fprintf(out,"/pix %d string def\n",wpad/8);
    fprintf(out,"\n%% dimensions of data\n");
    fprintf(out,"%d %d 1\n",ww,hh);
    fprintf(out,"\n%% mapping matrix\n");
    fprintf(out,"[%d 0 0 %d 0 %d]\n",ww,-hh,hh);
    fprintf(out,"\n{currentfile pix readhexstring pop}\n");
    fprintf(out,"image\n");
  }
  /* write greyscale image */
  else if ( key == 'G' ) {
    fprintf(out,"/pix %d string def\n",ww);
    fprintf(out,"\n%% dimensions of data\n");
    fprintf(out,"%d %d 8\n",ww,hh);
    fprintf(out,"\n%% mapping matrix\n");
    fprintf(out,"[%d 0 0 %d 0 %d]\n",ww,-hh,hh);
    fprintf(out,"\n{currentfile pix readhexstring pop}\n");
    fprintf(out,"image\n");
  }
  /* color image */
  else if ( key == 'C' ) {
    fprintf(out,"/pix %d string def\n",3*ww);

    fprintf(out,"\n%% dimensions of data\n");
    fprintf(out,"%d %d 8\n",ww,hh);
    fprintf(out,"\n%% mapping matrix\n");
    fprintf(out,"[%d 0 0 %d 0 %d]\n",ww,-hh,hh);

    fprintf(out,"\n%% define 'colorimage' if it isn't defined\n");
    fprintf(out,"/colorimage where   %% do we know about 'colorimage'?\n");
    fprintf(out,"  { pop }           %% yes: pop off the 'dict' returned\n");
    fprintf(out,"  {                 %% no:  define one\n");
    fprintf(out,"    /colortogray {  %% define an RGB->I function\n");
    fprintf(out,"      /rgbdata exch store    %% call input 'rgbdata'\n");
    fprintf(out,"      rgbdata length 3 idiv\n");
    fprintf(out,"      /npixls exch store\n");
    fprintf(out,"      /rgbindx 0 store\n");
    fprintf(out,"      0 1 npixls 1 sub {\n");
    fprintf(out,"        grays exch\n");
    fprintf(out,"        rgbdata rgbindx       get 20 mul    %% Red\n");
    fprintf(out,"        rgbdata rgbindx 1 add get 32 mul    %% Green\n");
    fprintf(out,"        rgbdata rgbindx 2 add get 12 mul    %% Blue\n");
    fprintf(out,"        add add 64 idiv     %% I = .5G + .31R + .18B\n");
    fprintf(out,"        put\n");
    fprintf(out,"       /rgbindx rgbindx 3 add store\n");
    fprintf(out,"      } for\n");
    fprintf(out,"      grays 0 npixls getinterval\n");
    fprintf(out,"    } bind def\n");

    fprintf(out,"\n    %% Utility procedure for colorimage operator.\n");
    fprintf(out,"    %% This procedure takes two procedures off the\n");
    fprintf(out,"    %% stack and merges them into a single procedure.\n");

    fprintf(out,"\n    /mergeprocs { %% def\n");
    fprintf(out,"      dup length\n");
    fprintf(out,"      3 -1 roll\n");
    fprintf(out,"      dup\n");
    fprintf(out,"      length\n");
    fprintf(out,"     dup\n");
    fprintf(out,"      5 1 roll\n");
    fprintf(out,"      3 -1 roll\n");
    fprintf(out,"      add\n");
    fprintf(out,"      array cvx\n");
    fprintf(out,"      dup\n");
    fprintf(out,"      3 -1 roll\n");
    fprintf(out,"      0 exch\n");
    fprintf(out,"      putinterval\n");
    fprintf(out,"      dup\n");
    fprintf(out,"      4 2 roll\n");
    fprintf(out,"      putinterval\n");
    fprintf(out,"    } bind def\n");

    fprintf(out,"    /colorimage { %% def\n");
    fprintf(out,"      pop pop     %% remove 'false 3' operands\n");
    fprintf(out,"      {colortogray} mergeprocs\n");
    fprintf(out,"      image\n");
    fprintf(out,"    } bind def\n");
    fprintf(out,"  } ifelse          %% end of 'false' case\n");
    fprintf(out,"\n{currentfile pix readhexstring pop}\n");
    fprintf(out,"false 3 colorimage\n\n");
  }
}

void writeEPStrailer(FILE *out) {
  fprintf(out,"\nshowpage\n");
  fprintf(out,"\n%% stop using temporary dictionary\n");
  fprintf(out,"end\n");
  fprintf(out,"\n%% restore original state\n");
  fprintf(out,"origstate restore\n");
  fprintf(out,"\n%%%%Trailer\n");
}

void writeEPSRow(FILE *out,char key,ubyte *buffer,int size,ubyte bckbyt) {
  int    c,k,l;
  ubyte  byte,bbyte;

  l = 0;
  switch (key) {
  case 'B':  /* black & white */
    byte = 0x0;
    c    = 0;
    for (k=0; k<3*size; k+=3) {
      bbyte = (ubyte)(0.30*buffer[k]   + 0.59*buffer[k+1] 
                    + 0.11*buffer[k+2] + 0.5);
      if ( bbyte == bckbyt )  byte |= (1 << (7-c));
      /*
      if ( bbyte > 253 )  byte |= (1 << (7-c));
      
      if ( bbyte && bbyte != 255 )  {
        printf("buffer %d %d %d\n",buffer[k],buffer[k+1],buffer[k+2]);
        printf("bbyte  %d  byte %d\n",bbyte,byte);
        exit(1);
      }
      */
      if ( ++c == 8 ) {
        fprintf(out,"%.2x",byte);
        if ( ++l == 36 ) {
          fprintf(out,"\n");
          l = 0;
        }
        byte = 0x0;
        c = 0;
      }
    }
    /* padding */
    if ( c ) {
      for (l=8; l>c; l--) byte |= (1 << l);
      fprintf(out,"%.2x",byte);
    }
    break;
  case 'G':  /* greyscale */
    for (k=0; k<3*size; k+=3) {
      byte = (ubyte)(0.30*buffer[k] + 0.59*buffer[k+1] + 0.11*buffer[k+2]);
      fprintf(out,"%.2x",byte);
      if ( ++l == 36 ) {
        fprintf(out,"\n");
        l = 0;
      }
    }
    break;
  case 'C':  /* color scale */
    for (k=0; k<3*size; k++) {
      fprintf(out,"%.2x",buffer[k]);
      if ( ++l == 36 ) {
        fprintf(out,"\n");
        l = 0;
      }
    }
    break;
  }
  fprintf(out,"\n");
}


#ifdef __cplusplus
}
#endif
