#ifndef IGL_TGA_H
#define IGL_TGA_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"
// See license in tga.cpp

/* tga.h - interface for TrueVision (TGA) image file loader */

#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "OpenGL_convenience.h"

namespace igl
{

typedef struct {

  GLsizei  width;
  GLsizei  height;
  GLint    components;
  GLenum   format;

  GLsizei  cmapEntries;
  GLenum   cmapFormat;
  GLubyte *cmap;

  GLubyte *pixels;
  
} gliGenericImage;

typedef struct {
  unsigned char idLength;
  unsigned char colorMapType;

  /* The image type. */
#define TGA_TYPE_MAPPED 1
#define TGA_TYPE_COLOR 2
#define TGA_TYPE_GRAY 3
#define TGA_TYPE_MAPPED_RLE 9
#define TGA_TYPE_COLOR_RLE 10
#define TGA_TYPE_GRAY_RLE 11
  unsigned char imageType;

  /* Color Map Specification. */
  /* We need to separately specify high and low bytes to avoid endianness
     and alignment problems. */
  unsigned char colorMapIndexLo, colorMapIndexHi;
  unsigned char colorMapLengthLo, colorMapLengthHi;
  unsigned char colorMapSize;

  /* Image Specification. */
  unsigned char xOriginLo, xOriginHi;
  unsigned char yOriginLo, yOriginHi;

  unsigned char widthLo, widthHi;
  unsigned char heightLo, heightHi;

  unsigned char bpp;

  /* Image descriptor.
     3-0: attribute bpp
     4:   left-to-right ordering
     5:   top-to-bottom ordering
     7-6: zero
     */
#define TGA_DESC_ABITS 0x0f
#define TGA_DESC_HORIZONTAL 0x10
#define TGA_DESC_VERTICAL 0x20
  unsigned char descriptor;

} TgaHeader;

typedef struct {
  unsigned int extensionAreaOffset;
  unsigned int developerDirectoryOffset;
#define TGA_SIGNATURE "TRUEVISION-XFILE"
  char signature[16];
  char dot;
  char null;
} TgaFooter;

IGL_INLINE extern gliGenericImage *gliReadTGA(FILE *fp, char *name, int hflip, int vflip);
IGL_INLINE int gli_verbose(int new_verbose);
IGL_INLINE extern int gliVerbose(int newVerbose);

IGL_INLINE void writeTGA( gliGenericImage* image, FILE *fp);



} // end of igl namespace

#ifndef IGL_STATIC_LIBRARY
#  include "tga.cpp"
#endif

#endif 
#endif 
