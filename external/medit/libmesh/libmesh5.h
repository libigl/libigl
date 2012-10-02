

/*----------------------------------------------------------*/
/*															*/
/*						LIBMESH V 5.0						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		handle .meshb file format I/O		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 16 2007							*/
/*	Last modification:	apr 10 2007							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define GmfStrSiz 1024
#define GmfMaxTyp 20
#define GmfMaxKwd 79
#define GmfMshVer 1
#define GmfRead 1
#define GmfWrite 2
#define GmfSca 1
#define GmfVec 2
#define GmfSymMat 3
#define GmfMat 4
#define GmfFloat 1
#define GmfDouble 2

enum GmfKwdCod
{
	GmfReserved1, \
	GmfVersionFormatted, \
	GmfReserved2, \
	GmfDimension, \
	GmfVertices, \
	GmfEdges, \
	GmfTriangles, \
	GmfQuadrilaterals, \
	GmfTetrahedra, \
	GmfPentahedra, \
	GmfHexahedra, \
	GmfReserved3, \
	GmfReserved4, \
	GmfCorners, \
	GmfRidges, \
	GmfRequiredVertices, \
	GmfRequiredEdges, \
	GmfRequiredTriangles, \
	GmfRequiredQuadrilaterals, \
	GmfTangentAtEdgeVertices, \
	GmfNormalAtVertices, \
	GmfNormalAtTriangleVertices, \
	GmfNormalAtQuadrilateralVertices, \
	GmfAngleOfCornerBound, \
	GmfReserved5, \
	GmfReserved6, \
	GmfReserved7, \
	GmfReserved8, \
	GmfReserved9, \
	GmfReserved10, \
	GmfReserved11, \
	GmfReserved12, \
	GmfReserved13, \
	GmfReserved14, \
	GmfReserved15, \
	GmfReserved16, \
	GmfReserved17, \
	GmfReserved18, \
	GmfReserved19, \
	GmfReserved20, \
	GmfReserved21, \
	GmfReserved22, \
	GmfReserved23, \
	GmfReserved24, \
	GmfReserved25, \
	GmfReserved26, \
	GmfReserved27, \
	GmfReserved28, \
	GmfReserved29, \
	GmfReserved30, \
	GmfBoundingBox, \
	GmfReserved31, \
	GmfReserved32, \
	GmfReserved33, \
	GmfEnd, \
	GmfReserved34, \
	GmfReserved35, \
	GmfReserved36, \
	GmfReserved37, \
	GmfTangents, \
	GmfNormals, \
	GmfTangentAtVertices, \
	GmfSolAtVertices, \
	GmfSolAtEdges, \
	GmfSolAtTriangles, \
	GmfSolAtQuadrilaterals, \
	GmfSolAtTetrahedra, \
	GmfSolAtPentahedra, \
	GmfSolAtHexahedra, \
	GmfDSolAtVertices, \
	GmfISolAtVertices, \
	GmfISolAtEdges, \
	GmfISolAtTriangles, \
	GmfISolAtQuadrilaterals, \
	GmfISolAtTetrahedra, \
	GmfISolAtPentahedra, \
	GmfISolAtHexahedra, \
	GmfIterations, \
	GmfTime, \
	GmfReserved38
};


/*----------------------------------------------------------*/
/* External procedures										*/
/*----------------------------------------------------------*/

extern int GmfOpenMesh(const char *, int, ...);
extern int GmfCloseMesh(int);
extern int GmfStatKwd(int, int, ...);
extern int GmfGotoKwd(int, int);
extern int GmfSetKwd(int, int, ...);
extern void GmfGetLin(int, int, ...);
extern void GmfSetLin(int, int, ...);


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif

int call(gmfopenmeshf77)(char *, int *, int *, int *, int);
int call(gmfclosemeshf77)(int *);
int call(gmfstatkwdf77)(int *, int *, int *, int *, int *);
int call(gmfgotokwdf77)(int *, int *);
int call(gmfsetkwdf77)(int *, int *, int *, int *, int *);
int call(gmfgetvertex2df77)(int *, float *, float *, int *);
int call(gmfgetvertex3df77)(int *, float *, float *, float *, int *);
int call(gmfsetvertex2df77)(int *, float *, float *, int *);
int call(gmfsetvertex3df77)(int *, float *, float *, float *, int *);
int call(gmfgettrianglef77)(int *, int *, int *, int *, int *);
int call(gmfsettrianglef77)(int *, int *, int *, int *, int *);
int call(gmfgettetrahedronf77)(int *, int *, int *, int *, int *, int *);
int call(gmfsettetrahedronf77)(int *, int *, int *, int *, int *, int *);
int call(gmfgetedgef77)(int *, int *, int *, int *);
int call(gmfsetedgef77)(int *, int *, int *, int *);
int call(gmfgetquadrilateralf77)(int *, int *, int *, int *, int *, int *);
int call(gmfsetquadrilateralf77)(int *, int *, int *, int *, int *, int *);
int call(gmfgethexahedronf77)(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
int call(gmfsethexahedronf77)(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
int call(gmfsetsolf77)(int *, int *, int *);
int call(gmfgetsolf77)(int *, int *, int *);


/*----------------------------------------------------------*/
/* Transmesh private API									*/
/*----------------------------------------------------------*/

#ifdef TRANSMESH

extern char *KwdFmt[ GmfMaxKwd + 1 ][3];
extern int GmfCpyLin(int, int, int);

#endif
