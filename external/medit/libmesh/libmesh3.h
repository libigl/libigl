

/*----------------------------------------------------------*/
/*															*/
/*						LIBMESH V 3.0						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		handle .meshb file format I/O		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		aug 02 2003							*/
/*	Last modification:	jan 25 2006							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define LM_NBKW 79
#define LM_MESH_VERSION 1
#define LM_INT 1
#define LM_REAL 2
#define LM_READ 1
#define LM_WRITE 2
#define LM_ASCII  1
#define LM_BINARY 2
#define LM_MESH   4
#define LM_SOL    8
#define LM_SCALAR 1
#define LM_VECTOR 2
#define LM_SYM_MATRIX 3
#define LM_MATRIX 4

enum LM_kw_tags
{
	LM_Reserved1, \
	LM_MeshVersionFormatted, \
	LM_Reserved2, \
	LM_Dimension, \
	LM_Vertices, \
	LM_Edges, \
	LM_Triangles, \
	LM_Quadrilaterals, \
	LM_Tetrahedra, \
	LM_Pentahedra, \
	LM_Hexahedra, \
	LM_SubDomainFromGeom, \
	LM_SubDomainFromMesh, \
	LM_Corners, \
	LM_Ridges, \
	LM_RequiredVertices, \
	LM_RequiredEdges, \
	LM_RequiredTriangles, \
	LM_RequiredQuadrilaterals, \
	LM_TangentAtEdgeVertices, \
	LM_NormalAtVertices, \
	LM_NormalAtTriangleVertices, \
	LM_NormalAtQuadrilateralVertices, \
	LM_AngleOfCornerBound, \
	LM_Geometry, \
	LM_VertexOnGeometricVertex, \
	LM_VertexOnGeometricEdge, \
	LM_VertexOnGeometricTriangle, \
	LM_VertexOnGeometricQuadrilateral, \
	LM_EdgeOnGeometricEdge, \
	LM_TriangleOnGeometricTriangle, \
	LM_TriangleOnGeometricQuadrilateral, \
	LM_QuadrilateralOnGeometricTriangle, \
	LM_QuadrilateralOnGeometricQuadrilateral, \
	LM_MeshSupportOfVertices, \
	LM_VertexOnSupportVertex, \
	LM_VertexOnSupportEdge, \
	LM_VertexOnSupportTriangle, \
	LM_VertexOnSupportQuadrilateral, \
	LM_VertexOnSupportTetrahedron, \
	LM_VertexOnSupportPentahedron, \
	LM_VertexOnSupportHexahedron, \
	LM_CrackedEdges, \
	LM_CrackedTriangles, \
	LM_CrackedQuadrilaterals, \
	LM_EquivalentEdges, \
	LM_EquivalentTriangles, \
	LM_EquivalentQuadrilaterals, \
	LM_PhysicsReference, \
	LM_IncludeFile, \
	LM_BoundingBox, \
	LM_Identifier, \
	LM_IdentityOfGeometry, \
	LM_IdentityOfMeshSupport, \
	LM_End, \
	LM_Reserved10, \
	LM_Reserved11, \
	LM_Reserved12, \
	LM_Reserved13, \
	LM_Tangents, \
	LM_Normals, \
	LM_TangentAtVertices, \
	LM_SolAtVertices, \
	LM_SolAtEdges, \
	LM_SolAtTriangles, \
	LM_SolAtQuadrilaterals, \
	LM_SolAtTetrahedra, \
	LM_SolAtPentahedra, \
	LM_SolAtHexahedra, \
	LM_DSolAtVertices, \
	LM_ISolAtVertices, \
	LM_ISolAtEdges, \
	LM_ISolAtTriangles, \
	LM_ISolAtQuadrilaterals, \
	LM_ISolAtTetrahedra, \
	LM_ISolAtPentahedra, \
	LM_ISolAtHexahedra, \
	LM_Iterations, \
	LM_Time, \
	LM_VertexHack
};


/*----------------------------------------------------------*/
/* Structures												*/
/*----------------------------------------------------------*/

typedef struct
{
	/* Public */

	int dimension;
	int kw_counters[ LM_NBKW + 1 ];
	int *sol_headers[ LM_NBKW + 1 ];

	/* Private */

	int version, mode, type, endian, current_kw;
	FILE *handle;
	char *filename;
	size_t kw_pos[ LM_NBKW + 1 ][3];
}LM_mesh_struct;


/*----------------------------------------------------------*/
/* External tables and procedures							*/
/*----------------------------------------------------------*/
#ifdef __cplusplus
extern "C" {
#endif
    
extern int LM_open_mesh(const char *, int, LM_mesh_struct *, ...);
extern int LM_read_field(LM_mesh_struct *, int, int, void *);
extern int LM_write_field(LM_mesh_struct *, int, int, void *, ...);
extern int LM_read_line(LM_mesh_struct *, int, ...);
extern int LM_write_line(LM_mesh_struct *, int, ...);
extern int LM_close_mesh(LM_mesh_struct *);
extern char *LM_kw_table[][3];
#ifdef __cplusplus
}
#endif
    
