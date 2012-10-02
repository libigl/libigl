
/* libmesh5 example : transform a quadrilateral mesh into a triangular one */

#include <stdio.h>
#include <stdlib.h>
#include <libmesh5.h>


int main()
{
	int i, NmbVer, NmbQad, InpMsh, OutMsh, ver, dim, *RefTab, (*QadTab)[5];
	float (*VerTab)[3];


	/*----------------------------------------------*/
	/* Open the quadrilateral mesh file for reading */
	/*----------------------------------------------*/

	if(!(InpMsh = GmfOpenMesh("quad.mesh", GmfRead, &ver, &dim)))
		return(1);

	printf("InpMsh : idx = %d, version = %d, dimension = %d\n", InpMsh, ver, dim);

	if( (ver != GmfFloat) || (dim != 3) )
		exit(1);

	/* Get some stats to allocate memory */

	NmbVer = GmfStatKwd(InpMsh, GmfVertices);
	printf("InpMsh : nmb vertices = %d\n", NmbVer);
	VerTab = malloc((NmbVer+1) * 3 * sizeof(float));
	RefTab = malloc((NmbVer+1) * sizeof(int));

	NmbQad = GmfStatKwd(InpMsh, GmfQuadrilaterals);
	printf("InpMsh : nmb quads = %d\n", NmbQad);
	QadTab = malloc((NmbQad+1) * 5 * sizeof(int));

	/* Read the vertices */

	GmfGotoKwd(InpMsh, GmfVertices);

	for(i=1;i<=NmbVer;i++)
		GmfGetLin(InpMsh, GmfVertices, &VerTab[i][0], &VerTab[i][1], &VerTab[i][2], &RefTab[i]);

	/* Read the quads */

	GmfGotoKwd(InpMsh, GmfQuadrilaterals);

	for(i=1;i<=NmbQad;i++)
		GmfGetLin(InpMsh, GmfQuadrilaterals, &QadTab[i][0], &QadTab[i][1], &QadTab[i][2], &QadTab[i][3], &QadTab[i][4]);

	/* Close the quadrilateral mesh */

	GmfCloseMesh(InpMsh);


	/*--------------------------*/
	/* Create a triangular mesh */
	/*--------------------------*/

	if(!(OutMsh = GmfOpenMesh("tri.mesh", GmfWrite, ver, dim)))
		return(1);

	/* Set the number of vertices */

	GmfSetKwd(OutMsh, GmfVertices, NmbVer);

	/* Then write them down */

	for(i=1;i<=NmbVer;i++)
		GmfSetLin(OutMsh, GmfVertices, VerTab[i][0], VerTab[i][1], VerTab[i][2], RefTab[i]);

	/* Build two triangles from each quad */

	GmfSetKwd(OutMsh, GmfTriangles, 2*NmbQad);
	printf("OutMsh : nmb triangles = %d\n", 2*NmbQad);

	for(i=1;i<=NmbQad;i++)
	{
		GmfSetLin(OutMsh, GmfTriangles, QadTab[i][0], QadTab[i][1], QadTab[i][2], QadTab[i][4]);
		GmfSetLin(OutMsh, GmfTriangles, QadTab[i][0], QadTab[i][2], QadTab[i][3], QadTab[i][4]);
	}

	/* Don't forget to close the file */

	GmfCloseMesh(OutMsh);

	free(QadTab);
	free(RefTab);
	free(VerTab);

	return(0);
}
