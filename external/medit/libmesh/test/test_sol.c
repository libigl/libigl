/* Extraction d'un champ de solution en utilisant la libmesh5 */

#include <stdio.h>
#include <libmesh5.h>

int main()
{
	int i, NmbSol, InpSol, OutSol, ver, dim, NmbTyp, SolSiz, TypTab[20];
	float InpSolTab[100], OutSolTab[100];

	/* Ouverture du fichier .sol a lire */

	if(!(InpSol = GmfOpenMesh("test.sol", GmfRead, &ver, &dim)))
		return(1);

	printf("\nInpSol : idx = %d, version = %d, dim = %d\n", InpSol, ver, dim);

	NmbSol = GmfStatKwd(InpSol, GmfSolAtVertices, &NmbTyp, &SolSiz, TypTab);

	printf("InpSol : NmbSol = %d, NmbTyp = %d, siz = %d floats\n\n", NmbSol,NmbTyp,SolSiz);

	if(!(OutSol = GmfOpenMesh("out.sol", GmfWrite, ver, dim)))
		return(1);

	printf("OutSol : idx = %d\n\n", OutSol);

	/* Lecture des infos sur les champs */

	for(i=0;i<NmbTyp;i++)
		if(TypTab[i] == GmfSca)
			printf("typ %d = scalaire\n", i);
		else if(TypTab[i] == GmfVec)
			printf("typ %d = vecteur\n", i);
		else if(TypTab[i] == GmfSymMat)
			printf("typ %d = matrice symetrique\n", i);
		else if(TypTab[i] == GmfMat)
			printf("typ %d = matrice pleine\n", i);

	/* Lecture et ecriture du 2eme champ (vecteur) */

	GmfGotoKwd(InpSol, GmfSolAtVertices);

	NmbTyp = 1;
	TypTab[0] = GmfVec;

	puts("\necriture du 2eme champ (vecteur)\n");

	GmfSetKwd(OutSol, GmfSolAtVertices, NmbSol, NmbTyp, TypTab);

	for(i=1;i<=NmbSol;i++)
	{
		GmfGetLin(InpSol, GmfSolAtVertices, InpSolTab);

		OutSolTab[0] = InpSolTab[1];
		OutSolTab[1] = InpSolTab[2];
		OutSolTab[2] = InpSolTab[3];

		GmfSetLin(OutSol, GmfSolAtVertices, OutSolTab);
	}

	/* Fermeture des fichiers */

	GmfCloseMesh(InpSol);
	GmfCloseMesh(OutSol);

	return(0);
}
