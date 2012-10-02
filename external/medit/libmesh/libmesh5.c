

/*----------------------------------------------------------*/
/*															*/
/*						LIBMESH V 5.0						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		handle .meshb file format I/O		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 16 2007							*/
/*	Last modification:	jul 12 2007							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "libmesh5.h"


/*----------------------------------------------------------*/
/* Structures												*/
/*----------------------------------------------------------*/

typedef struct
{
	int pos, typ, SolSiz, NmbLin, NmbTyp, TypTab[ GmfMaxTyp ];
	char fmt[ GmfMaxTyp ];
}KwdSct;

typedef struct
{
	int dim, ver, iter, mod, typ, cod, NexKwdPos;
	double angle, bbox[3][2], time;
	KwdSct KwdTab[ GmfMaxKwd + 1 ];
	FILE *hdl;
	char FilNam[ GmfStrSiz ];
}GmfMshSct;


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define Asc 1
#define Bin 2
#define MshFil 4
#define SolFil 8
#define MaxMsh 100
#define InfKwd 1
#define RegKwd 2
#define SolKwd 3
#define WrdSiz 4


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

int IniFlg=0;
GmfMshSct *MshTab[ MaxMsh + 1 ];
char *KwdFmt[ GmfMaxKwd + 1 ][3] = 
{	{"Reserved", "", ""},
	{"MeshVersionFormatted", "", "i"},
	{"Reserved", "", ""},
	{"Dimension", "", "i"},
	{"Vertices", "i", "dri"},
	{"Edges", "i", "iii"},
	{"Triangles", "i", "iiii"},
	{"Quadrilaterals", "i", "iiiii"},
	{"Tetrahedra", "i", "iiiii"},
	{"Pentahedra", "i", "iiiiiii"},
	{"Hexahedra", "i", "iiiiiiiii"},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Corners", "i", "i"},
	{"Ridges", "i", "i"},
	{"RequiredVertices", "i", "i"},
	{"RequiredEdges", "i", "i"},
	{"RequiredTriangles", "i", "i"},
	{"RequiredQuadrilaterals", "i", "i"},
	{"TangentAtEdgeVertices", "i", "iii"},
	{"NormalAtVertices", "i", "ii"},
	{"NormalAtTriangleVertices", "i", "iii"},
	{"NormalAtQuadrilateralVertices", "i", "iiii"},
	{"AngleOfCornerBound", "", "r"},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"BoundingBox", "", "drdr"},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"End", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Reserved", "", ""},
	{"Tangents", "i", "dr"},
	{"Normals", "i", "dr"},
	{"TangentAtVertices", "i", "ii"},
	{"SolAtVertices", "i", "sr"},
	{"SolAtEdges", "i", "sr"},
	{"SolAtTriangles", "i", "sr"},
	{"SolAtQuadrilaterals", "i", "sr"},
	{"SolAtTetrahedra", "i", "sr"},
	{"SolAtPentahedra", "i", "sr"},
	{"SolAtHexahedra", "i", "sr"},
	{"DSolAtVertices", "i", "sr"},
	{"ISolAtVertices", "i", "i"},
	{"ISolAtEdges", "i", "ii"},
	{"ISolAtTriangles", "i", "iii"},
	{"ISolAtQuadrilaterals", "i", "iiii"},
	{"ISolAtTetrahedra", "i", "iiii"},
	{"ISolAtPentahedra", "i", "iiiiii"},
	{"ISolAtHexahedra", "i", "iiiiiiii"},
	{"Iterations","","i"},
	{"Time","","r"},
	{"Reserved","",""}
 };


/*----------------------------------------------------------*/
/* Prototypes of local procedures							*/
/*----------------------------------------------------------*/

static void ScaWrd(GmfMshSct *, unsigned char *);
static void ScaDblWrd(GmfMshSct *, unsigned char *);
static void RecWrd(GmfMshSct *, unsigned char *);
static void RecDblWrd(GmfMshSct *, unsigned char *);
static int ScaKwdTab(GmfMshSct *);
static void ExpFmt(GmfMshSct *, int);
static void ScaKwdHdr(GmfMshSct *, int);


/*----------------------------------------------------------*/
/* Open a mesh file in read or write mod					*/
/*----------------------------------------------------------*/

int GmfOpenMesh(const char *FilNam, int mod, ...)
{
	int i, KwdCod, res, *PtrVer, *PtrDim, MshIdx=0;
	char str[ GmfStrSiz ];
	va_list par;
	GmfMshSct *msh;

	if(!IniFlg)
	{
	  for(i=0;i<MaxMsh;i++)
	    MshTab[i] = NULL;

	  IniFlg = 1;
	}

	/*---------------------*/
	/* MESH STRUCTURE INIT */
	/*---------------------*/

	for(i=1;i<MaxMsh;i++)
	  if(!MshTab[i])
	    {
	      MshIdx = i;
	      break;
	    }
	
	if( !MshIdx || !(msh = (GmfMshSct *) calloc(1, sizeof(GmfMshSct))) )
	  return(0);

	/* Copy the FilNam into the structure */

	if(strlen(FilNam) + 7 >= GmfStrSiz)
	  return(0);

	strcpy(msh->FilNam, FilNam);

	/* Store the opening mod (read or write) and guess the filetype (binary or ascii) depending on the extension */

	msh->mod = mod;

	if(strstr(msh->FilNam, ".meshb"))
		msh->typ |= (Bin | MshFil);
	else if(strstr(msh->FilNam, ".mesh"))
		msh->typ |= (Asc | MshFil);
	else if(strstr(msh->FilNam, ".solb"))
		msh->typ |= (Bin | SolFil);
	else if(strstr(msh->FilNam, ".sol"))
		msh->typ |= (Asc | SolFil);
	else
		return(0);

	/* Open the file in the required mod and initialyse the mesh structure */

	if(msh->mod == GmfRead)
	{

		/*-----------------------*/
		/* OPEN FILE FOR READING */
		/*-----------------------*/

		va_start(par, mod);
		PtrVer = va_arg(par, int *);
		PtrDim = va_arg(par, int *);
		va_end(par);

		/* Create the name string and open the file */

		if(!(msh->hdl = fopen(msh->FilNam, "rb")))
			return(0);

		/* Read the endian coding tag, the mesh version and the mesh dimension (mandatory kwd) */

		if(msh->typ & Bin)
		{
			fread((unsigned char *)&msh->cod, WrdSiz, 1, msh->hdl);

			if( (msh->cod != 1) && (msh->cod != 16777216) )
				return(0);

			ScaWrd(msh, (unsigned char *)&msh->ver);
			ScaWrd(msh, (unsigned char *)&KwdCod);

			if(KwdCod != GmfDimension)
				return(0);

			ScaWrd(msh, (unsigned char *)&KwdCod);
			ScaWrd(msh, (unsigned char *)&msh->dim);
		}
		else
		{
			do
			{
				res = fscanf(msh->hdl, "%s", str);
			}while( (res != EOF) && strcmp(str, "MeshVersionFormatted") );

			if(res == EOF)
				return(0);

			fscanf(msh->hdl, "%d", &msh->ver);

			do
			{
				res = fscanf(msh->hdl, "%s", str);
			}while( (res != EOF) && strcmp(str, "Dimension") );

			if(res == EOF)
				return(0);

			fscanf(msh->hdl, "%d", &msh->dim);
		}

		if( (msh->dim != 2) && (msh->dim != 3) )
			return(0);

		(*PtrVer) = msh->ver;
		(*PtrDim) = msh->dim;

		/*------------*/
		/* KW READING */
		/*------------*/

		/* Read the list of kw present in the file */

		if(!ScaKwdTab(msh))
			return(0);

		MshTab[ MshIdx ] = msh;

		return(MshIdx);
	}
	else if(msh->mod == GmfWrite)
	{

		/*-----------------------*/
		/* OPEN FILE FOR WRITING */
		/*-----------------------*/

		msh->cod = 1;

		/* Check if the user provided a valid version number and dimension */

		va_start(par, mod);
		msh->ver = va_arg(par, int);
		msh->dim = va_arg(par, int);
		va_end(par);

		if( (msh->ver != 1) && (msh->ver != 2) )
			return(0);
		if( (msh->dim != 2) && (msh->dim != 3) )
			return(0);

		/* Create the mesh file */

		if(!(msh->hdl = fopen(msh->FilNam, "wb")))
			return(0);

		MshTab[ MshIdx ] = msh;


		/*------------*/
		/* KW WRITING */
		/*------------*/

		/* Write the mesh version and dimension */

		if(msh->typ & Asc)
		{
			fprintf(msh->hdl, "%s %d\n\n", KwdFmt[ GmfVersionFormatted ][0], msh->ver);
			fprintf(msh->hdl, "%s %d\n", KwdFmt[ GmfDimension ][0], msh->dim);
		}
		else
		{
			RecWrd(msh, (unsigned char *)&msh->cod);
			RecWrd(msh, (unsigned char *)&msh->ver);
			GmfSetKwd(MshIdx, GmfDimension, 0);
			RecWrd(msh, (unsigned char *)&msh->dim);
		}

		return(MshIdx);
	}
	else
		return(0);
}


/*----------------------------------------------------------*/
/* Close a meshfile in the right way						*/
/*----------------------------------------------------------*/

int GmfCloseMesh(int MshIdx)
{
	int res = 1;
	GmfMshSct *msh;

	if( (MshIdx < 1) || (MshIdx > MaxMsh) )
	  return(0);
	
	fflush(stdout);
	msh = MshTab[ MshIdx ];

	/* In write down the "End" kw in write mode */

	if(msh->mod == GmfWrite)
		if(msh->typ & Asc)
			fprintf(msh->hdl, "\n%s\n", KwdFmt[ GmfEnd ][0]);
		else
			GmfSetKwd(MshIdx, GmfEnd, 0);

	/* Close the file and free the mesh structure */

	if(fclose(msh->hdl))
		res = 0;

	free(msh);
	MshTab[ MshIdx ] = NULL;
	return(res);
}


/*----------------------------------------------------------*/
/* Read the number of lines and set the position to this kwd*/
/*----------------------------------------------------------*/

int GmfStatKwd(int MshIdx, int KwdCod, ...)
{
	int i, *PtrNmbTyp, *PtrSolSiz, *TypTab;
	GmfMshSct *msh;
	KwdSct *kwd;
	va_list par;

	if( (MshIdx < 1) || (MshIdx > MaxMsh) )
		return(0);

	msh = MshTab[ MshIdx ];

	if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
		return(0);

	kwd = &msh->KwdTab[ KwdCod ];

	if(!kwd->NmbLin)
		return(0);

	/* Read further arguments if this kw is a sol */

	if(kwd->typ == SolKwd)
	{
		va_start(par, KwdCod);

		PtrNmbTyp = va_arg(par, int *);
		*PtrNmbTyp = kwd->NmbTyp;

		PtrSolSiz = va_arg(par, int *);
		*PtrSolSiz = kwd->SolSiz;

		TypTab = va_arg(par, int *);

		for(i=0;i<kwd->NmbTyp;i++)
			TypTab[i] = kwd->TypTab[i];

		va_end(par);
	}

	return(kwd->NmbLin);
}


/*----------------------------------------------------------*/
/* Set the current file position to a given kwd				*/
/*----------------------------------------------------------*/

int GmfGotoKwd(int MshIdx, int KwdCod)
{
	GmfMshSct *msh;
	KwdSct *kwd;

	if( (MshIdx < 1) || (MshIdx > MaxMsh) )
		return(0);

	msh = MshTab[ MshIdx ];

	if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
		return(0);

	kwd = &msh->KwdTab[ KwdCod ];

	if(!kwd->NmbLin)
		return(0);

	return(fseek(msh->hdl, kwd->pos, SEEK_SET));
}


/*----------------------------------------------------------*/
/* Write the kwd and set the number of lines				*/
/*----------------------------------------------------------*/

int GmfSetKwd(int MshIdx, int KwdCod, ...)
{
	int i, CurPos, NmbLin=0, NulPos=0, *TypTab;
	va_list par;
	GmfMshSct *msh;
	KwdSct *kwd;

	if( (MshIdx < 1) || (MshIdx > MaxMsh) )
		return(0);

	msh = MshTab[ MshIdx ];

	if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
		return(0);

	kwd = &msh->KwdTab[ KwdCod ];

	/* Read further arguments if this kw has a header */

	if(strlen(KwdFmt[ KwdCod ][1]))
	{
		va_start(par, KwdCod);
		NmbLin = va_arg(par, int);

		if(!strcmp(KwdFmt[ KwdCod ][2], "sr"))
		{
			kwd->NmbTyp = va_arg(par, int);
			TypTab = va_arg(par, int *);

			for(i=0;i<kwd->NmbTyp;i++)
				kwd->TypTab[i] = TypTab[i];
		}

		va_end(par);
	}

	/* Setup the kwd info */

	ExpFmt(msh, KwdCod);

	if(!kwd->typ)
		return(0);
	else if(kwd->typ == InfKwd)
		kwd->NmbLin = 1;
	else
		kwd->NmbLin = NmbLin;

	/* Store the next kwd position in binary file */

	if( (msh->typ & Bin) && msh->NexKwdPos )
	{
		CurPos = ftell(msh->hdl);
		fseek(msh->hdl, msh->NexKwdPos, SEEK_SET);
		RecWrd(msh, (unsigned char *)&CurPos);
		fseek(msh->hdl, CurPos, SEEK_SET);
	}

	/* Write the header */

	if(msh->typ & Asc)
	{
		fprintf(msh->hdl, "\n%s\n", KwdFmt[ KwdCod ][0]);

		if(kwd->typ != InfKwd)
			fprintf(msh->hdl, "%d\n", kwd->NmbLin);

		/* In case of solution field, write the extended header */

		if(kwd->typ == SolKwd)
		{
			fprintf(msh->hdl, "%d ", kwd->NmbTyp);

			for(i=0;i<kwd->NmbTyp;i++)
				fprintf(msh->hdl, "%d ", kwd->TypTab[i]);

			fprintf(msh->hdl, "\n\n");
		}
	}
	else
	{
		RecWrd(msh, (unsigned char *)&KwdCod);
		msh->NexKwdPos = ftell(msh->hdl);
		RecWrd(msh, (unsigned char *)&NulPos);

		if(kwd->typ != InfKwd)
			RecWrd(msh, (unsigned char *)&kwd->NmbLin);

		/* In case of solution field, write the extended header at once */

		if(kwd->typ == SolKwd)
		{
			RecWrd(msh, (unsigned char *)&kwd->NmbTyp);

			for(i=0;i<kwd->NmbTyp;i++)
				RecWrd(msh, (unsigned char *)&kwd->TypTab[i]);
		}
	}

	return(kwd->NmbLin);
}


/*----------------------------------------------------------*/
/* Read a full line from the current kwd					*/
/*----------------------------------------------------------*/

void GmfGetLin(int MshIdx, int KwdCod, ...)
{
	double *DblPtr, *DblSolTab;
	float *FltPtr, *FltSolTab;
	int i, j, *IntPtr;
	va_list par;
	GmfMshSct *msh = MshTab[ MshIdx ];
	KwdSct *kwd = &msh->KwdTab[ KwdCod ];

	/* Start decoding the arguments */

	va_start(par, KwdCod);

	if(kwd->typ != SolKwd)
	{
		if(msh->ver == 1)
		{
			if(msh->typ & Asc)
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						FltPtr = va_arg(par, float *);
						fscanf(msh->hdl, "%f", FltPtr);
					}
					else
					{
						IntPtr = va_arg(par, int *);
						fscanf(msh->hdl, "%d", IntPtr);
					}
				}
			}
			else
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						FltPtr = va_arg(par, float *);
						ScaWrd(msh, (unsigned char *)FltPtr);
					}
					else
					{
						IntPtr = va_arg(par, int *);
						ScaWrd(msh, (unsigned char *)IntPtr);
					}
				}
			}
		}
		else
		{
			if(msh->typ & Asc)
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						DblPtr = va_arg(par, double *);
						fscanf(msh->hdl, "%lf", DblPtr);
					}
					else
					{
						IntPtr = va_arg(par, int *);
						fscanf(msh->hdl, "%d", IntPtr);
					}
				}
			}
			else
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						DblPtr = va_arg(par, double *);
						ScaDblWrd(msh, (unsigned char *)DblPtr);
					}
					else
					{
						IntPtr = va_arg(par, int *);
						ScaWrd(msh, (unsigned char *)IntPtr);
					}
				}
			}
		}
	}
	else
	{
		if(msh->ver == 1)
		{
			FltSolTab = va_arg(par, float *);

			if(msh->typ & Asc)
				for(j=0;j<kwd->SolSiz;j++)
					fscanf(msh->hdl, "%f", &FltSolTab[j]);
			else
				for(j=0;j<kwd->SolSiz;j++)
					ScaWrd(msh, (unsigned char *)&FltSolTab[j]);
		}
		else if(msh->ver == 2)
		{
			DblSolTab = va_arg(par, double *);

			if(msh->typ & Asc)
				for(j=0;j<kwd->SolSiz;j++)
					fscanf(msh->hdl, "%lf", &DblSolTab[j]);
			else
				for(j=0;j<kwd->SolSiz;j++)
					ScaDblWrd(msh, (unsigned char *)&DblSolTab[j]);
		}
	}

	va_end(par);
}


/*----------------------------------------------------------*/
/* Write a full line from the current kwd					*/
/*----------------------------------------------------------*/

void GmfSetLin(int MshIdx, int KwdCod, ...)
{
	double d, *DblSolTab;
	float f, *FltSolTab;
	int i, j;
	va_list par;
	GmfMshSct *msh = MshTab[ MshIdx ];
	KwdSct *kwd = &msh->KwdTab[ KwdCod ];

	/* Start decoding the arguments */

	va_start(par, KwdCod);

	if(kwd->typ != SolKwd)
	{
		if(msh->ver == 1)
		{
			if(msh->typ & Asc)
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						d = va_arg(par, double);
						fprintf(msh->hdl, "%g ", (float)d);
					}
					else
					{
						j = va_arg(par, int);
						fprintf(msh->hdl, "%d ", j);
					}
				}
			}
			else
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						d = va_arg(par, double);
						f = d;
						RecWrd(msh, (unsigned char *)&f);
					}
					else
					{
						j = va_arg(par, int);
						RecWrd(msh, (unsigned char *)&j);
					}
				}
			}
		}
		else
		{
			if(msh->typ & Asc)
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						d = va_arg(par, double);
						fprintf(msh->hdl, "%.15lg ", d);
					}
					else
					{
						j = va_arg(par, int);
						fprintf(msh->hdl, "%d ", j);
					}
				}
			}
			else
			{
				for(i=0;i<kwd->SolSiz;i++)
				{
					if(kwd->fmt[i] == 'r')
					{
						d = va_arg(par, double);
						RecDblWrd(msh, (unsigned char *)&d);
					}
					else
					{
						j = va_arg(par, int);
						RecWrd(msh, (unsigned char *)&j);
					}
				}
			}
		}
	}
	else
	{
		if(msh->ver == 1)
		{
			FltSolTab = va_arg(par, float *);

			if(msh->typ & Asc)
				for(j=0;j<kwd->SolSiz;j++)
					fprintf(msh->hdl, "%g ", FltSolTab[j]);
			else
				for(j=0;j<kwd->SolSiz;j++)
					RecWrd(msh, (unsigned char *)&FltSolTab[j]);
		}
		else if(msh->ver == 2)
		{
			DblSolTab = va_arg(par, double *);

			if(msh->typ & Asc)
				for(j=0;j<kwd->SolSiz;j++)
					fprintf(msh->hdl, "%.15lg ", DblSolTab[j]);
			else
				for(j=0;j<kwd->SolSiz;j++)
					RecDblWrd(msh, (unsigned char *)&DblSolTab[j]);
		}
	}

	va_end(par);

	if(msh->typ & Asc)
		fprintf(msh->hdl, "\n");
}


/*----------------------------------------------------------*/
/* Private procedure for transmesh : copy a whole line		*/
/*----------------------------------------------------------*/

void GmfCpyLin(int InpIdx, int OutIdx, int KwdCod)
{
	double d;
	float f;
	int i, a;
	GmfMshSct *InpMsh = MshTab[ InpIdx ], *OutMsh = MshTab[ OutIdx ];
	KwdSct *kwd = &InpMsh->KwdTab[ KwdCod ];

	for(i=0;i<kwd->SolSiz;i++)
	{
		if(kwd->fmt[i] == 'r')
		{
			if(InpMsh->ver == 1)
			{
				if(InpMsh->typ & Asc)
					fscanf(InpMsh->hdl, "%f", &f);
				else
					ScaWrd(InpMsh, (unsigned char *)&f);

				d = f;
			}
			else
			{
				if(InpMsh->typ & Asc)
					fscanf(InpMsh->hdl, "%lf", &d);
				else
					ScaDblWrd(InpMsh, (unsigned char *)&d);

				f = (float)d;
			}

			if(OutMsh->ver == 1)
				if(OutMsh->typ & Asc)
					fprintf(OutMsh->hdl, "%g ", f);
				else
					RecWrd(OutMsh, (unsigned char *)&f);
			else
				if(OutMsh->typ & Asc)
					fprintf(OutMsh->hdl, "%.15g ", d);
				else
					RecDblWrd(OutMsh, (unsigned char *)&d);
		}
		else
		{
			if(InpMsh->typ & Asc)
				fscanf(InpMsh->hdl, "%d", &a);
			else
				ScaWrd(InpMsh, (unsigned char *)&a);

			if(OutMsh->typ & Asc)
				fprintf(OutMsh->hdl, "%d ", a);
			else
				RecWrd(OutMsh, (unsigned char *)&a);
		}
	}

	if(OutMsh->typ & Asc)
		fprintf(OutMsh->hdl, "\n");
}


/*----------------------------------------------------------*/
/* Find every kw present in a meshfile						*/
/*----------------------------------------------------------*/

static int ScaKwdTab(GmfMshSct *msh)
{
	int KwdCod, NexPos, CurPos, EndPos;
	char str[ GmfStrSiz ];

	if(msh->typ & Asc)
	{
		/* Scan each string in the file until the end */

		while(fscanf(msh->hdl, "%s", str) != EOF)
		{
			/* Fast test in order to reject quickly the numeric values */

			if(isalpha(str[0]))
			{
				/* Search which kwd code this string is associated with, 
					then get its header and save the curent position in file (just before the data) */

				for(KwdCod=1; KwdCod<= GmfMaxKwd; KwdCod++)
					if(!strcmp(str, KwdFmt[ KwdCod ][0]))
					{
						ScaKwdHdr(msh, KwdCod);
						break;
					}
			}
			else if(str[0] == '#')
				while(fgetc(msh->hdl) != '\n');
		}
	}
	else
	{
		/* Get file size */

		CurPos = ftell(msh->hdl);
		fseek(msh->hdl, 0, SEEK_END);
		EndPos = ftell(msh->hdl);
		fseek(msh->hdl, CurPos, SEEK_SET);

		/* Jump through kwd positions in the file */

		do
		{
			/* Get the kwd code and the next kwd position */

			ScaWrd(msh, (unsigned char *)&KwdCod);
			ScaWrd(msh, (unsigned char *)&NexPos);

			if(NexPos > EndPos)
				return(0);

			/* Check if this kwd belongs to this mesh version */

			if( (KwdCod >= 1) && (KwdCod <= GmfMaxKwd) )
				ScaKwdHdr(msh, KwdCod);

			/* Go to the next kwd */

			if(NexPos)
				fseek(msh->hdl, NexPos, SEEK_SET);
		}while(NexPos && (KwdCod != GmfEnd));
	}

	return(1);
}


/*----------------------------------------------------------*/
/* Read and setup the keyword's header						*/
/*----------------------------------------------------------*/

static void ScaKwdHdr(GmfMshSct *msh, int KwdCod)
{
	int i;
	KwdSct *kwd = &msh->KwdTab[ KwdCod ];

	if(!strcmp("i", KwdFmt[ KwdCod ][1]))
	{
		if(msh->typ & Asc)
			fscanf(msh->hdl, "%d", &kwd->NmbLin);
		else
			ScaWrd(msh, (unsigned char *)&kwd->NmbLin);
	}
	else
		kwd->NmbLin = 1;

	if(!strcmp("sr", KwdFmt[ KwdCod ][2]))
	{
		if(msh->typ & Asc)
		{
			fscanf(msh->hdl, "%d", &kwd->NmbTyp);

			for(i=0;i<kwd->NmbTyp;i++)
				fscanf(msh->hdl, "%d", &kwd->TypTab[i]);
		}
		else
		{
			ScaWrd(msh, (unsigned char *)&kwd->NmbTyp);

			for(i=0;i<kwd->NmbTyp;i++)
				ScaWrd(msh, (unsigned char *)&kwd->TypTab[i]);
		}
	}

	ExpFmt(msh, KwdCod);
	kwd->pos = ftell(msh->hdl);
}


/*----------------------------------------------------------*/
/* Expand the compacted format and compute the line size	*/
/*----------------------------------------------------------*/

static void ExpFmt(GmfMshSct *msh, int KwdCod)
{
	int i, j, TmpSiz=0;
	char chr, *InpFmt = KwdFmt[ KwdCod ][2];
	KwdSct *kwd = &msh->KwdTab[ KwdCod ];

	/* Set the kwd's type */

	if(!strlen(KwdFmt[ KwdCod ][1]))
		kwd->typ = InfKwd;
	else if(!strcmp(InpFmt, "sr"))
		kwd->typ = SolKwd;
	else
		kwd->typ = RegKwd;

	/* Get the solution-field's size */

	if(kwd->typ == SolKwd)
		for(i=0;i<kwd->NmbTyp;i++)
			switch(kwd->TypTab[i])
			{
				case GmfSca    : TmpSiz += 1; break;
				case GmfVec    : TmpSiz += msh->dim; break;
				case GmfSymMat : TmpSiz += (msh->dim * (msh->dim+1)) / 2; break;
				case GmfMat    : TmpSiz += msh->dim * msh->dim; break;
			}

	/* Scan each character from the format string */

	i = 0;

	while(i < (int) strlen(InpFmt))
	{
		chr = InpFmt[ i++ ];

		if(chr == 'd')
		{
			chr = InpFmt[i++];

			for(j=0;j<msh->dim;j++)
				kwd->fmt[ kwd->SolSiz++ ] = chr;
		}
		else if(chr == 's')
		{
			chr = InpFmt[i++];

			for(j=0;j<TmpSiz;j++)
				kwd->fmt[ kwd->SolSiz++ ] = chr;
		}
		else
			kwd->fmt[ kwd->SolSiz++ ] = chr;
	}
}


/*----------------------------------------------------------*/
/* Read a four bytes word in a mesh file					*/
/*----------------------------------------------------------*/

static void ScaWrd(GmfMshSct *msh, unsigned char *wrd)
{
	unsigned char swp;

	fread(wrd, WrdSiz, 1, msh->hdl);

	if(msh->cod == 1)
		return;

	swp = wrd[3];
	wrd[3] = wrd[0];
	wrd[0] = swp;

	swp = wrd[2];
	wrd[2] = wrd[1];
	wrd[1] = swp;
}


/*----------------------------------------------------------*/
/* Read an eight bytes word in a mesh file					*/
/*----------------------------------------------------------*/

static void ScaDblWrd(GmfMshSct *msh, unsigned char *wrd)
{
	int i;
	unsigned char swp;

	fread(wrd, WrdSiz, 2, msh->hdl);

	if(msh->cod == 1)
		return;

	for(i=0;i<4;i++)
	{
		swp = wrd[7-i];
		wrd[7-i] = wrd[i];
		wrd[i] = swp;
	}
}


/*----------------------------------------------------------*/
/* Write a four bytes word in a mesh file					*/
/*----------------------------------------------------------*/

static void RecWrd(GmfMshSct *msh, unsigned char *wrd)
{
	fwrite(wrd, WrdSiz, 1, msh->hdl);
}


/*----------------------------------------------------------*/
/* Write an eight bytes word in a mesh file					*/
/*----------------------------------------------------------*/

static void RecDblWrd(GmfMshSct *msh, unsigned char *wrd)
{
	fwrite(wrd, WrdSiz, 2, msh->hdl);
}


/*----------------------------------------------------------*/
/* Fortran 77 API											*/
/*----------------------------------------------------------*/

int call(gmfopenmeshf77)(char *FilNam, int *mod, int *ver, int *dim, int StrSiz)
{
	int i;
	char TmpNam[ GmfStrSiz ];

	for(i=0;i<StrSiz;i++)
		TmpNam[i] = FilNam[i];

	TmpNam[ StrSiz ] = 0;

	if(*mod == GmfRead)
		return(GmfOpenMesh(TmpNam, *mod, ver, dim));
	else
		return(GmfOpenMesh(TmpNam, *mod, *ver, *dim));
}

int call(gmfclosemeshf77)(int *idx)
{
	return(GmfCloseMesh(*idx));
}

int call(gmfstatkwdf77)(int *MshIdx, int *KwdIdx, int *NmbTyp, int *SolSiz, int *TypTab)
{
	if(!strcmp(KwdFmt[ *KwdIdx ][2], "sr"))
		return(GmfStatKwd(*MshIdx, *KwdIdx, NmbTyp, SolSiz, TypTab));
	else
		return(GmfStatKwd(*MshIdx, *KwdIdx));
}

int call(gmfgotokwdf77)(int *MshIdx, int *KwdIdx)
{
	return(GmfGotoKwd(*MshIdx, *KwdIdx));
}

int call(gmfsetkwdf77)(int *MshIdx, int *KwdIdx, int *NmbLin, int *NmbTyp, int *TypTab)
{
	if(!strcmp(KwdFmt[ *KwdIdx ][2], "sr"))
		return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin, *NmbTyp, TypTab));
	else if(strlen(KwdFmt[ *KwdIdx ][1]))
		return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin));
	else
		return(GmfSetKwd(*MshIdx, *KwdIdx));
}

int call(gmfgetvertex2df77)(int *MshIdx, float *x, float *y, int *ref)
{
	GmfGetLin(*MshIdx, GmfVertices, x, y, ref);
	return(1);
}

int call(gmfgetvertex3df77)(int *MshIdx, float *x, float *y, float *z, int *ref)
{
	GmfGetLin(*MshIdx, GmfVertices, x, y, z, ref);
	return(1);
}

int call(gmfsetvertex2df77)(int *MshIdx, float *x, float *y, int *ref)
{
	GmfSetLin(*MshIdx, GmfVertices, *x, *y, *ref);
	return(1);
}

int call(gmfsetvertex3df77)(int *MshIdx, float *x, float *y, float *z, int *ref)
{
	GmfSetLin(*MshIdx, GmfVertices, *x, *y, *z, *ref);
	return(1);
}

int call(gmfgettrianglef77)(int *MshIdx, int *p1, int *p2, int *p3, int *ref)
{
	GmfGetLin(*MshIdx, GmfTriangles, p1, p2, p3, ref);
	return(1);
}

int call(gmfsettrianglef77)(int *MshIdx, int *p1, int *p2, int *p3, int *ref)
{
	GmfSetLin(*MshIdx, GmfTriangles, *p1, *p2, *p3, *ref);
	return(1);
}

int call(gmfgettetrahedronf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *ref)
{
	GmfGetLin(*MshIdx, GmfTetrahedra, p1, p2, p3, p4, ref);
	return(1);
}

int call(gmfsettetrahedronf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *ref)
{
	GmfSetLin(*MshIdx, GmfTetrahedra, *p1, *p2, *p3, *p4, *ref);
	return(1);
}

int call(gmfgetedgef77)(int *MshIdx, int *p1, int *p2, int *ref)
{
	GmfGetLin(*MshIdx, GmfEdges, p1, p2, ref);
	return(1);
}

int call(gmfsetedgef77)(int *MshIdx, int *p1, int *p2, int *ref)
{
	GmfSetLin(*MshIdx, GmfEdges, *p1, *p2, *ref);
	return(1);
}

int call(gmfgetquadrilateralf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *ref)
{
	GmfGetLin(*MshIdx, GmfQuadrilaterals, p1, p2, p3, p4, ref);
	return(1);
}

int call(gmfsetquadrilateralf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *ref)
{
	GmfSetLin(*MshIdx, GmfQuadrilaterals, *p1, *p2, *p3, *p4, *ref);
	return(1);
}

int call(gmfgethexahedronf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *p5, int *p6, int *p7, int *p8, int *ref)
{
	GmfGetLin(*MshIdx, GmfHexahedra, p1, p2, p3, p4, p5, p6, p7, p8, ref);
	return(1);
}

int call(gmfsethexahedronf77)(int *MshIdx, int *p1, int *p2, int *p3, int *p4, int *p5, int *p6, int *p7, int *p8, int *ref)
{
	GmfSetLin(*MshIdx, GmfHexahedra, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *ref);
	return(1);
}

int call(gmfsetsolf77)(int *MshIdx, int *kwd, int *SolTab)
{
	GmfSetLin(*MshIdx, *kwd, SolTab);
	return(1);
}

int call(gmfgetsolf77)(int *MshIdx, int *kwd, int *SolTab)
{
	GmfGetLin(*MshIdx, *kwd, SolTab);
	return(1);
}
