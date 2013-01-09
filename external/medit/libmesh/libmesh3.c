

/*----------------------------------------------------------*/
/*															*/
/*						LIBMESH V 3.0						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		handle .meshb file format I/O		*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		aug  2 2003							*/
/*	Last modification:	jan 25 2006							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "libmesh3.h"


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

char *LM_kw_table[ LM_NBKW + 1 ][3] = 
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
	{"SubDomainFromGeom", "i", "iiii"},
	{"SubDomainFromMesh", "i", "iiii"},
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
	{"Geometry", "", "c"},
	{"VertexOnGeometricVertex", "i", "ii"},
	{"VertexOnGeometricEdge", "i", "iir"},
	{"VertexOnGeometricTriangle", "i", "iirr"},
	{"VertexOnGeometricQuadrilateral", "i", "iirr"},
	{"EdgeOnGeometricEdge", "i", "ii"},
	{"TriangleOnGeometricTriangle", "i", "ii"},
	{"TriangleOnGeometricQuadrilateral", "i", "ii"},
	{"QuadrilateralOnGeometricTriangle", "i", "ii"},
	{"QuadrilateralOnGeometricQuadrilateral", "i", "ii"},
	{"MeshSupportOfVertices", "", "c"},
	{"VertexOnSupportVertex", "i", "ii"},
	{"VertexOnSupportEdge", "i", "iir"},
	{"VertexOnSupportTriangle", "i", "iirr"},
	{"VertexOnSupportQuadrilateral", "i", "iirr"},
	{"VertexOnSupportTetrahedron", "i", "iirrr"},
	{"VertexOnSupportPentahedron", "i", "iirrr"},
	{"VertexOnSupportHexahedron", "i", "iirrr"},
	{"CrackedEdges", "i", "ii"},
	{"CrackedTriangles", "i", "ii"},
	{"CrackedQuadrilaterals", "i", "ii"},
	{"EquivalentEdges", "i", "ii"},
	{"EquivalentTriangles", "i", "ii"},
	{"EquivalentQuadrilaterals", "i", "ii"},
	{"PhysicsReference", "i", "ic"},
	{"IncludeFile", "", "c"},
	{"BoundingBox", "", "drdr"},
	{"Identifier", "", "c"},
	{"IdentityOfGeometry", "", "c"},
	{"IdentityOfMeshSupport", "", "c"},
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
	{"VertexHack","","drdr"}
 };


/*----------------------------------------------------------*/
/* Prototypes of local procedures							*/
/*----------------------------------------------------------*/

static void write_kw(LM_mesh_struct *, int);
static int read_int(LM_mesh_struct *);
static void write_int(LM_mesh_struct *, int);
static void file2kw_tab(LM_mesh_struct *);
static void kw_tab2file(LM_mesh_struct *);
static int expand_format(LM_mesh_struct *, int, char *);
static void swap_bytes(void *, void *, int);
static void read_sol_headers(LM_mesh_struct *);


/*----------------------------------------------------------*/
/* Open a mesh file in read or write mode					*/
/*----------------------------------------------------------*/

int LM_open_mesh(const char *filename, int mode, LM_mesh_struct *mesh, ...)
{
	int i;
	va_list pa;
	

	/*---------------------*/
	/* MESH STRUCTURE INIT */
	/*---------------------*/

	/* Init the kw table */

	for(i=0;i<=LM_NBKW;i++)
	{
		mesh->kw_counters[i] = 0;
		mesh->kw_pos[i][0] = 0;
		mesh->kw_pos[i][1] = 0;
		mesh->kw_pos[i][2] = 0;
		mesh->sol_headers[i] = NULL;
	}

	/* Allocate a string large enough to contain the full filename and path plus the ".meshb" extension. */

	mesh->filename = (char *)calloc((strlen(filename) + 7), sizeof(char));
	strcpy(mesh->filename, filename);

	/* Store the opening mode (read or write) and guess the filetype (binary or ascii) depending on the extension */

	mesh->mode = mode;
	mesh->current_kw = mesh->type = 0;
	mesh->endian = 1;

	if(strstr(mesh->filename, ".meshb"))
		mesh->type |= (LM_BINARY | LM_MESH);
	else if(strstr(mesh->filename, ".mesh"))
		mesh->type |= (LM_ASCII | LM_MESH);
	else if(strstr(mesh->filename, ".solb"))
		mesh->type |= (LM_BINARY | LM_SOL);
	else if(strstr(mesh->filename, ".sol"))
		mesh->type |= (LM_ASCII | LM_SOL);
	else
		return(0);

	/* Open the file in the required mode and initialyse the mesh structure */

	if(mesh->mode == LM_READ)
	{

		/*-----------------------*/
		/* OPEN FILE FOR READING */
		/*-----------------------*/

		/* Create the name string and open the file */

		if(!(mesh->handle = fopen(mesh->filename, "rb")))
			return(0);

		/* Read the endian tag and the mesh version in binary */

		if(mesh->type & LM_BINARY)
		{
			mesh->endian = read_int(mesh);

			if( (mesh->endian != 1) && (mesh->endian != 16777216) )
				return(0);

			mesh->version = read_int(mesh);
		}


		/*------------*/
		/* KW READING */
		/*------------*/

		/* Read the list of kw present in the file */

		file2kw_tab(mesh);

		/* Check the mesh dimension */

		if(!mesh->kw_counters[ LM_Dimension ])
			return(0);

		LM_read_field(mesh, LM_Dimension, 1, &mesh->dimension);

		if( (mesh->dimension != 2) && (mesh->dimension != 3) )
			return(0);

		/* Read the meshversion in ascii case */

		if(mesh->type & LM_ASCII)
			LM_read_field(mesh, LM_MeshVersionFormatted, 1, &mesh->version);

		/* And read the extended sol headers */

		read_sol_headers(mesh);
	}
	else if(mesh->mode == LM_WRITE)
	{

		/*-----------------------*/
		/* OPEN FILE FOR WRITING */
		/*-----------------------*/

		/* Check if the user provided a valid dimension */

		va_start(pa, mesh);
		mesh->dimension = va_arg(pa, int);
		va_end(pa);

		if( (mesh->dimension != 2) && (mesh->dimension != 3) )
			return(0);

		/* If no extension has been provided, create a binary file */

		if(!(mesh->handle = fopen(mesh->filename, "wb")))
			return(0);


		/*------------*/
		/* KW WRITING */
		/*------------*/

		/* Initialyse the required fields. The kw will be stored afterward. */

		mesh->version = LM_MESH_VERSION;
		mesh->endian = 1;

		/* Write the mesh version */

		if(mesh->type & LM_ASCII)
			LM_write_field(mesh, LM_MeshVersionFormatted, 1, &mesh->version);
		else
		{
			write_int(mesh, mesh->endian);
			write_int(mesh, mesh->version);
		}

		/* Write the mesh dimension */

		LM_write_field(mesh, LM_Dimension, 1, &mesh->dimension);
	}
	else
		return(0);

	return(1);
}


/*----------------------------------------------------------*/
/* Close a meshfile in the right way						*/
/*----------------------------------------------------------*/

int LM_close_mesh(LM_mesh_struct *mesh)
{
	if(mesh->mode == LM_WRITE)
	{
		/* Test if the user wrote the "End" kw */

		if(!mesh->kw_counters[ LM_End ])
			LM_write_field(mesh, LM_End, 0, NULL);

		/* Write down the number lines written in each field to the file */

		kw_tab2file(mesh);
	}

	if(fclose(mesh->handle))
		return(0);
	else
		return(1);
}


/*----------------------------------------------------------*/
/* Bufferized read of a whole orpart of a field				*/
/*----------------------------------------------------------*/

int LM_read_field(LM_mesh_struct *mesh, int kw_code, int nbl, void *buffer)
{
	int i, j, swaped, size, *int_buffer = (int *)buffer, string_size;
	float *flt_buffer = (float *)buffer;
	char format[256], letter, str_buf[256];

	/* Check if the kw code is valid */

	if( (kw_code < 1) || (kw_code > LM_NBKW) )
		return(0);

	/* Check if this kw has a format */

	if(!strlen(LM_kw_table[ kw_code ][2]))
		return(0);

	/* If this kw is only a header, the number of lines to be read is set to one */

	if(!strlen(LM_kw_table[ kw_code ][1]))
		nbl = 1;

	/* Check if the user is not asking more lines than the remaining lines in the file */

	if(nbl > mesh->kw_counters[ kw_code ] - mesh->kw_pos[ kw_code ][2])
		nbl = mesh->kw_counters[ kw_code ] - mesh->kw_pos[ kw_code ][2];

	if(!nbl)
		return(0);

	/* Set the curent position in file to the begining of this kw's data */

	fseek(mesh->handle, mesh->kw_pos[ kw_code ][1], SEEK_SET);

	/* Transform the internal format into a "c" format string for the scanf and compute the size of 
		field's line in order to compute the right adresses in the buffer */

	size = expand_format(mesh, kw_code, format);

	if(mesh->type & LM_ASCII)
	{
		for(i=0;i<nbl;i++)
			for(j=0;j<size;j++)
				if(format[j] == 'i')
					fscanf(mesh->handle, "%d", &int_buffer[ i * size + j ]);
				else if(format[j] == 'r')
					fscanf(mesh->handle, "%g", &flt_buffer[ i * size + j ]);
				else if(format[j] == 'c')
				{
					string_size = 0;

					do
					{
						fscanf(mesh->handle, "%c", &letter);
					}while(letter != '"');

					do
					{
						fscanf(mesh->handle, "%c", &letter);
						str_buf[ string_size++ ] = letter;
					}while( (letter != '"') && (string_size <= 256) );

					str_buf[ string_size-1 ] = 0;

					memset(&flt_buffer[ i * size + j ], 0, 256);
					strcpy((char *)&flt_buffer[ i * size + j ], str_buf);
				}
	}
	else
	{
		fread(buffer, nbl * size * 4, 1, mesh->handle);

		/* Swap the bytes in the whole buffer in case of different endian */

		if(mesh->endian != 1)
			for(i=0;i<nbl*size;i++)
			{
				swap_bytes((void *)&int_buffer[i], (void *)&swaped, 4);
				int_buffer[i] = swaped;
			}
	}

	/* Then store the curent position and the total number of read lines in case we didn't read the whole data */

	mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);
	mesh->kw_pos[ kw_code ][2] += nbl;

	return(nbl);
}


/*----------------------------------------------------------*/
/* Bufferized write of a whole field or part of it			*/
/*----------------------------------------------------------*/

int LM_write_field(LM_mesh_struct *mesh, int kw_code, int nbl, void *buffer, ...)
{
	int i, j, size, *int_buffer = (int *)buffer, nbsol;
	float *flt_buffer = (float *)buffer;
	char format[256];
	va_list pa;

	/* Check if the kw code is valid */

	if( (kw_code < 1) || (kw_code > LM_NBKW) )
		return(0);

	/* Read further arguments if this kw is solution field and the extra header was not provided by the user */

	if(!mesh->sol_headers[ kw_code ] && !strcmp(LM_kw_table[ kw_code ][2], "sr"))
	{
		va_start(pa, buffer);

		nbsol = va_arg(pa, int);

#ifdef IGL
		if(!(mesh->sol_headers[ kw_code ] = malloc((nbsol+2) * sizeof(int))))
#else
		if(!(mesh->sol_headers[ kw_code ] = malloc((nbsol+2) * sizeof(int))))
#endif
			return(0);

		mesh->sol_headers[ kw_code ][0] = nbsol;
		mesh->sol_headers[ kw_code ][1] = 0;

		for(i=1;i<=nbsol;i++)
		{
			mesh->sol_headers[ kw_code ][i+1] = va_arg(pa, int);

			switch(mesh->sol_headers[ kw_code ][i+1])
			{
				case 1 : mesh->sol_headers[ kw_code ][1] += 1; break;
				case 2 : mesh->sol_headers[ kw_code ][1] += mesh->dimension; break;
				case 3 : mesh->sol_headers[ kw_code ][1] += (mesh->dimension * (mesh->dimension+1)) / 2; break;
				case 4 : mesh->sol_headers[ kw_code ][1] += mesh->dimension * mesh->dimension; break;
			}
		}

		va_end(pa);
	}

	/* If this kw is only a header, the number of lines to be read is set to one */

	if(!strlen(LM_kw_table[ kw_code ][1]))
		nbl = 1;

	if(!mesh->kw_counters[ kw_code ])
		write_kw(mesh, kw_code);

	mesh->kw_counters[ kw_code ] += nbl;

	/* Check if this kw has a format */

	if(!strlen(LM_kw_table[ kw_code ][2]))
		return(0);

	size = expand_format(mesh, kw_code, format);

	if(mesh->type & LM_ASCII)
	{
		for(i=0;i<nbl;i++)
		{
			for(j=0;j<size;j++)
				if(format[j] == 'i')
					fprintf(mesh->handle, "%d ", int_buffer[ i * size + j ]);
				else if(format[j] == 'r')
					fprintf(mesh->handle, "%g ", flt_buffer[ i * size + j ]);
				else if(format[j] == 'c')
				{
					fputc('"', mesh->handle);
					fprintf(mesh->handle, "%s", (char *)&flt_buffer[ i * size + j ]);
					fputc('"', mesh->handle);
				}

			fprintf(mesh->handle, "\n");
		}
	}
	else
		fwrite(buffer, nbl * size * 4, 1, mesh->handle);

	mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);

	return(1);
}


/*----------------------------------------------------------*/
/* Single line read											*/
/*----------------------------------------------------------*/

int LM_read_line(LM_mesh_struct *mesh, int kw_code, ...)
{
	float buffer[10], *ptr_flt;
	int i, size;
	char format[256];
	va_list pa;

	/* Check wether this kw is not a simple header */

	if(!strlen(LM_kw_table[ kw_code ][2]))
		return(0);

	/* Get a one line buffer */

	LM_read_field(mesh, kw_code, 1, buffer);

	/* Start decoding the arguments */

	va_start(pa, kw_code);
	size = expand_format(mesh, kw_code, format);

	for(i=0;i<size;i++)
	{
		ptr_flt = va_arg(pa, float *);
		*ptr_flt = buffer[i];
	}

	va_end(pa);

	/* return the number of arguments filled */

	return(size);
}


/*----------------------------------------------------------*/
/* Single line write										*/
/*----------------------------------------------------------*/

int LM_write_line(LM_mesh_struct *mesh, int kw_code, ...)
{
	float buffer[10], *ptr_flt;
	int i, size;
	char format[256];
	va_list pa;

	/* Check wether this kw is not a simple header */

	if(!strlen(LM_kw_table[ kw_code ][2]))
		return(0);

	/* Start decoding the arguments */

	va_start(pa, kw_code);
	size = expand_format(mesh, kw_code, format);

	for(i=0;i<size;i++)
	{
		ptr_flt = va_arg(pa, float *);
		buffer[i] = *ptr_flt;
	}

	va_end(pa);

	/* Write a one line buffer */

	LM_write_field(mesh, kw_code, 1, buffer);

	/* return the number of arguments filled */

	return(size);
}


/*----------------------------------------------------------*/
/* Find every kw present in a meshfile						*/
/*----------------------------------------------------------*/

static void file2kw_tab(LM_mesh_struct *mesh)
{
	int kw_code, next_pos;
	char str[256];

	if(mesh->type & LM_ASCII)
	{
		/* Scan each string of the file until the end */

		while(fscanf(mesh->handle, "%s", str) != EOF)
		{
			/* Fast test in order to reject quickly the numeric values */

			if(isalpha(str[0]))
			{
				/* Search which kw code this string is associated with */

				for(kw_code=1; kw_code<= LM_NBKW; kw_code++)
					if(!strcmp(str, LM_kw_table[ kw_code ][0]))
					{
						/* The base position (0) in the file is set right after the kw */

						mesh->kw_pos[ kw_code ][0] = ftell(mesh->handle);

						/* If this kw has a header, read the number of lines */

						if(!strcmp(LM_kw_table[ kw_code ][1], "i"))
							mesh->kw_counters[ kw_code ] = read_int(mesh);
						else
							mesh->kw_counters[ kw_code ] = 1;

						/* The curent position (1) in the file is set right before the data */

						mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);
						break;
					}
			}
			else if(str[0] == '#')
				while(fgetc(mesh->handle) != '\n');
		}
	}
	else
	{
		/* Jump through kw positions in the file */

		do
		{
			kw_code = read_int(mesh);

			/* Check if this kw belongs to this mesh version */

			if( (kw_code >= 1) && (kw_code <= LM_NBKW) )
			{
				/* The base position (0) in the file is set right after the kw */

				mesh->kw_pos[ kw_code ][0] =  ftell(mesh->handle);

				/* Read the next kw position */

				next_pos =  read_int(mesh);

				/* If this kw has a header, read the number of lines */

				if(!strcmp(LM_kw_table[ kw_code ][1], "i"))
					mesh->kw_counters[ kw_code ] = read_int(mesh);
				else
					mesh->kw_counters[ kw_code ] = 1;

				/* The curent position (1) in the file is set right before the data */

				mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);
			}
			else
			{
				/* Otherwise, read the next kw position in order to skip these unknown kw */

				next_pos =  read_int(mesh);
			}

			/* Go to the next kw */

			if(next_pos)
				fseek(mesh->handle, next_pos, SEEK_SET);
		}while(next_pos);
	}
}


/*----------------------------------------------------------*/
/* Update the number of lines written for each kw			*/
/*----------------------------------------------------------*/

static void kw_tab2file(LM_mesh_struct *mesh)
{
	int i;

	for(i=1;i<=LM_NBKW;i++)
		if( mesh->kw_counters[i] && strlen(LM_kw_table[i][2]) )
			write_kw(mesh, i);
}


/*----------------------------------------------------------*/
/* Write the string associated with the kw code				*/
/*----------------------------------------------------------*/

static void write_kw(LM_mesh_struct *mesh, int kw_code)
{
	int i;

	if(mesh->type & LM_ASCII)
	{
		/* Test if it is the first time this kw is written */

		if(!mesh->kw_counters[ kw_code ])
		{
			/* If so, write the string and reserve some place afterward in order to store the 
				number of lines */

			fprintf(mesh->handle, "\n%s\n", LM_kw_table[ kw_code ][0]);

			mesh->kw_pos[ kw_code ][0] = ftell(mesh->handle);

			if(!strcmp("i", LM_kw_table[ kw_code ][1]))
				fprintf(mesh->handle, "          \n");

			/* In case of solution field, write the extended header at once */

			if(mesh->sol_headers[ kw_code ])
			{
				fprintf(mesh->handle, "%d ", mesh->sol_headers[ kw_code ][0]);

				for(i=1;i<=mesh->sol_headers[ kw_code ][0];i++)
					fprintf(mesh->handle, "%d ", mesh->sol_headers[ kw_code ][i+1]);

				fprintf(mesh->handle, "\n\n");
			}

			/* Store the positions right after the kw in order to write the number of lines at closing time.
				Store the position right before the data for the write_field process. */

			mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);
		}
		else
		{
			/* If this kw has already be written in the file and has a header, go to pos(0) and write
				down the final value */

			if(strcmp("i", LM_kw_table[ kw_code ][1]))
				return;

			fseek(mesh->handle, mesh->kw_pos[ kw_code][0], SEEK_SET);
			fprintf(mesh->handle, "%d\n", mesh->kw_counters[ kw_code ]);
		}
	}
	else
	{
		/* Test if it is the first time this kw is written */

		if(!mesh->kw_counters[ kw_code ])
		{
			/* If so, write the code, store the position afterward, reserve an int for the next kw position
				and write the number of lines */

			write_int(mesh, kw_code);
			mesh->kw_pos[ kw_code ][0] = ftell(mesh->handle);
			write_int(mesh, 0);

			/* Set the next kw pos in the previously written kw */

			if(mesh->current_kw)
			{
				fseek(mesh->handle, mesh->kw_pos[ mesh->current_kw ][0], SEEK_SET);
				write_int(mesh, mesh->kw_pos[ kw_code ][0] - 4);
				fseek(mesh->handle, mesh->kw_pos[ kw_code ][0] + 4, SEEK_SET);
			}

			mesh->current_kw = kw_code;

			if(!strcmp("i", LM_kw_table[ kw_code ][1]))
				write_int(mesh, 0);

			/* In case of solution field, write the extended header at once */

			if(mesh->sol_headers[ kw_code ])
			{
				write_int(mesh, mesh->sol_headers[ kw_code ][0]);

				for(i=1;i<=mesh->sol_headers[ kw_code ][0];i++)
					write_int(mesh, mesh->sol_headers[ kw_code ][i+1]);
			}

			mesh->kw_pos[ kw_code ][1] = ftell(mesh->handle);
		}
		else
		{
			/* Write the number of lines written at closing time */

			if(strcmp("i", LM_kw_table[ kw_code ][1]))
				return;

			fseek(mesh->handle, mesh->kw_pos[ kw_code ][0] + 4, SEEK_SET);
			write_int(mesh, mesh->kw_counters[ kw_code ]);
		}
	}
}


/*----------------------------------------------------------*/
/* Read an integer in a mesh file							*/
/*----------------------------------------------------------*/

static int read_int(LM_mesh_struct *mesh)
{
	int swaped, integer = 0;

	if(mesh->type & LM_ASCII)
		fscanf(mesh->handle, "%d", &integer);
	else
	{
		/* Read a 4 bytes block in the file */

		fread(&integer, 4, 1, mesh->handle);

		if(mesh->endian != 1)
		{
			swap_bytes((void *)&integer, (void *)&swaped, 4);
			integer = swaped;
		}
	}

	return(integer);
}


/*----------------------------------------------------------*/
/* Write an integer in a mesh file							*/
/*----------------------------------------------------------*/

static void write_int(LM_mesh_struct *mesh, int integer)
{
	if(mesh->type & LM_ASCII)
		fprintf(mesh->handle, "%d ", integer);
	else
		fwrite(&integer, 4, 1, mesh->handle);
}


/*----------------------------------------------------------*/
/* Convert little endian <-> big endian						*/
/*----------------------------------------------------------*/

static void swap_bytes(void *c1, void *c2, int nbytes)
{
  int   k;
  char *c11, *c22;

  c11 = (char*)c1;
  c22 = (char*)c2;

  for (k=0; k<nbytes; k++)
    c22[k] = c11[ nbytes-k-1 ];
}


/*----------------------------------------------------------*/
/* Expand the compacted format and compute the line size	*/
/*----------------------------------------------------------*/

static int expand_format(LM_mesh_struct *mesh, int kw_code, char *out_format)
{
	int i, j, duplicate_next = 0, size = 0;
	char *in_format = LM_kw_table[ kw_code ][2];

	out_format[0] = '\0';

	/* Scan each character of the format string */

	for(i=0;i<strlen(in_format);i++)
	{
		if( (in_format[i] == 'i') || (in_format[i] == 'r') )
		{
			if(!duplicate_next)
				duplicate_next = 1;

			for(j=1;j<=duplicate_next;j++)
			{
				strncat(out_format, &in_format[i], 1);
				size++;
			}

			duplicate_next = 0;
		}
		else if(in_format[i] == 'c')
		{
			strncat(out_format, &in_format[i], 1);
			size += 64;
		}
		else if(in_format[i] == 'd')
		{
			/* 'd' means duplicate the next character mesh->dimension times */

			duplicate_next = mesh->dimension;
		}
		else if(in_format[i] == 's')
		{
			/* 's' means duplicate the next character mesh->solsize times */

			duplicate_next = mesh->sol_headers[ kw_code ][1];
		}
	}

	/* Return the sum of item to be read */

	return(size);
}


/*----------------------------------------------------------*/
/* Read the extended sol headers							*/
/*----------------------------------------------------------*/

static void read_sol_headers(LM_mesh_struct *mesh)
{
	int i, j, nbsol, solsize;

	for(i=1;i<=LM_NBKW;i++)
	{
		/* If the kw is a special Sol field, read the extended header */

		if(!mesh->kw_counters[i] || strcmp(LM_kw_table[i][2], "sr"))
			continue;

		/* Set the curent position in file to the begining of this kw's data */

		fseek(mesh->handle, mesh->kw_pos[i][1], SEEK_SET);

		nbsol = read_int(mesh);

		if( (mesh->sol_headers[i] = malloc((nbsol+2) * sizeof(int))))
		{
			mesh->sol_headers[i][0] = nbsol;
			solsize = 0;

			for(j=1;j<=nbsol;j++)
			{
				mesh->sol_headers[i][j+1] = read_int(mesh);

				switch(mesh->sol_headers[i][j+1])
				{
					case LM_SCALAR     : solsize += 1; break;
					case LM_VECTOR     : solsize += mesh->dimension; break;
					case LM_SYM_MATRIX : solsize += (mesh->dimension * (mesh->dimension+1)) / 2; break;
					case LM_MATRIX     : solsize += mesh->dimension * mesh->dimension; break;
				}

				mesh->sol_headers[i][1] = solsize;
			}

			mesh->kw_pos[i][1] = ftell(mesh->handle);
		}
	}
}
