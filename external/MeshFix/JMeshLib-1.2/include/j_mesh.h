/****************************************************************************
* JMeshLib                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef _J_MESH_H
#define _J_MESH_H

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define JMESH_VERSION	"1.1"
#define JMESH_YEAR	2006

class JMesh
{
 public:

 static double acos_tolerance;
 static FILE *historyFile;
 static void (*display_message)(char *, int);
 static char *app_name;
 static char *app_version;
 static char *app_year;
 static char *app_authors;
 static char *app_url;
 static char *app_maillist;
 static bool quiet;

 static void init(double = 0.0001, FILE * =NULL, void (*)(char *, int) = NULL);

 static void info(const char *, ...);
 static void warning(const char *, ...);
 static void error(const char *, ...);
 static void begin_progress();
 static void report_progress(const char *, ...);
 static void end_progress();
};

#define DISPMSG_ACTION_SETWIDGET	1
#define DISPMSG_ACTION_PUTNEWLINE	2
#define DISPMSG_ACTION_PUTPROGRESS	3
#define DISPMSG_ACTION_PUTMESSAGE	4
#define DISPMSG_ACTION_ERRORDIALOG	5


typedef	double coord;
#define COORD_MAX	DBL_MAX
#define COORD_MIN	DBL_MIN

#ifndef _INC_WINDOWS
typedef unsigned char	UBYTE;
typedef   signed char	 BYTE;
typedef unsigned short UINT16;
typedef   signed short	INT16;
#endif

#ifdef IS64BITPLATFORM
typedef long int j_voidint;
#else
typedef int	 j_voidint;
#endif 

#define UBYTE_MAX	255
#define UINT16_MAX	65535

#define FABS(a) (((a)<0)?(-(a)):(a))
#define LOG2(a) (log(a)/log(2))
#define PI2	(M_PI/2.0)
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)(((a)>(b))?(a):(b))
#endif


//////// Swaps two pointers. ///////////////////////////////

inline void p_swap(void **a, void **b) {void *t = *a; *a = *b; *b = t;}

/////////////////////////////////////////////////////////////////////////////////////////////

#endif //_J_MESH_H

