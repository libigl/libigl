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

#include "jmesh.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

double JMesh::acos_tolerance = 0.0001;
FILE *JMesh::historyFile = NULL;
void (* JMesh::display_message)(char*, int) = NULL;

char *JMesh::app_name = NULL;
char *JMesh::app_version = NULL;
char *JMesh::app_year = NULL;
char *JMesh::app_authors = NULL;
char *JMesh::app_url = NULL;
char *JMesh::app_maillist = NULL;
bool JMesh::quiet = false;

void JMesh::init(double at, FILE *hf, void (*dm)(char *, int))
{
 acos_tolerance = at;
 historyFile = hf;
 display_message = dm;
 app_name = NULL;
 app_version = NULL;
 app_year = NULL;
 app_authors = NULL;
 app_url = NULL;
 app_maillist = NULL;
 quiet = false;
}


///////////// Prints a fatal error message and exits /////////////

void JMesh::error(const char *msg, ...)
{
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"\nERROR- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (historyFile != NULL)
 {
  fclose(historyFile);
  strcat(fms, "Try the '-rescue' command line option.\n");
 }

 if (display_message != NULL)
  display_message(fms, DISPMSG_ACTION_ERRORDIALOG);
 else
 {
  fprintf(stderr,fms);
  exit(-1);
 }
}

///////////// Prints a warning message /////////////

void JMesh::warning(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"WARNING- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (display_message != NULL) 
  display_message(fms, DISPMSG_ACTION_PUTMESSAGE);
 else
  fprintf(stderr,fms);

 va_end(ap);
}

///////////// Prints an information message /////////////

void JMesh::info(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"INFO- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (display_message != NULL) 
  display_message(fms, DISPMSG_ACTION_PUTMESSAGE);
 else
  printf(fms);

 va_end(ap);
}

///////// Reports progress status for a process //////////

void JMesh::begin_progress()
{
 if (quiet) return;
 if (display_message != NULL) 
  display_message("\n", DISPMSG_ACTION_PUTNEWLINE);
 else
  printf("\n");
}

void JMesh::report_progress(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048] = "\r";
 static char fms[4096];
 static char rotating_bar[5] = "-\\|/";
 static unsigned char wc=0;

 if (msg == NULL)
 {
  sprintf(fms,"%c",rotating_bar[wc++]); if (wc==4) wc=0;
  strcpy(fmt+1,fms);

  if (display_message != NULL) 
   display_message(fmt, DISPMSG_ACTION_PUTPROGRESS);
  else
  {
   printf("%s",fmt);
   fflush(stdout);
  }
 }
 else
 {
  va_list ap;
  va_start(ap, msg);
  strcpy(fmt+1,msg);
  vsprintf(fms,fmt,ap);

  if (display_message != NULL) 
   display_message(fms, DISPMSG_ACTION_PUTPROGRESS);
  else
  {
   printf("%s", fms);
   fflush(stdout);
  }
  va_end(ap);
 }
}

void JMesh::end_progress()
{
 if (quiet) return;
 if (display_message != NULL) 
  display_message("\n", DISPMSG_ACTION_PUTNEWLINE);
 else
  printf("\n");
}

