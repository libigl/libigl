#include "medit.h"
#include "extern.h"
#include "sproto.h"

extern ubyte      stereoMode;


int saveMeditFile(char *file,pScene sc) {
  FILE      *out;
  pMaterial  pm;
  time_t     timeptr;
  int        i,k,m;
  char      *ptr,data[128];

  strcpy(data,file);
  ptr = (char *)strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  if ( !strstr(data,".medit") )
    strcat(data,".medit");

  fprintf(stdout,"  Writing %s\n",data);
  out = fopen(data,"w");
  if ( !out ) {
    fprintf(stdout,"  ## Unnable to write file\n");
    return(0);
  }

  /* create standard parameter file */
  fprintf(out,"#  File created with Medit %s\n",ME_VER);
  fprintf(out,"#  Release %s\n",ME_REL);
  time(&timeptr);
  fprintf(out,"# Created: %s",ctime(&timeptr));

  fprintf(out,"\n#  Assign background color\n");
  fprintf(out,"BackgroundColor\n%f %f %f\n",
	  sc->par.back[0],sc->par.back[1],sc->par.back[2]);
  if ( sc->par.linc ) 
    fprintf(out,"LineColor\n%f %f %f\n",
	    sc->par.line[0],sc->par.line[1],sc->par.line[2]);

  if ( sc->mode == HIDDEN )
    fprintf(out,"\nRenderMode\nhidden\n");
  else if ( sc->mode == FILL )
    fprintf(out,"\nRenderMode\nfill\n");
  else if ( sc->mode == S_FILL+S_COLOR+S_MATERIAL )
    fprintf(out,"\nRenderMode\ncolorshading\n");
  else if ( sc->mode == SHADED )
    fprintf(out,"\nRenderMode\nshading\n");
  else if ( sc->mode == S_FILL+S_BDRY+S_COLOR+S_MATERIAL )
    fprintf(out,"\nRenderMode\ncolorshadingline\n");

  fprintf(out,"\n#  Window Size\n");
  fprintf(out,"WindowSize\n%d %d\n",sc->par.xs,sc->par.ys);
  if ( sc->par.sunp ) {
    fprintf(out,"\n#  Source Light\n");
    fprintf(out,"SunPosition\n%f %f %f\n",
            sc->par.sunpos[0]/(2.0*sc->dmax),
            sc->par.sunpos[1]/(2.0*sc->dmax),
            sc->par.sunpos[2]/(2.0*sc->dmax));
  }
  if ( sc->iso.palette ) {
    fprintf(out,"# Color palette\n"); 
    fprintf(out,"Palette\n");
    for (k=0; k<MAXISO; k++)
      fprintf(out,"%f ",sc->iso.val[k]);
    fprintf(out,"\n");
  }
  if ( sc->par.nbmat ) {
    m = 0;
    for (k=0; k<sc->par.nbmat; k++) {
      pm = &sc->material[k];
      for (i=LTria; i<=LHexa; i++)
        if ( pm->depmat[i] && !pm->flag )  m++;
    }
    fprintf(out,"\n# Subdomains colors\n");
    fprintf(out,"NbMaterials\n%d\n",m);
    for (k=0; k<sc->par.nbmat; k++) {
      pm = &sc->material[k];
      m = 0;
      for (i=LTria; i<=LHexa; i++)
        if ( pm->depmat[i] && !pm->flag )  break;
      if ( i > LHexa )  continue;
      fprintf(out,"# Material %s %d\n",pm->name,pm->ref);
      fprintf(out,"Material %s %d\n",pm->name,pm->ref);
      fprintf(out,"%.2f %.2f %.2f %.2f\n",
                  pm->amb[0],pm->amb[1],pm->amb[2],pm->amb[3]);
      fprintf(out,"%.2f %.2f %.2f %.2f\n",
                  pm->dif[0],pm->dif[1],pm->dif[2],pm->dif[3]);
      fprintf(out,"%.2f %.2f %.2f %.2f\n",
                  pm->spe[0],pm->spe[1],pm->spe[2],pm->spe[3]);
      fprintf(out,"%.2f %.2f %.2f %.2f\n",
                  pm->emi[0],pm->emi[1],pm->emi[2],pm->emi[3]);
      fprintf(out,"%.2f\n",pm->shininess);
    }
  }

  if ( stereoMode != MONO )
    fprintf(out,"EyeSep\n%f",sc->par.eyesep);

  fclose(out);
  return(1);
}

/* setup default values */
void iniopt(pScene sc,pMesh mesh) {
  GLfloat  sunpos[4]   = { 0.0, 0.0, 1.0, 1.0};
  int      k;

  sc->par.back[0] = sc->par.back[1] = sc->par.back[2] = 0.0;
  sc->par.back[3] = 1.0;
  sc->par.line[0] = sc->par.line[1] = sc->par.line[2] = 1.0;
  sc->par.line[3] = 1.0;
  sc->par.edge[0] = 1.0;
  sc->par.edge[1] = 0.7;
  sc->par.edge[2] = 0.1;
  sc->par.edge[3] = 1.0;

  memcpy(sc->par.sunpos,sunpos,4*sizeof(GLfloat));
  sc->par.nbmat = -1;
  sc->par.sunp  =  0;
  sc->par.linc  =  0;
  sc->par.linewidth = mesh->dim==3 ? 3 : 2;
  sc->par.pointsize = 3.0;

  for (k=0; k<MAXISO; k++) {
    sc->iso.val[k] = 0.0;
    sc->iso.col[k] = 0.0;
  }
  sc->iso.palette = 0;

  /* postscript */
  sc->par.dpi = 300.;
  sc->par.cm  = 10.;
  sc->par.coeff = 0.0;
  
  /* streamlines */
  sc->par.maxtime = FLT_MAX;
  sc->par.dt      = FLT_MAX;
  sc->par.nbpart  = 1;

  /* clip plane */
  sc->par.clip[0] = 0.0;
  sc->par.clip[1] = 0.0;
  sc->par.clip[2] = 0.0;

  sc->par.clip[3] = -1.0;
  sc->par.clip[4] =  0.0;
  sc->par.clip[5] =  0.0;
  
  sc->par.eyesep  = -1.0;
}


/* parse the program options */
int parsop(pScene sc,pMesh mesh) {
  FILE      *in;
  pMaterial  pm;
  double     dd;
  float      r,g,b,ca,cb,cc,na,nb,nc;
  int        k,i,m,n,xs,ys,ret,ref,nbmat;
  char      *ptr,ub,data[128],key[256],buf[256],pscol[32];
  
  /* check if user-specified parameters */
  iniopt(sc,mesh);

  strcpy(data,mesh->name);
  ptr = (char *)strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  in = 0;
  if ( !strstr(data,".medit") ) {
    strcat(data,".medit");
    in = fopen(data,"r"); 
  }
  if ( !in ) {
    sprintf(data,"%s",DEFAULT_FILE);
    in = fopen(data,"r");
    if ( !in ) {
      fprintf(stdout,"   Loading default options\n");
      sc->par.nbmat = MAX_MATERIAL;
      matInit(sc);
      return(1);
    }
  }
  if ( !quiet )  fprintf(stdout,"   Reading %s\n",data);
  
  m = n = 0;
  while ( !feof(in) ) {
    fscanf(in,"%s",key);
    if ( feof(in) )  break;
    for (i=0; i<strlen(key); i++) key[i] = tolower(key[i]);

    if ( key[0] == '#' ) {
      fgets(key,255,in);
    }
    else if ( !strcmp(key,"backgroundcolor") ) {
      fscanf(in,"%f %f %f",&r,&g,&b);
      sc->par.back[0] = r;
      sc->par.back[1] = g;
      sc->par.back[2] = b;
      sc->par.back[3] = 1.0f;
    }
    else if ( !strcmp(key,"boundingbox") ) {
      fscanf(in,"%f %f %f %f %f %f",&ca,&cb,&cc,&na,&nb,&nc);
      mesh->xmin = ca;
      mesh->ymin = cb;
      mesh->zmin = cc;
      mesh->xmax = na;
      mesh->ymax = nb;
      mesh->zmax = nc;
    }
    else if ( !strcmp(key,"clipplane") ) {
      fscanf(in,"%f %f %f %f %f %f",&ca,&cb,&cc,&na,&nb,&nc);
      sc->par.clip[0] = ca - mesh->xtra;
      sc->par.clip[1] = cb - mesh->ytra;
      sc->par.clip[2] = cc - mesh->ztra;
      dd = sqrt(na*na + nb*nb + nc*nc);
      if ( dd > EPS ) {
        sc->par.clip[3] = na / dd;
        sc->par.clip[4] = nb / dd;
        sc->par.clip[5] = nc / dd;
      }
    }
    else if ( !strcmp(key,"linecolor") ) {
      fscanf(in,"%f %f %f",&r,&g,&b);
      sc->par.line[0] = r;
      sc->par.line[1] = g;
      sc->par.line[2] = b;
      sc->par.linc = 1;
    }
    else if ( !strcmp(key,"linewidth") ) {
      fscanf(in,"%f",&r);
      sc->par.linewidth = MEDIT_MAX(1.0,MEDIT_MIN(10.0,r));
      sc->par.linc = 1;
    }
    else if ( !strcmp(key,"pointsize") ) {
      fscanf(in,"%f",&r);
      sc->par.pointsize = MEDIT_MAX(1.0,MEDIT_MIN(10.0,r));
      sc->par.linc = 1;
    }
    else if ( !strcmp(key,"edgecolor") ) {
      fscanf(in,"%f %f %f",&r,&g,&b);
      sc->par.edge[0] = r;  sc->par.edge[1] = g;  sc->par.edge[2] = b;
      sc->par.linc = 1;
    }
    else if ( !strcmp(key,"sunposition") ) {
      fscanf(in,"%f %f %f",&r,&g,&b);
      sc->dmax = mesh->xmax - mesh->xmin;
      sc->dmax = MEDIT_MAX(sc->dmax,mesh->ymax - mesh->ymin);
      sc->dmax = MEDIT_MAX(sc->dmax,mesh->zmax - mesh->zmin);
      sc->dmax = fabs(sc->dmax);
      sc->par.sunpos[0] = 2.0*sc->dmax*r;  
      sc->par.sunpos[1] = 2.0*sc->dmax*g;
      sc->par.sunpos[2] = 2.0*sc->dmax*b;
      sc->par.sunp = 1;
    }
    else if ( !strcmp(key,"windowsize") ) {
      fscanf(in,"%d %d",&xs,&ys);
      sc->par.xs = (short)xs;  sc->par.ys = (short)ys;
    }
    else if ( !strcmp(key,"windowposition") ) {
      fscanf(in,"%d %d",&xs,&ys);
      sc->par.pxs = (short)xs;  sc->par.pys = (short)ys;
    }
    else if ( !strcmp(key,"rendermode") ) {
      fscanf(in,"%s",buf);
      for (i=0; i<strlen(buf); i++) buf[i] = tolower(buf[i]);
      if ( strstr(buf,"hidden") )
        sc->mode = HIDDEN;
      else if ( strstr(buf,"fill") )
        sc->mode = FILL;
      else if ( strstr(buf,"colorshadingline") )
        sc->mode = SHADED + S_MATERIAL;
      else if ( strstr(buf,"colorshading") )
        sc->mode = S_FILL + S_COLOR + S_MATERIAL;
      else if ( strstr(buf,"shading") )
        sc->mode = SHADED; 
    }
    else if ( strstr(key,"palette") ) {
      sc->iso.palette = 1;
      if ( !strcmp(key,"palettet") ) 
        sc->iso.palette = 1;
      else if ( !strcmp(key,"paletteb") )
        sc->iso.palette = 2;
      else if ( !strcmp(key,"palettel") )
        sc->iso.palette = 3;
      else if ( !strcmp(key,"paletter") )
        sc->iso.palette = 4;

      for (k=0; k<MAXISO; k++)
        ret = fscanf(in,"%f",&sc->iso.val[k]);
      if ( sc->iso.val[MAXISO-1] < sc->iso.val[0] ) sc->iso.palette = 0;
    }
    else if ( !strcmp(key,"postscript") ) {
      fscanf(in,"%f %f %s %s",&sc->par.cm,&sc->par.dpi,
                              buf,pscol);
      strncpy(sc->par.pscolor,pscol,10);
      sc->par.coeff = atof(buf);
      if ( sc->par.coeff < 0.0f ) sc->par.coeff = 0.0f;
      if ( sc->par.coeff > 1.0f ) sc->par.coeff = 1.0f;
    }
    else if ( !strcmp(key,"time") ) {
      ret = fscanf(in,"%f %f %f",&sc->par.maxtime,&sc->par.pertime,&sc->par.dt);
      if ( !EatSpace(in) ) {
        fscanf(in,"%c",&ub);
        sc->par.nbpart = MEDIT_MAX(atoi(&ub),1);
      }
    }
    else if ( !strcmp(key,"nbmaterial") ) {
      fscanf(in,"%d",&nbmat);
      sc->par.nbmat = MEDIT_MAX(2,nbmat);
      matInit(sc);
    }
    else if ( !strcmp(key,"material") ) {
      if ( sc->par.nbmat == -1 ) {
        sc->par.nbmat = MAX_MATERIAL;
        matInit(sc);
      }
      fgets(buf,255,in);
      if ( n > sc->par.nbmat ) continue;
      ret = sscanf(buf,"%s %d",buf,&ref);
      ptr = strstr(buf,"DEFAULT");
      pm  = ptr ? &sc->material[DEFAULT_MAT] : &sc->material[++n];
      strcpy(pm->name,buf);
      if ( ret < 2 )  ref = 0;
      pm->ref = ref ? ref : n;
      fscanf(in,"%f %f %f %f",&pm->amb[0],&pm->amb[1],&pm->amb[2],&pm->amb[3]);
      fscanf(in,"%f %f %f %f",&pm->dif[0],&pm->dif[1],&pm->dif[2],&pm->dif[3]);
      fscanf(in,"%f %f %f %f",&pm->spe[0],&pm->spe[1],&pm->spe[2],&pm->spe[3]);
      fscanf(in,"%f %f %f %f",&pm->emi[0],&pm->emi[1],&pm->emi[2],&pm->emi[3]);
      fscanf(in,"%f",&pm->shininess);

      if ( pm->amb[3] == 0.0 )  pm->amb[3] = 1.0;
      if ( pm->dif[3] == 0.0 )  pm->dif[3] = 1.0;
      if ( pm->spe[3] == 0.0 )  pm->spe[3] = 1.0;
      if ( pm->emi[3] == 0.0 )  pm->emi[3] = 1.0;
      pm->shininess = MEDIT_MIN(fabs(pm->shininess),128.0f);
      pm->shininess = MEDIT_MAX(pm->shininess,3.0f);
      ++m;
    }
    
    /* stereo mode */
    else if ( !strcmp(key,"eyesep") ) {
      fscanf(in,"%f",&sc->par.eyesep);
    }
    //// Alec: send key commands
    //else if ( !strcmp(key,"keyscene") )
    //{
    //  char key[256];
    //  if(fscanf(in,"%1s",key) == 1)
    //  {
    //    keyScene(key[0],0,0);
    //  }
    //}
    else if(!strcmp(key,"mode"))
    {
      char mode[255];
      fscanf(in,"%s",mode);
      if(!strcmp(mode,"S_MAP"))
      {
        // Alec: Q: Why doesn't this work?
        //sc->mode |= S_MAP;
      }
    }

  }
  fclose(in);

  if ( sc->par.nbmat < 0 ) {
   sc->par.nbmat = MAX_MATERIAL;
   matInit(sc);
  }
  else if ( m == n )  sc->par.nbmat++;
  if ( !sc->par.linc ) {
    sc->par.line[0] = 1.0 - sc->par.back[0];
    sc->par.line[1] = 1.0 - sc->par.back[1];
    sc->par.line[2] = 1.0 - sc->par.back[2];
  }

  if ( ddebug ) fprintf(stdout,"    Materials %8d\n",sc->par.nbmat);
  return(1);
}
