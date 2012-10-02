#include "medit.h"
#include "extern.h"
#include "sproto.h"

int EatLine(FILE  *in) {
  int    k,c;  
  
  k = 0;
  while ( (c=fgetc(in)) != EOF ) {
    k++;
    if ( c==10 || c== 13)  return(1);
  }
  return(0);
}


int EatSpace(FILE  *in) {
  int   k,c,ret=0;  
  
  k=0;
  while (isspace(c=fgetc(in)) ) {
    k++;
    if ( c==10 || c== 13 ) ret=1;
  }
  if ( c==EOF) return -1;
  
  /*  fprintf(stdout,"EatSpace %d %d last char=%d\n",ret,k,c); */
  ungetc(c,in);
  
  return ret;
}


int bbfile(pMesh mesh) {
  FILE      *in;
  pSolution  ps;
  double     a,b,c,lambda[3],eigv[3][3],m[6],vp[2][2];
  float      dummy;
  int        j,k,l,dim,np,nfield,nf,i1,i2,i3,i4,err,iord;
  char      *ptr,data[128],tmp[128];
  ubyte      bigbb;

  /* default */
  strcpy(tmp,mesh->name);
  ptr = (char *)strstr(tmp,".mesh");
  if ( ptr ) *ptr = '\0';

  sprintf(data,"%s.bb",tmp);
  in = fopen(data,"r");
  bigbb = 0;
  if ( !in ) {
    sprintf(data,"%s.pbb",tmp);
    in = fopen(data,"r");
  }  
  if ( !in ) {
    bigbb = 1;
    sprintf(data,"%s.BB",tmp);
    in = fopen(data,"r");  
    if ( !in ) { /* hack FH pour le mac */
      sprintf(data,"%s.gbb",tmp);
      in = fopen(data,"r");
    }
    
  }
  if ( !in )    
    return(0);
  
  /* if ( !quiet )  fprintf(stdout,"  Reading %s\n",data); */
  i1=i2=i3=i4=-1;
  /* read file format */
  err=0;
  fscanf(in,"%d",&dim);
  if(EatSpace(in)) err++;
  fscanf(in,"%d",&i1);
  if(EatSpace(in)) err++;
  fscanf(in,"%d",&i2);
  if(EatSpace(in)) err++;
  fscanf(in,"%d",&i3);
  bigbb=  (EatSpace(in)==0); /* not nl after the 4 integer => BB */

  if ( !quiet )
    if(bigbb)  fprintf(stdout,"  Reading BB file %s\n",data);
    else  fprintf(stdout,"  Reading bb file %s\n",data);

  if ( dim < 2 || dim > 3 || err ) {
    fprintf(stderr,"  %%%% Wrong file (dim=%d) (err=%d). Ignored\n",dim,err);
    return(0);
  }
  /* read number of field(s) */
  nf = 0;

  if ( bigbb ) {
    /* get only 1st field */
    /* fscanf(in,"%d",&nfield);*/
    nfield=i1;
    /*fscanf(in,"%d",&mesh->nfield);*/
    mesh->nfield = i2;
    if (nfield>1) 
      {
	nf += i3;
	for (k=1; k<nfield-1; k++) {
	  fscanf(in,"%d",&np);
	  nf += np;
	}
	fscanf(in,"%d",&np);
      }
    else 
      np = i3;
    /* read file type */
    fscanf(in,"%d",&mesh->typage);
    printf(" np= %d, type= %d\n",np,mesh->typage);
  }
  else {
   /* fscanf(in,"%d",&mesh->nfield);
      fscanf(in,"%d",&np);*/
    /* read file type */
    /* fscanf(in,"%d",&mesh->typage);*/
    mesh->nfield=i1;
    np=i2;
    mesh->typage=i3;
  }  
  

  if ( mesh->typage == 2 ) {
    if ( np < mesh->np ) {
      fprintf(stderr,"  %%%% Wrong solution number (%d , %d). Ignored\n",np,mesh->np);
      fclose(in);
      return(0);
    }
    mesh->nbb = mesh->np;
  }
  else if ( mesh->typage == 1 ) {
    if ( np < mesh->ne ) {
      fprintf(stderr,"  %%%% Wrong solution number (%d , %d). Ignored\n",np,mesh->ne);
      fclose(in);
      return(0);
    }
    mesh->nbb = mesh->ne;
  }
  else {
    fprintf(stderr,"  %%%% Wrong typage (%d). Ignored\n",mesh->typage);
    fclose(in);
    return(0);
  }

  /* read solutions */
  mesh->bbmin  =  1.e10;
  mesh->bbmax  = -1.e10;

  /* allocate memory */
  if ( !zaldy2(mesh) ) {
    mesh->nbb = 0;
    fclose(in);
    return(0);
  }

  /* scalar field */
  if ( mesh->nfield == 1 ) {
    if ( ddebug )  printf("   scalar (isotropic) field\n");
    for (k=1; k<=mesh->nbb; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0;
      if ( fscanf(in,"%s",data) != 1 )  continue;
      if ( ptr = strpbrk(data,"dD") ) *ptr = 'E';
      sscanf(data,"%f",&ps->bb);
      if ( ps->bb < mesh->bbmin )  mesh->bbmin = ps->bb;
      if ( ps->bb > mesh->bbmax )  mesh->bbmax = ps->bb;
      for (j=1; j<=nf; j++)  fscanf(in,"%f",&dummy);
    }
  }

  /* vector field */
  else if ( mesh->nfield == mesh->dim ) {
    if ( ddebug )  fprintf(stdout,"   vector field \n");
    for (k=1; k<=mesh->nbb; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0;
      for (l=0; l<mesh->dim; l++) {
        if ( fscanf(in,"%s",data) != 1 )  continue;
        if ( ptr = strpbrk(data,"dD") ) *ptr = 'E';
        sscanf(data,"%f",&ps->m[l]);
        ps->bb += ps->m[l]*ps->m[l];
      }
      ps->bb = sqrt(ps->bb);
      if ( ps->bb < mesh->bbmin )  mesh->bbmin = ps->bb;
      if ( ps->bb > mesh->bbmax )
        mesh->bbmax = ps->bb;
      for (j=1; j<nf; j++)  fscanf(in,"%f",&dummy);
    }
    fclose(in);
    return(0);
  }
  else if ( dim == 2 && mesh->nfield == 3 ) {
    if ( ddebug )  fprintf(stdout,"   2D metric field\n");
    for (k=1; k<=mesh->np; k++) {
      ps = &mesh->sol[k];
      fscanf(in,"%lf %lf %lf",&a,&b,&c);
      ps->m[0] = a;
      ps->m[1] = b;
      ps->m[2] = c;
      m[0] = a;
      m[1] = b;
      m[2] = c;
      eigen2(m,lambda,vp);
      ps->bb = min(lambda[0],lambda[1]);
      if ( ps->bb < mesh->bbmin )  mesh->bbmin = ps->bb;
      if ( ps->bb > mesh->bbmax )  mesh->bbmax = ps->bb;
      for (j=1; j<nf; j++)  fscanf(in,"%f",&dummy);
    }
  }
  else if ( dim == 3 && mesh->nfield == 6 ) {
    if ( ddebug )  fprintf(stdout,"   3D metric field\n");
    for (k=1; k<=mesh->np; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0f;
      for (l=0; l<6; l++) {
        if ( fscanf(in,"%s",data) != 1 )  continue;
        if ( ptr = strpbrk(data,"dD") ) *ptr = 'E';
        sscanf(data,"%f",&dummy);
        m[l] = dummy;
      }
      ps->m[0] = m[0];
      ps->m[1] = m[1];
      ps->m[2] = m[3];
      ps->m[3] = m[2];
      ps->m[4] = m[4];
      ps->m[5] = m[5];

      m[2] = ps->m[2];
      m[3] = ps->m[3];
      iord = eigenv(1,m,lambda,eigv);
      if ( iord ) {
        ps->bb = lambda[0];
        ps->bb = max(ps->bb,lambda[1]);
        ps->bb = max(ps->bb,lambda[2]);
        if ( ps->bb < mesh->bbmin )  mesh->bbmin = ps->bb;
        if ( ps->bb > mesh->bbmax )  mesh->bbmax = ps->bb;
      }
      else {
        fprintf(stdout,"  ## Eigenvalue problem.\n");
      }
      
      for (j=1; j<nf; j++)  fscanf(in,"%f",&dummy);
    }
  }
  else {
    fprintf(stderr," %%%% Solution not suitable. Ignored\n");
    mesh->nbb = 0;
  }

  fclose(in);
  return(np);
}
