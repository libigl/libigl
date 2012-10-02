#include "medit.h"
#include "extern.h"
#include "sproto.h"


int inmsh2(pMesh mesh) {
  FILE      *inp,*inf;
  pPoint     ppt,pp0,pp1,pp2,pp3;
  pTriangle  pt1;
  pQuad      pq1;
  pEdge      pr;
  int        k,disc,ret,degree,dum,ref,tag;
  char      *ptr,data[256],sx[128],sy[128],sz[128];

  /* check for .points */
  strcpy(data,mesh->name);
  strcat(data,".points");
  inp = fopen(data,"r");
  if ( !inp ) return(0);

  /* check for .faces */
  strcpy(data,mesh->name);
  strcat(data,".faces");
  inf = fopen(data,"r");
  if ( !inf ) return(0);
  if ( !quiet )  fprintf(stdout,"  Reading %s.{points,.faces}\n",mesh->name);

  /* get number of elements */
  /*fgets(data,255,inp);
    sscanf(data,"%d",&mesh->np);*/
  fscanf(inp,"%d",&mesh->np);
  EatLine(inp);
  
  /*fgets(data,255,inf);
  sscanf(data,"%d",&mesh->ne);*/
  fscanf(inf,"%d",&mesh->ne);
  EatLine(inf);
  if ( !mesh->np ) { /*|| (mesh->dim == 3 && !mesh->ne) ) {*/
	fprintf(stdout,"  ## No vertex.\n");
    fclose(inp);
    return(-1);
  }
  mesh->dim = 3;
  mesh->nt = mesh->nq = mesh->ntet = mesh->nhex = mesh->nvn = 0;

  /* first pass get number of faces */
  for (k=1; k<=mesh->ne; k++) {
    fscanf(inf,"%d",&degree);
    if ( degree < 2 || degree > 4 ) {
      fprintf(stdout,"  ## Wrong degree\n");
      fclose(inp);
      fclose(inf);
      return(0);
    }
    else if ( degree == 2 )
      mesh->na++;
    else if ( degree == 3 )
      mesh->nt++;
    else if ( degree == 4 )
      mesh->nq++;
    /*fgets(data,80,inf);*/
    EatLine(inf);
  }
  
  /* check if vertices and elements found */
  if ( !mesh->np ) { /*|| mesh->ne == 0 ) {*/
    fclose(inp);
    fclose(inf);
    return(0);
  }

  /* memory allocation for mesh */
  if ( zaldy1(mesh) != TRUE ) {
    fclose(inp);
    fclose(inf);
    return(0);
  }
  
  /* read mesh vertices */
  for(k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    /* parse coordinates into strings */
    ret = fscanf(inp,"%s %s %s %d",sx,sy,sz,&ref);
    if ( ret != 4 ) {
      fclose(inp);
      fclose(inf);
      return(0);
    }
    if ( ptr = strpbrk(sx,"dD") )
      *ptr = 'E';
    if ( ptr = strpbrk(sy,"dD") )
      *ptr = 'E';
    if ( ptr = strpbrk(sz,"dD") )
      *ptr = 'E';
    /*
    sscanf(sx,"%f",&ppt->c[0]);
    sscanf(sy,"%f",&ppt->c[1]);
    sscanf(sz,"%f",&ppt->c[2]);
    */
    ppt->c[0] = atof(sx);
    ppt->c[1] = atof(sy);
    ppt->c[2] = atof(sz);

    ppt->ref = ref;
    ppt->tag = M_UNUSED;
  }
  fclose(inp);

  /* allocate memory for mesh edges */
  if ( mesh->na > 0 ) {
    mesh->edge = (pEdge)M_calloc(mesh->na+1,sizeof(Edge),"inmsh2.edge");
    if ( !mesh->edge ) {
      fprintf(stderr,"  ## WARN 0004, INMESH, %d\n",mesh->na);
      fclose(inf);
      return(1);
    }
  }

  /* read mesh faces */
  rewind(inf);
  /*fgets(data,255,inf);
  sscanf(data,"%d",&mesh->ne);*/
  fscanf(inf,"%d",&mesh->ne);
  EatLine(inf);
  mesh->nt = 0;
  mesh->nq = 0;
  mesh->na = 0;
  disc     = 0;
  for (k=1; k<=mesh->ne; k++) {
    fscanf(inf,"%d",&degree);
    
    if ( degree == 2 ) {
      pr = &mesh->edge[++mesh->na];
      fscanf(inf,"%d %d %d %d %d\n",&pr->v[0],&pr->v[1],&tag,&dum,&dum);
      pr->tag = tag == 0 ? M_NOTAG : M_TAG; 
      pp0 = &mesh->point[pr->v[0]];
      pp1 = &mesh->point[pr->v[1]];
      pp0->tag = M_NOTAG;
      pp1->tag = M_NOTAG;
    }
    else if ( degree == 3 ) {
      pt1 = &mesh->tria[++mesh->nt];
      fscanf(inf,"%d %d %d %d %d %d %d\n",&pt1->v[0],&pt1->v[1],&pt1->v[2],
	     &ref,&dum,&dum,&dum);
      if ( pt1->v[0] <= 0 || pt1->v[0] > mesh->np ||
           pt1->v[1] <= 0 || pt1->v[1] > mesh->np ||
           pt1->v[2] <= 0 || pt1->v[2] > mesh->np ) {
	    fprintf(stdout,"  ## Wrong index\n");
	    disc++;
	    pt1->v[0] = 0;
	    continue;
      }
      pt1->ref = fabs(ref);
      pp0 = &mesh->point[pt1->v[0]];
      pp1 = &mesh->point[pt1->v[1]];
      pp2 = &mesh->point[pt1->v[2]];
      pp0->tag = M_NOTAG;
      pp1->tag = M_NOTAG;
      pp2->tag = M_NOTAG;
    }
    else if ( degree == 4 ) {
      pq1 = &mesh->quad[++mesh->nq];
      fscanf(inf,"%d %d %d %d",&pq1->v[0],&pq1->v[1],&pq1->v[2],&pq1->v[3]);
      fscanf(inf,"%d %d %d %d %d",&ref,&dum,&dum,&dum,&dum);
      if ( pq1->v[0] <= 0 || pq1->v[0] > mesh->np ||
           pq1->v[1] <= 0 || pq1->v[1] > mesh->np ||
           pq1->v[2] <= 0 || pq1->v[2] > mesh->np ||
           pq1->v[3] <= 0 || pq1->v[3] > mesh->np ) {
		fprintf(stdout,"  ## Wrong index\n");
	    disc++;
	    pq1->v[0] = 0;
	    continue;
      }

      pq1->ref = fabs(ref);
      pp0 = &mesh->point[ pq1->v[0] ];
      pp1 = &mesh->point[ pq1->v[1] ];
      pp2 = &mesh->point[ pq1->v[2] ];
      pp3 = &mesh->point[ pq1->v[3] ];
      pp0->tag = M_NOTAG;
      pp1->tag = M_NOTAG;
      pp2->tag = M_NOTAG;
      pp3->tag = M_NOTAG;
    }
  }

  fclose(inf);
  if ( disc > 0 ) {
    fprintf(stdout," ## %d entities discarded\n",disc);
  }
  
  return(1);
}
