#include "exttrimesh.h"
#include <string.h>
#include <stdlib.h>
#include <jrs_predicates.h>

const char *input_filename;
double epsilon_angle = 0.0;

// Simulates the ASCII rounding error
void asciiAlign(ExtTriMesh& tin)
{
 char outname[2048];
 Vertex *v;
 Node *m;
 float a;
 FOREACHVVVERTEX((&(tin.V)), v, m)
 {
  sprintf(outname,"%f",v->x); sscanf(outname,"%f",&a); v->x = a;
  sprintf(outname,"%f",v->y); sscanf(outname,"%f",&a); v->y = a;
  sprintf(outname,"%f",v->z); sscanf(outname,"%f",&a); v->z = a;
 }
}


// Return TRUE if the triangle is exactly degenerate

inline bool isDegenerateEdge(Edge *e)
{
 return ((*(e->v1))==(*(e->v2)));
}

bool isDegenerateTriangle(Triangle *t)
{
 double xy1[2], xy2[2], xy3[2];
 xy1[0] = t->v1()->x; xy1[1] = t->v1()->y; 
 xy2[0] = t->v2()->x; xy2[1] = t->v2()->y; 
 xy3[0] = t->v3()->x; xy3[1] = t->v3()->y; 
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 xy1[0] = t->v1()->y; xy1[1] = t->v1()->z; 
 xy2[0] = t->v2()->y; xy2[1] = t->v2()->z; 
 xy3[0] = t->v3()->y; xy3[1] = t->v3()->z; 
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 xy1[0] = t->v1()->z; xy1[1] = t->v1()->x; 
 xy2[0] = t->v2()->z; xy2[1] = t->v2()->x; 
 xy3[0] = t->v3()->z; xy3[1] = t->v3()->x; 
 if (orient2d(xy1, xy2, xy3)!=0.0) return false;
 return true;
}


Edge *getLongestEdge(Triangle *t)
{
 double l1 = t->e1->squaredLength();
 double l2 = t->e2->squaredLength();
 double l3 = t->e3->squaredLength();
 if (l1>=l2 && l1>=l3) return t->e1;
 if (l2>=l1 && l2>=l3) return t->e2;
 return t->e3;
}


// Iterate on all the selected triangles as long as possible.
// Keep the selection only on the degeneracies that could not be removed.
// Return the number of degeneracies that could not be removed
int swap_and_collapse(ExtTriMesh *tin)
{
 Node *n;
 Triangle *t;

 if (epsilon_angle != 0.0)
 {
  FOREACHVTTRIANGLE((&(tin->T)), t, n) UNMARK_VISIT(t);
  JMesh::quiet = true; tin->removeDegenerateTriangles(); JMesh::quiet = false; 
  int failed = 0;
  FOREACHVTTRIANGLE((&(tin->T)), t, n) if (IS_VISITED(t)) failed++;
  return failed;
 }

 List triangles;
 Edge *e;
 const int MAX_ATTEMPTS = 10;

 FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info=0;

 // VISIT2 means that the triangle is in the list
 FOREACHVTTRIANGLE((&(tin->T)), t, n) if (IS_VISITED(t))
 {
  UNMARK_VISIT(t);
  if (isDegenerateTriangle(t)) {triangles.appendTail(t); MARK_VISIT2(t);}
 }

 while ((t=(Triangle *)triangles.popHead())!=NULL)
 {
  UNMARK_VISIT2(t);
  if (t->isLinked())
  {
   if (isDegenerateEdge(t->e1)) t->e1->collapse();
   else if (isDegenerateEdge(t->e2)) t->e2->collapse();
   else if (isDegenerateEdge(t->e3)) t->e3->collapse();
   else if ((e=getLongestEdge(t))!=NULL)
   {
    if (e->swap())
	{
	 t=e->t1;
           // Alec: replaced "int" with "j_voidint"
	 if (isDegenerateTriangle(t) && !IS_VISITED2(t) && ((j_voidint)t->info < MAX_ATTEMPTS))
	  {triangles.appendTail(t); MARK_VISIT2(t); t->info = (void *)(((j_voidint)t->info)+1);}
	 t=e->t2;
	 if (isDegenerateTriangle(t) && !IS_VISITED2(t) && ((j_voidint)t->info < MAX_ATTEMPTS))
	  {triangles.appendTail(t); MARK_VISIT2(t); t->info = (void *)(((j_voidint)t->info)+1);}
	}
   }
  }
 }

 tin->removeUnlinkedElements();

 int failed=0;
 // This should check only on actually processed triangles
 FOREACHVTTRIANGLE((&(tin->T)), t, n) if (isDegenerateTriangle(t)) {failed++; MARK_VISIT(t);}

 JMesh::info("%d degeneracies selected\n",failed);
 return failed;
}

// returns true on success

bool removeDegenerateTriangles(ExtTriMesh& tin, int max_iters)
{
 int n, iter_count = 0;

 printf("Removing degeneracies...\n");
 while ((++iter_count) <= max_iters && swap_and_collapse(&tin))
 {
  for (n=1; n<iter_count; n++) tin.growSelection();
  tin.removeSelectedTriangles();
  tin.removeSmallestComponents();
  JMesh::quiet = true; tin.fillSmallBoundaries(tin.E.numels()); JMesh::quiet = false;
  asciiAlign(tin);
 }

 if (iter_count > max_iters) return false;
 return true;
}









bool appendCubeToList(Triangle *t0, List& l)
{
 if (!IS_VISITED(t0) || IS_VISITED2(t0)) return false;

 Triangle *t, *s;
 Vertex *v;
 List triList(t0);
 MARK_VISIT2(t0);
 double minx=DBL_MAX, maxx=-DBL_MAX, miny=DBL_MAX, maxy=-DBL_MAX, minz=DBL_MAX, maxz=-DBL_MAX;

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  v = t->v1();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v2();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v3();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  if ((s = t->t1()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
 }

 l.appendTail(new Point(minx, miny, minz));
 l.appendTail(new Point(maxx, maxy, maxz));
 return true;
}

bool isVertexInCube(Vertex *v, List& loc)
{
 Node *n;
 Point *p1, *p2;
 FOREACHNODE(loc, n)
 {
  p1 = (Point *)n->data; n=n->next(); p2 = (Point *)n->data;
  if (!(v->x < p1->x || v->y < p1->y || v->z < p1->z ||
      v->x > p2->x || v->y > p2->y || v->z > p2->z)) return true;
 }

 return false;
}

void selectTrianglesInCubes(ExtTriMesh& tin)
{
 Triangle *t;
 Vertex *v;
 Node *n;
 List loc;
 FOREACHVTTRIANGLE((&(tin.T)), t, n) appendCubeToList(t, loc);
 FOREACHVVVERTEX((&(tin.V)), v, n) if (isVertexInCube(v, loc)) MARK_VISIT(v);
 FOREACHVTTRIANGLE((&(tin.T)), t, n)
 {
  UNMARK_VISIT2(t);
  if (IS_VISITED(t->v1()) || IS_VISITED(t->v2()) || IS_VISITED(t->v3())) MARK_VISIT(t);
 }
 FOREACHVVVERTEX((&(tin.V)), v, n) UNMARK_VISIT(v);
 loc.freeNodes();
}







// returns true on success

bool removeSelfIntersections(ExtTriMesh& tin, int max_iters)
{
 int n, iter_count = 0;

 printf("Removing self-intersections...\n");
 while ((++iter_count) <= max_iters && tin.selectIntersectingTriangles())
 {
  for (n=1; n<iter_count; n++) tin.growSelection();
  tin.removeSelectedTriangles();
  tin.removeSmallestComponents();
  JMesh::quiet = true; tin.fillSmallBoundaries(tin.E.numels()); JMesh::quiet = false;
  asciiAlign(tin);
  selectTrianglesInCubes(tin);
 }

 if (iter_count > max_iters) return false;
 return true;
}


bool isDegeneracyFree(ExtTriMesh& tin)
{
 Node *n;
 Triangle *t;

 if (epsilon_angle != 0.0)
 {FOREACHVTTRIANGLE((&(tin.T)), t, n) if (t->isDegenerate()) return false;}
 else
 {FOREACHVTTRIANGLE((&(tin.T)), t, n) if (isDegenerateTriangle(t)) return false;}

 return true;
}


// returns true on success

bool meshclean(ExtTriMesh& tin, int max_iters = 10, int inner_loops = 3)
{
 bool ni, nd;

 tin.deselectTriangles();
 tin.invertSelection();

 for (int n=0; n<max_iters; n++)
 {
  printf("********* ITERATION %d *********\n",n);
  nd=removeDegenerateTriangles(tin, inner_loops);
  tin.deselectTriangles(); tin.invertSelection();
  ni=removeSelfIntersections(tin, inner_loops);
  if (ni && nd && isDegeneracyFree(tin)) return true;
 }

 return false;
}



double closestPair(List *bl1, List *bl2, Vertex **closest_on_bl1, Vertex **closest_on_bl2)
{
 Node *n, *m;
 Vertex *v,*w;
 double adist, mindist = DBL_MAX;

 FOREACHVVVERTEX(bl1, v, n)
  FOREACHVVVERTEX(bl2, w, m)
   if ((adist = w->squaredDistance(v))<mindist)
   {
	mindist=adist;
	*closest_on_bl1 = v;
	*closest_on_bl2 = w;
   }

 return mindist;
}

bool joinClosestComponents(ExtTriMesh *tin)
{
  Vertex *v,*w, *gv, *gw;
  Triangle *t, *s;
  Node *n;
  List triList, boundary_loops, *one_loop;
  List **bloops_array;
  int i, j, numloops;

  i=0;
  FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
  FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info == NULL)
  {
   i++;
   triList.appendHead(t);
   t->info = (void *)i;

   while(triList.numels())
   {
    t = (Triangle *)triList.popHead();
    if ((s = t->t1()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
    if ((s = t->t2()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
    if ((s = t->t3()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)i;}
   }
  }

  if (i<2)
  {
   FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
//   JMesh::info("Mesh is a single component. Nothing done.");
   return false;
  }

  FOREACHVTTRIANGLE((&(tin->T)), t, n)
  {
   t->v1()->info = t->v2()->info = t->v3()->info = t->info;
  }

  FOREACHVVVERTEX((&(tin->V)), v, n) if (!IS_VISITED2(v) && v->isOnBoundary())
  {
   w = v;
   one_loop = new List;
   do
   {
    one_loop->appendHead(w); MARK_VISIT2(w);
    w = w->nextOnBoundary();
   } while (w != v);
   boundary_loops.appendHead(one_loop);
  }
  FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_VISIT2(v);

  bloops_array = (List **)boundary_loops.toArray();
  numloops = boundary_loops.numels();

  int numtris = tin->T.numels();
  double adist, mindist=DBL_MAX;

  gv=NULL;
  for (i=0; i<numloops; i++)
   for (j=0; j<numloops; j++)
	if (((Vertex *)bloops_array[i]->head()->data)->info != ((Vertex *)bloops_array[j]->head()->data)->info)
	{
	 adist = closestPair(bloops_array[i], bloops_array[j], &v, &w);
	 if (adist<mindist) {mindist=adist; gv=v; gw=w;}
	}

  if (gv!=NULL) tin->joinBoundaryLoops(gv, gw, 1, 0, 0);

  FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
  FOREACHVVVERTEX((&(tin->V)), v, n) v->info = NULL;

  free(bloops_array);
  while ((one_loop=(List *)boundary_loops.popHead())!=NULL) delete one_loop;

  return (gv!=NULL);
}



//#define DISCLAIMER

void usage()
{
 printf("\nMeshFix V1.0 - by Marco Attene\n------\n");
 printf("Usage: MeshFix meshfile [-a epsilon_angle] [-w] [-n]\n");
 printf("  Processes 'meshfile' and saves the result to 'meshfile_fixed.off'\n");
 printf("  By default, epsilon_angle is 0.\n  If specified, it must be in the range (0 - 2) degrees.\n");
 printf("  With '-w', the result is saved in VRML1.0 format instead of OFF.\n");
 printf("  With '-n', only the biggest input component is kept.\n");
 printf("  Accepted input formats are OFF, PLY and STL.\n  Other formats are supported only partially.\n");
 printf("  See http://jmeshlib.sourceforge.net for details on supported formats.\n");
 printf("\nIf MeshFix is used for research purposes, please cite the following paper:\n");
 printf("\n   M. Attene.\n   A lightweight approach to repairing digitized polygon meshes.\n   The Visual Computer, 2010. (c) Springer.\n");
 printf("\nHIT ENTER TO EXIT.\n");
 getchar();
 exit(0);
}

char *createFilename(const char *iname, const char *subext, const char *newextension)
{
 static char tname[2048];
 char *oname = (char *)malloc(strlen(iname)+strlen(subext)+strlen(newextension)+1);
 strcpy(tname, iname);
 for (int n=strlen(tname)-1; n>0; n--) if (tname[n]=='.') {tname[n] = '\0'; break;}
 sprintf(oname,"%s%s%s",tname,subext,newextension);
 return oname;
}

int main(int argc, char *argv[])
{
 char subext[128]="_fixed";
 JMesh::init();
 JMesh::app_name = "MeshFix";
 JMesh::app_version = "1.0";
 JMesh::app_year = "2010";
 JMesh::app_authors = "Marco Attene";
 JMesh::app_maillist = "attene@ge.imati.cnr.it";

 ExtTriMesh tin;

#ifdef DISCLAIMER
 printf("\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
 printf("This software can be used ONLY with an explicit authorization of the author.\n");
 printf("If you do not have such an authorization, you must delete this software.\n");
 printf("In no event this version of MeshFix can be redistributed.\n");
 printf("\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
#endif

 if (argc < 2) usage();

 bool keep_all_components = true;
 bool save_vrml = false;
 float par;
 for (int i=2; i<argc; i++)
 {
  if (i<argc-1) par = (float)atof(argv[i+1]); else par = 0;
  if      (!strcmp(argv[i], "-a"))
  {
   if (par < 0) JMesh::error("Epsilon angle must be > 0.\n");
   if (par > 2) JMesh::error("Epsilon angle must be < 2 degrees.\n");
   epsilon_angle = par;
   if (epsilon_angle)
   {
	JMesh::acos_tolerance = asin((M_PI*epsilon_angle)/180.0);
	printf("Fixing asin tolerance to %e\n",JMesh::acos_tolerance);
   }
  }
  else if (!strcmp(argv[i], "-n")) keep_all_components = false;
  else if (!strcmp(argv[i], "-w")) save_vrml = true;
  else if (argv[i][0] == '-') JMesh::warning("%s - Unknown operation.\n",argv[i]);

  if (par) i++;
 }

 // The loader performs the conversion to a set of oriented manifolds
 if (tin.load(argv[1]) != 0) JMesh::error("Can't open file.\n");
 input_filename = argv[1];

 if (keep_all_components)
 {
  printf("\nJoining input components ...\n");
  JMesh::begin_progress();
  while (joinClosestComponents(&tin)) JMesh::report_progress("Num. components: %d       ",tin.shells());
  JMesh::end_progress();
  tin.deselectTriangles();
 }

 // Keep only the biggest component
 int sc = tin.removeSmallestComponents();
 if (sc) JMesh::warning("Removed %d small components\n",sc);

 // Fill holes by taking into account both sampling density and normal field continuity
 tin.fillSmallBoundaries(tin.E.numels(), true, true);

 // Run geometry correction
 if (tin.boundaries() || !meshclean(tin))
 {
  fprintf(stderr,"MeshFix failed!\n");
  fprintf(stderr,"Please try manually using ReMESH v1.2 or later (http://remesh.sourceforge.net).\n");
  FILE *fp = fopen("meshfix_log.txt","a");
  fprintf(fp,"MeshFix failed on %s\n",input_filename);
  fclose(fp);
 }

 char *fname = createFilename(argv[1], subext, (save_vrml)?(".wrl"):(".off"));
 printf("Saving output mesh to '%s'\n",fname);
 if (save_vrml) tin.saveVRML1(fname); else tin.saveOFF(fname);

 return 0;
}

