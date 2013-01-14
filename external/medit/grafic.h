#ifndef _GRAFIC_H
#define _GRAFIC_H

#ifdef IGL
#  include "IGLParams.h"
#  include <igl/ReAntTweakBar.h>
#endif

#define MAX_LIST  4
#define MAXISO    5

#define C_ON     (1 << 0)
#define C_EDIT   (1 << 1)   
#define C_VOL    (1 << 2)
#define C_UPDATE (1 << 3)
#define C_FREEZE (1 << 4)
#define C_CAP    (1 << 5)
#define C_REDO   (1 << 6)
#define C_HIDE   (1 << 7)

/* view modes */
#define S_BDRY     (1<<0)
#define S_FILL     (1<<1)
#define S_COLOR    (1<<2)
#define S_MAP      (1<<3)
#define S_MATERIAL (1<<4)
#define S_ALTITUDE (1<<5)
#define S_DISPL    (1<<6)

#define MONO     0
#define LEFT     1
#define RIGHT    2

enum {WIRE   = S_BDRY,
      HIDDEN = S_BDRY + S_FILL,
      DEPTH  = S_BDRY + S_COLOR,
      FILL   = S_FILL + S_COLOR,
      SHADED = S_BDRY + S_FILL + S_COLOR,
      SIZEMAP= S_BDRY + S_FILL + S_MAP
     };
enum {LTria, LQuad, LTets, LHexa, LEdges, LPoint};
enum {PERSPECTIVE, CAMERA, ORTHO};
enum {X_AXIS=0, Y_AXIS, Z_AXIS};
enum {VECTOR,CONE};

/* items */
#define S_AXIS    (1<<0)
#define S_BOX     (1<<1)
#define S_GRID    (1<<2)
#define S_GEOM    (1<<3)
#define S_ISO     (1<<4)
#define S_PALETTE (1<<5)
#define S_NUMP    (1<<6)
#define S_NUMF    (1<<7)

/* type */
#define S_FLAT     (1<<0)     /* render with facet normals  */
#define S_SCISSOR  (1<<1)     /* scissoring mode            */
#define S_FOLLOW   (1<<2)     /* picking mode               */
#define S_NORMAL   (1<<3)
#define S_OPPOS    (1<<4)
#define S_DECO     (1<<5)
#define S_PATH     (1<<6)
#define S_RESET    (1<<7)

/* iso-values */
#define MAX_ISO    5
#define S_ISOLINE  (1<<1)
#define S_ISOSURF  (1<<2)
#define S_STREAML  (1<<3)
#define S_STREAMR  (1<<4)
#define S_VECTOR   (1<<5)
#define S_CRITP    (1<<6)
#define S_PARTICLE (1<<7)

    
typedef struct sperspective {
  float      fovy,depth;
  float      matrix[16],alpha,gamma;
  int        rubix,rubfx,rubiy,rubfy;
  ubyte      pmode,rubber;
} Persp;
typedef Persp * pPersp;

typedef struct triangle {
  float a[3],b[3],c[3];
  float va,vb,vc;
  float na[3],nb[3],nc[3];
} triangle;

typedef struct material {
  float   amb[4],emi[4],dif[4],spe[4],shininess;
  float   ext[6];
  GLint   list;
  int     depmat[MAX_LIST];
  int     ref,sort;
  char    name[128];
  ubyte   flag;
} Material ;
typedef Material * pMaterial;

typedef struct _cell {
  int   id;
  int   x,y;
  float min,max;
  float value;
  float step;
  char* info;
  char* format;
} cell;

typedef struct transform {
  float    pos[3];                /* current mouse position */
  float    angle,axis[3];         /* rotation angle + axis  */ // Alec: really this is delta-angle+axis
  float    panx,pany,opanx,opany; /* screen translation     */ // Alec: also this is delta
  float    matrix[16],oldmat[16]; /* transformation matrix  */
  float    rot[16],tra[16];
  int      mstate,mbutton,manim;
} Transform;
typedef Transform * pTransform;

typedef struct cube {
  pTransform   cubetr;
  float        cmi[3],cma[3];
  ubyte        active;
} Cube;
typedef Cube *pCube;

typedef struct clip {
  pTransform   cliptr;
  double       eqn[4];
  ubyte        active;
} Clip;
typedef Clip *pClip;


typedef struct camera {
  GLfloat  eye[3];                /* Position of the camera */
  GLfloat  speed[3],spmod,altinc;
  int      vecup;
} Camera;
typedef Camera * pCamera;

/* scene parameters */
typedef struct sparam {
  float     back[4],line[4],edge[4],sunpos[4],clip[6];
  float     cm,dpi,coeff,cumtim,cumpertime,maxtime,pertime,dt;
  float     eyesep,linewidth,pointsize;
  short     xi,yi,xs,ys,pxs,pys;
  int       nbmat;
  char      pscolor[10];
  ubyte     sunp,linc,advtim,nbpart;
} Param;

/* trajectoire */
typedef struct straj {
  int      np;
  float   *pt;
  float   *tg;
  float    sec;
  GLuint   tlist;
} Trajet;

/* streamlines */
typedef struct sstream {
  double   size,norm;
  float    xmin,xmax,ymin,ymax,zmin,zmax;
  float   *listp;
  float    stpt[4][3],stcol[4];
  int      stnp,stiso[4];
  short    nbstl;
  ubyte    typtrack;
} Stream;
typedef Stream * pStream;

typedef struct strgrd {
  pTransform   strtr;
  GLuint       grid;
  ubyte        active;
} Strgrd;
typedef Strgrd * pStrgrd;

typedef struct siso {
  float  val[MAXISO+2];
  // Alec: ranges from 0 to 240
  float  col[MAXISO+2];
  ubyte  palette,ptyp;
} Iso;

typedef struct scene {
  pTransform view;
  pClip      clip;
  pCube      cube;
  pPersp     persp;
  pCamera    camera;
  pMaterial  material;
  pStream    stream;
  /*pStrgrd    stg;*/
  Param      par;
  Trajet     path;
  Iso        iso;

  float      dmin,dmax;       /* scene size    */
  float      shrink;          /* shrink value  */
  float      cx,cy,cz;        /* center of scene */
  
  GLuint     dlist[MAX_LIST];    /* display lists  */
  GLuint     mlist[MAX_LIST];    /* metric lists   */
  GLuint     ilist[MAX_LIST];    /* iso-surfaces   */
  GLuint     clist[MAX_LIST];
  GLuint     cmlist[MAX_LIST];   /* clipped elts   */
  GLuint     vlist[MAX_LIST];    /* vector list    */
  GLuint    *slist,cplist;       /* streamlines    */
  GLuint     glist,nlist;        /* geometry lists */
  GLuint     grid;
  GLuint     picklist;

  int       *matsort;
  short      idwin,idmesh;    /* window, mesh id */
  short      master,slave;

  ubyte      item;            /* display items */
  ubyte      mode;            /* render mode   */
  ubyte      type;
  ubyte      isotyp;
  ubyte      picked;
#ifdef IGL
  igl::ReTwBar rebar;
  // Pointer so recompiling is easier
  IGLParams * igl_params; 
#endif
} Scene;
typedef Scene * pScene;


#endif
