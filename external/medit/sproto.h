/* animat.c */
int    loadNextMesh(pMesh ,int ,int );
int    animat();
int    playAnim(pScene ,pMesh ,int ,int );
void   glutIdle(void);
int animParticle(pScene sc,pMesh mesh);

/* camera.c */
double  Azimuth(pCamera );
double  Elevation(pCamera );
void    updateSun(pScene ,pCamera );
void    updateCamera(pScene ,pCamera ,double ,double );
pCamera initCamera(pScene ,int );

/* clip.c */
void   updateClip(pClip ,pMesh );
void   clipVertices(pMesh ,pScene ,pClip );
void   invertClip(pScene sc,pClip );
void   drawClip(pScene ,pClip ,pMesh ,GLboolean );
void   copyClip(pClip );
int    pasteClip(pClip );
void tiltClip(pScene sc,pClip clip);
void   resetClip(pScene ,pClip ,pMesh );
pClip  createClip(pScene ,pMesh );
void resetCube(pScene sc,pCube cube,pMesh mesh);
void dumpCube(pScene sc,pMesh mesh,pCube cube);

void   updateCube(pCube ,pMesh );
pCube  createCube(pScene ,pMesh );

/* clipvol.c */
GLuint capTetra(pMesh );
GLuint capTetraMap(pMesh );
GLuint capTetraIso(pMesh );

/* critip.c */
GLuint listCritPoint(pScene ,pMesh );

/* dlists.c */
GLuint  listTria(pScene ,pMesh );
GLuint  listQuad(pScene ,pMesh );
GLuint  listTetra(pScene ,pMesh ,ubyte );
GLuint  listHexa(pScene ,pMesh ,ubyte );

/* geometry.c */
GLuint geomList(pScene ,pMesh );

/* gisfil.c */
int loadGIS(pMesh );

/* hash.c */
int hashTria(pMesh );
int hashTetra(pMesh );
int hashHexa(pMesh mesh);

/* image.c */
PPMimage *loadPPM(const char *imgname,int *type);
int  savePPM(const char *imgname,pPPMimage img,int typimg);
void saveEPS(pScene sc,const char *imgname,pPPMimage pixels);
int  imgHard(pScene sc,char *data,char key);

/* input/output */
int    EatLine(FILE  *in);
int    inmsh2(pMesh mesh);
int    bbfile(pMesh );
int    loadMesh(pMesh );
int    saveMesh(pScene ,pMesh ,char *,ubyte );
int    loadSol(pMesh mesh,char *filename,int numsol);
int    loadMesh_popen(pMesh); //pour popen

/* ilists.c */
GLuint listTriaIso(pScene ,pMesh );
GLuint listQuadIso(pScene ,pMesh );
GLuint listTetraIso(pScene ,pMesh );
int tetraIsoPOVray(pScene sc,pMesh mesh);

/* items.c */
void   drawAxis(pScene ,int );
void   drawBox(pScene ,pMesh ,int );
void   drawCube(pScene ,pMesh );
void   drawGrid(pScene ,pMesh );
void   rubberBand(pPersp );
void   drawBase(pScene ,pMesh );
void   drawIso(pScene );

/* keyboard.c */
void   specCamera(pScene ,int );
void   special(int ,int ,int );
void   keyScene(unsigned char ,int ,int );

/* listnum.c */
void   updatePoints(pScene sc,pMesh mesh,int refmat);
void   listNum(pScene sc,pMesh mesh);

/* material.c */
void   matInit(pScene );
void   matSort(pScene );
int    matRef(pScene ,int );
void   matReshape(int ,int );
void   matsubReshape(int ,int );
void   matsubDisplay();
void   matDisplay();
void   matMouse(int ,int ,int ,int );
void   matKeyboard(unsigned char ,int ,int );
void   matEdit(pScene );

/* medit.c */
int    medit0();
int    medit1();

/* menus.c */
void   doLists(pScene ,pMesh );
void   doMapLists(pScene ,pMesh ,int );
void   doIsoLists(pScene ,pMesh ,int );
void   keyFile(unsigned char ,int ,int );
void   menuFile(int );
void   keyItem(unsigned char ,int ,int );
void   menuItem(int );
void   keyAnim(unsigned char ,int ,int );
void   menuAnim(int );
void   keyTrajet(unsigned char ,int ,int );
void   menuTrajet(int );
void   keyMode(unsigned char ,int ,int );
void   menuMode(int );
void   menuScene(int ); 
void   keyView(unsigned char ,int ,int );
void   menuView(int );
void   keyColor(unsigned char ,int ,int );
void   menuColor(int );
void   keyClip(unsigned char ,int ,int );
void   menuClip(int );
void keyCube(unsigned char key,int x,int y);
void   keyFeature(unsigned char ,int ,int );
void   menuFeature(int );
void   menuImage(int );
void   keyMetric(unsigned char key,int x,int y);
int    createMenus(pScene ,pMesh );
 
/* mesh.c */
void   meshInfo(pMesh );
int    meshSurf(pMesh );
void   meshCoord(pMesh ,int );
void   meshBox(pMesh mesh,int bb);
void   meshRef(pScene sc,pMesh mesh);
int    meshUpdate(pScene sc,pMesh mesh);

/* mlists.c */
GLuint  listTriaMap(pScene ,pMesh );
GLuint  listQuadMap(pScene ,pMesh );
GLuint  listTetraMap(pScene ,pMesh ,ubyte );
GLuint  listHexaMap(pScene ,pMesh ,ubyte );
void    cutTriangle(pScene ,triangle );
GLuint  alt2dList(pScene ,pMesh ,int ,float ,float );
void    setupPalette(pScene ,pMesh );
GLuint  drawPalette(pScene );

/* morphing */
int    morphMesh(pScene sc,pMesh mesh1);
int    modeMorphing();
 
/* mouse.c */
void   mouse(int button,int state,int x,int y);
void   passive_motion(int x,int y);
void   motion(int x,int y);
void   redrawOverlay(int stretchX,int stretchY);
void   motionCamera(int x,int y);
void   mouseCamera(int button,int state,int x,int y);
void   animateCamera();

/* normal.c */
GLuint drawNormals(pMesh mesh,pScene sc);

/* outmsh.c */
int    outmsh(pScene ,pMesh ,char *name,ubyte clipon);

/* parsar.c */
int    parsar(int argc,char *argv[]);

/* parsop.c */
int    saveMeditFile(char *,pScene ); 
void   iniopt(pScene ,pMesh );
int    parsop(pScene ,pMesh );

/* particle.c */
int createParticle(pScene ,pMesh );
#ifdef IGL
int advectParticle(pScene sc,pMesh mesh);
int displayParticle(pScene sc,pMesh mesh);
#endif


/* path.c */
int    pathAdd(pScene ,int, int);
GLuint pathList(pScene );
int    pathLoad(char *data,pScene );
int    pathSave(char *file,pScene );
void   pathFollow(pScene );

/* persp.c */
void   setPersp(pScene ,pPersp ,int );
pPersp initPersp(pPersp ,float );

/* picking.c */
GLuint pickingList(pScene ,int ,int );
GLuint pickingPoint(pScene sc,int x,int y);
GLuint pickItem(pMesh ,pScene ,int );
GLuint pickingScene(pScene sc,int x,int y,int ident);
 
/* prierr.c */
void   prierr(int typerr,int indice);

/* psfile.c */
void writeEPSheader(FILE *,char *,char ,int ,int ,float ,float);
void writeEPStrailer(FILE *);
void writeEPSRow(FILE *,char ,ubyte *,int ,ubyte );

/* scene.c */
int    currentScene();
void   checkErrors(void);
void   oglerr(GLenum error);
void   farclip(GLboolean );
void   reshapeScene(int width,int height);
void   setupView(pScene sc);
void   drawBackTex(pScene sc);
void   drawModel(pScene sc);  
void   drawScene(pScene );
void   redrawScene();
void   deleteScene(pScene sc);
void   initGrafix(pScene sc,pMesh mesh);
int    createScene(pScene sc,int idmesh);
void   streamIdle();

/* scissor.c */
void   scissorScene();

/* status.c */
void reshapeStatusBar(pScene sc,int width,int height);
void redrawStatusBar(pScene sc);
void mouseStatus(int button,int state,int x,int y);

/* stream.c */
double  sizeTetra(pMesh ,int );
double  sizeHexa(pMesh ,int );
double  sizeTria(pMesh ,int );
int     locateTria(pMesh mesh,int nsdep,int base,float *p,double *cb);
int     locateTetra(pMesh mesh,int nsdep,int base,float *p,double *cb);
int     inTria(pMesh ,int ,float *,double *);
int     listTetraStream(pScene ,pMesh ,float *,int );
int     listHexaStream(pScene ,pMesh ,float *,int );
int     listTriaStream(pScene ,pMesh ,float *);
pStream createStream(pScene ,pMesh );
int     streamRefTria(pScene sc,pMesh mesh);
int     streamRefQuad(pScene sc,pMesh mesh);
int     streamRefPoint(pScene sc,pMesh mesh);
int     listSaddleStream(pScene sc,pMesh mesh,int depart,
                         float *pp,float *vv,double lambda);
int     nxtPoint3D(pMesh mesh,int nsdep,float *p,float step,double *v);
int     nxtPoint2D(pMesh mesh,int nsdep,float *p,float step,double *v);
double  field2DInterp(pMesh mesh,int iel,double *cb,double *v);
double  vector3DInterp(pMesh mesh,pPoint pt[4],double *cb,double *v);
double  field3DInterp(pMesh mesh,int iel,double *cb,double *v);
double  sizeTria(pMesh mesh,int k);
double  sizeQuad(pMesh mesh,int k);
double  sizeTetra(pMesh mesh,int k);
int streamIsoPoint(pScene sc,pMesh mesh);
int inTetra(pMesh mesh,int nsdep,float *p,double *cb);


/* tensor.c */
GLuint listPointVector(pMesh ,ubyte );

/* texture.c */
pPPMimage texDistortion(pPPMimage );

/* tiles.c */
int   imgTiling(pScene sc,char *data,char key);

/* transform.c */
void   resetTransform(pTransform );
pTransform createTransform();

/* util.c */
void   setFont(char* name,int size);
void   drwstr(GLuint x,GLuint y,char* format, ...);
void   output2(GLfloat x,GLfloat y,char *format,...); 
void   output3(GLfloat x,GLfloat y,GLfloat z,char *format,...);
#ifdef IGL
template <typename T>
void hsvrgb(T *hsv,double *rgb);
#else
void   hsvrgb(double *hsv,double *rgb);
#endif
void   transformPoint(double u[4],float v[4],float m[16]);
void transformPoint2(double u[4],float v[4],float m[16]) ;
void   transformPointd(double u[4],double v[4],double m[16]);
void   transformVector(float u[4],float v[4],float m[16]); 
void   multMatrix(GLfloat *p,GLfloat *a,GLfloat *b);
void   rotateMatrix(GLfloat angle,GLfloat x,GLfloat y,GLfloat z,GLfloat rm[16]);
int    invertMatrix(float src[16],float inverse[16]);
int    filnum(char *data,int numdep,char *ext);

/* vector.c */
void   drawVector2D(float p[2],double u[2],double scal);
void   drawVector3D(float p[3],double u[3],double scal);
GLuint listTria2dVector(pMesh mesh);
GLuint listTria3dVector(pMesh mesh);
GLuint listClipTetraVector(pMesh mesh);
GLuint listClipHexaVector(pMesh mesh);
GLuint listQuad2dVector(pMesh mesh);

/* view.c */
void   copyView(pTransform view,pCamera cam,pPersp persp);
int    pasteView(pTransform view,pCamera cam,pPersp persp);
int    linkView(pScene sc1);
void   unlinkView(pScene sc1);

/* zaldy.c */
int    zaldy1(pMesh mesh);
int    zaldy2(pMesh mesh);

#ifdef IGL
/* cenrad.c */
int cenrad(pMesh mesh,int iel,double *c,double *rad);
/* param.c */
void parEdit(pScene sc);
/* inout_popenbinaire.c */
int loadMesh_popen_bin(pMesh mesh);
int loadSol_popen_bin(pMesh mesh,char *filename,int numsol);
/* inout_morice.c */
int loadSol_popen(pMesh mesh,char *filename,int numsol);
/* sftcpy.c */
int sftcpy(pScene sc,pMesh mesh);
/* ellipse.c */
void drawEllipsoid(pScene sc,pMesh mesh,int typel,int k);
void drawEllipse(pScene sc,pMesh mesh,int typel,int k);
GLuint drawAllEllipse(pScene sc,pMesh mesh);
void circumSphere(pScene sc,pMesh mesh,int typel,int k);
/* bbfile.c */
int EatSpace(FILE  *in);
#endif
