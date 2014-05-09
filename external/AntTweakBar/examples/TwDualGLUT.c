//  ---------------------------------------------------------------------------
//
//  @file       TwDualGLUT.c
//  @brief      This example illustrates the use of AntTweakBar with multiple
//              windows. Here we are using GLUT to create a main window having
//              two sub-windows, each sub-window has a tweak bar.
//              This example extends the TwSimpleGLUT example.
//
//              AntTweakBar: http://anttweakbar.sourceforge.net/doc
//              OpenGL:      http://www.opengl.org
//              GLUT:        http://opengl.org/resources/libraries/glut
//  
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
//  MiniGLUT.h is provided to avoid the need of having GLUT installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual GLUT library SDK.
#   define USE_MINI_GLUT
#endif

#if defined(USE_MINI_GLUT)
#   include "../src/MiniGLUT.h"
#elif defined(_MACOSX)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif


// This example displays one of the following shapes in each sub-window
typedef enum { SHAPE_TEAPOT=1, SHAPE_TORUS, SHAPE_CONE } Shape;
#define NUM_SHAPES 3

typedef struct
{
	int   WinID;
	TwBar *Bar;
	Shape ObjectShape;
	float Zoom;
	float Rotation[4];
	int   AutoRotate;
	int   RotateTime;
	float RotateStart[4];
	float MatAmbient[4];
	float MatDiffuse[4];
	float LightMultiplier;
	float LightDirection[3];
} SubWindowData;

SubWindowData g_SubWindowData[2];


// Routine to set a quaternion from a rotation axis and angle
// ( input axis = float[3] angle = float  output: quat = float[4] )
void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
{
    float sina2, norm;
    sina2 = (float)sin(0.5f * angle);
    norm = (float)sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    quat[0] = sina2 * axis[0] / norm;
    quat[1] = sina2 * axis[1] / norm;
    quat[2] = sina2 * axis[2] / norm;
    quat[3] = (float)cos(0.5f * angle);
}


// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void ConvertQuaternionToMatrix(const float *quat, float *mat)
{
    float yy2 = 2.0f * quat[1] * quat[1];
    float xy2 = 2.0f * quat[0] * quat[1];
    float xz2 = 2.0f * quat[0] * quat[2];
    float yz2 = 2.0f * quat[1] * quat[2];
    float zz2 = 2.0f * quat[2] * quat[2];
    float wz2 = 2.0f * quat[3] * quat[2];
    float wy2 = 2.0f * quat[3] * quat[1];
    float wx2 = 2.0f * quat[3] * quat[0];
    float xx2 = 2.0f * quat[0] * quat[0];
    mat[0*4+0] = - yy2 - zz2 + 1.0f;
    mat[0*4+1] = xy2 + wz2;
    mat[0*4+2] = xz2 - wy2;
    mat[0*4+3] = 0;
    mat[1*4+0] = xy2 - wz2;
    mat[1*4+1] = - xx2 - zz2 + 1.0f;
    mat[1*4+2] = yz2 + wx2;
    mat[1*4+3] = 0;
    mat[2*4+0] = xz2 + wy2;
    mat[2*4+1] = yz2 - wx2;
    mat[2*4+2] = - xx2 - yy2 + 1.0f;
    mat[2*4+3] = 0;
    mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
    mat[3*4+3] = 1;
}


// Routine to multiply 2 quaternions (ie, compose rotations)
// ( input q1 = float[4] q2 = float[4]  output: qout = float[4] )
void MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
    float qr[4];
	qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}


// Return elapsed time in milliseconds
int GetTimeMs()
{
#if !defined(_WIN32)
    return glutGet(GLUT_ELAPSED_TIME);
#else
    // glutGet(GLUT_ELAPSED_TIME) seems buggy on Windows
    return (int)GetTickCount(); 
#endif
}


// Find the current WindowData
SubWindowData *GetCurrentSubWindowData()
{
    int i, currWinID;
    currWinID = glutGetWindow();
    for (i = 0; i < 2; i++)
        if (g_SubWindowData[i].WinID == currWinID)
            return &g_SubWindowData[i];
    return NULL;
}


// Callback function called by GLUT to render the main window content
void DisplayMainWindow(void)
{
    // Clear frame buffer
    glClearColor(1, 0.8f, 0.4f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Present frame buffer
    glutSwapBuffers();
}


// Callback function called by GLUT when the main window size has changed
void ReshapeMainWindow(int width, int height)
{
    if (width > 32 && height > 32) 
    {
	    glutSetWindow(g_SubWindowData[0].WinID);
	    glutPositionWindow(8, 8);
	    glutReshapeWindow(width/2-16, height-16);

	    glutSetWindow(g_SubWindowData[1].WinID);
	    glutPositionWindow(width/2+8, 8);
	    glutReshapeWindow(width/2-16, height-16);
    }
}


// Callback function called by GLUT to render sub-window content
void DisplaySubWindow(void)
{
    float v[4]; // will be used to set light parameters
    float mat[4*4]; // rotation matrix
    SubWindowData *win;

    win = GetCurrentSubWindowData();
    if (win == NULL) return;

    // Tell AntTweakBar which is the current window
    TwSetCurrentWindow(win->WinID);

    // Clear frame buffer
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    // Set light
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    v[0] = v[1] = v[2] = win->LightMultiplier*0.4f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
    v[0] = v[1] = v[2] = win->LightMultiplier*0.8f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
    v[0] = -win->LightDirection[0]; v[1] = -win->LightDirection[1]; v[2] = -win->LightDirection[2]; v[3] = 0.0f;
    glLightfv(GL_LIGHT0, GL_POSITION, v);

    // Set material
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, win->MatAmbient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, win->MatDiffuse);

    // Rotate and draw shape
    glPushMatrix();
    glTranslatef(0.5f, -0.3f, 0.0f);
    if( win->AutoRotate ) 
    {
        float axis[3] = { 0, 1, 0 };
        float angle = (float)(GetTimeMs()-win->RotateTime)/1000.0f;
        float quat[4];
        SetQuaternionFromAxisAngle(axis, angle, quat);
        MultiplyQuaternions(win->RotateStart, quat, win->Rotation);
    }
    ConvertQuaternionToMatrix(win->Rotation, mat);
    glMultMatrixf(mat);
    glScalef(win->Zoom, win->Zoom, win->Zoom);
    glCallList(win->ObjectShape);
    glPopMatrix();

    // Draw tweak bars
    TwDraw();

    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();
}


// Callback function called by GLUT when sub-window size has changed
void ReshapeSubWindow(int width, int height)
{
	SubWindowData *win;
	
    win = GetCurrentSubWindowData();
    if (win == NULL) return;

    // Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, (double)width/height, 1, 10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5, 0,0,0, 0,1,0);
    glTranslatef(0, 0.6f, -1);

    // Send the new window size to AntTweakBar
    TwSetCurrentWindow(win->WinID);
    TwWindowSize(width, height);
}


// Function called at exit
void Terminate(void)
{ 
    int i;
	for (i=0; i<2; i++) {
		glutSetWindow(g_SubWindowData[i].WinID);	
    	glDeleteLists(SHAPE_TEAPOT, NUM_SHAPES);
	}
    TwTerminate();
}


//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetAutoRotateCB(const void *value, void *clientData)
{
	SubWindowData *win;
	
	win = (SubWindowData *)clientData;
    win->AutoRotate = *(const int *)value; // copy value to win->AutoRotate

    if (win->AutoRotate != 0) 
    {
        // init rotation
        win->RotateTime = GetTimeMs();
        win->RotateStart[0] = win->Rotation[0];
        win->RotateStart[1] = win->Rotation[1];
        win->RotateStart[2] = win->Rotation[2];
        win->RotateStart[3] = win->Rotation[3];
    }

    // make Rotation variable read-only or read-write
    TwSetCurrentWindow(win->WinID);
    TwSetParam(win->Bar, "ObjRotation", "readonly", TW_PARAM_INT32, 1, &win->AutoRotate);
}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetAutoRotateCB(void *value, void *clientData)
{
	SubWindowData *win;
	
	win = (SubWindowData *)clientData;
    *(int *)value = win->AutoRotate; // copy win->AutoRotate to value
}


// Mouse Button event callbacks
int MouseButtonCB(int glutButton, int glutState, int mouseX, int mouseY) 
{
	TwSetCurrentWindow(glutGetWindow());
	return TwEventMouseButtonGLUT(glutButton,glutState,mouseX,mouseY);
}


// Mouse Motion event callbacks
int MouseMotionCB(int mouseX, int mouseY) 
{
	TwSetCurrentWindow(glutGetWindow());
	return TwEventMouseMotionGLUT(mouseX,mouseY);
}


// Keyboard event callbacks
int KeyboardCB(unsigned char glutKey, int mouseX, int mouseY) 
{
	TwSetCurrentWindow(glutGetWindow());
	return TwEventKeyboardGLUT(glutKey,mouseX,mouseY);	
}


// Special key event callbacks
int SpecialKeyCB(int glutKey, int mouseX, int mouseY) 
{
	TwSetCurrentWindow(glutGetWindow());
	return TwEventSpecialGLUT(glutKey,mouseX,mouseY);	
}


// Setup new sub-window
void SetupSubWindow(int subWinIdx) 
{
    float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
    float angle = 0.8f;
    SubWindowData *win;
    
    win = &g_SubWindowData[subWinIdx];
    win->ObjectShape = (subWinIdx == 0) ? SHAPE_TEAPOT : SHAPE_TORUS;
	win->Zoom = 1;
	win->AutoRotate = (subWinIdx == 0);
    win->MatAmbient[0] = (subWinIdx == 1) ? 0.0f : 0.5f;; win->MatAmbient[1] = win->MatAmbient[2] = 0.2f; win->MatAmbient[3] = 1;
    win->MatDiffuse[0] = (subWinIdx == 1) ? 0.0f : 1.0f; win->MatDiffuse[1] = 1; win->MatDiffuse[2] = 0; win->MatDiffuse[3] = 1;
	win->LightMultiplier = 1;
    win->LightDirection[0] = win->LightDirection[1] = win->LightDirection[2] = -0.57735f;
    win->RotateTime = GetTimeMs();
    SetQuaternionFromAxisAngle(axis, angle, win->Rotation);
    SetQuaternionFromAxisAngle(axis, angle, win->RotateStart);
	
	glutSetWindow(win->WinID);
    // Set GLUT callbacks
    glutDisplayFunc(DisplaySubWindow);
    glutReshapeFunc(ReshapeSubWindow);
    // Set GLUT event callbacks
    // - Register mouse button events callback
    glutMouseFunc((GLUTmousebuttonfun)MouseButtonCB);
    // - Register mouse motion events callback
    glutMotionFunc((GLUTmousemotionfun)MouseMotionCB);
    // - Register mouse "passive" motion events (same as Motion)
    glutPassiveMotionFunc((GLUTmousemotionfun)MouseMotionCB);
    // - Register keyboard events callback
    glutKeyboardFunc((GLUTkeyboardfun)KeyboardCB);
    // - Register special key events callback
    glutSpecialFunc((GLUTspecialfun)SpecialKeyCB);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);

    // Create some 3D objects (stored in display lists)
    glNewList(SHAPE_TEAPOT, GL_COMPILE);
    glutSolidTeapot(1.0);
    glEndList();
    glNewList(SHAPE_TORUS, GL_COMPILE);
    glutSolidTorus(0.3, 1.0, 16, 32);
    glEndList();
    glNewList(SHAPE_CONE, GL_COMPILE);
    glutSolidCone(1.0, 1.5, 64, 4);
    glEndList();

    // Declare this window as current for AntTweakBar.
    // Here, this will create a new AntTweakBar context for this window,
    // which will be identified by the number WinID.
    TwSetCurrentWindow(win->WinID);

    // Create a tweak bar
    win->Bar = TwNewBar("TweakBar");
    TwDefine(" GLOBAL help='This example shows how to use AntTweakBar with multiple windows.' "); // Message added to the help bar.
    TwDefine(" TweakBar size='200 400' color='96 216 224' "); // change default tweak bar size and color

    // Add 'win->Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
    TwAddVarRW(win->Bar, "Zoom", TW_TYPE_FLOAT, &win->Zoom, 
               " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

    // Add 'win->Rotation' to 'bar': this is a variable of type TW_TYPE_QUAT4F which defines the object's orientation
    TwAddVarRW(win->Bar, "ObjRotation", TW_TYPE_QUAT4F, &win->Rotation, 
               " label='Object rotation' open help='Change the object orientation.' ");

    // Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(win->Bar, "AutoRotate", TW_TYPE_BOOL32, SetAutoRotateCB, GetAutoRotateCB, win, 
               " label='Auto-rotate' key=space help='Toggle auto-rotate mode.' ");

    // Add 'win->LightMultiplier' to 'bar': this is a variable of type TW_TYPE_FLOAT. Its key shortcuts are [+] and [-].
    TwAddVarRW(win->Bar, "Multiplier", TW_TYPE_FLOAT, &win->LightMultiplier, 
               " label='Light booster' min=0.1 max=4 step=0.02 keyIncr='+' keyDecr='-' help='Increase/decrease the light power.' ");

    // Add 'win->LightDirection' to 'bar': this is a variable of type TW_TYPE_DIR3F which defines the light direction
    TwAddVarRW(win->Bar, "LightDir", TW_TYPE_DIR3F, &win->LightDirection, 
               " label='Light direction' open help='Change the light direction.' ");

    // Add 'win->MatAmbient' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
    // and is inserted into a group named 'Material'.
    TwAddVarRW(win->Bar, "Ambient", TW_TYPE_COLOR3F, &win->MatAmbient, " group='Material' ");

    // Add 'win->MatDiffuse' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
    // and is inserted into group 'Material'.
    TwAddVarRW(win->Bar, "Diffuse", TW_TYPE_COLOR3F, &win->MatDiffuse, " group='Material' ");

    // Add the enum variable 'win->ObjectShape' to 'bar'
    // (before adding an enum variable, its enum type must be declared to AntTweakBar as follow)
    {
        // ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
        TwEnumVal shapeEV[NUM_SHAPES] = { {SHAPE_TEAPOT, "Teapot"}, {SHAPE_TORUS, "Torus"}, {SHAPE_CONE, "Cone"} };
        // Create a type for the enum shapeEV
        TwType shapeType = TwDefineEnum("ShapeType", shapeEV, NUM_SHAPES);
        // add 'win->CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
        TwAddVarRW(win->Bar, "Shape", shapeType, &win->ObjectShape, " keyIncr='<' keyDecr='>' help='Change object shape.' ");
    }
}


// Main
int main(int argc, char *argv[])
{	
    int mainWinID;

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(960, 480);
    mainWinID = glutCreateWindow("AntTweakBar multi-window example");
    glutCreateMenu(NULL);
	glutDisplayFunc(DisplayMainWindow);
	glutReshapeFunc(ReshapeMainWindow);

    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);

	// Create two sub-windows
    g_SubWindowData[0].WinID = glutCreateSubWindow(mainWinID, 8, 8, 480, 464);
	SetupSubWindow(0);
    g_SubWindowData[1].WinID = glutCreateSubWindow(mainWinID, 488, 8, 464, 464);
    SetupSubWindow(1);
    
    atexit(Terminate);  // Called after glutMainLoop ends

    // Call the GLUT main loop
    glutMainLoop();

    return 0;
}

