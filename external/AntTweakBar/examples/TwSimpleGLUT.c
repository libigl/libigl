//  ---------------------------------------------------------------------------
//
//  @file       TwSimpleGLUT.c
//  @brief      A simple example that uses AntTweakBar with OpenGL and GLUT.
//
//              AntTweakBar: http://www.antisphere.com/Wiki/tools:anttweakbar
//              OpenGL:      http://www.opengl.org
//              GLUT:        http://opengl.org/resources/libraries/glut
//  
//  @author     Philippe Decaudin - http://www.antisphere.com
//  @date       2006/05/20
//
//  Compilation:
//  http://www.antisphere.com/Wiki/tools:anttweakbar:examples#twsimpleglut
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

// This example displays one of the following shapes
typedef enum { SHAPE_TEAPOT=1, SHAPE_TORUS, SHAPE_CONE } Shape;
#define NUM_SHAPES 3
Shape g_CurrentShape = SHAPE_TORUS;
// Shapes scale
float g_Zoom = 1.0f;
// Shape orientation (stored as a quaternion)
float g_Rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };
// Auto rotate
int g_AutoRotate = 0;
int g_RotateTime = 0;
float g_RotateStart[] = { 0.0f, 0.0f, 0.0f, 1.0f };
// Shapes material
float g_MatAmbient[] = { 0.5f, 0.0f, 0.0f, 1.0f };
float g_MatDiffuse[] = { 1.0f, 1.0f, 0.0f, 1.0f };
// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { -0.57735f, -0.57735f, -0.57735f };


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


// Callback function called by GLUT to render screen
void Display(void)
{
    float v[4]; // will be used to set light parameters
    float mat[4*4]; // rotation matrix

    // Clear frame buffer
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    // Set light
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.4f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.8f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
    v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
    glLightfv(GL_LIGHT0, GL_POSITION, v);

    // Set material
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, g_MatAmbient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, g_MatDiffuse);

    // Rotate and draw shape
    glPushMatrix();
    glTranslatef(0.5f, -0.3f, 0.0f);
    if( g_AutoRotate ) 
    {
        float axis[3] = { 0, 1, 0 };
        float angle = (float)(GetTimeMs()-g_RotateTime)/1000.0f;
        float quat[4];
        SetQuaternionFromAxisAngle(axis, angle, quat);
        MultiplyQuaternions(g_RotateStart, quat, g_Rotation);
    }
    ConvertQuaternionToMatrix(g_Rotation, mat);
    glMultMatrixf(mat);
    glScalef(g_Zoom, g_Zoom, g_Zoom);
    glCallList(g_CurrentShape);
    glPopMatrix();

    // Draw tweak bars
    TwDraw();

    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();
}


// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
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
    TwWindowSize(width, height);
}


// Function called at exit
void Terminate(void)
{ 
    glDeleteLists(SHAPE_TEAPOT, NUM_SHAPES);

    TwTerminate();
}


//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetAutoRotateCB(const void *value, void *clientData)
{
    (void)clientData; // unused

    g_AutoRotate = *(const int *)value; // copy value to g_AutoRotate
    if( g_AutoRotate!=0 ) 
    {
        // init rotation
        g_RotateTime = GetTimeMs();
        g_RotateStart[0] = g_Rotation[0];
        g_RotateStart[1] = g_Rotation[1];
        g_RotateStart[2] = g_Rotation[2];
        g_RotateStart[3] = g_Rotation[3];

        // make Rotation variable read-only
        TwDefine(" TweakBar/ObjRotation readonly ");
    }
    else
        // make Rotation variable read-write
        TwDefine(" TweakBar/ObjRotation readwrite ");
}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetAutoRotateCB(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = g_AutoRotate; // copy g_AutoRotate to value
}

 
// Main
int main(int argc, char *argv[])
{
    TwBar *bar; // Pointer to the tweak bar
    float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
    float angle = 0.8f;

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutCreateWindow("AntTweakBar simple example using GLUT");
    glutCreateMenu(NULL);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    atexit(Terminate);  // Called after glutMainLoop ends

    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);

    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
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

    // Create a tweak bar
    bar = TwNewBar("TweakBar");
    TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLUT and OpenGL.' "); // Message added to the help bar.
    TwDefine(" TweakBar size='200 400' color='96 216 224' "); // change default tweak bar size and color

    // Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
    TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &g_Zoom, 
               " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

    // Add 'g_Rotation' to 'bar': this is a variable of type TW_TYPE_QUAT4F which defines the object's orientation
    TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &g_Rotation, 
               " label='Object rotation' opened=true help='Change the object orientation.' ");

    // Add callback to toggle auto-rotate mode (callback functions are defined above).
    TwAddVarCB(bar, "AutoRotate", TW_TYPE_BOOL32, SetAutoRotateCB, GetAutoRotateCB, NULL, 
               " label='Auto-rotate' key=space help='Toggle auto-rotate mode.' ");

    // Add 'g_LightMultiplier' to 'bar': this is a variable of type TW_TYPE_FLOAT. Its key shortcuts are [+] and [-].
    TwAddVarRW(bar, "Multiplier", TW_TYPE_FLOAT, &g_LightMultiplier, 
               " label='Light booster' min=0.1 max=4 step=0.02 keyIncr='+' keyDecr='-' help='Increase/decrease the light power.' ");

    // Add 'g_LightDirection' to 'bar': this is a variable of type TW_TYPE_DIR3F which defines the light direction
    TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection, 
               " label='Light direction' opened=true help='Change the light direction.' ");

    // Add 'g_MatAmbient' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
    // and is inserted into a group named 'Material'.
    TwAddVarRW(bar, "Ambient", TW_TYPE_COLOR3F, &g_MatAmbient, " group='Material' ");

    // Add 'g_MatDiffuse' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
    // and is inserted into group 'Material'.
    TwAddVarRW(bar, "Diffuse", TW_TYPE_COLOR3F, &g_MatDiffuse, " group='Material' ");

    // Add the enum variable 'g_CurrentShape' to 'bar'
    // (before adding an enum variable, its enum type must be declared to AntTweakBar as follow)
    {
        // ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
        TwEnumVal shapeEV[NUM_SHAPES] = { {SHAPE_TEAPOT, "Teapot"}, {SHAPE_TORUS, "Torus"}, {SHAPE_CONE, "Cone"} };
        // Create a type for the enum shapeEV
        TwType shapeType = TwDefineEnum("ShapeType", shapeEV, NUM_SHAPES);
        // add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
        TwAddVarRW(bar, "Shape", shapeType, &g_CurrentShape, " keyIncr='<' keyDecr='>' help='Change object shape.' ");
    }

    // Store time
    g_RotateTime = GetTimeMs();
    // Init rotation
    SetQuaternionFromAxisAngle(axis, angle, g_Rotation);
    SetQuaternionFromAxisAngle(axis, angle, g_RotateStart);

    // Call the GLUT main loop
    glutMainLoop();

    return 0;
}

