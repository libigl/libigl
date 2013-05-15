//  ---------------------------------------------------------------------------
//
//  @file       TwGLCoreSDL.c
//  @brief      An example that uses AntTweakBar with OpenGL Core Profile 
//              and SDL 1.3.
//
//              AntTweakBar: http://anttweakbar.sourceforge.net/doc
//              OpenGL:      http://www.opengl.org
//              SDL:         http://www.libsdl.org
//  
//  @author     Philippe Decaudin
//
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#ifdef _WIN32
//  MiniSDL13.h is provided to avoid the need of having SDL installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual SDL library SDK.
#   define USE_MINI_SDL
#endif

//#define GL3_PROTOTYPES 1 ////
//#include <GL3/gl3.h>     ////

#ifdef USE_MINI_SDL
#   include "../src/MiniSDL13.h"
#else
#   include <SDL/SDL.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _WIN32
#   include <windows.h> // required by gl.h
#endif
#include <GL/gl.h>
#include <GL/glu.h>


// In this example, we draw a simple rotating square using the OpenGL core profile
// (which requires much more code than with the compatibility profile).
// A tweak bar is created to allow the user to change the color of the 
// rectangle and see its rotation.

// Part of OpenGL core interface is not directly accessible from the common
// OpenGL header and library (at least on windows) so we have to retrieve the
// core functions using glGetProcAddress. These functions are prefixed by
// underscore to avoid possible confict if a modified gl.h has been installed.

#ifdef _WIN32
#   define glGetProcAddress wglGetProcAddress
#else
#   define GLX_GLXEXT_LEGACY
#   include <GL/glx.h>
#   define glGetProcAddress glXGetProcAddressARB
#endif
#ifndef APIENTRY
#   define APIENTRY
#endif
typedef GLuint (APIENTRY *PFNGLCreateShader)(GLenum type);
typedef void (APIENTRY *PFNGLDeleteShader)(GLuint shader);
typedef void (APIENTRY *PFNGLShaderSource)(GLuint shader, GLsizei count, const char* *str, const GLint *length);
typedef void (APIENTRY *PFNGLCompileShader)(GLuint shader);
typedef void (APIENTRY *PFNGLAttachShader)(GLuint program, GLuint shader);
typedef GLuint (APIENTRY *PFNGLCreateProgram)(void);
typedef void (APIENTRY *PFNGLLinkProgram)(GLuint program);
typedef void (APIENTRY *PFNGLUseProgram)(GLuint program);
typedef void (APIENTRY *PFNGLDeleteProgram)(GLuint program);
typedef void (APIENTRY *PFNGLGenBuffers)(GLsizei n, GLuint *buffers);
typedef void (APIENTRY *PFNGLBindBuffer)(GLenum target, GLuint buffer);
typedef void (APIENTRY *PFNGLVertexAttribPointer)(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid *pointer);
typedef void (APIENTRY *PFNGLEnableVertexAttribArray)(GLuint index);
typedef void (APIENTRY *PFNGLGenVertexArrays)(GLsizei n, GLuint *arrays);
typedef void (APIENTRY *PFNGLBindVertexArray)(GLuint array);
typedef void (APIENTRY *PFNGLDeleteVertexArrays)(GLsizei n, const GLuint *arrays);
typedef GLint (APIENTRY *PFNGLGetAttribLocation)(GLuint program, const char *name);
typedef GLint (APIENTRY *PFNGLGetUniformLocation)(GLuint program, const char *name);
typedef void (APIENTRY *PFNGLUniform1f)(GLint location, GLfloat v0);
typedef void (APIENTRY *PFNGLUniform3f)(GLint location, GLfloat v0, GLfloat v1, GLfloat v2);
typedef void (APIENTRY *PFNGLBufferData)(GLenum target, ptrdiff_t size, const GLvoid *data, GLenum usage);
typedef void (APIENTRY *PFNGLDeleteBuffers)(GLsizei n, const GLuint *buffers);
PFNGLCreateShader _glCreateShader;
PFNGLDeleteShader _glDeleteShader;
PFNGLShaderSource _glShaderSource;
PFNGLCompileShader _glCompileShader;
PFNGLAttachShader _glAttachShader;
PFNGLCreateProgram _glCreateProgram;
PFNGLLinkProgram _glLinkProgram;
PFNGLUseProgram _glUseProgram;
PFNGLDeleteProgram _glDeleteProgram;
PFNGLGenBuffers _glGenBuffers;
PFNGLBindBuffer _glBindBuffer;
PFNGLVertexAttribPointer _glVertexAttribPointer;
PFNGLEnableVertexAttribArray _glEnableVertexAttribArray;
PFNGLGenVertexArrays _glGenVertexArrays;
PFNGLBindVertexArray _glBindVertexArray;
PFNGLDeleteVertexArrays _glDeleteVertexArrays;
PFNGLGetAttribLocation _glGetAttribLocation;
PFNGLGetUniformLocation _glGetUniformLocation;
PFNGLUniform1f _glUniform1f;
PFNGLUniform3f _glUniform3f;
PFNGLBufferData _glBufferData;
PFNGLDeleteBuffers _glDeleteBuffers;
#ifndef GL_ARRAY_BUFFER
#   define GL_ARRAY_BUFFER      0x8892
#endif
#ifndef GL_STATIC_DRAW
#   define GL_STATIC_DRAW       0x88E4
#endif
#ifndef GL_VERTEX_SHADER
#   define GL_VERTEX_SHADER     0x8B31
#endif
#ifndef GL_FRAGMENT_SHADER
#   define GL_FRAGMENT_SHADER   0x8B30
#endif

int LoadGLCoreFunctions()
{
    _glCreateShader = (PFNGLCreateShader)glGetProcAddress("glCreateShader");
    _glDeleteShader = (PFNGLDeleteShader)glGetProcAddress("glDeleteShader");
    _glShaderSource = (PFNGLShaderSource)glGetProcAddress("glShaderSource"); 
    _glCompileShader = (PFNGLCompileShader)glGetProcAddress("glCompileShader");
    _glAttachShader = (PFNGLAttachShader)glGetProcAddress("glAttachShader"); 
    _glCreateProgram = (PFNGLCreateProgram)glGetProcAddress("glCreateProgram");
    _glLinkProgram = (PFNGLLinkProgram)glGetProcAddress("glLinkProgram");
    _glUseProgram = (PFNGLUseProgram)glGetProcAddress("glUseProgram");
    _glDeleteProgram = (PFNGLDeleteProgram)glGetProcAddress("glDeleteProgram");
    _glGenBuffers = (PFNGLGenBuffers)glGetProcAddress("glGenBuffers");
    _glBindBuffer = (PFNGLBindBuffer)glGetProcAddress("glBindBuffer");
    _glVertexAttribPointer = (PFNGLVertexAttribPointer)glGetProcAddress("glVertexAttribPointer");
    _glEnableVertexAttribArray = (PFNGLEnableVertexAttribArray)glGetProcAddress("glEnableVertexAttribArray");
    _glGenVertexArrays = (PFNGLGenVertexArrays)glGetProcAddress("glGenVertexArrays");
    _glBindVertexArray = (PFNGLBindVertexArray)glGetProcAddress("glBindVertexArray");
    _glDeleteVertexArrays = (PFNGLDeleteVertexArrays)glGetProcAddress("glDeleteVertexArrays");
    _glGetAttribLocation = (PFNGLGetAttribLocation)glGetProcAddress("glGetAttribLocation");
    _glGetUniformLocation = (PFNGLGetUniformLocation)glGetProcAddress("glGetUniformLocation");
    _glUniform1f = (PFNGLUniform1f)glGetProcAddress("glUniform1f");
    _glUniform3f = (PFNGLUniform3f)glGetProcAddress("glUniform3f");
    _glBufferData = (PFNGLBufferData)glGetProcAddress("glBufferData");
    _glDeleteBuffers = (PFNGLDeleteBuffers)glGetProcAddress("glDeleteBuffers");

    if (_glCreateShader == NULL || _glDeleteShader == NULL || _glShaderSource == NULL || _glCompileShader == NULL
        || _glAttachShader == NULL || _glCreateProgram == NULL || _glLinkProgram == NULL || _glUseProgram  == NULL
        || _glDeleteProgram == NULL || _glGenBuffers == NULL || _glBindBuffer == NULL || _glVertexAttribPointer == NULL
        || _glEnableVertexAttribArray == NULL || _glGenVertexArrays == NULL || _glBindVertexArray == NULL
        || _glDeleteVertexArrays == NULL || _glGetAttribLocation == NULL || _glGetUniformLocation == NULL
        || _glUniform1f == NULL || _glUniform3f == NULL || _glBufferData == NULL || _glDeleteBuffers == NULL)
        return 0;
    else 
        return 1;
}


// Shaders globals
GLuint vshader, fshader, program, varray, buffer;
GLint cosa, sina, colorloc;
float angle = 0, quat[4];
float color[] = {0.8f, 1.0f, 0.2f};
float FLOAT_PI = 3.14159265f;

void InitRender()
{
    // Vertex shader
    char *vsource[] = {
        "#version 150 core\n"
        "uniform float cosa, sina;"
        "in vec3 vertex;"
        "void main() { gl_Position = vec4(cosa*vertex.x-sina*vertex.y, sina*vertex.x+cosa*vertex.y, 0, 1); }"
    };
    // Fragment shader
    char *fsource[] = {
        "#version 150 core\n"
        "uniform vec3 color;"
        "out vec4 fcolor;"
        "void main() { fcolor = vec4(color, 1); }"
    };
    // Geometry vertex array
    GLfloat vertices[] = { 
        -0.5f,  0.5f, 0,
        -0.5f, -0.5f, 0, 
         0.5f,  0.5f, 0,
         0.5f, -0.5f, 0
    };
    GLint vlocation;

    // Compile and bind shaders
    vshader = _glCreateShader(GL_VERTEX_SHADER);
    fshader = _glCreateShader(GL_FRAGMENT_SHADER);
    program = _glCreateProgram();
    _glShaderSource(vshader, 1, vsource, NULL);
    _glShaderSource(fshader, 1, fsource, NULL);
    _glCompileShader(vshader);
    _glCompileShader(fshader);
    _glAttachShader(program, vshader);
    _glAttachShader(program, fshader);
    _glLinkProgram(program);
    _glUseProgram(program);

    cosa = _glGetUniformLocation(program, "cosa");
    sina = _glGetUniformLocation(program, "sina");
    colorloc = _glGetUniformLocation(program, "color");

    // Create and bind vertex buffer
    _glGenVertexArrays(1, &varray);
    _glBindVertexArray(varray);
    _glGenBuffers(1, &buffer);
    _glBindBuffer(GL_ARRAY_BUFFER, buffer);
    _glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    vlocation = _glGetAttribLocation(program, "vertex");
    _glVertexAttribPointer(vlocation, 3, GL_FLOAT, GL_TRUE, 0, NULL);
    _glEnableVertexAttribArray(vlocation);

    // GL states
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
}

void UninitRender()
{
    _glDeleteProgram(program);
    _glDeleteShader(fshader);
    _glDeleteShader(vshader);
    _glDeleteBuffers(1, &buffer);
    _glDeleteVertexArrays(1, &varray);
}

void Render()
{
    _glUniform1f(cosa, (float)cos(angle));
    _glUniform1f(sina, (float)sin(angle));
    _glUniform3f(colorloc, color[0], color[1], color[2]);
    _glBindVertexArray(varray);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}


TwBar *CreateTweakBar()
{
    TwBar *bar;

    // Create a tweak bar
    bar = TwNewBar("TweakBar");
    TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with SDL and OpenGL Core Profile.\n' "); // Message added to the help bar.

    // Add variables
    TwAddVarRW(bar, "Rotation", TW_TYPE_QUAT4F, &quat, " opened=true help='Rectangle rotation' ");
    TwAddVarRW(bar, "Color", TW_TYPE_COLOR3F, &color, " opened=true help='Rectangle color' ");
 
    return bar;
}


// SDL redefines main
#ifdef main
#   undef main
#endif

int main()
{
    const SDL_VideoInfo* video = NULL;
    int width  = 480, height = 480;
    int bpp, flags;
    int quit = 0;

    // Initialize SDL, then get the current video mode and use it to create a SDL window.
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        fprintf(stderr, "Video initialization failed: %s\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }
    video = SDL_GetVideoInfo();
    if (!video) 
    {
        fprintf(stderr, "Video query failed: %s\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }
    // Request GL context to be OpenGL 3.2 Core Profile
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    // Other GL attributes
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    bpp = video->vfmt->BitsPerPixel;
    flags = SDL_OPENGL | SDL_HWSURFACE;
    if (!SDL_SetVideoMode(width, height, bpp, flags))
    {
        fprintf(stderr, "Video mode set failed: %s\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }
    SDL_WM_SetCaption("AntTweakBar example using OpenGL Core Profile and SDL", "AntTweakBar+GLCore+SDL");

    // Enable SDL unicode and key-repeat
    SDL_EnableUNICODE(1);
    SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY, SDL_DEFAULT_REPEAT_INTERVAL);

    // Load some OpenGL core functions
    if (!LoadGLCoreFunctions())
    {
        fprintf(stderr, "OpenGL 3.2 not supported.\n");
        SDL_Quit();
        exit(1);
    }

    // Initialize AntTweakBar
    if (!TwInit(TW_OPENGL_CORE, NULL)) {
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        SDL_Quit();
        exit(1);
    }
    // Tell the window size to AntTweakBar
    TwWindowSize(width, height);
    // Create a tweak bar
    CreateTweakBar();

    // Set OpenGL viewport
    glViewport(0, 0, width, height);

    // Prepare GL shaders and programs for drawing
    InitRender();

    // Main loop:
    // - Draw scene
    // - Process events
    while (!quit)
    {
        SDL_Event event;
        int handled;
        GLenum error;

        // Clear screen
        glClearColor(0.5f, 0.75f, 0.8f, 1);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        // Update angle and draw geometry
        angle = (float)SDL_GetTicks()/25.0f * (FLOAT_PI/180.0f);
        quat[0] = quat[1] = 0;
        quat[2] = (float)sin(angle/2.0f);
        quat[3] = (float)cos(angle/2.0f);
        Render();

        // Draw tweak bars
        TwDraw();

        // Present frame buffer
        SDL_GL_SwapBuffers();

        // Process incoming events
        while (SDL_PollEvent(&event)) 
        {
            // Send event to AntTweakBar
            handled = TwEventSDL(&event, SDL_MAJOR_VERSION, SDL_MINOR_VERSION);

            // If event has not been handled by AntTweakBar, process it
            if (!handled)
            {
                switch (event.type)
                {
                case SDL_QUIT:  // Window is closed
                    quit = 1;
                    break;

                case SDL_VIDEORESIZE:   // Window size has changed
                    // Resize SDL video mode
                    width = event.resize.w;
                    height = event.resize.h;
                    if (!SDL_SetVideoMode(width, height, bpp, flags))
                        fprintf(stderr, "WARNING: Video mode set failed: %s\n", SDL_GetError());

                    // Resize OpenGL viewport
                    glViewport(0, 0, width, height);
                    
                    // Restore OpenGL states
                    InitRender();
                    
                    // TwWindowSize has been called by TwEventSDL, 
                    // so it is not necessary to call it again here.

                    break;
                }
            }
        }

        while ((error = glGetError()) != GL_NO_ERROR) 
            fprintf(stderr, "GL error detected: 0x%04X\n", error);

    } // End of main loop

    // Terminate AntTweakBar
    TwTerminate();

    // Delete GL shaders and buffer
    UninitRender();

    // Terminate SDL
    SDL_Quit();

    return 0;
}  

