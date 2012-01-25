//  ---------------------------------------------------------------------------
//
//  @file       TwGLCoreSDL.c
//  @brief      An example that uses AntTweakBar with OpenGL Core Profile 
//              and SDL 1.3.
//
//              AntTweakBar: http://www.antisphere.com/Wiki/tools:anttweakbar
//              OpenGL:      http://www.opengl.org
//              SDL:         http://www.libsdl.org
//  
//  @author     Philippe Decaudin - http://www.antisphere.com
//
//  Compilation:
//  http://www.antisphere.com/Wiki/tools:anttweakbar:examples#twsimplesdl
//
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#if defined(_WIN32) || defined(_WIN64)
//  MiniSDL13.h is provided to avoid the need of having SDL installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual SDL library SDK.
#   define USE_MINI_SDL
#endif

#define GL3_PROTOTYPES 1 ////
#include <GL3/gl3.h> ////

#ifdef USE_MINI_SDL
#   include "../src/MiniSDL13.h"
#else
#   include <SDL/SDL.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#if defined(_WIN32) || defined(_WIN64)
#   include <windows.h> // required by gl.h
#endif
#include <GL/gl.h>
#include <GL/glu.h>
*/

// In this example, we draw a simple rectangle using the OpenGL core profile
// (which requires much more code than with the compatibility profile).
// A tweak bar is created to allow the user to change the size and color
// of the rectangle.

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
typedef void (APIENTRY *PFNGLShaderSource)(GLuint shader, GLsizei count, const GLchar* *string, const GLint *length);
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
typedef GLint (APIENTRY *PFNGLGetAttribLocation)(GLuint program, const GLchar *name);
typedef void (APIENTRY *PFNGLBufferData)(GLenum target, GLsizeiptr size, const GLvoid *data, GLenum usage);
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
PFNGLBufferData _glBufferData;
PFNGLDeleteBuffers _glDeleteBuffers;

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
    _glBufferData = (PFNGLBufferData)glGetProcAddress("glBufferData");
    _glDeleteBuffers = (PFNGLDeleteBuffers)glGetProcAddress("glDeleteBuffers");

    if (_glCreateShader == NULL || _glDeleteShader == NULL || _glShaderSource == NULL || _glCompileShader == NULL
        || _glAttachShader == NULL || _glCreateProgram == NULL || _glLinkProgram == NULL || _glUseProgram  == NULL
        || _glDeleteProgram == NULL || _glGenBuffers == NULL || _glBindBuffer == NULL || _glVertexAttribPointer == NULL
        || _glEnableVertexAttribArray == NULL || _glGenVertexArrays == NULL || _glBindVertexArray == NULL
        || _glDeleteVertexArrays == NULL || _glGetAttribLocation == NULL || _glBufferData == NULL || _glDeleteBuffers == NULL)
        return 0;
    else 
        return 1;
}

GLuint vshader, fshader, program, varray, buffer;

void InitRender()
{
    // vertex shader
    GLchar *vsource[] = {
        "#version 150 core\n"
        "in vec3 vertex;"
        "void main() { gl_Position = vec4(vertex, 1); }"
    };
    // fragment shader
    GLchar *fsource[] = {
        "#version 150 core\n"
        "out vec4 color;"
        "void main() { color = vec4(0, 0, 1, 0); }"
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
    glEnable(GL_DEPTH_TEST);
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

TwBar *CreateTweakBar()
{
    TwBar *bar;

    // Create a tweak bar
    bar = TwNewBar("TweakBar");
    TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with SDL and OpenGL Core Profile.\nPress [Space] to toggle fullscreen.' "); // Message added to the help bar.
/*
    // Add 'numCurves' to 'bar': this is a modifiable variable of type TW_TYPE_INT32. Its shortcuts are [c] and [C].
    TwAddVarRW(bar, "NumCubes", TW_TYPE_INT32, &numCubes, 
               " label='Number of cubes' min=1 max=100 keyIncr=c keyDecr=C help='Defines the number of cubes in the scene.' ");

    // Add 'ka', 'kb and 'kc' to 'bar': they are modifiable variables of type TW_TYPE_DOUBLE
    TwAddVarRW(bar, "ka", TW_TYPE_DOUBLE, &ka, 
               " label='X path coeff' keyIncr=1 keyDecr=CTRL+1 min=-10 max=10 step=0.01 ");
    TwAddVarRW(bar, "kb", TW_TYPE_DOUBLE, &kb, 
               " label='Y path coeff' keyIncr=2 keyDecr=CTRL+2 min=-10 max=10 step=0.01 ");
    TwAddVarRW(bar, "kc", TW_TYPE_DOUBLE, &kc, 
               " label='Z path coeff' keyIncr=3 keyDecr=CTRL+3 min=-10 max=10 step=0.01 ");

    // Add 'color0' and 'color1' to 'bar': they are modifable variables of type TW_TYPE_COLOR3F (3 floats color)
    TwAddVarRW(bar, "color0", TW_TYPE_COLOR3F, &color0, 
               " label='Start color' help='Color of the first cube.' ");
    TwAddVarRW(bar, "color1", TW_TYPE_COLOR3F, &color1, 
               " label='End color' help='Color of the last cube. Cube colors are interpolated between the Start and End colors.' ");

    // Add 'quit' to 'bar': this is a modifiable (RW) variable of type TW_TYPE_BOOL32 
    // (a boolean stored in a 32 bits integer). Its shortcut is [ESC].
    TwAddVarRW(bar, "Quit", TW_TYPE_BOOL32, &quit, 
               " label='Quit?' true='+' false='-' key='ESC' help='Quit program.' ");
*/
 
    return bar;
}


// SDL redefines main
#ifdef main
#   undef main
#endif

int main()
{
    const SDL_VideoInfo* video = NULL;
    int width  = 640, height = 480;
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
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    bpp = video->vfmt->BitsPerPixel;
    flags = SDL_OPENGL | SDL_HWSURFACE | SDL_RESIZABLE;
    //flags |= SDL_FULLSCREEN;
    if (!SDL_SetVideoMode(width, height, bpp, flags))
    {
        fprintf(stderr, "Video mode set failed: %s\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }
    SDL_WM_SetCaption("AntTweakBar example using SDL and OpenGL Core Profile", "AntTweakBar+SDL+GLCore");

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

        // Draw geometry
//        _glBindVertexArray(varray);
//        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

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

                case SDL_KEYDOWN:
                    if (event.key.keysym.sym==' ') // toggle fullscreen if Space key is pressed
                    {
                        flags ^= SDL_FULLSCREEN;
                        SDL_SetVideoMode(800, 600, bpp, flags);

                        // Push a resize event because SDL does not do it for us
                        event.type = SDL_VIDEORESIZE;
                        event.resize.w = 800;
                        event.resize.h = 600;
                        SDL_PushEvent(&event);
                    }
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

