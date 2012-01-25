//  ---------------------------------------------------------------------------
//
//  @file       TwSimpleSFML.cpp
//  @brief      A simple example that uses AntTweakBar with OpenGL and SFML.
//              This example draws moving cubic particles with some
//              interactive control on particles generation.
//
//              AntTweakBar: http://www.antisphere.com/Wiki/tools:anttweakbar
//              OpenGL:      http://www.opengl.org
//              SFML:        http://www.sfml-dev.org
//  
//  @author     Philippe Decaudin - http://www.antisphere.com
//
//  Compilation:
//  http://www.antisphere.com/Wiki/tools:anttweakbar:examples#twsimplesfml
//
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#if defined(_WIN32)
//  MiniSFML16.h is provided to avoid the need of having SFML installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual SFML library SDK.
#   define USE_MINI_SFML
#endif

#ifdef USE_MINI_SFML
#   include "../src/MiniSFML16.h"
#else
#   include <SFML/Graphics.hpp>
#endif

#if defined(_WIN32)
#   include <windows.h> // required by gl.h
#endif
#include <GL/gl.h>
#include <GL/glu.h>

#include <list>
#include <cstdlib>
#include <cmath>


// Pseudo-random value between -1 and 1
float Random() 
{
    return 2.0f * ((float)rand() / RAND_MAX) - 1.0f;
}

// Particle randomly initialized
struct Particle
{
    float Size;
    float Position[3];      // [px, py, pz]
    float Speed[3];         // [vx, vy, vz]
    float RotationAxis[3];  // [rx, ry, rz]
    float RotationAngle;    // in degree
    float RotationSpeed;
    float Color[3];         // [r, g, b]
    float Age;
    Particle(float size, float speedDir[3], float speedNorm, float color[3]) // Constructor
    {
        Size = size * (1.0f + 0.2f * Random());
        Position[0] = Position[1] = Position[2] = 0;
        Speed[0] = speedNorm * (speedDir[0] + 0.1f * Random());
        Speed[1] = speedNorm * (speedDir[1] + 0.1f * Random());
        Speed[2] = speedNorm * (speedDir[2] + 0.1f * Random());
        RotationAxis[0] = Random();
        RotationAxis[1] = Random();
        RotationAxis[2] = Random();
        RotationAngle = 360.0f * Random();
        RotationSpeed = 360.0f * Random();
        Color[0] = color[0] + 0.2f * Random();
        Color[1] = color[1] + 0.2f * Random();
        Color[2] = color[2] + 0.2f * Random();
        Age = 0;
    }
    void Update(float dt) // Apply one animation step
    {
        Position[0] += dt * Speed[0];
        Position[1] += dt * Speed[1];
        Position[2] += dt * Speed[2];
        Speed[1] -= dt * 9.81f; // gravity
        RotationAngle += dt * RotationSpeed;
        Age += dt;
    }
};


// Main 
int main()
{
    // Create main window
    sf::RenderWindow app(sf::VideoMode(800, 600), "AntTweakBar simple example using SFML");
    app.PreserveOpenGLStates(true);

    // Particules
    std::list<Particle> particles;
    std::list<Particle>::iterator p;
    float birthCount = 0;
    float birthRate = 20;               // number of particles generated per second
    float maxAge = 3.0f;                // particles life time
    float speedDir[3] = {0, 1, 0};      // initial particles speed direction
    float speedNorm = 7.0f;             // initial particles speed amplitude
    float size = 0.1f;                  // particles size
    float color[3] = {0.8f, 0.6f, 0};   // particles color
    float bgColor[3] = {0, 0.6f, 0.6f}; // background color


    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);

    // Tell the window size to AntTweakBar
    TwWindowSize(app.GetWidth(), app.GetHeight());

    // Create a tweak bar
    TwBar *bar = TwNewBar("Particles");
    TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with SFML and OpenGL.' "); // Message added to the help bar.
    
    // Change bar position
    int barPos[2] = {16, 240};
    TwSetParam(bar, NULL, "position", TW_PARAM_INT32, 2, &barPos);
    
    // Add 'birthRate' to 'bar': this is a modifiable variable of type TW_TYPE_FLOAT in range [0.1, 100]. Its shortcuts are [+] and [-].
    TwAddVarRW(bar, "Birth rate", TW_TYPE_FLOAT, &birthRate, " min=0.1 max=100 step=0.1 keyIncr='+' keyDecr='-' ");

    // Add 'speedNorm' to 'bar': this is a modifiable variable of type TW_TYPE_FLOAT in range [0.1, 10]. Its shortcuts are [s] and [S].
    TwAddVarRW(bar, "Speed", TW_TYPE_FLOAT, &speedNorm, " min=0.1 max=10 step=0.1 keyIncr='s' keyDecr='S' ");

    // Add 'speedDir' to 'bar': this is a modifiable variable of type TW_TYPE_DIR3F. Just displaying the arrow widget
    TwAddVarRW(bar, "Direction", TW_TYPE_DIR3F, &speedDir, " opened=true showval=false ");

    // Add 'color' to 'bar': this is a modifiable variable of type TW_TYPE_COLOR3F. Switched to HLS
    TwAddVarRW(bar, "Color", TW_TYPE_COLOR3F, &color, " colorMode=hls opened=true ");

    // Add 'bgColor' to 'bar': this is a modifiable variable of type TW_TYPE_COLOR3F. Switched to HLS
    TwAddVarRW(bar, "Background color", TW_TYPE_COLOR3F, &bgColor, " colorMode=hls opened=true ");

    // Initialize OpenGL states
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.f, (float)app.GetWidth()/app.GetHeight(), 0.1f, 100.f);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

    // Init time
    sf::Clock clock;
    float time = clock.GetElapsedTime();

    // Main loop
    while (app.IsOpened())
    {
        // Process events
        sf::Event event;
        while (app.GetEvent(event))
        {
            // Send event to AntTweakBar
            int handled = TwEventSFML(&event, 1, 6); // Assume SFML version 1.6 here

            // If event has not been handled by AntTweakBar, process it
            if( !handled )
            {
                // Close or Escape
                if (event.Type == sf::Event::Closed
                    || (event.Type == sf::Event::KeyPressed && event.Key.Code == sf::Key::Escape))
                    app.Close();

                // Resize
                if (event.Type == sf::Event::Resized)
                {
                    glViewport(0, 0, event.Size.Width, event.Size.Height);
                    glMatrixMode(GL_PROJECTION);
                    glLoadIdentity();
                    gluPerspective(90.f, (float)event.Size.Width/event.Size.Height, 1.f, 500.f);
                    glMatrixMode(GL_MODELVIEW);

                    // TwWindowSize has been called by TwEventSFML, 
                    // so it is not necessary to call it again here.
                }
            }
        }

        if (!app.IsOpened())
            continue;

        // Update time
        float dt = clock.GetElapsedTime() - time;
        if (dt < 0) dt = 0;
        time += dt;

        // Update particles
        p = particles.begin(); 
        while (p != particles.end())
        {
            p->Update(dt);
            if (p->Age >= maxAge) 
                p = particles.erase(p); // Die!
            else
                ++p;
        }

        // Generate new particles
        birthCount += dt * birthRate;
        while (birthCount >= 1.0f) 
        {
            particles.push_back(Particle(size, speedDir, speedNorm, color));
            birthCount--;
        }

        // Clear depth buffer
        glClearColor(bgColor[0], bgColor[1], bgColor[2], 1);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        // Draw particles
        for (p = particles.begin(); p != particles.end(); ++p)
        {
            glColor4fv(p->Color);
            glLoadIdentity();
            glTranslatef(0.0f, -1.0f, -3.0f); // Camera position
            glTranslatef(p->Position[0], p->Position[1], p->Position[2]);
            glScalef(p->Size, p->Size, p->Size);
            glRotatef(p->RotationAngle, p->RotationAxis[0], p->RotationAxis[1], p->RotationAxis[2]);

            // Draw a cube
            glBegin(GL_QUADS);
                glNormal3f(0,0,-1); glVertex3f(0,0,0); glVertex3f(0,1,0); glVertex3f(1,1,0); glVertex3f(1,0,0); // front face
                glNormal3f(0,0,+1); glVertex3f(0,0,1); glVertex3f(1,0,1); glVertex3f(1,1,1); glVertex3f(0,1,1); // back face
                glNormal3f(-1,0,0); glVertex3f(0,0,0); glVertex3f(0,0,1); glVertex3f(0,1,1); glVertex3f(0,1,0); // left face
                glNormal3f(+1,0,0); glVertex3f(1,0,0); glVertex3f(1,1,0); glVertex3f(1,1,1); glVertex3f(1,0,1); // right face
                glNormal3f(0,-1,0); glVertex3f(0,0,0); glVertex3f(1,0,0); glVertex3f(1,0,1); glVertex3f(0,0,1); // bottom face  
                glNormal3f(0,+1,0); glVertex3f(0,1,0); glVertex3f(0,1,1); glVertex3f(1,1,1); glVertex3f(1,1,0); // top face
            glEnd();
        }

        TwDraw();

        // Finally, display the rendered frame on screen
        app.Display();
    }

    // Un-initialize AntTweakBar 
    TwTerminate();

    return EXIT_SUCCESS;
}
