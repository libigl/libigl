#include <GL/glut.h>
#include "YImage.hpp"

#include <iostream>
#include <cstdlib>

void displayFunc( void ) ;
void reshapeFunc( int width, int height ) ;

void keyboardFunc( unsigned char key, int x, int y ) ;
void menu( int id ) ;

void cleanup( void ) ;

static YImage* gImg ;

int main( int argc, char* argv[] )
{
    glutInit( &argc, argv ) ;
    
    std::string toLoad ;
    if( argc == 1 )
    {
        toLoad = "sakura.png" ;
    }
    else if( argc == 2 )
    {
        toLoad = argv[1] ;
    }
    else
    {
        std::cerr << "Usage: " << argv[0] << " image.png" << std::endl ;
        exit(1) ;
    }
    
    gImg = new YImage() ;
    if( !gImg->load( toLoad.c_str() ) )
    {
        std::cerr << "Error loading PNG file: " << toLoad << std::endl ;
        exit(1) ;
    }
    gImg->flip() ;
    
    atexit( cleanup ) ;
    
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE ) ;
    glutInitWindowSize( gImg->width(), gImg->height() ) ;
    glutCreateWindow( toLoad.c_str() ) ;
    
    glutDisplayFunc( displayFunc ) ;
    glutReshapeFunc( reshapeFunc ) ;
    // call reshapeFunc once to set up my projection matrix
    reshapeFunc( glutGet( GLUT_WINDOW_WIDTH ), glutGet( GLUT_WINDOW_HEIGHT ) ) ;
    
    glutKeyboardFunc( keyboardFunc ) ;
    glutCreateMenu( menu ) ;
    glutAddMenuEntry( "Flip 'f'", 'f' ) ;
    glutAddMenuEntry( "Mirror 'm'", 'm' ) ;
    glutAddMenuEntry( "Greyscale 'g'", 'g' ) ;
    glutAddMenuEntry( "", 0 ) ;
    glutAddMenuEntry( "Save image to \"output.png\" 's'", 's' ) ;
    glutAddMenuEntry( "Capture screen to \"screen-capture.png\" 'c'", 'c' ) ;
    glutAddMenuEntry( "", 0 ) ;
    glutAddMenuEntry( "Quit 'q'", 'q' ) ;
    glutAttachMenu( GLUT_LEFT_BUTTON ) ;
    
    glutMainLoop() ;
    
    exit(0) ;
}

void keyboardFunc( unsigned char key, int x, int y )
{
    menu( key ) ;
}

void menu( int id )
{
    switch( id )
    {
        case 'f':
            gImg->flip() ;
            break ;
        
        case 'm':
            gImg->mirror() ;
            break ;
        
        case 'g':
            gImg->greyscale() ;
            break ;
        
        case 's':
            gImg->flip() ;
            gImg->save( "output.png" ) ;
            gImg->flip() ;
            break ;
        
        case 'c':
        {
            GLint view[4] ;
            glGetIntegerv( GL_VIEWPORT, view ) ;
            
            YImage img ;
            img.resize( view[2], view[3] ) ;
            glReadBuffer( GL_FRONT );
            glReadPixels( view[0], view[1], view[2], view[3], GL_RGBA, GL_UNSIGNED_BYTE, img.data() ) ;
            glReadBuffer( GL_BACK );
            
            img.flip() ;
            img.save( "screen-capture.png" ) ;
        }
        break ;
        
        case 'q':
        case 'Q':
            exit(0) ;
            break ;
        
        default:
            break ;
    }
    
    glutPostRedisplay() ;
}

void displayFunc( void )
{
    glClearColor( 0,1,0,1 ) ;
    glClear( GL_COLOR_BUFFER_BIT ) ;
    
    glEnable( GL_BLEND ) ;
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA ) ;
    
    glRasterPos2f( 0,0 ) ;
    glDrawPixels( gImg->width(), gImg->height(), GL_RGBA, GL_UNSIGNED_BYTE, gImg->data() ) ;
    
    glutSwapBuffers() ;
}

void reshapeFunc( int width, int height )
{
    // make each pixel a unit, origin at lower-left
    
    glMatrixMode( GL_PROJECTION ) ;
    glLoadIdentity() ;
    glOrtho( 0,width, 0,height, -10000,10000 ) ;
    
    glViewport( 0,0, width,height ) ;
}

void cleanup( void )
{
    delete gImg ;
    gImg = NULL ;
}
