//  ---------------------------------------------------------------------------
//
//  @file       TwString.cpp
//  @brief      This example illustrates the use of the different types of
//              AntTweakBar string variables.
//              The graphic window is created by GLUT.
//
//              AntTweakBar: http://anttweakbar.sourceforge.net/doc
//              OpenGL:      http://www.opengl.org
//              GLUT:        http://opengl.org/resources/libraries/glut
//  
//  @author     Philippe Decaudin
//
//  ---------------------------------------------------------------------------

#include <AntTweakBar.h>

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

#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#if !defined _MSC_VER
#   define _snprintf snprintf
#endif


// ---------------------------------------------------------------------------
// 1) Callback functions for std::string variables
// ---------------------------------------------------------------------------

// Function called to copy the content of a std::string (souceString) handled 
// by the AntTweakBar library to destinationClientString handled by our application
void TW_CALL CopyStdStringToClient(std::string& destinationClientString, const std::string& sourceLibraryString)
{
    destinationClientString = sourceLibraryString;
}

std::vector<std::string> g_BarTitles;

// Callback function called by AntTweakBar to set the "EditTitle" std::string variable
void TW_CALL SetBarTitleCB(const void *value, void *clientData)
{
    const std::string *newTitle = (const std::string *)(value);
    int barIndex = *(int *)(&clientData);   // clientData stores the bar index

    // Stores the new bar title
    g_BarTitles[barIndex] = *newTitle;

    // Create the def command to change the bar label (ie., its title)
    std::stringstream def;
    def << "bar_" << barIndex << " label=`" << g_BarTitles[barIndex] << "`";
    // Execute the command
    TwDefine(def.str().c_str());
}

// Callback function called by AntTweakBar to get the "EditTitle" std::string variable
void TW_CALL GetBarTitleCB(void *value, void *clientData)
{
    std::string *destStringPtr = (std::string *)(value);
    int barIndex = *(int *)(&clientData);   // clientData stores the bar index
    std::string title = g_BarTitles[barIndex];

    // Do not assign destStringPtr directly (see TwCopyStdStringToLibrary doc for explanation):
    // Use TwCopyStdStringToLibrary to copy the bar title string to AntTweakBar
    TwCopyStdStringToLibrary(*destStringPtr, title);
}

// Callback function to create a bar with a given title
void TW_CALL CreateBarCB(void *clientData)
{
    const std::string *title = (const std::string *)(clientData);

    // Create a unique bar name
    int barIndex = (int)g_BarTitles.size();
    std::stringstream name;
    name << "bar_" << barIndex;

    g_BarTitles.push_back(*title);

    // Create bar
    TwBar *bar = TwNewBar(name.str().c_str());
    TwAddButton(bar, "Info", NULL, NULL, " label='std::string variable:' ");

    // Set bar label (ie. the title)
    std::stringstream def;
    def << name.str() << " label=`" << *title << "` ";
    TwDefine(def.str().c_str());

    // Cast barNum as void* to use it as clientData 
    // (doing so it could be directly sent through the get/set callbacks)
    void *barIndexAsVoidPtr = *(void **)&barIndex;

    // Add a std::string variable to the bar to edit its title,
    // The variable will be accessed through callbacks
    TwAddVarCB(bar, "EditTitle", TW_TYPE_STDSTRING, SetBarTitleCB, GetBarTitleCB, barIndexAsVoidPtr, 
               " label='Edit bar title' help='Edit this string to change the tweak bar title.' ");
}


// ---------------------------------------------------------------------------
// 2) Callback functions for C-Dynamic string variables
// ---------------------------------------------------------------------------

// Function called to copy the content of a C-Dynamic String (src) handled by
// the AntTweakBar library to a C-Dynamic string (*destPtr) handled by our application
void TW_CALL CopyCDStringToClient(char **destPtr, const char *src)
{
    size_t srcLen = (src!=NULL) ? strlen(src) : 0;
    size_t destLen = (*destPtr!=NULL) ? strlen(*destPtr) : 0;

    // Alloc or realloc dest memory block if needed
    if( *destPtr==NULL )
        *destPtr = (char *)malloc(srcLen+1);
    else if( srcLen>destLen )
        *destPtr = (char *)realloc(*destPtr, srcLen+1);

    // Copy src
    if( srcLen>0 )
        strncpy(*destPtr, src, srcLen);
    (*destPtr)[srcLen] = '\0'; // null-terminated string
}

// Callback function called by AntTweakBar to set the "TextLine" CDString variable
void TW_CALL SetTextLineCB(const void *value, void *clientData)
{
    const char *src = *(const char **)value;
    char **destPtr = (char **)clientData;
    
    // Copies src to *destPtr (destPtr might be reallocated)
    CopyCDStringToClient(destPtr, src);

    // Change the label of the "Echo" inactive button
    size_t srcLen = strlen(src);
    if( srcLen>0 )
    {
        char *def = (char *)malloc(128+srcLen);
        _snprintf(def, 128+srcLen, " Main/Echo label=`%s` ", src);
        TwDefine(def);
        free(def);
    }
    else
        TwDefine(" Main/Echo label=` ` ");
}

// Callback function called by AntTweakBar to get the "TextLine" CDString variable
void TW_CALL GetTextLineCB(void *value, void *clientData)
{
    char **destPtr = (char **)value;
    char *src = *(char **)clientData;

    // Do not assign destPtr directly:
    // Use TwCopyCDStringToLibrary to copy TextLine to AntTweakBar
    TwCopyCDStringToLibrary(destPtr, src);
}


// ---------------------------------------------------------------------------
// 3) Callback functions for C-Static sized string variables
// ---------------------------------------------------------------------------

// A static sized string
char g_CapStr[17] = "16 chars max"; // 17 = 16 + the null termination char

// A utility function: Convert a C string to lower or upper case
void CaseCopy(char *dest, const char *src, size_t maxLength, int capCase)
{
    size_t i;
    if( capCase==0 ) // lower case
        for( i=0; i<maxLength-1 && src[i]!='\0'; ++i )
            dest[i] = (char)tolower(src[i]);
    else // upper case
        for( i=0; i<maxLength-1 && src[i]!='\0'; ++i )
            dest[i] = (char)toupper(src[i]);
    dest[i] = '\0'; // ensure that dest is null-terminated
}

// Callback function called by AntTweakBar to set the "CapStr" CSString variable
void TW_CALL SetCapStrCB(const void *value, void *clientData)
{
    const char *src = (const char *)value;
    int capCase = *(int *)clientData;
    CaseCopy(g_CapStr, src, sizeof(g_CapStr), capCase);
}

// Callback function called by AntTweakBar to get the "CapStr" CSString variable
void TW_CALL GetCapStrCB(void *value, void *clientData)
{
    char *dest = (char *)value;
    int capCase = *(int *)clientData;
    CaseCopy(dest, g_CapStr, sizeof(g_CapStr), capCase);
}


// ---------------------------------------------------------------------------
// GLUT callbacks
// ---------------------------------------------------------------------------

// Callback function called by GLUT to render screen
void OnDisplay(void)
{
    // Clear frame buffer
    glClearColor(0.5f, 0.5f, 0.6f, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    // App drawing here
    // ...

    // Draw tweak bars
    TwDraw();

    // Present frame buffer
    glutSwapBuffers();
}

// Callback function called by GLUT when window size changes
void OnReshape(int width, int height)
{
    // Set OpenGL viewport
    glViewport(0, 0, width, height);

    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);
}

// Function called at exit
void OnTerminate(void)
{ 
    // terminate AntTweakBar
    TwTerminate();
}

// Event callbacks
void OnMouseButton(int glutButton, int glutState, int mouseX, int mouseY)
{
    // send event to AntTweakBar
    if (TwEventMouseButtonGLUT(glutButton, glutState, mouseX, mouseY))
        glutPostRedisplay(); // request redraw if event has been handled
}

void OnMouseMotion(int mouseX, int mouseY)
{
    // send event to AntTweakBar
    if (TwEventMouseMotionGLUT(mouseX, mouseY))
        glutPostRedisplay(); // request redraw if event has been handled
}

void OnKeyboard(unsigned char glutKey, int mouseX, int mouseY)
{
    // send event to AntTweakBar
    if (TwEventKeyboardGLUT(glutKey, mouseX, mouseY))
        glutPostRedisplay(); // request redraw if event has been handled
}

void OnSpecial(int glutKey, int mouseX, int mouseY)
{
    // send event to AntTweakBar
    if (TwEventSpecialGLUT(glutKey, mouseX, mouseY))
        glutPostRedisplay(); // request redraw if event has been handled
}


// ---------------------------------------------------------------------------
// Main function (application based on GLUT)
// ---------------------------------------------------------------------------

int main(int argc, char *argv[]) 
{
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 480);
    glutCreateWindow("AntTweakBar string example");
    glutCreateMenu(NULL);

    // Set GLUT callbacks
    glutDisplayFunc(OnDisplay);
    glutReshapeFunc(OnReshape);
    atexit(OnTerminate);  // Called after glutMainLoop ends

    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);

    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc(OnMouseButton);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc(OnMouseMotion);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc(OnMouseMotion);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc(OnKeyboard);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc(OnSpecial);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);


    // Create a tweak bar
    TwBar *bar = TwNewBar("Main");
    TwDefine(" Main label='~ String variable examples ~' fontSize=3 position='180 16' size='270 440' valuesWidth=100 ");


    //
    // 1) C++ std::string variable example
    //

    TwAddButton(bar, "Info1.1", NULL, NULL, " label='1) This example uses' ");
    TwAddButton(bar, "Info1.2", NULL, NULL, " label='std::string variables' ");
    
    // Define the required callback function to copy a std::string (see TwCopyStdStringToClientFunc documentation)
    TwCopyStdStringToClientFunc(CopyStdStringToClient);
    
    // Adding a std::string variable
    std::string newBarTitle = "a title";
    TwAddVarRW(bar, "NewBarTitle", TW_TYPE_STDSTRING, &newBarTitle, 
               " label='Bar title' group=StdString help='Define a title for the new tweak bar.' ");
    
    // Add a button to create a new bar using the title
    TwAddButton(bar, "NewBarCreate", CreateBarCB, &newBarTitle, 
                " label='--> Create' group=StdString key=c help='Create a new tweak bar.' ");
    
    // Set the group label & separator
    TwDefine(" Main/StdString label='Create a new tweak bar' help='This example demonstates different use of std::string variables.' ");
    TwAddSeparator(bar, "Sep1", "");
    TwAddButton(bar, "Blank1", NULL, NULL, " label=' ' ");


    //
    // 2) C-Dynamic string variable example
    //

    TwAddButton(bar, "Info2.1", NULL, NULL, "label='2) This example uses' ");
    TwAddButton(bar, "Info2.2", NULL, NULL, "label='C-Dynamic string variables' ");
    
    // Define the required callback function to copy a CDString (see TwCopyCDStringToClientFunc documentation)
    TwCopyCDStringToClientFunc(CopyCDStringToClient);
    
    // Add a CDString variable
    char *someText = NULL;
    TwAddVarRW(bar, "Input", TW_TYPE_CDSTRING, &someText, 
               " label='Text input' group=CDString help=`The text to be copied to 'Text output'.` ");
    TwAddVarRO(bar, "Output", TW_TYPE_CDSTRING, &someText, 
               " label='Text output' group=CDString help=`Carbon copy of the text entered in 'Text input'.` ");
    
    // Add a line of text (we will use the label of a inactive button)
    #define TEXTLINE "a line of text"
    TwAddButton(bar, "Echo", NULL, NULL, 
                " label=`" TEXTLINE "` group=CDString help='Echo of the text entered in the next field' ");
    
    // Add a CDString variable accessed through callbacks
    char *textLine = (char *)malloc(sizeof(TEXTLINE)+1);
    strncpy(textLine, TEXTLINE, sizeof(TEXTLINE));
    TwAddVarCB(bar, "TextLine", TW_TYPE_CDSTRING, SetTextLineCB, GetTextLineCB, &textLine, 
               " label='Change text above' group=CDString help='The text to be echoed.' ");

    // Set the group label & separator
    TwDefine(" Main/CDString label='Echo some text' help='This example demonstates different use of C-Dynamic string variables.' ");
    TwAddSeparator(bar, "Sep2", "");
    TwAddButton(bar, "Blank2", NULL, NULL, " label=' ' ");


    //
    // 3) C-Static string variable example
    //

    TwAddButton(bar, "Info3.1", NULL, NULL, "label='3) This example uses' ");
    TwAddButton(bar, "Info3.2", NULL, NULL, "label='C strings of fixed size' ");

    // Add a CSString
    char tenStr[] = "0123456789"; // 10 characters + null_termination_char -> size = 11
    TwAddVarRW(bar, "Ten", TW_TYPE_CSSTRING(sizeof(tenStr)), tenStr, 
               " label='10 chars max' group=CSString help='A string with a length of 10 characters max.' ");

    // Add a CSString accessed through callbacks. The callbacks will convert the string characters to upper or lower case
    int capCase = 1; // O: lower-case, 1: upper-case
    TwAddVarCB(bar, "Capitalize", TW_TYPE_CSSTRING(sizeof(g_CapStr)), SetCapStrCB, GetCapStrCB, &capCase, 
               " group=CSString help='A string of fixed size to be converted to upper or lower case.' ");

    // Add a bool variable
    TwAddVarRW(bar, "Case", TW_TYPE_BOOL32, &capCase, 
               " false=lower true=UPPER group=CSString key=Space help=`Changes the characters case of the 'Capitalize' string.` ");

    // Set the group label & separator
    TwDefine(" Main/CSString label='Character capitalization' help='This example demonstates different use of C-Static sized variables.' ");
    TwAddSeparator(bar, "Sep3", "");


    // Call the GLUT main loop
    glutMainLoop();

    return 0;
}

