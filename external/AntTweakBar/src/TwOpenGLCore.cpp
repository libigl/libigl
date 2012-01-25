//  ---------------------------------------------------------------------------
//
//  @file       TwOpenGLCore.cpp
//  @author     Philippe Decaudin - http://www.antisphere.com
//  @license    This file is part of the AntTweakBar library.
//              For conditions of distribution and use, see License.txt
//
//  note:       Work In Progress
//
//  ---------------------------------------------------------------------------

#pragma warning GL3
#define GL3_PROTOTYPES 1 ////
#include <GL3/gl3.h> ////
#define ANT_OGL_HEADER_INCLUDED ////
#include "TwPrecomp.h"
#include "LoadOGLCore.h"
#include "TwOpenGLCore.h"
#include "TwMgr.h"

using namespace std;

extern const char *g_ErrCantLoadOGL;
extern const char *g_ErrCantUnloadOGL;

//  ---------------------------------------------------------------------------

#ifdef _DEBUG
    static void CheckGLCoreError(const char *file, int line, const char *func)
    {
        int err=0;
        char msg[256];
        while( (err=_glGetError())!=0 )
        {
            sprintf(msg, "%s(%d) : [%s] GL_CORE_ERROR=0x%x\n", file, line, func, err);
            #ifdef ANT_WINDOWS
                OutputDebugString(msg);
            #endif
            fprintf(stderr, msg);
        }
    }
#   ifdef __FUNCTION__
#       define CHECK_GL_ERROR CheckGLCoreError(__FILE__, __LINE__, __FUNCTION__)
#   else
#       define CHECK_GL_ERROR CheckGLCoreError(__FILE__, __LINE__, "")
#   endif
#else
#   define CHECK_GL_ERROR ((void)(0))
#endif

//  ---------------------------------------------------------------------------

static GLuint BindFont(const CTexFont *_Font)
{
    GLuint TexID = 0;
/*
    _glGenTextures(1, &TexID);
    _glBindTexture(GL_TEXTURE_2D, TexID);
    _glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
    _glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
    _glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    _glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    _glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    _glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    _glPixelTransferf(GL_ALPHA_SCALE, 1);
    _glPixelTransferf(GL_ALPHA_BIAS, 0);
    _glPixelTransferf(GL_RED_BIAS, 1);
    _glPixelTransferf(GL_GREEN_BIAS, 1);
    _glPixelTransferf(GL_BLUE_BIAS, 1);
    _glTexImage2D(GL_TEXTURE_2D, 0, 4, _Font->m_TexWidth, _Font->m_TexHeight, 0, GL_ALPHA, GL_UNSIGNED_BYTE, _Font->m_TexBytes);
    _glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    _glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    _glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    _glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    _glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    _glBindTexture(GL_TEXTURE_2D, 0);
    _glPixelTransferf(GL_ALPHA_BIAS, 0);
    _glPixelTransferf(GL_RED_BIAS, 0);
    _glPixelTransferf(GL_GREEN_BIAS, 0);
    _glPixelTransferf(GL_BLUE_BIAS, 0);
*/
    return TexID;
}

static void UnbindFont(GLuint _FontTexID)
{
/*
    if( _FontTexID>0 )
        _glDeleteTextures(1, &_FontTexID);
*/
}

//  ---------------------------------------------------------------------------

int CTwGraphOpenGLCore::Init()
{
    m_Drawing = false;
    m_FontTexID = 0;
    m_FontTex = NULL;

    if( LoadOpenGLCore()==0 )
    {
        g_TwMgr->SetLastError(g_ErrCantLoadOGL);
        return 0;
    }

/*
    m_MaxClipPlanes = -1;

    // Get extensions
    _glBindBufferARB = reinterpret_cast<PFNGLBindBufferARB>(_glGetProcAddress("glBindBufferARB"));
    _glBindProgramARB = reinterpret_cast<PFNGLBindProgramARB>(_glGetProcAddress("glBindProgramARB"));
    _glGetHandleARB = reinterpret_cast<PFNGLGetHandleARB>(_glGetProcAddress("glGetHandleARB"));
    _glUseProgramObjectARB = reinterpret_cast<PFNGLUseProgramObjectARB>(_glGetProcAddress("glUseProgramObjectARB"));
    _glTexImage3D = reinterpret_cast<PFNGLTexImage3D>(_glGetProcAddress("glTexImage3D"));
    _glActiveTextureARB = reinterpret_cast<PFNGLActiveTextureARB>(_glGetProcAddress("glActiveTextureARB"));
    _glClientActiveTextureARB = reinterpret_cast<PFNGLClientActiveTextureARB>(_glGetProcAddress("glClientActiveTextureARB"));
    _glBlendEquation = reinterpret_cast<PFNGLBlendEquation>(_glGetProcAddress("glBlendEquation"));
    _glBlendEquationSeparate = reinterpret_cast<PFNGLBlendEquationSeparate>(_glGetProcAddress("glBlendEquationSeparate"));
    _glBlendFuncSeparate = reinterpret_cast<PFNGLBlendFuncSeparate>(_glGetProcAddress("glBlendFuncSeparate"));

#if !defined(ANT_OSX)
    const char *ext = (const char *)_glGetString(GL_EXTENSIONS);
    if( ext!=0 && strlen(ext)>0 )
        m_SupportTexRect = (strstr(ext, "GL_ARB_texture_rectangle")!=NULL);
    else
#endif
        m_SupportTexRect = false;
*/

    // Create shaders
    const GLchar *lineRectVS[] = {
        "#version 150 core\n"
        "in vec3 vertex;"
        "void main() { gl_Position = vec4(vertex, 1); }"
    };
    m_LineRectVS = _glCreateShader(GL_VERTEX_SHADER);
    _glShaderSource(m_LineRectVS, 1, lineRectVS, NULL);
    _glCompileShader(m_LineRectVS);

    const GLchar *lineRectFS[] = {
        "#version 150 core\n"
        "out vec4 color;"
        "void main() { color = vec4(1, 0, 1, 1); }"
    };
    m_LineRectFS = _glCreateShader(GL_FRAGMENT_SHADER);
    _glShaderSource(m_LineRectFS, 1, lineRectFS, NULL);
    _glCompileShader(m_LineRectFS);

    m_LineRectProgram = _glCreateProgram();
    _glAttachShader(m_LineRectProgram, m_LineRectVS);
    _glAttachShader(m_LineRectProgram, m_LineRectFS);
    _glLinkProgram(m_LineRectProgram);

    // Create line/rect vertex buffer
    const GLfloat lineRectInitBuffer[] = { 0,0,0, 0,0,0, 0,0,0, 0,0,0 };
    _glGenVertexArrays(1, &m_LineRectVArray);
    _glBindVertexArray(m_LineRectVArray);
    _glGenBuffers(1, &m_LineRectBuffer);
    _glBindBuffer(GL_ARRAY_BUFFER, m_LineRectBuffer);
    _glBufferData(GL_ARRAY_BUFFER, sizeof(lineRectInitBuffer), lineRectInitBuffer, GL_DYNAMIC_DRAW);

    CHECK_GL_ERROR;
    return 1;
}

//  ---------------------------------------------------------------------------

int CTwGraphOpenGLCore::Shut()
{
    assert(m_Drawing==false);

    UnbindFont(m_FontTexID);

    int Res = 1;
    if( UnloadOpenGLCore()==0 )
    {
        g_TwMgr->SetLastError(g_ErrCantUnloadOGL);
        Res = 0;
    }

    return Res;
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::BeginDraw(int _WndWidth, int _WndHeight)
{
    assert(m_Drawing==false && _WndWidth>0 && _WndHeight>0);
    m_Drawing = true;
    m_WndWidth = _WndWidth;
    m_WndHeight = _WndHeight;
    m_OffsetX = 0;
    m_OffsetY = 0;

    CHECK_GL_ERROR;

    //_glPushAttrib(GL_ALL_ATTRIB_BITS);
    //_glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);
/*
    _glGetIntegerv(GL_ACTIVE_TEXTURE, &m_PrevActiveTexture);
    GLint maxTexUnits = 1;
    _glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &maxTexUnits);
    if( maxTexUnits<1 ) 
        maxTexUnits = 1;
    else if( maxTexUnits > 32 )
        maxTexUnits = 32;
    for( GLint i=0; i<maxTexUnits; ++i )
    {
        _glActiveTexture(GL_TEXTURE0+i);
        m_PrevActiveTexture1D[i] = _glIsEnabled(GL_TEXTURE_1D);
        m_PrevActiveTexture2D[i] = _glIsEnabled(GL_TEXTURE_2D);
        m_PrevActiveTexture3D[i] = _glIsEnabled(GL_TEXTURE_3D);
        _glDisable(GL_TEXTURE_1D);
        _glDisable(GL_TEXTURE_2D);
        _glDisable(GL_TEXTURE_3D);
    }
    _glActiveTexture(GL_TEXTURE0);
    CHECK_GL_ERROR;
*/
    GLint Vp[4];
    _glGetIntegerv(GL_VIEWPORT, Vp);

    if( _WndWidth>0 && _WndHeight>0 )
    {
        Vp[0] = 0;
        Vp[1] = 0;
        Vp[2] = _WndWidth-1;
        Vp[3] = _WndHeight-1;
        _glViewport(Vp[0], Vp[1], Vp[2], Vp[3]);
    }
    _glGetIntegerv(GL_VIEWPORT, m_ViewportInit);

    _glGetFloatv(GL_LINE_WIDTH, &m_PrevLineWidth);
    _glLineWidth(1);
    _glDisable(GL_LINE_SMOOTH);
    _glDisable(GL_CULL_FACE);
    _glDisable(GL_DEPTH_TEST);
    _glEnable(GL_BLEND);
    _glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    _glDisable(GL_SCISSOR_TEST);
    m_PrevTexture = 0;
    _glGetIntegerv(GL_TEXTURE_BINDING_2D, &m_PrevTexture);

/*
    _glDisableClientState(GL_VERTEX_ARRAY);
    _glDisableClientState(GL_NORMAL_ARRAY);
    _glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    _glDisableClientState(GL_INDEX_ARRAY);
    _glDisableClientState(GL_COLOR_ARRAY);
    _glDisableClientState(GL_EDGE_FLAG_ARRAY);

    if( _glBindBuffer!=NULL )
    {
        m_PrevArrayBufferARB = m_PrevElementArrayBufferARB = 0;
        _glGetIntegerv(GL_ARRAY_BUFFER_BINDING_ARB, &m_PrevArrayBufferARB);
        _glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING_ARB, &m_PrevElementArrayBufferARB);
        _glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
        _glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
    }
    if( _glBindProgramARB!=NULL )
    {
        m_PrevVertexProgramARB = _glIsEnabled(GL_VERTEX_PROGRAM_ARB);
        m_PrevFragmentProgramARB = _glIsEnabled(GL_FRAGMENT_PROGRAM_ARB);
        _glDisable(GL_VERTEX_PROGRAM_ARB);
        _glDisable(GL_FRAGMENT_PROGRAM_ARB);
    }
    if( _glGetHandleARB!=NULL && _glUseProgramObjectARB!=NULL )
    {
        m_PrevProgramObjectARB = _glGetHandleARB(GL_PROGRAM_OBJECT_ARB);
        _glUseProgramObjectARB(0);
    }
*/
/*
    _glDisable(GL_TEXTURE_1D);
    _glDisable(GL_TEXTURE_2D);
    m_PrevTexture3D = _glIsEnabled(GL_TEXTURE_3D);
    _glDisable(GL_TEXTURE_3D);
    m_PrevTexRect = _glIsEnabled(GL_TEXTURE_RECTANGLE);
    _glDisable(GL_TEXTURE_RECTANGLE);
    if( _glBlendEquationSeparate!=NULL )
    {
        _glGetIntegerv(GL_BLEND_EQUATION_RGB, &m_PrevBlendEquationRGB);
        _glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &m_PrevBlendEquationAlpha);
        _glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    }
    if( _glBlendFuncSeparate!=NULL )
    {
        _glGetIntegerv(GL_BLEND_SRC_RGB, &m_PrevBlendSrcRGB);
        _glGetIntegerv(GL_BLEND_DST_RGB, &m_PrevBlendDstRGB);
        _glGetIntegerv(GL_BLEND_SRC_ALPHA, &m_PrevBlendSrcAlpha);
        _glGetIntegerv(GL_BLEND_DST_ALPHA, &m_PrevBlendDstAlpha);
        _glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    if( _glBlendEquation!=NULL )
    {
        _glGetIntegerv(GL_BLEND_EQUATION, &m_PrevBlendEquation);
        _glBlendEquation(GL_FUNC_ADD);
    }
*/
    CHECK_GL_ERROR;

//    _glUseProgram(m_LineRectProgram);
//    GLint projLoc = _glGetUniformLocation(m_LineRectProgram, "proj");
//    _glUniformMatrix4fv(projLoc, 1, false, proj);

    CHECK_GL_ERROR;
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::EndDraw()
{
    assert(m_Drawing==true);
    m_Drawing = false;
/*
    _glBindTexture(GL_TEXTURE_2D, m_PrevTexture);
    if( _glBindBufferARB!=NULL )
    {
        _glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_PrevArrayBufferARB);
        _glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, m_PrevElementArrayBufferARB);
    }
    if( _glBindProgramARB!=NULL )
    {
        if( m_PrevVertexProgramARB )
            _glEnable(GL_VERTEX_PROGRAM_ARB);
        if( m_PrevFragmentProgramARB )
            _glEnable(GL_FRAGMENT_PROGRAM_ARB);
    }
    if( _glGetHandleARB!=NULL && _glUseProgramObjectARB!=NULL )
        _glUseProgramObjectARB(m_PrevProgramObjectARB);
    if( _glTexImage3D!=NULL && m_PrevTexture3D )
        _glEnable(GL_TEXTURE_3D);
    if( m_SupportTexRect && m_PrevTexRectARB )
        _glEnable(GL_TEXTURE_RECTANGLE_ARB);
    if( _glBlendEquation!=NULL )
        _glBlendEquation(m_PrevBlendEquation);
    if( _glBlendEquationSeparate!=NULL )
        _glBlendEquationSeparate(m_PrevBlendEquationRGB, m_PrevBlendEquationAlpha);
    if( _glBlendFuncSeparate!=NULL )
        _glBlendFuncSeparate(m_PrevBlendSrcRGB, m_PrevBlendDstRGB, m_PrevBlendSrcAlpha, m_PrevBlendDstAlpha);
    
    _glPolygonMode(GL_FRONT, m_PrevPolygonMode[0]);
    _glPolygonMode(GL_BACK, m_PrevPolygonMode[1]);
    _glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, m_PrevTexEnv);
    _glLineWidth(m_PrevLineWidth);
    _glMatrixMode(GL_PROJECTION);
    _glPopMatrix();
    _glMatrixMode(GL_MODELVIEW);
    _glPopMatrix();
    _glMatrixMode(GL_TEXTURE);
    _glPopMatrix();
    _glPopClientAttrib();
    _glPopAttrib();

    if( _glActiveTextureARB )
    {
        GLint maxTexUnits = 1;
        _glGetIntegerv(GL_MAX_TEXTURE_UNITS_ARB, &maxTexUnits);
        if( maxTexUnits<1 ) 
            maxTexUnits = 1;
        else if( maxTexUnits > 32 )
            maxTexUnits = 32;
        for( GLint i=0; i<maxTexUnits; ++i )
        {
            _glActiveTextureARB(GL_TEXTURE0_ARB+i);
            if( m_PrevActiveTexture1D[i] )
                _glEnable(GL_TEXTURE_1D);
            if( m_PrevActiveTexture2D[i] )
                _glEnable(GL_TEXTURE_2D);
            if( m_PrevActiveTexture3D[i] )
                _glEnable(GL_TEXTURE_3D);
        }
        _glActiveTextureARB(m_PrevActiveTextureARB);
    }
*/
    CHECK_GL_ERROR;
}

//  ---------------------------------------------------------------------------

bool CTwGraphOpenGLCore::IsDrawing()
{
    return m_Drawing;
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::Restore()
{
    UnbindFont(m_FontTexID);
    m_FontTexID = 0;
    m_FontTex = NULL;
}

//  ---------------------------------------------------------------------------

static inline float ToNormScreenX(int x, int wndWidth)
{
    return 2.0f*((float)x-0.5f)/wndWidth - 1.0f;
}

static inline float ToNormScreenY(int y, int wndHeight)
{
    return 1.0f - 2.0f*((float)y-0.5f)/wndHeight;
}

static inline color32 ToR8G8B8A8(color32 col)
{
    return (col & 0xff00ff00) | ((col>>16) & 0xff) | ((col<<16) & 0xff0000);
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::DrawLine(int _X0, int _Y0, int _X1, int _Y1, color32 _Color0, color32 _Color1, bool _AntiAliased)
{
    assert(m_Drawing==true);
    /*
    //const GLfloat dx = +0.0f;
    const GLfloat dx = +0.5f;
    //GLfloat dy = -0.2f;
    const GLfloat dy = -0.5f;
    if( _AntiAliased )
        _glEnable(GL_LINE_SMOOTH);
    else
        _glDisable(GL_LINE_SMOOTH);
    _glDisable(GL_TEXTURE_2D);
    _glMatrixMode(GL_MODELVIEW);
    _glLoadIdentity();
    _glBegin(GL_LINES);
        _glColor4ub(GLubyte(_Color0>>16), GLubyte(_Color0>>8), GLubyte(_Color0), GLubyte(_Color0>>24));
        _glVertex2f((GLfloat)_X0+dx, (GLfloat)_Y0+dy);
        _glColor4ub(GLubyte(_Color1>>16), GLubyte(_Color1>>8), GLubyte(_Color1), GLubyte(_Color1>>24));
        _glVertex2f((GLfloat)_X1+dx, (GLfloat)_Y1+dy);
        //_glVertex2i(_X0, _Y0);
        //_glVertex2i(_X1, _Y1);
    _glEnd();
    _glDisable(GL_LINE_SMOOTH);
    */

    GLfloat x0 = ToNormScreenX(_X0 + m_OffsetX, m_WndWidth);
    GLfloat y0 = ToNormScreenY(_Y0 + m_OffsetY, m_WndHeight);
    GLfloat x1 = ToNormScreenX(_X1 + m_OffsetX, m_WndWidth);
    GLfloat y1 = ToNormScreenY(_Y1 + m_OffsetY, m_WndHeight);
    GLfloat vertices[] = { x0, y0, 0,  x1, y1, 0 };
    _glBindBuffer(GL_ARRAY_BUFFER, m_LineRectBuffer);
    _glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);

    _glUseProgram(m_LineRectProgram);
    _glBindVertexArray(m_LineRectVArray);
    GLint vlocation = _glGetAttribLocation(m_LineRectProgram, "vertex");
    _glVertexAttribPointer(vlocation, 3, GL_FLOAT, GL_TRUE, 0, NULL);
    _glEnableVertexAttribArray(vlocation);
    _glDrawArrays(GL_LINES, 0, 2);

    CHECK_GL_ERROR;
}
  
//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::DrawRect(int _X0, int _Y0, int _X1, int _Y1, color32 _Color00, color32 _Color10, color32 _Color01, color32 _Color11)
{
    assert(m_Drawing==true);

    /*
    // border adjustment
    if(_X0<_X1)
        ++_X1;
    else if(_X0>_X1)
        ++_X0;
    if(_Y0<_Y1)
        --_Y0;
    else if(_Y0>_Y1)
        --_Y1;
    const GLfloat dx = +0.0f;
    const GLfloat dy = +0.0f;

    _glDisable(GL_TEXTURE_2D);
    _glMatrixMode(GL_MODELVIEW);
    _glLoadIdentity();
    //GLubyte r = GLubyte(_Color>>16);
    //GLubyte g = GLubyte(_Color>>8);
    //GLubyte b = GLubyte(_Color);
    //GLubyte a = GLubyte(_Color>>24);
    //_glColor4ub(GLubyte(_Color>>16), GLubyte(_Color>>8), GLubyte(_Color), GLubyte(_Color>>24));
    //_glColor4ub(r, g, b, a);
    _glBegin(GL_QUADS);
        _glColor4ub(GLubyte(_Color00>>16), GLubyte(_Color00>>8), GLubyte(_Color00), GLubyte(_Color00>>24));
        _glVertex2f((GLfloat)_X0+dx, (GLfloat)_Y0+dy);
        _glColor4ub(GLubyte(_Color10>>16), GLubyte(_Color10>>8), GLubyte(_Color10), GLubyte(_Color10>>24));
        _glVertex2f((GLfloat)_X1+dx, (GLfloat)_Y0+dy);
        _glColor4ub(GLubyte(_Color11>>16), GLubyte(_Color11>>8), GLubyte(_Color11), GLubyte(_Color11>>24));
        _glVertex2f((GLfloat)_X1+dx, (GLfloat)_Y1+dy);
        _glColor4ub(GLubyte(_Color01>>16), GLubyte(_Color01>>8), GLubyte(_Color01), GLubyte(_Color01>>24));
        _glVertex2f((GLfloat)_X0+dx, (GLfloat)_Y1+dy);
    _glEnd();
    */
}

//  ---------------------------------------------------------------------------

void *CTwGraphOpenGLCore::NewTextObj()
{
    return new CTextObj;
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::DeleteTextObj(void *_TextObj)
{
    assert(_TextObj!=NULL);
    delete static_cast<CTextObj *>(_TextObj);
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::BuildText(void *_TextObj, const std::string *_TextLines, color32 *_LineColors, color32 *_LineBgColors, int _NbLines, const CTexFont *_Font, int _Sep, int _BgWidth)
{
    assert(m_Drawing==true);
    assert(_TextObj!=NULL);
    assert(_Font!=NULL);

    if( _Font != m_FontTex )
    {
        UnbindFont(m_FontTexID);
        m_FontTexID = BindFont(_Font);
        m_FontTex = _Font;
    }
    CTextObj *TextObj = static_cast<CTextObj *>(_TextObj);
    TextObj->m_TextVerts.resize(0);
    TextObj->m_TextUVs.resize(0);
    TextObj->m_BgVerts.resize(0);
    TextObj->m_Colors.resize(0);
    TextObj->m_BgColors.resize(0);

    int x, x1, y, y1, i, Len;
    unsigned char ch;
    const unsigned char *Text;
    color32 LineColor = COLOR32_RED;
    for( int Line=0; Line<_NbLines; ++Line )
    {
        x = 0;
        y = Line * (_Font->m_CharHeight+_Sep);
        y1 = y+_Font->m_CharHeight;
        Len = (int)_TextLines[Line].length();
        Text = (const unsigned char *)(_TextLines[Line].c_str());
        if( _LineColors!=NULL )
            LineColor = (_LineColors[Line]&0xff00ff00) | GLubyte(_LineColors[Line]>>16) | (GLubyte(_LineColors[Line])<<16);

        for( i=0; i<Len; ++i )
        {
            ch = Text[i];
            x1 = x + _Font->m_CharWidth[ch];

            TextObj->m_TextVerts.push_back(Vec2(x , y ));
            TextObj->m_TextVerts.push_back(Vec2(x1, y ));
            TextObj->m_TextVerts.push_back(Vec2(x , y1));
            TextObj->m_TextVerts.push_back(Vec2(x1, y ));
            TextObj->m_TextVerts.push_back(Vec2(x1, y1));
            TextObj->m_TextVerts.push_back(Vec2(x , y1));

            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU0[ch], _Font->m_CharV0[ch]));
            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU1[ch], _Font->m_CharV0[ch]));
            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU0[ch], _Font->m_CharV1[ch]));
            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU1[ch], _Font->m_CharV0[ch]));
            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU1[ch], _Font->m_CharV1[ch]));
            TextObj->m_TextUVs.push_back(Vec2(_Font->m_CharU0[ch], _Font->m_CharV1[ch]));

            if( _LineColors!=NULL )
            {
                TextObj->m_Colors.push_back(LineColor);
                TextObj->m_Colors.push_back(LineColor);
                TextObj->m_Colors.push_back(LineColor);
                TextObj->m_Colors.push_back(LineColor);
                TextObj->m_Colors.push_back(LineColor);
                TextObj->m_Colors.push_back(LineColor);
            }

            x = x1;
        }
        if( _BgWidth>0 )
        {
            TextObj->m_BgVerts.push_back(Vec2(-1        , y ));
            TextObj->m_BgVerts.push_back(Vec2(_BgWidth+1, y ));
            TextObj->m_BgVerts.push_back(Vec2(-1        , y1));
            TextObj->m_BgVerts.push_back(Vec2(_BgWidth+1, y ));
            TextObj->m_BgVerts.push_back(Vec2(_BgWidth+1, y1));
            TextObj->m_BgVerts.push_back(Vec2(-1        , y1));

            if( _LineBgColors!=NULL )
            {
                color32 LineBgColor = (_LineBgColors[Line]&0xff00ff00) | GLubyte(_LineBgColors[Line]>>16) | (GLubyte(_LineBgColors[Line])<<16);
                TextObj->m_BgColors.push_back(LineBgColor);
                TextObj->m_BgColors.push_back(LineBgColor);
                TextObj->m_BgColors.push_back(LineBgColor);
                TextObj->m_BgColors.push_back(LineBgColor);
                TextObj->m_BgColors.push_back(LineBgColor);
                TextObj->m_BgColors.push_back(LineBgColor);
            }
        }
    }
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::DrawText(void *_TextObj, int _X, int _Y, color32 _Color, color32 _BgColor)
{
    assert(m_Drawing==true);
    assert(_TextObj!=NULL);
    CTextObj *TextObj = static_cast<CTextObj *>(_TextObj);

    if( TextObj->m_TextVerts.size()<4 && TextObj->m_BgVerts.size()<4 )
        return; // nothing to draw
/*
    _glMatrixMode(GL_MODELVIEW);
    _glLoadIdentity();
    _glTranslatef((GLfloat)_X, (GLfloat)_Y, 0);
    _glEnableClientState(GL_VERTEX_ARRAY);
    if( (_BgColor!=0 || TextObj->m_BgColors.size()==TextObj->m_BgVerts.size()) && TextObj->m_BgVerts.size()>=4 )
    {
        _glDisable(GL_TEXTURE_2D);
        _glVertexPointer(2, GL_FLOAT, 0, &(TextObj->m_BgVerts[0]));
        if( TextObj->m_BgColors.size()==TextObj->m_BgVerts.size() && _BgColor==0 )
        {
            _glEnableClientState(GL_COLOR_ARRAY);
            _glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(TextObj->m_BgColors[0]));
        }
        else
        {
            _glDisableClientState(GL_COLOR_ARRAY);
            _glColor4ub(GLubyte(_BgColor>>16), GLubyte(_BgColor>>8), GLubyte(_BgColor), GLubyte(_BgColor>>24));
        }
        _glDrawArrays(GL_TRIANGLES, 0, (int)TextObj->m_BgVerts.size());
    }
    _glEnable(GL_TEXTURE_2D);
    _glBindTexture(GL_TEXTURE_2D, m_FontTexID);
    _glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    if( TextObj->m_TextVerts.size()>=4 )
    {
        _glVertexPointer(2, GL_FLOAT, 0, &(TextObj->m_TextVerts[0]));
        _glTexCoordPointer(2, GL_FLOAT, 0, &(TextObj->m_TextUVs[0]));
        if( TextObj->m_Colors.size()==TextObj->m_TextVerts.size() && _Color==0 )
        {
            _glEnableClientState(GL_COLOR_ARRAY);
            _glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(TextObj->m_Colors[0]));
        }
        else
        {
            _glDisableClientState(GL_COLOR_ARRAY);
            _glColor4ub(GLubyte(_Color>>16), GLubyte(_Color>>8), GLubyte(_Color), GLubyte(_Color>>24));
        }

        _glDrawArrays(GL_TRIANGLES, 0, (int)TextObj->m_TextVerts.size());
    }
    
    _glDisableClientState(GL_VERTEX_ARRAY);
    _glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    _glDisableClientState(GL_COLOR_ARRAY);
*/
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::ChangeViewport(int _X0, int _Y0, int _Width, int _Height, int _OffsetX, int _OffsetY)
{
/*
    if( _Width>0 && _Height>0 )
    {
        GLint vp[4];
        vp[0] = _X0;
        vp[1] = _Y0;
        vp[2] = _Width-1;
        vp[3] = _Height-1;
        _glViewport(vp[0], m_WndHeight-vp[1]-vp[3], vp[2], vp[3]);

        GLint matrixMode = 0;
        _glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
        _glMatrixMode(GL_PROJECTION);
        _glLoadIdentity();
        _glOrtho(_OffsetX, _OffsetX+vp[2], vp[3]-_OffsetY, -_OffsetY, -1, 1);
        _glMatrixMode(matrixMode);
    }
*/
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::RestoreViewport()
{
/*
    _glViewport(m_ViewportInit[0], m_ViewportInit[1], m_ViewportInit[2], m_ViewportInit[3]);

    GLint matrixMode = 0;
    _glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
    _glMatrixMode(GL_PROJECTION);
    _glLoadMatrixf(m_ProjMatrixInit);
    _glMatrixMode(matrixMode);
*/
}

//  ---------------------------------------------------------------------------

void CTwGraphOpenGLCore::DrawTriangles(int _NumTriangles, int *_Vertices, color32 *_Colors, Cull _CullMode)
{
    assert(m_Drawing==true);

    const GLfloat dx = +0.0f;
    const GLfloat dy = +0.0f;
/*
    GLint prevCullFaceMode, prevFrontFace;
    _glGetIntegerv(GL_CULL_FACE_MODE, &prevCullFaceMode);
    _glGetIntegerv(GL_FRONT_FACE, &prevFrontFace);
    GLboolean prevCullEnable = _glIsEnabled(GL_CULL_FACE);
    _glCullFace(GL_BACK);
    _glEnable(GL_CULL_FACE);
    if( _CullMode==CULL_CW )
        _glFrontFace(GL_CCW);
    else if( _CullMode==CULL_CCW )
        _glFrontFace(GL_CW);
    else
        _glDisable(GL_CULL_FACE);

    _glDisable(GL_TEXTURE_2D);
    _glMatrixMode(GL_MODELVIEW);
    _glLoadIdentity();
    _glBegin(GL_TRIANGLES);
    for(int i=0; i<3*_NumTriangles; ++i)
    {
        color32 col = _Colors[i];
        _glColor4ub(GLubyte(col>>16), GLubyte(col>>8), GLubyte(col), GLubyte(col>>24));
        _glVertex2f((GLfloat)_Vertices[2*i+0]+dx, (GLfloat)_Vertices[2*i+1]+dy);
    }
    _glEnd();

    _glCullFace(prevCullFaceMode);
    _glFrontFace(prevFrontFace);
    if( prevCullEnable )
        _glEnable(GL_CULL_FACE);
    else
        _glDisable(GL_CULL_FACE);
*/
}

//  ---------------------------------------------------------------------------
