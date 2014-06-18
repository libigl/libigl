//  ---------------------------------------------------------------------------
//
//  @file       TwFonts.h
//  @brief      Bitmaps fonts
//  @author     Philippe Decaudin
//  @license    This file is part of the AntTweakBar library.
//              For conditions of distribution and use, see License.txt
//
//  note:       Private header
//
//  ---------------------------------------------------------------------------


#if !defined ANT_TW_FONTS_INCLUDED
#define ANT_TW_FONTS_INCLUDED

//#include <AntTweakBar.h>

/*
A source bitmap includes 224 characters starting from ascii char 32 (i.e. space) 
to ascii char 255 (extended ASCII Latin1/CP1252):
  
 !"#$%&'()*+,-./0123456789:;<=>?
@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_
`abcdefghijklmnopqrstuvwxyz{|}~
€‚ƒ„…†‡ˆ‰Š‹Œ‘’“”•–—˜™š›œŸ
 ¡¢£¤¥¦§¨©ª«¬­®¯°±²³´µ¶·¸¹º»¼½¾¿
ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏĞÑÒÓÔÕÖ×ØÙÚÛÜİŞß
àáâãäåæçèéêëìíîïğñòóôõö÷øùúûüışÿ

First pixel column of a source bitmap is a delimiter with color=zero at the end of each line of characters.
Last pixel row of a line of characters is a delimiter with color=zero at the last pixel of each character.

*/

// ----------------------------------------------------------------------------
//  OS specific definitions
// ----------------------------------------------------------------------------

#if (defined(_WIN32) || defined(_WIN64)) && !defined(TW_STATIC)
#   define TW_CALL          __stdcall
#   define TW_CDECL_CALL    __cdecl
#   define TW_EXPORT_API    __declspec(dllexport)
#   define TW_IMPORT_API    __declspec(dllimport)
#else
#   define TW_CALL
#   define TW_CDECL_CALL
#   define TW_EXPORT_API
#   define TW_IMPORT_API
#endif

#if defined TW_EXPORTS
#   define TW_API TW_EXPORT_API
#elif defined TW_STATIC
#   define TW_API
#   if defined(_MSC_VER) && !defined(TW_NO_LIB_PRAGMA)
#       ifdef _WIN64
#           pragma comment(lib, "AntTweakBarStatic64")
#       else
#           pragma comment(lib, "AntTweakBarStatic")
#       endif
#   endif
#else
#   define TW_API TW_IMPORT_API
#   if defined(_MSC_VER) && !defined(TW_NO_LIB_PRAGMA)
#       ifdef _WIN64
#           pragma comment(lib, "AntTweakBar64")
#       else
#           pragma comment(lib, "AntTweakBar")
#       endif
#   endif
#endif


struct CTexFont
{
    unsigned char * m_TexBytes;
    int             m_TexWidth;     // power of 2
    int             m_TexHeight;    // power of 2
    float           m_CharU0[256];
    float           m_CharV0[256];
    float           m_CharU1[256];
    float           m_CharV1[256];
    int             m_CharWidth[256];
    int             m_CharHeight;
    int             m_NbCharRead;

    CTexFont();
    ~CTexFont();
};


CTexFont *TwGenerateFont(const unsigned char *_Bitmap, int _BmWidth, int _BmHeight, float _Scaling=1.0f);


TW_API extern CTexFont *g_DefaultSmallFont;
TW_API extern CTexFont *g_DefaultNormalFont;
TW_API extern CTexFont *g_DefaultLargeFont;
TW_API extern CTexFont *g_DefaultFixed1Font;
TW_API extern CTexFont *g_DefaultFixedRuFont;

void TwGenerateDefaultFonts(float _Scaling=1.0f);
void TwDeleteDefaultFonts();


#endif  // !defined ANT_TW_FONTS_INCLUDED
