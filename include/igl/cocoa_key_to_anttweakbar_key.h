#ifndef IGL_COCOA_KEY_TO_ANTTWEAKBAR_KEY_H
#define IGL_COCOA_KEY_TO_ANTTWEAKBAR_KEY_H
#ifndef IGL_NO_ANTTWEAKBAR
#include "igl_inline.h"


namespace igl
{
  // Convert an unsigned char (like that from Cocoa apps) to AntTweakBar key
  // code.
  // See also: TranslateKey() in TwMgr.cpp in AntTweakBar source
  // Inputs:
  //   key  unsigned char key from keyboard
  // Returns int of new key code 
  IGL_INLINE int cocoa_key_to_anttweakbar_key(int key);
}

#ifdef IGL_HEADER_ONLY
#  include "cocoa_key_to_anttweakbar_key.cpp"
#endif

#endif
#endif
