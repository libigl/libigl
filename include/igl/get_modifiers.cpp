#include "get_modifiers.h"

/* glutGetModifiers return mask. */
#ifndef GLUT_ACTIVE_SHIFT
#  define GLUT_ACTIVE_SHIFT 1
#endif
#ifndef GLUT_ACTIVE_CTRL
#  define GLUT_ACTIVE_CTRL 2
#endif
#ifndef GLUT_ACTIVE_ALT
#  define GLUT_ACTIVE_ALT 4
#endif
#ifndef GLUT_ACTIVE_COMMAND
#  define GLUT_ACTIVE_COMMAND 8
#endif

#ifdef __APPLE__
#include <Carbon/Carbon.h>
#endif

IGL_INLINE int igl::get_modifiers()
{
  int mod = 0;
#ifdef __APPLE__
  // http://stackoverflow.com/a/18082326/148668
  KeyMap keyStates;
  const auto & carbon_is_keydown = [&keyStates]( uint16_t vKey )->bool
  {
    uint8_t index = vKey / 32 ;
    uint8_t shift = vKey % 32 ;
    return keyStates[index].bigEndianValue & (1 << shift) ;
  };
  GetKeys(keyStates) ;
  mod |= (carbon_is_keydown(kVK_Command)?GLUT_ACTIVE_COMMAND:0);
  mod |= (carbon_is_keydown(kVK_Shift)?GLUT_ACTIVE_SHIFT:0);
  mod |= (carbon_is_keydown(kVK_Option)?GLUT_ACTIVE_ALT:0);
  mod |= (carbon_is_keydown(kVK_Control)?GLUT_ACTIVE_CTRL:0);
#else
#  error "Not supported."
#endif
  return mod;
}
