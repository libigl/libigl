#ifndef GET_MODIFIERS_H
#define GET_MODIFIERS_H
//#include "igl_inline.h"
namespace igl
{
  enum Modifier
  {
    MODIFIER_OPTION = 1,
    MODIFIER_SHIFT = 2,
    MODIFIER_CONTROL = 4,
    MODIFIER_COMMAND = 8,
    NUM_MODIFIERS = 4,
  };
  // Retrieve current modifier constellation. 
  //
  // Returns int that's an "or" of the active modifiers above.
  //
  // FORCED INLINE
  inline int get_modifiers();
}

// Implementation 

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
//#include <Carbon/HIToolbox/Events.h>
#include <Carbon/Carbon.h>
#endif

#warning "igl::get_modifiers is deprecated. If using GLUT, try Alec's glut patch www.alecjacobson.com/weblog/?p=3659 and use glutGetModifiers"

// FORCED INLINE
inline int igl::get_modifiers()
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
#  warning "igl::get_modifiers not supported on your OS, some demos may not work correctly."
#endif
  return mod;
}

#endif
