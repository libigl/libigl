#ifndef GET_MODIFIERS_H
#define GET_MODIFIERS_H
#include "igl_inline.h"
namespace igl
{
  enum Modifier
  {
    MODIFIER_OPTION = 1,
    MODIFIER_SHIFT = 3,
    MODIFIER_CONTROL = 5,
    MODIFIER_COMMAND = 9,
    NUM_MODIFIERS = 4,
  };
  // Retrieve current modifier constellation. 
  //
  // Returns int that's an "or" of the active modifiers above.
  IGL_INLINE int get_modifiers();
}
#ifndef IGL_STATIC_LIBRARY
#include "get_modifiers.cpp"
#endif
#endif
