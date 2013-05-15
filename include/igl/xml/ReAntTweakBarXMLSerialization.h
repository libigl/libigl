#ifndef IGL_REANTTWEAKBAR_XML_SERIALIZATION_H
#define IGL_REANTTWEAKBAR_XML_SERIALIZATION_H
#include "../igl_inline.h"

namespace igl
{  
  IGL_INLINE bool save_ReAntTweakBar(igl::ReTwBar* bar, const char* file_name);
  IGL_INLINE bool save_ReAntTweakBar(igl::ReTwBar* bar, tinyxml2::XMLDocument* doc);
  IGL_INLINE bool load_ReAntTweakBar(igl::ReTwBar* bar, const char *file_name);
  IGL_INLINE bool load_ReAntTweakBar(igl::ReTwBar* bar, tinyxml2::XMLDocument* doc);
}

#ifdef IGL_HEADER_ONLY
#  include "ReAntTweakBarXMLSerialization.cpp"
#endif

#endif
