// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifndef __EMBREE_PARMS_H__
#define __EMBREE_PARMS_H__

#include "variant.h"

#include <map>

namespace embree
{
  /*! Parameter container. Implements parameter container as a mapping
   *  from a string to variant values. This container is used to pass
   *  parameters for constructing objects from the API to the
   *  constructors of that objects. All the extraction functions
   *  return a default values in case the parameter is not found. */
  class Parms
  {
  public:

    /*! clears the parameter container */
    void clear() {
      m.clear();
    }

    /*! Extracts a named boolean out of the container. */
    bool getBool(const char* name, bool def = false) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::BOOL1) return def;
      return (*i).second.getBool();
    }

    /*! Extracts a named integer out of the container. */
    int getInt(const char* name, int def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::INT1) return def;
      return (*i).second.getInt();
    }

    /*! Extracts a named float out of the container. */
    float getFloat(const char* name, float def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT1) return def;
      return (*i).second.getFloat();
    }

    /*! Extracts a named Vec2f out of the container. */
    Vec2f getVec2f(const char* name, const Vec2f& def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT2) return def;
      return (*i).second.getVec2f();
    }

    /*! Extracts a named Vector3f out of the container. */
    Vector3f getVector3f(const char* name, const Vector3f& def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT3) return def;
      return (*i).second.getVector3f();
    }

    /*! Extracts a named color out of the container. */
    Color getColor(const char* name, const Color& def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT3) return def;
      return (*i).second.getColor();
    }

    /*! Extracts a named string out of the container. */
    std::string getString(const char* name, std::string def = "") const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::STRING) return def;
      return (*i).second.getString();
    }

    /*! Extracts a named image reference out of the container. */
    Ref<Image> getImage(const char* name, Ref<Image> def = null) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::IMAGE) return def;
      return (*i).second.getImage();
    }

    /*! Extracts a named texture reference out of the container. */
    Ref<Texture> getTexture(const char* name, Ref<Texture> def = null) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::TEXTURE) return def;
      return (*i).second.getTexture();
    }

    /*! Extracts a named transformation out of the container. */
    AffineSpace3f getTransform(const char* name, const AffineSpace3f& def = one) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::TRANSFORM) return def;
      return (*i).second.getTransform();
    }

    /*! Extracts a named data stream out of the container. */
    Variant getData(const char* name) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end()) return Variant();
      return (*i).second;
    }

    /*! Adds a new named element to the container. */
    void add(const std::string& name, Variant data) {
      m[name] = data;
    }

  private:

    /*! Implementation of the container as an STL map. */
    std::map<std::string,Variant> m;
  };
}

#endif
