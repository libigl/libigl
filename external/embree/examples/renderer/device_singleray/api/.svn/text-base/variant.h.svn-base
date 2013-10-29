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

#ifndef __EMBREE_VARIANT_H__
#define __EMBREE_VARIANT_H__

#include "data.h"
#include "datastream.h"
#include "../common/image/image.h"
#include "../textures/texture.h"

namespace embree
{
  /*! Variant type. We use the variant type to pass values from the
   *  API to the constructors of different objects. */
  struct Variant
  {
    ALIGNED_CLASS
  public:

    /*! Determines which kind of value is stored in the variant. */
    enum Type {
      EMPTY,     /*!< variant is empty                 */
      BOOL1,     /*!< variant stores bool value        */
      BOOL2,     /*!< variant stores bool2 value       */
      BOOL3,     /*!< variant stores bool3 value       */
      BOOL4,     /*!< variant stores bool4 value       */
      INT1,      /*!< variant stores int value         */
      INT2,      /*!< variant stores int2 value        */
      INT3,      /*!< variant stores int3 value        */
      INT4,      /*!< variant stores int4 value        */
      FLOAT1,    /*!< variant stores float value       */
      FLOAT2,    /*!< variant stores float2 value      */
      FLOAT3,    /*!< variant stores float3 value      */
      FLOAT4,    /*!< variant stores float4 value      */
      STRING,    /*!< variant stores string value      */
      IMAGE,     /*!< variant stores image reference   */
      TEXTURE,   /*!< variant stores texture reference */
      TRANSFORM  /*!< variant stores transformation    */
    };

    /*! Constructs an empty variant object. */
    Variant ( ) : type(EMPTY) { }

    /*! Constructs a variant object holding a bool value. */
    Variant (bool b0                           ) : type(BOOL1) { b[0] = b0; }

    /*! Constructs a variant object holding a bool2 value. */
    Variant (bool b0, bool b1                  ) : type(BOOL2) { b[0] = b0; b[1] = b1; }

    /*! Constructs a variant object holding a bool3 value. */
    Variant (bool b0, bool b1, bool b2         ) : type(BOOL3) { b[0] = b0; b[1] = b1; b[2] = b2; }

    /*! Constructs a variant object holding a bool4 value. */
    Variant (bool b0, bool b1, bool b2, bool b3) : type(BOOL4) { b[0] = b0; b[1] = b1; b[2] = b2; b[3] = b3; }

    /*! Constructs a variant object holding an int value. */
    Variant (int i0                            ) : type(INT1) { i[0] = i0; }

    /*! Constructs a variant object holding an int2 value. */
    Variant (int i0, int i1                    ) : type(INT2) { i[0] = i0; i[1] = i1; }

    /*! Constructs a variant object holding an int3 value. */
    Variant (int i0, int i1, int i2            ) : type(INT3) { i[0] = i0; i[1] = i1; i[2] = i2; }

    /*! Constructs a variant object holding an int4 value. */
    Variant (int i0, int i1, int i2, int i3    ) : type(INT4) { i[0] = i0; i[1] = i1; i[2] = i2; i[3] = i3; }

    /*! Constructs a variant object holding a float value. */
    Variant (float f0                              ) : type(FLOAT1) { f[0] = f0; }

    /*! Constructs a variant object holding a float2 value. */
    Variant (float f0, float f1                    ) : type(FLOAT2) { f[0] = f0; f[1] = f1; }

    /*! Constructs a variant object holding a float3 value. */
    Variant (float f0, float f1, float f2          ) : type(FLOAT3) { f[0] = f0; f[1] = f1; f[2] = f2; }

    /*! Constructs a variant object holding a float4 value. */
    Variant (float f0, float f1, float f2, float f3) : type(FLOAT4) { f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3; }

    /*! Constructs a variant object holding a pointer to a typed array. */
    Variant (const Ref<Data>& ptr, Type type, size_t size, size_t stride, size_t ofs) : type(type), data(new DataStream(ptr,size,stride,ofs)) { }

    /*! Constructs a variant object holding a string value. */
    Variant (const char* str) : type(STRING), str(str) {}

    /*! Constructs a variant object holding a string value. */
    Variant (const std::string& str) : type(STRING), str(str) {}

    /*! Constructs a variant object holding a texture reference. */
    Variant (const Ref<Texture>& texture) : type(TEXTURE), texture(texture) {}

    /*! Constructs a variant object holding an image reference. */
    Variant (const Ref<Image>& image) : type(IMAGE), image(image) {}

    /*! Constructs a variant object holding a transformation. */
    Variant (const AffineSpace3f& xfm) : type(TRANSFORM)
    {
      f[0] = xfm.l.vx.x; f[1] = xfm.l.vx.y; f[2] = xfm.l.vx.z;
      f[3] = xfm.l.vy.x; f[4] = xfm.l.vy.y; f[5] = xfm.l.vy.z;
      f[6] = xfm.l.vz.x; f[7] = xfm.l.vz.y; f[8] = xfm.l.vz.z;
      f[9] = xfm.p   .x; f[10]= xfm.p   .y; f[11]= xfm.p   .z;
    }

    /*! Extracts a boolean from the variant type. */
    bool  getBool () const { return b[0]; }

    /*! Extracts an integer from the variant type. */
    int   getInt  () const { return i[0]; }

    /*! Extracts a float from the variant type. */
    float getFloat() const { return f[0]; }

    /*! Extracts a Vec2f from the variant type. */
    Vec2f getVec2f() const { return Vec2f(f[0],f[1]); }

    /*! Extracts a Vector3f from the variant type. */
    Vector3f getVector3f() const { return Vector3f(f[0],f[1],f[2]); }

    /*! Extracts a color from the variant type. */
    Color getColor() const { return Color(f[0],f[1],f[2]); }

    /*! Extracts a string from the variant type. */
    std::string getString() const { return str;   }

    /*! Extracts an image reference from the variant type. */
    Ref<Image> getImage  () const { return image; }

    /*! Extracts a texture reference from the variant type. */
    Ref<Texture> getTexture() const { return texture; }

    /*! Extracts a transformation from the variant type. */
    AffineSpace3f getTransform() const {
      return AffineSpace3f(Vector3f(f[0],f[1],f[2]),Vector3f(f[3],f[4],f[5]),Vector3f(f[6],f[7],f[8]),Vector3f(f[9],f[10],f[11]));
    }

    operator bool() const {
      return type != EMPTY;
    }

  public:
    Type type;            //!< Type of the data contained in the variant.
    union {
      bool b[4];          //!< Storage for single bool,bool2,bool3, and bool4 values.
      int i[4];           //!< Storage for single int,int2,int3, and int4 values.
      float f[12];         //!< Storage for single float,float2,float3, float4, and AffineSpace3f values.
    };
    std::string str;      //!< Storage for string values.
    Ref<DataStream> data; //!< Data stream
    Ref<Image> image;     //!< Storage for image references.
    Ref<Texture> texture; //!< Storage for texture references.
  };
}

#endif
