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

#include "xml_loader.h"
#include "xml_parser.h"
#include "scene.h"
#include "math/affinespace.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include <stack>

namespace embree
{
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
    Variant (int i0) : type(INT1) { i[0] = i0; }

    /*! Constructs a variant object holding an int2 value. */
    Variant (Vec2i v) : type(INT2) { i[0] = v.x; i[1] = v.y; }

    /*! Constructs a variant object holding an int3 value. */
    Variant (Vec3i v) : type(INT3) { i[0] = v.x; i[1] = v.y; i[2] = v.z; }

    /*! Constructs a variant object holding an int4 value. */
    Variant (Vec4i v) : type(INT4) { i[0] = v.x; i[1] = v.y; i[2] = v.z; i[3] = v.w; }

    /*! Constructs a variant object holding a float value. */
    Variant (float f0) : type(FLOAT1) { f[0] = f0; }

    /*! Constructs a variant object holding a float2 value. */
    Variant (Vec2f v) : type(FLOAT2) { f[0] = v.x; f[1] = v.y; }

    /*! Constructs a variant object holding a float3 value. */
    Variant (Vec3f v) : type(FLOAT3) { f[0] = v.x; f[1] = v.y; f[2] = v.z; }

    /*! Constructs a variant object holding a float4 value. */
    Variant (Vec4f v) : type(FLOAT4) { f[0] = v.x; f[1] = v.y; f[2] = v.z; f[3] = v.w; }

    /*! Constructs a variant object holding a string value. */
    Variant (const char* str) : type(STRING), str(str) {}

    /*! Constructs a variant object holding a string value. */
    Variant (const std::string& str) : type(STRING), str(str) {}

    /*! Extracts a boolean from the variant type. */
    bool  getBool () const { return b[0]; }

    /*! Extracts an integer from the variant type. */
    int   getInt  () const { return i[0]; }

    /*! Extracts a float from the variant type. */
    float getFloat() const { return f[0]; }

    /*! Extracts a Vec2f from the variant type. */
    Vec2f getVec2f() const { return Vec2f(f[0],f[1]); }

    /*! Extracts a Vec3f from the variant type. */
    Vec3f getVec3f() const { return Vec3f(f[0],f[1],f[2]); }

    /*! Extracts a Vec3fa from the variant type. */
    Vec3f getVec3fa() const { return Vec3fa(f[0],f[1],f[2]); }

    /*! Extracts a string from the variant type. */
    std::string getString() const { return str;   }

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
  };

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

    /*! Extracts a named Vec3f out of the container. */
    Vec3f getVec3f(const char* name, const Vec3f& def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT3) return def;
      return (*i).second.getVec3f();
    }

    /*! Extracts a named Vec3f out of the container. */
    Vec3fa getVec3fa(const char* name, const Vec3fa& def = zero) const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::FLOAT3) return def;
      return (*i).second.getVec3fa();
    }

    /*! Extracts a named string out of the container. */
    std::string getString(const char* name, std::string def = "") const {
      std::map<std::string,Variant>::const_iterator i = m.find(name);
      if (i == m.end() || (*i).second.type != Variant::STRING) return def;
      return (*i).second.getString();
    }

    /*! Adds a new named element to the container. */
    void add(const std::string& name, Variant data) {
      m[name] = data;
    }

  private:

    /*! Implementation of the container as an STL map. */
    std::map<std::string,Variant> m;
  };

  class XMLLoader
  {
  public:

    XMLLoader(const FileName& fileName, const AffineSpace3f& space, OBJScene& scene);
   ~XMLLoader();

  public:
    void loadPointLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadSpotLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadDirectionalLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadDistantLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadAmbientLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadTriangleLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadHDRILight(const Ref<XML>& xml, const AffineSpace3f& space);
    Parms loadMaterialParms(const Ref<XML>& parms);
    int loadMaterial(const Ref<XML>& xml, std::string* name = NULL);
    void loadTriangleMesh(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadSubdivMesh(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadSphere(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadDisk(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadQuadLight(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadScene(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadTransformNode(const Ref<XML>& xml, const AffineSpace3f& space);
    void loadGroupNode(const Ref<XML>& xml, const AffineSpace3f& space);

  private:
    template<typename T> T load(const Ref<XML>& xml) { return T(zero); }
    template<typename T> T load(const Ref<XML>& xml, const T& opt) { return T(zero); }
    char* loadBinary(const Ref<XML>& xml, size_t eltSize, size_t& size);

    std::vector<float> loadFloatArray(const Ref<XML>& xml);
    std::vector<Vec2f> loadVec2fArray(const Ref<XML>& xml);
    std::vector<Vec3f> loadVec3fArray(const Ref<XML>& xml);
    std::vector<int>   loadIntArray(const Ref<XML>& xml);
    std::vector<Vec2i> loadVec2iArray(const Ref<XML>& xml);
    std::vector<Vec3i> loadVec3iArray(const Ref<XML>& xml);

  private:
    FileName path;         //!< path to XML file
    FILE* binFile;         //!< .bin file for reading binary data
    FileName binFileName;  //!< name of the .bin file

  private:
    std::map<std::string,int> materialMap;              //!< named materials
    std::map<Ref<XML>, int> materialCache;              //!< map for detecting repeated materials

  public:
    OBJScene& scene;
  };

  //////////////////////////////////////////////////////////////////////////////
  //// Loading standard types from an XML node
  //////////////////////////////////////////////////////////////////////////////

  template<> std::string XMLLoader::load<std::string>(const Ref<XML>& xml) {
    if (xml->body.size() < 1) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong string body");
    return xml->body[0].String();
  }

  template<> bool XMLLoader::load<bool>(const Ref<XML>& xml, const bool& opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 1) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong bool body");
    return xml->body[0].Int() != 0;
  }

  template<> int XMLLoader::load<int>(const Ref<XML>& xml) {
    if (xml->body.size() != 1) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong int body");
    return xml->body[0].Int();
  }

  template<> Vec2i XMLLoader::load<Vec2i>(const Ref<XML>& xml) {
    if (xml->body.size() != 2) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong int2 body");
    return Vec2i(xml->body[0].Int(),xml->body[1].Int());
  }

  template<> Vec3i XMLLoader::load<Vec3i>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong int3 body");
    return Vec3i(xml->body[0].Int(),xml->body[1].Int(),xml->body[2].Int());
  }

  template<> Vec4i XMLLoader::load<Vec4i>(const Ref<XML>& xml) {
    if (xml->body.size() != 4) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong int4 body");
    return Vec4i(xml->body[0].Int(),xml->body[1].Int(),xml->body[2].Int(),xml->body[3].Int());
  }

  template<> float XMLLoader::load<float>(const Ref<XML>& xml) {
    if (xml->body.size() != 1) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float body");
    return xml->body[0].Float();
  }

  template<> float XMLLoader::load<float>(const Ref<XML>& xml, const float& opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 1) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float body");
    return xml->body[0].Float();
  }

  template<> Vec2f XMLLoader::load<Vec2f>(const Ref<XML>& xml) {
    if (xml->body.size() != 2) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float2 body");
    return Vec2f(xml->body[0].Float(),xml->body[1].Float());
  }

  template<> Vec3f XMLLoader::load<Vec3f>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float3 body");
    return Vec3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }

  template<> Vec3fa XMLLoader::load<Vec3fa>(const Ref<XML>& xml, const Vec3fa& opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 3) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float3 body");
    return Vec3fa(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }

  template<> Vec4f XMLLoader::load<Vec4f>(const Ref<XML>& xml) {
    if (xml->body.size() != 4) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong float4 body");
    return Vec4f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float(),xml->body[3].Float());
  }

  template<> AffineSpace3f XMLLoader::load<AffineSpace3f>(const Ref<XML>& xml) 
  {
    if (xml->parm("translate") != "") {
      float x,y,z; sscanf(xml->parm("translate").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::translate(Vec3f(x,y,z));
    } else if (xml->parm("scale") != "") {
      float x,y,z; sscanf(xml->parm("scale").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::scale(Vec3f(x,y,z));
    } else if (xml->parm("rotate_x") != "") {
      float degrees; sscanf(xml->parm("rotate_x").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vec3f(1,0,0),deg2rad(degrees));
    } else if (xml->parm("rotate_y") != "") {
      float degrees; sscanf(xml->parm("rotate_y").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vec3f(0,1,0),deg2rad(degrees));
    } else if (xml->parm("rotate_z") != "") {
      float degrees; sscanf(xml->parm("rotate_z").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vec3f(0,0,1),deg2rad(degrees));
    } else if (xml->parm("rotate") != "" && xml->parm("axis") != "") {
      float degrees; sscanf(xml->parm("rotate").c_str(),"%f",&degrees);
      float x,y,z; sscanf(xml->parm("axis").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::rotate(Vec3f(x,y,z),deg2rad(degrees));
    } else {
      if (xml->body.size() != 12) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong AffineSpace body");
      return AffineSpace3f(LinearSpace3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[ 2].Float(),
					 xml->body[4].Float(),xml->body[5].Float(),xml->body[ 6].Float(),
					 xml->body[8].Float(),xml->body[9].Float(),xml->body[10].Float()),
			   Vec3f(xml->body[3].Float(),xml->body[7].Float(),xml->body[11].Float()));
    }
  }

  char* XMLLoader::loadBinary(const Ref<XML>& xml, size_t eltSize, size_t& size)
  {
    if (!binFile) 
      THROW_RUNTIME_ERROR("cannot open file "+binFileName.str()+" for reading");

    size_t ofs = atol(xml->parm("ofs").c_str());
    fseek(binFile,long(ofs),SEEK_SET);

    size = atol(xml->parm("size").c_str());
    char* data = (char*) alignedMalloc(size*eltSize);

    if (size != fread(data, eltSize, size, binFile)) 
      THROW_RUNTIME_ERROR("error reading from binary file: "+binFileName.str());

    return data;
  }

  std::vector<float> XMLLoader::loadFloatArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<float>(); }

    size_t size = 0;
    float* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (float*)loadBinary(xml,sizeof(float),size);
    } else {
      size_t elts = xml->body.size();
      size = elts;
      data = (float*) alignedMalloc(size*sizeof(float));
      for (size_t i=0; i<size; i++) 
        data[i] = xml->body[i].Float();
    }
    std::vector<float> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  std::vector<Vec2f> XMLLoader::loadVec2fArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<Vec2f>(); }

    size_t size = 0;
    Vec2f* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec2f*)loadBinary(xml,2*sizeof(float),size);
    } else {
      size_t elts = xml->body.size();
      if (elts % 2 != 0) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong vector<float2> body");
      size = elts/2;
      data = (Vec2f*) alignedMalloc(size*sizeof(Vec2f));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec2f(xml->body[2*i+0].Float(),xml->body[2*i+1].Float());
    }
    std::vector<Vec2f> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  std::vector<Vec3f> XMLLoader::loadVec3fArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<Vec3f>(); }

    size_t size = 0;
    Vec3f* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec3f*) loadBinary(xml,3*sizeof(float),size);
    }
    else {
      size_t elts = xml->body.size();
      if (elts % 3 != 0) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong vector<float3> body");
      size = elts/3;
      data = (Vec3f*) alignedMalloc(size*sizeof(Vec3f));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec3f(xml->body[3*i+0].Float(),xml->body[3*i+1].Float(),xml->body[3*i+2].Float());
    }
    std::vector<Vec3f> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  std::vector<int> XMLLoader::loadIntArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<int>(); }

    size_t size = 0;
    int* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (int*)loadBinary(xml,sizeof(int),size);
    } else {
      size_t elts = xml->body.size();
      size = elts;
      data = (int*) alignedMalloc(size*sizeof(int));
      for (size_t i=0; i<size; i++) 
        data[i] = xml->body[i].Int();
    }
    std::vector<int> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  std::vector<Vec2i> XMLLoader::loadVec2iArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<Vec2i>(); }
    
    size_t size = 0;
    Vec2i* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec2i*) loadBinary(xml,2*sizeof(int),size);
    }
    else {
      size_t elts = xml->body.size();
      if (elts % 2 != 0) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong vector<int2> body");
      size = elts/2;
      data = (Vec2i*) alignedMalloc(size*sizeof(Vec2i));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec2i(xml->body[2*i+0].Int(),xml->body[2*i+1].Int());
    }
    std::vector<Vec2i> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  std::vector<Vec3i> XMLLoader::loadVec3iArray(const Ref<XML>& xml)
  {
    /*! do not fail of array does not exist */
    if (!xml) { return std::vector<Vec3i>(); }
    
    size_t size = 0;
    Vec3i* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec3i*) loadBinary(xml,3*sizeof(int),size);
    }
    else {
      size_t elts = xml->body.size();
      if (elts % 3 != 0) THROW_RUNTIME_ERROR(xml->loc.str()+": wrong vector<int3> body");
      size = elts/3;
      data = (Vec3i*) alignedMalloc(size*sizeof(Vec3i));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec3i(xml->body[3*i+0].Int(),xml->body[3*i+1].Int(),xml->body[3*i+2].Int());
    }
    std::vector<Vec3i> res;
    for (size_t i=0; i<size; i++) res.push_back(data[i]);
    alignedFree(data);
    return res;
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Loading of objects from XML file
  //////////////////////////////////////////////////////////////////////////////

  void XMLLoader::loadPointLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa I = load<Vec3f>(xml->child("I"));
    Vec3fa P = space.p;
    scene.pointLights.push_back(OBJScene::PointLight(P,I));
  }

  void XMLLoader::loadSpotLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa I = load<Vec3fa>(xml->child("I"));
    Vec3fa P = space.p;
    Vec3fa D = space.l.vz;
    float angleMin = load<float>(xml->child("angleMin"));
    float angleMax = load<float>(xml->child("angleMax"));
  }

  void XMLLoader::loadDirectionalLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa E = load<Vec3fa>(xml->child("E"));
    Vec3fa D = space.l.vz;
    scene.directionalLights.push_back(OBJScene::DirectionalLight(D,E));
  }

  void XMLLoader::loadDistantLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa L = load<Vec3fa>(xml->child("L"));
    Vec3fa D = space.l.vz;
    float halfAngle = load<float>(xml->child("halfAngle"));
    scene.distantLights.push_back(OBJScene::DistantLight(D,L,halfAngle));
  }

  void XMLLoader::loadAmbientLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    Vec3fa L = load<Vec3fa>(xml->child("L"));
    scene.ambientLights.push_back(OBJScene::AmbientLight(L));
  }

  void XMLLoader::loadTriangleLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa L = load<Vec3fa>(xml->child("L"));
    Vec3fa v0 = xfmPoint(space, Vec3fa(1, 0, 0));
    Vec3fa v1 = xfmPoint(space, Vec3fa(0, 1, 0));
    Vec3fa v2 = xfmPoint(space, Vec3fa(0, 0, 0));
  }

  void XMLLoader::loadQuadLight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa L = load<Vec3fa>(xml->child("L"));
    Vec3fa v0 = xfmPoint(space, Vec3fa(0, 0, 0));
    Vec3fa v1 = xfmPoint(space, Vec3fa(0, 1, 0));
    Vec3fa v2 = xfmPoint(space, Vec3fa(1, 1, 0));
    Vec3fa v3 = xfmPoint(space, Vec3fa(1, 0, 0));
  }

  void XMLLoader::loadHDRILight(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->child("AffineSpace"));
    Vec3fa L = load<Vec3fa>(xml->child("L"));
    //image =  rtLoadImage(path + load<std::string>(xml->child("image"))));
  }

  Parms XMLLoader::loadMaterialParms(const Ref<XML>& parms)
  {
    Parms material;
    for (size_t i=0; i<parms->children.size(); i++) 
    {
      Ref<XML> entry = parms->children[i];
      std::string name = entry->parm("name");
      if      (entry->name == "int"    ) { material.add(name,load<int>  (entry)); }
      else if (entry->name == "int2"   ) { material.add(name, load<Vec2i>(entry)); }
      else if (entry->name == "int3"   ) { material.add(name, load<Vec3i>(entry)); }
      else if (entry->name == "int4"   ) { material.add(name, load<Vec4i>(entry)); }
      else if (entry->name == "float"  ) { material.add(name, load<float>(entry)); }
      else if (entry->name == "float2" ) { material.add(name, load<Vec2f>(entry)); }
      else if (entry->name == "float3" ) { material.add(name, load<Vec3f>(entry)); }
      else if (entry->name == "float4" ) { material.add(name, load<Vec4f>(entry)); }
      else if (entry->name == "texture") { material.add(name, (path + load<std::string>(entry)).str()); }
      else THROW_RUNTIME_ERROR(entry->loc.str()+": invalid type: "+entry->name);
    }
    return material;
  }

  int XMLLoader::loadMaterial(const Ref<XML>& xml, std::string* name) 
  {
    if (xml->parm("id") != "") {
      if (name) *name = xml->parm("id");
      return materialMap[xml->parm("id")];
    }

    Ref<XML> parameters = xml->child("parameters");
    if (materialCache.find(parameters) != materialCache.end()) {
      return materialCache[parameters];
    }

    std::string type = load<std::string>(xml->child("code")).c_str();
    Parms parms = loadMaterialParms(parameters);

    OBJScene::Material material;
    if (type == "Matte")
    {
      const Vec3fa reflectance = parms.getVec3fa("reflectance",one);
      new (&material) OBJScene::MatteMaterial(reflectance);
    }
    else if (type == "Mirror")
    {
      const Vec3fa reflectance = parms.getVec3fa("reflectance",one);
      new (&material) OBJScene::MirrorMaterial(reflectance);
    }
    else if (type == "OBJ") 
    {
      //map_d = parms.getTexture("map_d");  
      const float d = parms.getFloat("d", 1.0f);
      //map_Kd = parms.getTexture("map_Kd");  
      const Vec3fa Kd = parms.getVec3fa("Kd", one);
      //map_Ks = parms.getTexture("map_Ks");  
      const Vec3fa Ks = parms.getVec3fa("Ks", zero);
      //map_Ns = parms.getTexture("map_Ns");  
      const float Ns = parms.getFloat("Ns", 10.0f);
      //map_Bump = parms.getTexture("map_Bump");
      new (&material) OBJScene::OBJMaterial(d,Kd,Ks,Ns);
    }
    else if (type == "ThinDielectric" || type == "ThinGlass")
    {
      const Vec3fa transmission = parms.getVec3fa("transmission",one);
      const float eta          = parms.getFloat("eta",1.4f);
      const float thickness    = parms.getFloat("thickness",0.1f);
      new (&material) OBJScene::ThinDielectricMaterial(transmission,eta,thickness);
    }
    else if (type == "Plastic")
    {
      const Vec3fa pigmentColor = parms.getVec3fa("pigmentColor",one);
      const float eta          = parms.getFloat("eta",1.4f);
      const float roughness    = parms.getFloat("roughness",0.01f);
      new (&material) OBJScene::MetallicPaintMaterial(pigmentColor,pigmentColor,roughness,eta);
    }
    else if (type == "Metal")
    {
      const Vec3fa reflectance  = parms.getVec3fa("reflectance",one);
      const Vec3fa eta          = parms.getVec3fa("eta",Vec3fa(1.4f));
      const Vec3fa k            = parms.getVec3fa("k",Vec3fa(0.0f));
      const float roughness     = parms.getFloat("roughness",0.01f);
      if (roughness == 0.0f)
        new (&material) OBJScene::MetalMaterial(reflectance,eta,k);
      else 
        new (&material) OBJScene::MetalMaterial(reflectance,eta,k,roughness);
    }
    else if (type == "Velvet")
    {
      const Vec3fa reflectance = parms.getVec3fa("reflectance",one);
      const float backScattering = parms.getFloat("backScattering",zero);
      const Vec3fa horizonScatteringColor = parms.getVec3fa("horizonScatteringColor",one);
      const float horizonScatteringFallOff = parms.getFloat("horizonScatteringFallOff",zero);
      new (&material) OBJScene::VelvetMaterial(reflectance,backScattering,horizonScatteringColor,horizonScatteringFallOff);
    }
    else if (type == "Dielectric")
    {
      const Vec3fa transmissionOutside = parms.getVec3fa("transmissionOutside",one);
      const Vec3fa transmissionInside  = parms.getVec3fa("transmission",one);
      const float etaOutside = parms.getFloat("etaOutside",1.0f);
      const float etaInside  = parms.getFloat("etaInside",1.4f);
      new (&material) OBJScene::DielectricMaterial(transmissionOutside,transmissionInside,etaOutside,etaInside);
    }
    else if (type == "MetallicPaint")
    {
      const Vec3fa shadeColor    = parms.getVec3fa("shadeColor",one);
      const Vec3fa glitterColor  = parms.getVec3fa("glitterColor",zero);
      const float glitterSpread = parms.getFloat("glitterSpread",1.0f);
      const float eta           = parms.getFloat("eta",1.4f);
      new (&material) OBJScene::MetallicPaintMaterial(shadeColor,glitterColor,glitterSpread,eta);
    }
    else {
      std::cout << "Warning: unsupported material " << type << std::endl;
      new (&material) OBJScene::OBJMaterial(1.0f,0.5f,0.0f,0.0f);
    }
    int materialID = scene.materials.size();
    scene.materials.push_back(material);
    materialCache[parameters] = materialID;
    return materialID;
  }

  void XMLLoader::loadSubdivMesh(const Ref<XML>& xml, const AffineSpace3f& space) 
  {
    std::string materialName;
    int materialID = loadMaterial(xml->child("material"),&materialName);

    OBJScene::SubdivMesh* mesh = new OBJScene::SubdivMesh;
    std::vector<Vec3f> positions = loadVec3fArray(xml->childOpt("positions"));
    for (size_t i=0; i<positions.size(); i++) mesh->positions.push_back(xfmPoint(space,positions[i]));
    std::vector<Vec3f> normals = loadVec3fArray(xml->childOpt("normals"));
    for (size_t i=0; i<normals.size(); i++) mesh->normals.push_back(xfmNormal(space,normals[i]));
    mesh->texcoords = loadVec2fArray(xml->childOpt("texcoords"));
    mesh->position_indices = loadIntArray(xml->childOpt("position_indices"));
    mesh->normal_indices   = loadIntArray(xml->childOpt("normal_indices"));
    mesh->texcoord_indices = loadIntArray(xml->childOpt("texcoord_indices"));
    mesh->verticesPerFace  = loadIntArray(xml->childOpt("faces"));
    mesh->holes            = loadIntArray(xml->childOpt("holes"));
    mesh->edge_creases     = loadVec2iArray(xml->childOpt("edge_creases"));
    mesh->edge_crease_weights = loadFloatArray(xml->childOpt("edge_crease_weights"));
    mesh->vertex_creases      = loadIntArray(xml->childOpt("vertex_creases"));
    mesh->vertex_crease_weights = loadFloatArray(xml->childOpt("vertex_crease_weights"));
    mesh->materialID = materialID;
    scene.subdiv.push_back(mesh);
  }

  void XMLLoader::loadTriangleMesh(const Ref<XML>& xml, const AffineSpace3f& space) 
  {
    std::string materialName;
    int materialID = loadMaterial(xml->child("material"),&materialName);
    std::vector<Vec3f> positions = loadVec3fArray(xml->childOpt("positions"));
    std::vector<Vec3f> motions   = loadVec3fArray(xml->childOpt("motions"  ));
    std::vector<Vec3f> normals   = loadVec3fArray(xml->childOpt("normals"  ));
    std::vector<Vec2f> texcoords = loadVec2fArray(xml->childOpt("texcoords"));
    std::vector<Vec3i> triangles = loadVec3iArray(xml->childOpt("triangles"));

    OBJScene::Mesh* mesh = new OBJScene::Mesh;
    for (size_t i=0; i<positions.size(); i++)
      mesh->v.push_back(xfmPoint(space,positions[i]));
    for (size_t i=0; i<normals.size(); i++)
      mesh->vn.push_back(xfmVector(space,normals[i]));     
    for (size_t i=0; i<texcoords.size(); i++)
      mesh->vt.push_back(texcoords[i]);
    for (size_t i=0; i<triangles.size(); i++)
      mesh->triangles.push_back(OBJScene::Triangle(triangles[i].x,triangles[i].y,triangles[i].z,materialID));

    scene.meshes.push_back(mesh);
  }

  void XMLLoader::loadSphere(const Ref<XML>& xml, const AffineSpace3f& space) {
    std::cout << "Warning: ignoring sphere" << std::endl;
  }

  void XMLLoader::loadDisk(const Ref<XML>& xml, const AffineSpace3f& space) {
    std::cout << "Warning: ignoring disk" << std::endl;
  }

  void XMLLoader::loadTransformNode(const Ref<XML>& xml, const AffineSpace3f& space_in) 
  {
    AffineSpace3f space = space_in*load<AffineSpace3f>(xml->children[0]);
    for (size_t i=1; i<xml->children.size(); i++)
      loadScene(xml->children[i],space);
  }

  void XMLLoader::loadGroupNode(const Ref<XML>& xml, const AffineSpace3f& space) 
  {
    for (size_t i=0; i<xml->children.size(); i++)
      loadScene(xml->children[i],space);
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Loading of scene graph node from XML file
  //////////////////////////////////////////////////////////////////////////////
  
  void XMLLoader::loadScene(const Ref<XML>& xml, const AffineSpace3f& space)
  {
    if (xml->name == "assign") 
    {
      if (xml->parm("type") == "material")
        materialMap[xml->parm("id")] = loadMaterial(xml->child(0));
      //else if (xml->parm("type") == "scene")
      //sceneMap[xml->parm("id")] = loadScene(xml->child(0));
      else 
        THROW_RUNTIME_ERROR(xml->loc.str()+": unknown type: "+xml->parm("type"));
    }
    else 
    {
      if (xml->name == "xml") {
        loadXML(path + xml->parm("src"),space,scene);
      }
      else if (xml->name == "obj") {
        loadOBJ(path + xml->parm("src"),space,scene);
      }
      else if (xml->name == "extern") {
        FileName fname = path + xml->parm("src");
        if      (fname.ext() == "xml") loadXML(path + xml->parm("src"),space,scene);
        else if (fname.ext() == "obj") loadOBJ(path + xml->parm("src"),space,scene);
        else THROW_RUNTIME_ERROR("unknown file type:" + fname.str());
      }
      //else if (xml->name == "ref") {
      //  prims = sceneMap[xml->parm("id")];
      //  for (size_t i=0; i<prims.size(); i++)
      //    prims[i] = g_device->rtTransformPrimitive(prims[i],copyToArray(transforms.top()));
      //}
      
      else if (xml->name == "PointLight"      ) loadPointLight      (xml,space);
      else if (xml->name == "SpotLight"       ) loadSpotLight       (xml,space);
      else if (xml->name == "DirectionalLight") loadDirectionalLight(xml,space);
      else if (xml->name == "DistantLight"    ) loadDistantLight    (xml,space);
      else if (xml->name == "AmbientLight"    ) loadAmbientLight    (xml,space);
      else if (xml->name == "TriangleLight"   ) loadTriangleLight   (xml,space);
      else if (xml->name == "QuadLight"       ) loadQuadLight       (xml,space);
      else if (xml->name == "HDRILight"       ) loadHDRILight       (xml,space);

      else if (xml->name == "TriangleMesh"    ) loadTriangleMesh    (xml,space);
      else if (xml->name == "SubdivisionMesh" ) loadSubdivMesh      (xml,space);
      else if (xml->name == "Sphere"          ) loadSphere          (xml,space);
      else if (xml->name == "Disk"            ) loadDisk            (xml,space);
      else if (xml->name == "Group"           ) loadGroupNode       (xml,space);
      else if (xml->name == "Transform"       ) loadTransformNode   (xml,space);
      
      else THROW_RUNTIME_ERROR(xml->loc.str()+": unknown tag: "+xml->name);
    }
  }

  XMLLoader::XMLLoader(const FileName& fileName, const AffineSpace3f& space, OBJScene& scene) : binFile(NULL), scene(scene)
  {
    path = fileName.path();
    binFileName = fileName.setExt(".bin");
    binFile = fopen(binFileName.c_str(),"rb");

    Ref<XML> xml = parseXML(fileName);
    if (xml->name != "scene") THROW_RUNTIME_ERROR(xml->loc.str()+": invalid scene tag");
    for (size_t i=0; i<xml->children.size(); i++) {
      loadScene(xml->children[i],space);
    }
  }

  XMLLoader::~XMLLoader() {
    if (binFile) fclose(binFile);
  }

  /*! read from disk */
  void loadXML(const FileName& fileName, const AffineSpace3f& space, OBJScene& scene) {
    XMLLoader loader(fileName,space,scene);
  }
}
