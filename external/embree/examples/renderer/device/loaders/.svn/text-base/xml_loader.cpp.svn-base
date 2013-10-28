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

#include "loaders.h"
#include "xml_parser.h"
#include "obj_loader.h"
#include "image/image.h"
#include "math/affinespace.h"
#include "math/color.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include <stack>

namespace embree
{
  class XMLLoader
  {
  public:

    XMLLoader(const FileName& fileName);
   ~XMLLoader();

  public:
    Handle<Device::RTPrimitive> loadPointLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadSpotLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadDirectionalLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadDistantLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadAmbientLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadTriangleLight(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadHDRILight(const Ref<XML>& xml);
    void                loadMaterialParms(Handle<Device::RTMaterial> material, const Ref<XML>& parms);
    Handle<Device::RTMaterial>  loadMaterial(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadTriangleMesh(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadSphere(const Ref<XML>& xml);
    Handle<Device::RTPrimitive> loadDisk(const Ref<XML>& xml);
    std::vector<Handle<Device::RTPrimitive> > loadQuadLight(const Ref<XML>& xml);
    std::vector<Handle<Device::RTPrimitive> > loadScene(const Ref<XML>& xml);
    std::vector<Handle<Device::RTPrimitive> > loadTransformNode(const Ref<XML>& xml);
    std::vector<Handle<Device::RTPrimitive> > loadGroupNode(const Ref<XML>& xml);

  private:
    template<typename T> T load(const Ref<XML>& xml) { return T(zero); }
    template<typename T> T load(const Ref<XML>& xml, T opt) { return T(zero); }
    char* loadBinary(const Ref<XML>& xml, size_t eltSize, size_t& size);
    Handle<Device::RTData> loadVec2fArray(const Ref<XML>& xml, size_t& size);
    Handle<Device::RTData> loadVec3fArray(const Ref<XML>& xml, size_t& size);
    Handle<Device::RTData> loadVector3iArray(const Ref<XML>& xml, size_t& size);

  private:
    FileName path;         //!< path to XML file
    FILE* binFile;         //!< .bin file for reading binary data
    FileName binFileName;  //!< name of the .bin file

  private:
    std::stack<AffineSpace3f> transforms;                              //!< stack of transformations
    std::map<std::string,Handle<Device::RTMaterial> > materialMap;              //!< named materials
    std::map<Ref<XML>, Handle<Device::RTMaterial> > materialCache;              //!< map for detecting repeated materials
    std::map<std::string,std::vector<Handle<Device::RTPrimitive> > > sceneMap;  //!< named parts of the scene

  public:
    std::vector<Handle<Device::RTPrimitive> > model;   //!< stores the output scene
  };

  //////////////////////////////////////////////////////////////////////////////
  //// Loading standard types from an XML node
  //////////////////////////////////////////////////////////////////////////////

  template<> std::string XMLLoader::load<std::string>(const Ref<XML>& xml) {
    if (xml->body.size() < 1) throw std::runtime_error(xml->loc.str()+": wrong string body");
    return xml->body[0].String();
  }

  template<> bool XMLLoader::load<bool>(const Ref<XML>& xml, bool opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 1) throw std::runtime_error(xml->loc.str()+": wrong bool body");
    return xml->body[0].Int() != 0;
  }

  template<> int XMLLoader::load<int>(const Ref<XML>& xml) {
    if (xml->body.size() != 1) throw std::runtime_error(xml->loc.str()+": wrong int body");
    return xml->body[0].Int();
  }

  template<> Vec2i XMLLoader::load<Vec2i>(const Ref<XML>& xml) {
    if (xml->body.size() != 2) throw std::runtime_error(xml->loc.str()+": wrong int2 body");
    return Vec2i(xml->body[0].Int(),xml->body[1].Int());
  }

  template<> Vector3i XMLLoader::load<Vector3i>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) throw std::runtime_error(xml->loc.str()+": wrong int3 body");
    return Vector3i(xml->body[0].Int(),xml->body[1].Int(),xml->body[2].Int());
  }

  template<> Vec4i XMLLoader::load<Vec4i>(const Ref<XML>& xml) {
    if (xml->body.size() != 4) throw std::runtime_error(xml->loc.str()+": wrong int4 body");
    return Vec4i(xml->body[0].Int(),xml->body[1].Int(),xml->body[2].Int(),xml->body[3].Int());
  }

  template<> float XMLLoader::load<float>(const Ref<XML>& xml) {
    if (xml->body.size() != 1) throw std::runtime_error(xml->loc.str()+": wrong float body");
    return xml->body[0].Float();
  }

  template<> float XMLLoader::load<float>(const Ref<XML>& xml, float opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 1) throw std::runtime_error(xml->loc.str()+": wrong float body");
    return xml->body[0].Float();
  }

  template<> Vec2f XMLLoader::load<Vec2f>(const Ref<XML>& xml) {
    if (xml->body.size() != 2) throw std::runtime_error(xml->loc.str()+": wrong float2 body");
    return Vec2f(xml->body[0].Float(),xml->body[1].Float());
  }

  template<> Vec3f XMLLoader::load<Vec3f>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) throw std::runtime_error(xml->loc.str()+": wrong float3 body");
    return Vec3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }

  template<> Vector3f XMLLoader::load<Vector3f>(const Ref<XML>& xml, Vector3f opt) {
    if (xml == null) return opt;
    if (xml->body.size() != 3) throw std::runtime_error(xml->loc.str()+": wrong float3 body");
    return Vector3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }

#if defined(__X86_64__)
  template<> Vector3f XMLLoader::load<Vector3f>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) throw std::runtime_error(xml->loc.str()+": wrong float3 body");
    return Vector3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }
#endif

  template<> Vec4f XMLLoader::load<Vec4f>(const Ref<XML>& xml) {
    if (xml->body.size() != 4) throw std::runtime_error(xml->loc.str()+": wrong float4 body");
    return Vec4f(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float(),xml->body[3].Float());
  }

  template<> Color XMLLoader::load<Color>(const Ref<XML>& xml) {
    if (xml->body.size() != 3) throw std::runtime_error(xml->loc.str()+": wrong color body");
    return Color(xml->body[0].Float(),xml->body[1].Float(),xml->body[2].Float());
  }

  template<> AffineSpace3f XMLLoader::load<AffineSpace3f>(const Ref<XML>& xml) 
  {
    if (xml->parm("translate") != "") {
      float x,y,z; sscanf(xml->parm("translate").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::translate(Vector3f(x,y,z));
    } else if (xml->parm("scale") != "") {
      float x,y,z; sscanf(xml->parm("scale").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::scale(Vector3f(x,y,z));
    } else if (xml->parm("rotate_x") != "") {
      float degrees; sscanf(xml->parm("rotate_x").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vector3f(1,0,0),deg2rad(degrees));
    } else if (xml->parm("rotate_y") != "") {
      float degrees; sscanf(xml->parm("rotate_y").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vector3f(0,1,0),deg2rad(degrees));
    } else if (xml->parm("rotate_z") != "") {
      float degrees; sscanf(xml->parm("rotate_z").c_str(),"%f",&degrees);
      return AffineSpace3f::rotate(Vector3f(0,0,1),deg2rad(degrees));
    } else if (xml->parm("rotate") != "" && xml->parm("axis") != "") {
      float degrees; sscanf(xml->parm("rotate").c_str(),"%f",&degrees);
      float x,y,z; sscanf(xml->parm("axis").c_str(),"%f %f %f",&x,&y,&z);
      return AffineSpace3f::rotate(Vector3f(x,y,z),deg2rad(degrees));
    } else {
      if (xml->body.size() != 12) throw std::runtime_error(xml->loc.str()+": wrong AffineSpace body");
      return AffineSpace3f(LinearSpace3f(xml->body[0].Float(),xml->body[1].Float(),xml->body[ 2].Float(),
					 xml->body[4].Float(),xml->body[5].Float(),xml->body[ 6].Float(),
					 xml->body[8].Float(),xml->body[9].Float(),xml->body[10].Float()),
			   Vector3f(xml->body[3].Float(),xml->body[7].Float(),xml->body[11].Float()));
    }
  }

  char* XMLLoader::loadBinary(const Ref<XML>& xml, size_t eltSize, size_t& size)
  {
    if (!binFile) 
      throw std::runtime_error("cannot open file "+binFileName.str()+" for reading");

    size_t ofs = atol(xml->parm("ofs").c_str());
    fseek(binFile,long(ofs),SEEK_SET);

    size = atol(xml->parm("size").c_str());
    char* data = (char*) alignedMalloc(size*eltSize);

    if (size != fread(data, eltSize, size, binFile)) 
      throw std::runtime_error("error reading from binary file: "+binFileName.str());

    return data;
  }

  Handle<Device::RTData> XMLLoader::loadVec2fArray(const Ref<XML>& xml, size_t& size)
  {
    /*! do not fail of array does not exist */
    if (!xml) { size = 0; return NULL; }

    Vec2f* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec2f*)loadBinary(xml,2*sizeof(float),size);
    } else {
      size_t elts = xml->body.size();
      if (elts % 2 != 0) throw std::runtime_error(xml->loc.str()+": wrong vector<float2> body");
      size = elts/2;
      data = (Vec2f*) alignedMalloc(size*sizeof(Vec2f));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec2f(xml->body[2*i+0].Float(),xml->body[2*i+1].Float());
    }
    return g_device->rtNewData("immutable_managed",size*sizeof(Vec2f),data);
  }

  Handle<Device::RTData> XMLLoader::loadVec3fArray(const Ref<XML>& xml, size_t& size)
  {
    /*! do not fail of array does not exist */
    if (!xml) { size = 0; return NULL; }

    Vec3f* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec3f*) loadBinary(xml,3*sizeof(float),size);
    }
    else {
      size_t elts = xml->body.size();
      if (elts % 3 != 0) throw std::runtime_error(xml->loc.str()+": wrong vector<float3> body");
      size = elts/3;
      data = (Vec3f*) alignedMalloc(size*sizeof(Vec3f));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec3f(xml->body[3*i+0].Float(),xml->body[3*i+1].Float(),xml->body[3*i+2].Float());
    }
    return g_device->rtNewData("immutable_managed",size*sizeof(Vec3f),data);
  }

  Handle<Device::RTData> XMLLoader::loadVector3iArray(const Ref<XML>& xml, size_t& size)
  {
    /*! do not fail of array does not exist */
    if (!xml) { size = 0; return NULL; }

    Vec3i* data = NULL;
    if (xml->parm("ofs") != "") {
      data = (Vec3i*) loadBinary(xml,3*sizeof(int),size);
    }
    else {
      size_t elts = xml->body.size();
      if (elts % 3 != 0) throw std::runtime_error(xml->loc.str()+": wrong vector<int3> body");
      size = elts/3;
      data = (Vec3i*) alignedMalloc(size*sizeof(Vec3i));
      for (size_t i=0; i<size; i++) 
        data[i] = Vec3i(xml->body[3*i+0].Int(),xml->body[3*i+1].Int(),xml->body[3*i+2].Int());
    }
    return g_device->rtNewData("immutable_managed",size*sizeof(Vec3i),data);
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Loading of objects from XML file
  //////////////////////////////////////////////////////////////////////////////

  Handle<Device::RTPrimitive> XMLLoader::loadPointLight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color I = load<Color>(xml->child("I"));

    Handle<Device::RTLight> light = g_device->rtNewLight("pointlight");
    g_device->rtSetFloat3(light, "P", space.p.x, space.p.y, space.p.z);
    g_device->rtSetFloat3(light, "I", I.r, I.g, I.b);
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadSpotLight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color I = load<Color>(xml->child("I"));

    Handle<Device::RTLight> light = g_device->rtNewLight("spotlight");
    g_device->rtSetFloat3(light, "P", space.p.x, space.p.y, space.p.z);
    g_device->rtSetFloat3(light, "D", space.l.vz.x, space.l.vz.y, space.l.vz.z);
    g_device->rtSetFloat3(light, "I", I.r, I.g, I.b);
    g_device->rtSetFloat1(light, "angleMin", load<float>(xml->child("angleMin")));
    g_device->rtSetFloat1(light, "angleMax", load<float>(xml->child("angleMax")));
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadDirectionalLight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color E = load<Color>(xml->child("E"));

    Handle<Device::RTLight> light = g_device->rtNewLight("directionallight");
    g_device->rtSetFloat3(light, "D", space.l.vz.x, space.l.vz.y, space.l.vz.z);
    g_device->rtSetFloat3(light, "E", E.r, E.g, E.b);
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadDistantLight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color L = load<Color>(xml->child("L"));

    Handle<Device::RTLight> light = g_device->rtNewLight("distantlight");
    g_device->rtSetFloat3(light, "D", space.l.vz.x, space.l.vz.y, space.l.vz.z);
    g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
    g_device->rtSetFloat1(light, "halfAngle", load<float>(xml->child("halfAngle")));
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadAmbientLight(const Ref<XML>& xml) 
  {
    Color L = load<Color>(xml->child("L"));

    Handle<Device::RTLight> light = g_device->rtNewLight("ambientlight");
    g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadTriangleLight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color L = load<Color>(xml->child("L"));
    Vector3f v0 = xfmPoint(space, Vector3f(1, 0, 0));
    Vector3f v1 = xfmPoint(space, Vector3f(0, 1, 0));
    Vector3f v2 = xfmPoint(space, Vector3f(0, 0, 0));

    Handle<Device::RTLight> light = g_device->rtNewLight("trianglelight");
    g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
    g_device->rtSetFloat3(light, "v0", v0.x, v0.y, v0.z);
    g_device->rtSetFloat3(light, "v1", v1.x, v1.y, v1.z);
    g_device->rtSetFloat3(light, "v2", v2.x, v2.y, v2.z);
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  std::vector<Handle<Device::RTPrimitive> > XMLLoader::loadQuadLight(const Ref<XML>& xml) 
  {
    std::vector<Handle<Device::RTPrimitive> > lights;

    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color L = load<Color>(xml->child("L"));
    Vector3f v0 = xfmPoint(space, Vector3f(0, 0, 0));
    Vector3f v1 = xfmPoint(space, Vector3f(0, 1, 0));
    Vector3f v2 = xfmPoint(space, Vector3f(1, 1, 0));
    Vector3f v3 = xfmPoint(space, Vector3f(1, 0, 0));

    Handle<Device::RTLight> light0 = g_device->rtNewLight("trianglelight");
    g_device->rtSetFloat3(light0, "L", L.r, L.g, L.b);
    g_device->rtSetFloat3(light0, "v0", v1.x, v1.y, v1.z);
    g_device->rtSetFloat3(light0, "v1", v3.x, v3.y, v3.z);
    g_device->rtSetFloat3(light0, "v2", v0.x, v0.y, v0.z);
    g_device->rtCommit(light0);
    lights.push_back(g_device->rtNewLightPrimitive(light0, NULL, copyToArray(transforms.top())));

    Handle<Device::RTLight> light1 = g_device->rtNewLight("trianglelight");
    g_device->rtSetFloat3(light1, "L", L.r, L.g, L.b);
    g_device->rtSetFloat3(light1, "v0", v2.x, v2.y, v2.z);
    g_device->rtSetFloat3(light1, "v1", v3.x, v3.y, v3.z);
    g_device->rtSetFloat3(light1, "v2", v1.x, v1.y, v1.z);
    g_device->rtCommit(light1);
    lights.push_back(g_device->rtNewLightPrimitive(light1, NULL, copyToArray(transforms.top())));

    return lights;
  }

  Handle<Device::RTPrimitive> XMLLoader::loadHDRILight(const Ref<XML>& xml) 
  {
    AffineSpace3f space = load<AffineSpace3f>(xml->child("AffineSpace"));
    Color L = load<Color>(xml->child("L"));

    Handle<Device::RTLight> light = g_device->rtNewLight("hdrilight");
    g_device->rtSetTransform(light, "local2world", copyToArray(space));
    g_device->rtSetFloat3(light, "L", L.r, L.g, L.b);
    g_device->rtSetImage(light, "image", rtLoadImage(path + load<std::string>(xml->child("image"))));
    g_device->rtCommit(light);
    return g_device->rtNewLightPrimitive(light, NULL, copyToArray(transforms.top()));
  }

  void XMLLoader::loadMaterialParms(Handle<Device::RTMaterial> material, const Ref<XML>& parms)
  {
    for (size_t i=0; i<parms->children.size(); i++) 
    {
      Ref<XML> entry = parms->children[i];
      std::string name = entry->parm("name");
      if      (entry->name == "int"    ) { int   value = load<int>  (entry);  g_device->rtSetInt1  (material, name.c_str(), value); }
      else if (entry->name == "int2"   ) { Vec2i value = load<Vec2i>(entry);  g_device->rtSetInt2  (material, name.c_str(), value.x, value.y); }
      else if (entry->name == "int3"   ) { Vector3i value = load<Vector3i>(entry);  g_device->rtSetInt3  (material, name.c_str(), value.x, value.y, value.z); }
      else if (entry->name == "int4"   ) { Vec4i value = load<Vec4i>(entry);  g_device->rtSetInt4  (material, name.c_str(), value.x, value.y, value.z, value.w); }
      else if (entry->name == "float"  ) { float value = load<float>(entry);  g_device->rtSetFloat1(material, name.c_str(), value); }
      else if (entry->name == "float2" ) { Vec2f value = load<Vec2f>(entry);  g_device->rtSetFloat2(material, name.c_str(), value.x, value.y); }
      else if (entry->name == "float3" ) { Vec3f value = load<Vec3f>(entry);  g_device->rtSetFloat3(material, name.c_str(), value.x, value.y, value.z); }
      else if (entry->name == "float4" ) { Vec4f value = load<Vec4f>(entry);  g_device->rtSetFloat4(material, name.c_str(), value.x, value.y, value.z, value.w); }
      else if (entry->name == "texture") g_device->rtSetTexture(material, name.c_str(), rtLoadTexture(path + load<std::string>(entry)));
      else throw std::runtime_error(entry->loc.str()+": invalid type: "+entry->name);
    }

    g_device->rtCommit(material);
  }

  Handle<Device::RTMaterial> XMLLoader::loadMaterial(const Ref<XML>& xml) 
  {
    if (xml->parm("file") != "") 
    {
      FileName fileName = path+xml->parm("file");
      FileName oldPath = path;
      path = fileName.path();
      Ref<XML> xml1 = parseXML(fileName);
      if (xml1->name != "material") throw std::runtime_error(xml->loc.str()+": invalid material tag");
      Handle<Device::RTMaterial> material = loadMaterial(xml1);
      loadMaterialParms(material,xml);
      path = oldPath;
      return material;
    }
    
    if (xml->parm("id") != "")
      return materialMap[xml->parm("id")];

    Ref<XML> parms = xml->child("parameters");
    if (materialCache.find(parms) != materialCache.end()) {
      return materialCache[parms];
    }

    Handle<Device::RTMaterial> material = g_device->rtNewMaterial(load<std::string>(xml->child("code")).c_str());
    loadMaterialParms(material,parms);
    materialCache[parms] = material;
    return material;
  }

  Handle<Device::RTPrimitive> XMLLoader::loadTriangleMesh(const Ref<XML>& xml) 
  {
    Handle<Device::RTMaterial> material = loadMaterial(xml->child("material"));
    size_t numPositions = 0;  Handle<Device::RTData> positions = loadVec3fArray(xml->childOpt("positions"), numPositions);
    size_t numMotions   = 0;  Handle<Device::RTData> motions   = loadVec3fArray(xml->childOpt("motions"  ), numMotions);
    size_t numNormals   = 0;  Handle<Device::RTData> normals   = loadVec3fArray(xml->childOpt("normals"  ), numNormals);
    size_t numTexCoords = 0;  Handle<Device::RTData> texcoords = loadVec2fArray(xml->childOpt("texcoords"), numTexCoords);
    size_t numTriangles = 0;  Handle<Device::RTData> triangles = loadVector3iArray(xml->childOpt("triangles"), numTriangles);

    Handle<Device::RTShape> mesh = g_device->rtNewShape("trianglemesh");
    if (numPositions) g_device->rtSetArray(mesh, "positions", "float3", positions, numPositions, sizeof(Vec3f), 0);
    if (numMotions  ) g_device->rtSetArray(mesh, "motions"  , "float3", motions  , numMotions  , sizeof(Vec3f), 0);
    if (numNormals  ) g_device->rtSetArray(mesh, "normals"  , "float3", normals  , numNormals  , sizeof(Vec3f), 0);
    if (numTexCoords) g_device->rtSetArray(mesh, "texcoords", "float2", texcoords, numTexCoords, sizeof(Vec2f ), 0);
    if (numTriangles) g_device->rtSetArray(mesh, "indices"  , "int3"  , triangles, numTriangles, sizeof(Vec3i ), 0);
    g_device->rtSetString(mesh,"accel",g_mesh_accel.c_str());
    g_device->rtSetString(mesh,"builder",g_mesh_builder.c_str());
    g_device->rtSetString(mesh,"traverser",g_mesh_traverser.c_str());
    g_device->rtCommit(mesh);
    g_device->rtClear(mesh);

    return g_device->rtNewShapePrimitive(mesh, material, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadSphere(const Ref<XML>& xml) 
  {
    Handle<Device::RTMaterial> material  = loadMaterial(xml->child("material"));
    Vector3f P = load<Vector3f>(xml->child("position"));
    Vector3f dPdt = load<Vector3f>(xml->childOpt("motion"),Vector3f(0,0,0));

    Handle<Device::RTShape> sphere = g_device->rtNewShape("sphere");
    g_device->rtSetFloat3(sphere, "P", P.x, P.y, P.z);
    g_device->rtSetFloat3(sphere, "dPdt", dPdt.x, dPdt.y, dPdt.z);
    g_device->rtSetFloat1(sphere, "r", load<float>(xml->child("radius")));
    g_device->rtSetInt1(sphere, "numTheta", load<int>(xml->child("numTheta")));
    g_device->rtSetInt1(sphere, "numPhi"  , load<int>(xml->child("numPhi"  )));
    g_device->rtSetString(sphere,"accel",g_mesh_accel.c_str());
    g_device->rtSetString(sphere,"builder",g_mesh_builder.c_str());
    g_device->rtSetString(sphere,"traverser",g_mesh_traverser.c_str());
    g_device->rtCommit(sphere);

    return g_device->rtNewShapePrimitive(sphere, material, copyToArray(transforms.top()));
  }

  Handle<Device::RTPrimitive> XMLLoader::loadDisk(const Ref<XML>& xml) 
  {
    Handle<Device::RTMaterial> material  = loadMaterial(xml->child("material"));
    Vector3f P = load<Vector3f>(xml->child("position"));

    Handle<Device::RTShape> disk = g_device->rtNewShape("disk");
    g_device->rtSetFloat3(disk, "P", P.x, P.y, P.z);
    g_device->rtSetFloat1(disk, "r", load<float>(xml->child("radius")));
    g_device->rtSetFloat1(disk, "h", load<float>(xml->childOpt("height"),0.0f));
    g_device->rtSetInt1(disk, "numTriangles", load<int>(xml->child("numTriangles")));
    g_device->rtSetString(disk,"accel",g_mesh_accel.c_str());
    g_device->rtSetString(disk,"builder",g_mesh_builder.c_str());
    g_device->rtSetString(disk,"traverser",g_mesh_traverser.c_str());
    g_device->rtCommit(disk);

    return g_device->rtNewShapePrimitive(disk, material, copyToArray(transforms.top()));
  }

  std::vector<Handle<Device::RTPrimitive> > XMLLoader::loadTransformNode(const Ref<XML>& xml) 
  {
    std::vector<Handle<Device::RTPrimitive> > group;

    /*! Pre-multiply the incoming transform with the top of the transform stack. */
    transforms.push(transforms.top()*load<AffineSpace3f>(xml->children[0]));

    /*! The transform at the top of the stack will be applied to new lights and shapes. */
    for (size_t i=1; i<xml->children.size(); i++) {
      std::vector<Handle<Device::RTPrimitive> > child = loadScene(xml->children[i]);
      group.insert(group.end(), child.begin(), child.end());
    }

    /*! Remove the transform associated with this node from the transform stack. */
    transforms.pop();  
    return group;
  }

  std::vector<Handle<Device::RTPrimitive> > XMLLoader::loadGroupNode(const Ref<XML>& xml) 
  {
    std::vector<Handle<Device::RTPrimitive> > group;
    
    for (size_t i=0; i<xml->children.size(); i++) {
      std::vector<Handle<Device::RTPrimitive> > child = loadScene(xml->children[i]);
      group.insert(group.end(), child.begin(), child.end());
    }

    return group;
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Loading of scene graph node from XML file
  //////////////////////////////////////////////////////////////////////////////
  
  std::vector<Handle<Device::RTPrimitive> > XMLLoader::loadScene(const Ref<XML>& xml)
  {
    std::vector<Handle<Device::RTPrimitive> > prims;
    
    if (xml->name == "assign") 
    {
      if (xml->parm("type") == "material")
        materialMap[xml->parm("id")] = loadMaterial(xml->child(0));
      else if (xml->parm("type") == "scene")
        sceneMap[xml->parm("id")] = loadScene(xml->child(0));
      else 
        throw std::runtime_error(xml->loc.str()+": unknown type: "+xml->parm("type"));
    }
    else 
    {
      if (xml->name == "xml") {
        prims = loadXML(path + xml->parm("src"));
        for (size_t i=0; i<prims.size(); i++)
          prims[i] = g_device->rtTransformPrimitive(prims[i],copyToArray(transforms.top()));
      }
      else if (xml->name == "obj") {
        prims = loadOBJ(path + xml->parm("src"));
        for (size_t i=0; i<prims.size(); i++)
          prims[i] = g_device->rtTransformPrimitive(prims[i],copyToArray(transforms.top()));
      }
      else if (xml->name == "extern") {
        prims = rtLoadScene(path + xml->parm("src"));
        for (size_t i=0; i<prims.size(); i++)
          prims[i] = g_device->rtTransformPrimitive(prims[i],copyToArray(transforms.top()));
      }
      else if (xml->name == "ref") {
        prims = sceneMap[xml->parm("id")];
        for (size_t i=0; i<prims.size(); i++)
          prims[i] = g_device->rtTransformPrimitive(prims[i],copyToArray(transforms.top()));
      }
      
      else if (xml->name == "PointLight"      ) prims.push_back(loadPointLight      (xml));
      else if (xml->name == "SpotLight"       ) prims.push_back(loadSpotLight       (xml));
      else if (xml->name == "DirectionalLight") prims.push_back(loadDirectionalLight(xml));
      else if (xml->name == "DistantLight"    ) prims.push_back(loadDistantLight    (xml));
      else if (xml->name == "AmbientLight"    ) prims.push_back(loadAmbientLight    (xml));
      else if (xml->name == "TriangleLight"   ) prims.push_back(loadTriangleLight   (xml));
      else if (xml->name == "QuadLight"       ) prims =         loadQuadLight       (xml);
      else if (xml->name == "HDRILight"       ) prims.push_back(loadHDRILight       (xml));

      else if (xml->name == "TriangleMesh"    ) prims.push_back(loadTriangleMesh    (xml));
      else if (xml->name == "Sphere"          ) prims.push_back(loadSphere          (xml));
      else if (xml->name == "Disk"            ) prims.push_back(loadDisk            (xml));
      else if (xml->name == "Group"           ) prims =         loadGroupNode       (xml);
      else if (xml->name == "Transform"       ) prims =         loadTransformNode   (xml);
      
      else throw std::runtime_error(xml->loc.str()+": unknown tag: "+xml->name);
    }

    return prims;
  }

  XMLLoader::XMLLoader(const FileName& fileName) : binFile(NULL)
  {
    path = fileName.path();
    binFileName = fileName.setExt(".bin");
    binFile = fopen(binFileName.c_str(),"rb");

    transforms.push(AffineSpace3f(one));

    Ref<XML> xml = parseXML(fileName);
    if (xml->name != "scene") throw std::runtime_error(xml->loc.str()+": invalid scene tag");
    for (size_t i=0; i<xml->children.size(); i++) {
      std::vector<Handle<Device::RTPrimitive> > prims = loadScene(xml->children[i]);
      model.insert(model.end(), prims.begin(), prims.end());
    }
  }

  XMLLoader::~XMLLoader() {
    if (transforms.size()) transforms.pop();
    if (binFile) fclose(binFile);
    rtClearImageCache();
    rtClearTextureCache();
  }

  std::vector<Handle<Device::RTPrimitive> > loadXML(const FileName &fileName) {
    XMLLoader loader(fileName); return loader.model;
  }
}
