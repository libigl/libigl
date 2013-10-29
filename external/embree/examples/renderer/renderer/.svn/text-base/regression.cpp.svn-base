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

#include "regression.h"
#include "sys/intrinsics.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include "math/affinespace.h"
#include "image/image.h"
#include <vector>

namespace embree
{
  extern std::string g_accel;
  extern std::string g_builder;
  extern std::string g_traverser;

  Handle<Device::RTImage> createRandomImage(Device *device, size_t width, size_t height)
  {
    if (random<bool>())
    {
      char* data = new char[3 * width * height];
      for (size_t y=0; y < height; y++) {
        for (size_t x=0; x < width; x++) {
          size_t ofs = y * width + x;
          data[3 * ofs + 0] = char(x * y);
          data[3 * ofs + 1] = char(y * x);
          data[3 * ofs + 2] = char(x + y);
        }
      }
      Handle<Device::RTImage> image = device->rtNewImage("RGB8", width, height, data);
      delete[] data;
      return(image);
    }
    else {
      Col3f* data = new Col3f[3 * width * height];
      for (size_t y=0; y < height; y++) {
        for (size_t x=0; x < width; x++) {
          size_t ofs = y * width + x;
          data[ofs].r = char(x * y) / 255.0f;
          data[ofs].g = char(y * x) / 255.0f;
          data[ofs].b = char(x + y) / 255.0f;
        }
      }
      Handle<Device::RTImage> image = device->rtNewImage("RGB_FLOAT32", width, height, data);
      delete[] data;
      return(image);
    }
  }

  Handle<Device::RTTexture> createRandomTexture(Device *device, size_t width, size_t height)
  {
    Handle<Device::RTTexture> texture = device->rtNewTexture("image");
    device->rtSetImage(texture, "image", createRandomImage(device, width, height));
    device->rtCommit(texture);
    return(texture);
  }

  Handle<Device::RTMaterial> createRandomMaterial(Device *device)
  {
    switch (random<int>() % 8)
    {
    case 0: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Matte");
      device->rtSetFloat3(material, "reflectance", random<float>(), random<float>(), random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 1: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Plastic");
      device->rtSetFloat3(material, "pigmentColor", random<float>(), random<float>(), random<float>());
      device->rtSetFloat1(material, "eta", 1.0f + random<float>());
      device->rtSetFloat1(material, "roughness", 0.1f * random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 2: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Dielectric");
      device->rtSetFloat3(material, "transmission", 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f);
      device->rtSetFloat1(material, "etaOutside", 1.0f);
      device->rtSetFloat1(material, "etaInside", 1.0f + random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 3: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("ThinDielectric");
      device->rtSetFloat3(material, "transmission", 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f);
      device->rtSetFloat1(material, "eta", 1.0f + random<float>());
      device->rtSetFloat1(material, "thickness", 0.5f * random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 4: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Mirror");
      device->rtSetFloat3(material, "reflectance", 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f);
      device->rtCommit(material);
      return(material);
    }
    case 5: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Metal");
      device->rtSetFloat3(material, "reflectance", 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f);
      device->rtSetFloat3(material, "eta", 1.0f + random<float>(), 1.0f + random<float>(), 1.0f + random<float>());
      device->rtSetFloat3(material, "k", 0.3f * random<float>(), 0.3f * random<float>(), 0.3f * random<float>());
      device->rtSetFloat1(material, "roughness", 0.3f * random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 6: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("MetallicPaint");
      device->rtSetFloat3(material, "shadeColor", 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f, 0.5f * random<float>() + 0.5f);
      device->rtSetFloat3(material, "glitterColor", random<float>(), random<float>(), random<float>());
      device->rtSetFloat1(material, "glitterSpread", 0.5f + random<float>());
      device->rtSetFloat1(material, "eta", 1.0f + random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 7: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("MatteTextured");
      device->rtSetTexture(material, "Kd", createRandomTexture(device, 32, 32));
      device->rtSetFloat2(material, "s0", random<float>(), random<float>());
      device->rtSetFloat2(material, "ds", 5.0f * random<float>(), 5.0f * random<float>());
      device->rtCommit(material);
      return(material);
    }
    case 8: {
      Handle<Device::RTMaterial> material = device->rtNewMaterial("Obj");
      device->rtSetFloat1(material, "d", random<float>());
      if (random<bool>()) device->rtSetTexture(material, "map_d", createRandomTexture(device, 32, 32));
      device->rtSetFloat3(material, "Tf", random<float>(), random<float>(), random<float>());
      device->rtSetFloat3(material, "Kd", random<float>(), random<float>(), random<float>());
      if (random<bool>()) device->rtSetTexture(material, "map_Kd", createRandomTexture(device, 32, 32));
      device->rtSetFloat3(material, "Ks", random<float>(), random<float>(), random<float>());
      if (random<bool>()) device->rtSetTexture(material, "map_Ks", createRandomTexture(device, 32, 32));
      device->rtSetFloat1(material, "Ns", 10.0f * random<float>());
      if (random<bool>()) device->rtSetTexture(material, "map_Ns", createRandomTexture(device, 32, 32));
      device->rtCommit(material);
      return(material);
    }
    }
    return(NULL);
  }

  Handle<Device::RTLight> createRandomLight(Device *device)
  {
    Handle<Device::RTLight> light = device->rtNewLight("ambientlight");
    device->rtSetFloat3(light, "L", 1.0f, 1.0f, 1.0f);
    device->rtCommit(light);
    return(light);
  }

  Handle<Device::RTShape> createRandomShape(Device *device, size_t numTriangles)
  {
    if (numTriangles < 20)
    {
      std::vector<Vec3fa> positions;
      std::vector<Vec2f> texcoords;
      std::vector<Vec3i>   indices;

      Vector3f pos = 2.0f * Vector3f(random<float>(), random<float>(), random<float>()) - Vector3f(1.0f);
      for (size_t i=0; i < numTriangles; i++) {
        positions.push_back(pos + 0.3f * Vec3fa(random<float>(), random<float>(), random<float>()));
        texcoords.push_back(Vec2f(random<float>(), random<float>()));
        indices  .push_back(Vec3i(random<int>() % numTriangles, random<int>() % numTriangles, random<int>() % numTriangles));
      }

      Handle<Device::RTData> dataPositions = device->rtNewData("immutable", positions.size() * sizeof(Vec3fa), positions.size() ? &positions[0] : NULL);
      Handle<Device::RTData> dataTexCoords = device->rtNewData("immutable", texcoords.size() * sizeof(Vec2f), texcoords.size() ? &texcoords[0] : NULL);
      Handle<Device::RTData> dataIndices   = device->rtNewData("immutable", indices.size()   * sizeof(Vec3i), indices.size()   ? &indices[0]   : NULL);

      Handle<Device::RTShape> shape = device->rtNewShape("trianglemesh");
      device->rtSetArray(shape, "positions", "float3", dataPositions, positions.size(), sizeof(Vec3fa), 0);
      device->rtSetArray(shape, "texcoords", "float2", dataTexCoords, texcoords.size(), sizeof(Vec2f), 0);
      device->rtSetArray(shape, "indices"  , "int3"  , dataIndices,   indices.size(),   sizeof(Vec3i), 0);
      device->rtCommit(shape);
      return(shape);
    }
    else
    {
      Handle<Device::RTShape> shape = device->rtNewShape("sphere");
      device->rtSetFloat3(shape, "P", 2.0f * random<float>() - 1.0f, 2.0f * random<float>() - 1.0f, 2.0f * random<float>() - 1.0f);
      device->rtSetFloat1(shape, "r", 0.2f * random<float>());
      device->rtSetInt1(shape, "numTheta", (int) numTriangles / 20);
      device->rtSetInt1(shape, "numPhi", 20);
      device->rtCommit(shape);
      return(shape);
    }
  }

  Handle<Device::RTScene> createRandomScene(Device *device, size_t numLights, size_t numObjects, size_t numTriangles)
  {
    std::vector<Device::RTPrimitive> prims;

    //for (size_t i=0; i < numLights; i++)
    prims.push_back(device->rtNewLightPrimitive(createRandomLight(device), NULL, copyToArray(AffineSpace3f(one))));

    for (size_t i=0; i < numObjects; i++) {
      size_t s = numTriangles ? random<int>() % numTriangles : 0;
      prims.push_back(device->rtNewShapePrimitive(createRandomShape(device, s), createRandomMaterial(device), copyToArray(AffineSpace3f(one))));
    }

    Handle<Device::RTScene> scene = device->rtNewScene("default");
    g_device->rtSetString(scene,"accel",g_accel.c_str());
    g_device->rtSetString(scene,"builder",g_builder.c_str());
    g_device->rtSetString(scene,"traverser",g_traverser.c_str());
    for (size_t i=0; i<prims.size(); i++) g_device->rtSetPrimitive(scene,i,prims[i]);
    g_device->rtCommit(scene);
    
    for (size_t i=0; i<prims.size(); i++)
      device->rtDecRef(prims[i]);

    return scene;
  }
}

