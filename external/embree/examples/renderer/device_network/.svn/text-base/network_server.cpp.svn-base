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

#include "default.h"
#include "network_common.h"
#include "network_server.h"
#include "sys/sysinfo.h"

namespace embree 
{
  NetworkServer::NetworkServer(network::socket_t socket, Device* device, int encoding, bool verbose) 
    : device(device), socket(socket), serverID(0), serverCount(1), encoding(encoding), verbose(verbose)
  {
    try {
      while (true) 
        receive();
    } 
    catch (const network::Disconnect&) {
    }
    catch (const std::exception &e) {
      std::cout << "Error: " << e.what() << std::endl;
    } 
    catch (...) {
      std::cout << "Error: Unknown exception caught." << std::endl;
    }
    delete this;
  }

  NetworkServer::~NetworkServer() 
  {
    /* free all objects of the device */
    for (size_t i=0; i<numRefs.size(); i++) 
      for (int j=0; j<numRefs[i]; j++)
        device->rtDecRef(get<Device::RTHandle>(i));

    /* clear all swapchains */
    swapChains.clear();

    /* delete the device */
    delete device; device = NULL;
  }

  size_t encodeRGBFloat32_to_RGB8(unsigned char *buffer, Col3f *pixels, size_t width, size_t height) 
  {
    for (size_t i=0, j=0 ; i < width * height ; i++) {
      buffer[j++] = (unsigned char) (255.0f * clamp(pixels[i].r));
      buffer[j++] = (unsigned char) (255.0f * clamp(pixels[i].g));
      buffer[j++] = (unsigned char) (255.0f * clamp(pixels[i].b));
    }   
    return(width * height * 3);
  }

  Vec4f encodeRGBE8(const Col3f &pixel) 
  {
    float largest = max(pixel.r, pixel.g, pixel.b);
    int   bound   = cast_f2i(largest) & 0x7F800000;
    int   exp     = (bound >> 23) - 7;
    int   nexp    = (-exp + 254) << 23;
    float scale   = cast_i2f(nexp);
    
    return(Vec4f(scale * pixel.r, scale * pixel.g, scale * pixel.b, float(exp)));
  }
  
  size_t encodeRGBFloat32_to_RGBE8(unsigned char *buffer, Col3f *pixels, size_t width, size_t height) 
  {
    for (size_t i=0, j=0 ; i < width * height ; i++) 
    {
      Vec4f pixel = encodeRGBE8(pixels[i]);
      buffer[j++] = (unsigned char) (pixel.x);
      buffer[j++] = (unsigned char) (pixel.y);
      buffer[j++] = (unsigned char) (pixel.z);
      buffer[j++] = (unsigned char) (pixel.w);
    }   
    return(width * height * 4);
  }

  void NetworkServer::receive() 
  {
    /*! read the magick number */
    if (network::read_int(socket) != magick) 
      throw std::runtime_error("received invalid command block");

    /*! read the command */
    int command = network::read_int(socket);

    /*! decode the command */
    switch (command) 
    {

      /*******************************************************************
                      creation of objects
      *******************************************************************/

    case EMBREE_NEW_CAMERA: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);

      if (verbose) printf("handle %06d = rtNewCamera(%s)\n", id, type.c_str());
      set(id, device->rtNewCamera(type.c_str()));
      break;
    }

    case EMBREE_NEW_DATA: 
    {
      int id = network::read_int(socket);
      size_t bytes = network::read_int(socket);
      char *data = (char*) alignedMalloc(bytes);  
      network::read(socket, data, bytes);
      
      if (verbose) printf("handle %06d = rtNewData(%zu)\n", id, bytes);
      set(id, device->rtNewData("immutable_managed", bytes, data));
      break;
    }

    case EMBREE_NEW_DATA_FROM_FILE: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      std::string fileName = network::read_string(socket);
      size_t offset = network::read_int(socket);
      size_t bytes = network::read_int(socket);
      
      if (verbose) printf("handle %06d = rtNewDataFromFile(%s, %s, %zu, %zu)\n", id, type.c_str(), fileName.c_str(), offset, bytes);
      set(id, device->rtNewDataFromFile(type.c_str(), fileName.c_str(), offset, bytes));
      break;
    } 

    case EMBREE_NEW_FRAMEBUFFER: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      size_t width = network::read_int(socket);
      size_t height = network::read_int(socket);
      size_t depth = network::read_int(socket);
      
      if (verbose) printf("handle %06d = rtNewFrameBuffer(%s, %zu, %zu, %zu)\n", id, type.c_str(), width, height, depth);
      Device::RTFrameBuffer fb = device->rtNewFrameBuffer(type.c_str(), width, height, depth);
      set(id, fb, new SwapChain(this,type,fb,id,width,height,depth,serverID,serverCount));
      break;
    }

    case EMBREE_NEW_IMAGE: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      size_t width = network::read_int(socket);
      size_t height = network::read_int(socket);
      
      size_t bytes = 0;
      if      (!strcmp(type.c_str(), "RGB8"        )) bytes = width * height * 3 * sizeof(char);
      else if (!strcmp(type.c_str(), "RGBA8"       )) bytes = width * height * 4 * sizeof(char);
      else if (!strcmp(type.c_str(), "RGB_FLOAT32" )) bytes = width * height * 3 * sizeof(float);
      else if (!strcmp(type.c_str(), "RGBA_FLOAT32")) bytes = width * height * 4 * sizeof(float);
      else throw std::runtime_error("unknown image type: " + std::string(type));
      char* data = (char*) malloc(bytes);  
      network::read(socket, data, bytes);
      
      if (verbose) printf("handle %06d = rtNewImage(%s, %zu, %zu)\n", id, type.c_str(), width, height);
      set(id, device->rtNewImage(type.c_str(), width, height, data, true));  
	  free(data);
      break;
    }
    
    case EMBREE_NEW_IMAGE_FROM_FILE: 
    {
      int id = network::read_int(socket);
      std::string fileName = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewImageFromFile(%s)\n", id, fileName.c_str());
      set(id, device->rtNewImageFromFile(fileName.c_str()));
      break;
    }
      
    case EMBREE_NEW_LIGHT: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewLight(%s)\n", id, type.c_str());
      set(id,device->rtNewLight(type.c_str()));
      break;
    } 

    case EMBREE_NEW_LIGHT_PRIMITIVE: 
    {
      int id = network::read_int(socket);
      int lightID = network::read_int(socket);
      int materialID = network::read_int(socket);
      float transform[12]; network::read(socket, transform, 12 * sizeof(float));
      
      if (verbose) printf("handle %06d = rtNewLightPrimitive(%06d, %06d, %s)\n", id, lightID, materialID, std::stringOf(copyFromArray(transform)).c_str());
      set(id, device->rtNewLightPrimitive(get<Device::RTLight>(lightID), get<Device::RTMaterial>(materialID), (float*) transform));
      break;
    } 

    case EMBREE_NEW_MATERIAL: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewMaterial(%s)\n", id, type.c_str());
      set(id, device->rtNewMaterial(type.c_str()));
      break;
    } 

    case EMBREE_NEW_RENDERER: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewRenderer(%s)\n", id, type.c_str());
      set(id, device->rtNewRenderer(type.c_str()));
      break;
    } 

    case EMBREE_NEW_SCENE: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      if (verbose) printf("handle %06d = rtNewScene(%s)\n", id, type.c_str());
      set(id, device->rtNewScene(type.c_str()));
      break;
    } 

    case EMBREE_SET_SCENE_PRIMITIVE:
    {
      int sceneID = network::read_int(socket);
      int slot    = network::read_int(socket);
      int primID  = network::read_int(socket);
      if (verbose) printf("rtSetScenePrimitive(%06d, %06d, %06d)\n", sceneID, slot, primID);
      device->rtSetPrimitive(get<Device::RTScene>(sceneID),slot,get<Device::RTPrimitive>(primID));
      break;
    }

    case EMBREE_NEW_SHAPE: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewShape(%s)\n", id, type.c_str());
      set(id, device->rtNewShape(type.c_str()));
      break;
    } 

    case EMBREE_NEW_SHAPE_PRIMITIVE: 
    {
      int id = network::read_int(socket);
      int shapeID = network::read_int(socket);
      int materialID = network::read_int(socket);
      float transform[12];  network::read(socket, transform, 12 * sizeof(float));
      if (verbose) printf("handle %06d = rtNewShapePrimitive(%06d, %06d, %s)\n", id, shapeID, materialID, std::stringOf(copyFromArray(transform)).c_str());
      set(id, device->rtNewShapePrimitive(get<Device::RTShape>(shapeID), get<Device::RTMaterial>(materialID), (float *) transform));
      break;
    } 

    case EMBREE_NEW_TEXTURE: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewTexture(%s)\n", id, type.c_str());
      set(id, device->rtNewTexture(type.c_str()));
      break;
    }

    case EMBREE_NEW_TONEMAPPER: 
    {
      int id = network::read_int(socket);
      std::string type = network::read_string(socket);
      
      if (verbose) printf("handle %06d = rtNewToneMapper(%s)\n", id, type.c_str());
      set(id, device->rtNewToneMapper(type.c_str()));
      break;
    }

    case EMBREE_TRANSFORM_PRIMITIVE: 
    {
      int id = network::read_int(socket);
      int primID = network::read_int(socket);
      float transform[12]; network::read(socket, transform, 12 * sizeof(float));
      
      if (verbose) printf("handle %06d = rtTransformPrimitive(%06d, %s)\n", id, primID, std::stringOf(copyFromArray(transform)).c_str());
      set(id, device->rtTransformPrimitive(get<Device::RTPrimitive>(primID), (float*) transform));
      break;
    }
    
    case EMBREE_CLEAR: 
    {
      int id = network::read_int(socket);
      if (verbose) printf("rtClear(%06d)\n", id);
      device->rtClear(get<Device::RTHandle>(id));
      break;
    } 

    case EMBREE_COMMIT: 
    {
      int id = network::read_int(socket);
      if (verbose) printf("rtCommit(%06d)\n", id);
      device->rtCommit(get<Device::RTHandle>(id));
      break;
    } 

    case EMBREE_DECREF: 
    {
      int id = network::read_int(socket);
      if (verbose) printf("rtDecRef(%06d)\n", id);
      device->rtDecRef(get<Device::RTHandle>(id));
      if (--numRefs[id] == 0) {
        handles[id] = NULL;
        swapChains[id] = NULL;
      }
      if (numRefs[id] < 0) 
        throw std::runtime_error("reference counter for object got negative");
      break;
    } 

    case EMBREE_INCREF: 
    {
      int id = network::read_int(socket);
      if (verbose) printf("rtIncRef(%06d)\n", id);
      device->rtIncRef(get<Device::RTHandle>(id));
      numRefs[id]++;
      break;
    }

  /*******************************************************************
                    setting of parameters
  *******************************************************************/

    case EMBREE_SET_BOOL1: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      bool x = network::read_bool(socket);
      
      if (verbose) printf("rtSetBool1(%06d, %s, [%d])\n", id, property.c_str(), x);
      device->rtSetBool1(get<Device::RTHandle>(id), property.c_str(), x);
      break;
    } 

    case EMBREE_SET_BOOL2: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      bool x = network::read_bool(socket);
      bool y = network::read_bool(socket);
      
      if (verbose) printf("rtSetBool2(%06d, %s, [%d, %d])\n", id, property.c_str(), x, y);
      device->rtSetBool2(get<Device::RTHandle>(id), property.c_str(), x, y);
      break;
    } 

    case EMBREE_SET_BOOL3: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      bool x = network::read_bool(socket);
      bool y = network::read_bool(socket);
      bool z = network::read_bool(socket);
      
      if (verbose) printf("rtSetBool3(%06d, %s, [%d, %d, %d])\n", id, property.c_str(), x, y, z);
      device->rtSetBool3(get<Device::RTHandle>(id), property.c_str(), x, y, z);
      break;
    }
     
    case EMBREE_SET_BOOL4: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      bool x = network::read_bool(socket);
      bool y = network::read_bool(socket);
      bool z = network::read_bool(socket);
      bool w = network::read_bool(socket);
      
      if (verbose) printf("rtSetBool4(%06d, %s, [%d, %d, %d, %d])\n", id, property.c_str(), x, y, z, w);
      device->rtSetBool4(get<Device::RTHandle>(id), property.c_str(), x, y, z, w);
      break;
    } 

    case EMBREE_SET_FLOAT1: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      float x = network::read_float(socket);
      
      if (verbose) printf("rtSetFloat1(%06d, %s, [%f])\n", id, property.c_str(), x);
      device->rtSetFloat1(get<Device::RTHandle>(id), property.c_str(), x);
      break;
    } 

    case EMBREE_SET_FLOAT2: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      float x = network::read_float(socket);
      float y = network::read_float(socket);
      
      if (verbose) printf("rtSetFloat2(%06d, %s, [%f, %f])\n", id, property.c_str(), x, y);
      device->rtSetFloat2(get<Device::RTHandle>(id), property.c_str(), x, y);
      break;
    } 

    case EMBREE_SET_FLOAT3: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      float x = network::read_float(socket);
      float y = network::read_float(socket);
      float z = network::read_float(socket);
      
      if (verbose) printf("rtSetFloat3(%06d, %s, [%f, %f, %f])\n", id, property.c_str(), x, y, z);
      device->rtSetFloat3(get<Device::RTHandle>(id), property.c_str(), x, y, z);
      break;
    } 

    case EMBREE_SET_FLOAT4: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      float x = network::read_float(socket);
      float y = network::read_float(socket);
      float z = network::read_float(socket);
      float w = network::read_float(socket);
      
      if (verbose) printf("rtSetFloat4(%06d, %s, [%f, %f, %f, %f])\n", id, property.c_str(), x, y, z, w);
      device->rtSetFloat4(get<Device::RTHandle>(id), property.c_str(), x, y, z, w);
      break;
    } 

    case EMBREE_SET_INT1: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int x = network::read_int(socket);

      if (id == 0) {
        if      (property == "serverID"   ) serverID = x;
        else if (property == "serverCount") serverCount = x;
      }
      
      if (verbose) printf("rtSetInt1(%06d, %s, [%d])\n", id, property.c_str(), x);
      device->rtSetInt1(get<Device::RTHandle>(id), property.c_str(), x);
      break;
      
    } 

    case EMBREE_SET_INT2: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int x = network::read_int(socket);
      int y = network::read_int(socket);
      
      if (verbose) printf("rtSetInt2(%06d, %s, [%d, %d])\n", id, property.c_str(), x, y);
      device->rtSetInt2(get<Device::RTHandle>(id), property.c_str(), x, y);
      break;
    } 

    case EMBREE_SET_INT3: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int x = network::read_int(socket);
      int y = network::read_int(socket);
      int z = network::read_int(socket);
      
      if (verbose) printf("rtSetInt3(%06d, %s, [%d, %d, %d])\n", id, property.c_str(), x, y, z);
      device->rtSetInt3(get<Device::RTHandle>(id), property.c_str(), x, y, z);
      break;
    } 

    case EMBREE_SET_INT4: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int x = network::read_int(socket);
      int y = network::read_int(socket);
      int z = network::read_int(socket);
      int w = network::read_int(socket);
      
      if (verbose) printf("rtSetInt4(%06d, %s, [%d, %d, %d, %d])\n", id, property.c_str(), x, y, z, w);
      device->rtSetInt4(get<Device::RTHandle>(id), property.c_str(), x, y, z, w);
      break;
    } 

    case EMBREE_SET_ARRAY: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      std::string type = network::read_string(socket);
      int dataID = network::read_int(socket);
      size_t size = network::read_int(socket);
      size_t stride = network::read_int(socket);
      size_t ofs = network::read_int(socket);

      if (verbose) printf("rtSetArray(%06d, %s, %s, %06d, %zu, %zu, %zu)\n", id, property.c_str(), type.c_str(), dataID, size, stride, ofs);
      device->rtSetArray(get<Device::RTHandle>(id), property.c_str(), type.c_str(), get<Device::RTData>(dataID), size, stride, ofs);
      break;
    }
    
    case EMBREE_SET_STRING: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      std::string string = network::read_string(socket);
      
      if (verbose) printf("rtSetString(%06d, %s, %s)\n", id, property.c_str(), string.c_str());
      device->rtSetString(get<Device::RTHandle>(id), property.c_str(), string.c_str());
      break;
    } 

    case EMBREE_SET_IMAGE:
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int imageID = network::read_int(socket);
      
      if (verbose) printf("rtSetImage(%06d, %s, %06d)\n", id, property.c_str(), imageID);
      device->rtSetImage(get<Device::RTHandle>(id), property.c_str(), get<Device::RTImage>(imageID));
      break;
    } 

    case EMBREE_SET_TEXTURE: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      int textureID = network::read_int(socket);
      
      if (verbose) printf("rtSetTexture(%06d, %s, %06d)\n", id, property.c_str(), textureID);
      device->rtSetTexture(get<Device::RTHandle>(id), property.c_str(), get<Device::RTTexture>(textureID));
      break;
    } 

    case EMBREE_SET_TRANSFORM: 
    {
      int id = network::read_int(socket);
      std::string property = network::read_string(socket);
      float transform[12];  network::read(socket, transform, 12 * sizeof(float));
      
      if (verbose) printf("rtSetTransform(%06d, %s, %s)\n", id, property.c_str(), std::stringOf(copyFromArray(transform)).c_str());
      device->rtSetTransform(get<Device::RTHandle>(id), property.c_str(), (float *) transform);
      break;
    }

  /*******************************************************************
                          render calls
  *******************************************************************/

    case EMBREE_RENDER_FRAME: 
    {
      int rendererID = network::read_int(socket);
      int cameraID = network::read_int(socket);
      int sceneID = network::read_int(socket);
      int toneMapperID = network::read_int(socket);
      int frameBufferID = network::read_int(socket);
      int accumulate = network::read_int(socket);
      
      /* verbosity */
      if (verbose) printf("rtRenderFrame(%06d, %06d, %06d, %06d, %06d, %d)\n", rendererID, cameraID, sceneID, toneMapperID, frameBufferID, accumulate);

      /* render the frame */
      device->rtRenderFrame(get<Device::RTRenderer>(rendererID),
                            get<Device::RTCamera>(cameraID),
                            get<Device::RTScene>(sceneID),
                            get<Device::RTToneMapper>(toneMapperID),
                            get<Device::RTFrameBuffer>(frameBufferID),
                            accumulate);
      break;
    } 

    case EMBREE_PICK: 
    {
      int cameraID = network::read_int(socket);
      float x = network::read_float(socket);
      float y = network::read_float(socket);
      int sceneID = network::read_int(socket);
      if (verbose) printf("rtPick([%f, %f], %06d, %06d)\n", x, y, cameraID, sceneID);
      
      /*! find the hit point */
      float px, py, pz;
      bool isHit = device->rtPick(get<Device::RTCamera>(cameraID), x, y, get<Device::RTScene>(sceneID), px, py, pz);
      if (verbose) printf("  returning %s, [%f, %f, %f]\n", (isHit ? "TRUE" : "FALSE"), px, py, pz);
      
      /*! return a packet with the hit information */
      network::write(socket, (int) magick);
      network::write(socket, (int) EMBREE_PICK_RESULT);
      network::write(socket, (bool) isHit);
      network::write(socket, (float) px);
      network::write(socket, (float) py);
      network::write(socket, (float) pz);
      network::flush(socket);
      break;
    } 

    case EMBREE_SWAP_BUFFERS: 
    {
      int frameBufferID = network::read_int(socket);
      if (verbose) printf("rtSwapBuffers(%06d)\n", frameBufferID);
      device->rtSwapBuffers(get<Device::RTFrameBuffer>(frameBufferID));
      swapChains[frameBufferID]->send();
      swapChains[frameBufferID]->nextWriteBuffer();
      break;
    } 

    default: { 
      throw std::runtime_error("unknown network command: " + std::stringOf(command)); 
    }
    }
  }

  void NetworkServer::SwapChain::nextWriteBuffer() {
    writeID++;
  }

  void NetworkServer::SwapChain::nextReadBuffer() {
    readID++;
  }

  void NetworkServer::SwapChain::send()
  {
    /*! map framebuffer data */
    int swapID = readID%numBuffers;
    void* data = server->device->rtMapFrameBuffer(frameBuffer, swapID);

    /*! count active rows */
    int height1 = 0;
    for (size_t y=0; y<height; y++)
      if (((y>>2)+serverID) % serverCount == 0)
        height1++;

    /*! encode image */
    size_t bytes = 0;    
    void* sendbuffer = NULL;
    switch (server->encoding) 
    {
    case EMBREE_FRAME_DATA_NATIVE: {
      sendbuffer = data;
      if      (type == "RGB_FLOAT32") bytes = 3*width*height1*sizeof(float);
      else if (type == "RGBA8"      ) bytes = 4*width*height1;       
      else if (type == "RGB8"       ) bytes = ((3*width+3)/4*4)*height1;
      else throw std::runtime_error("unsupported framebuffer format: "+type);
      break;
    }

    case EMBREE_FRAME_DATA_RGB8: {
      sendbuffer = encoded;
      if (type == "RGB_FLOAT32") 
        bytes = encodeRGBFloat32_to_RGB8 (encoded, (Col3f*)data, width, height1); 
      else if (type == "RGBA8"      ) {
        char* src = (char*)data;
        char* dst = (char*)encoded;
        for (size_t i=0; i<width*height1; i++) {
          dst[0] = src[0];
          dst[1] = src[1];
          dst[2] = src[2];
          dst+=3; src+=4;
        }
        bytes = 3*width*height1;
      }
      else throw std::runtime_error("unsupported framebuffer format: "+type);
      break;
    }

    case EMBREE_FRAME_DATA_RGBE8: {
      sendbuffer = encoded;
      if (type == "RGB_FLOAT32") bytes = encodeRGBFloat32_to_RGBE8(encoded, (Col3f*)data, width, height1);
      else throw std::runtime_error("unsupported framebuffer format: "+type);
      break;
    }
    default: throw std::runtime_error("invalid encoding");
    }

    /*! send the encoded frame buffer */
    network::write(server->socket, (int) magick);
    network::write(server->socket, (int) server->encoding);
    network::write(server->socket, (int) frameBufferID);
    network::write(server->socket, (int) swapID);
    network::write(server->socket, (int) width);  
    network::write(server->socket, (int) height1);  
    network::write(server->socket, sendbuffer, bytes);  
    network::flush(server->socket);  

    /*! unmap framebuffer data */
    server->device->rtUnmapFrameBuffer(frameBuffer, swapID);

    /*! goto next read buffer */
    nextReadBuffer();
  }

  /*! NetworkServer bookkeeping methods ================================== */
  
  template<typename T> T NetworkServer::get(size_t id) {
    if (id == 0) return (T) NULL;
    return((T) handles[id]);
  }

  void NetworkServer::set(int id, Device::RTHandle handle, SwapChain* swapChain) 
  {
    if (handles.size() <= (size_t) id) { 
      handles.resize(id + 1);  
      swapChains.resize(id + 1);
      numRefs.resize(id + 1); 
    }
    handles[id] = handle;
    swapChains[id] = swapChain;
    numRefs[id] = 1;
  }
}

