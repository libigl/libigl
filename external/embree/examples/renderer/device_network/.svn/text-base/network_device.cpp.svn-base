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

#include "network_device.h"
#include "image/image.h"

#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  #define strcasecmp lstrcmpiA
#endif

namespace embree
{
  /*******************************************************************
                  creation of device
  *******************************************************************/
  
#if !defined(NETWORK_NO_CREATE)

  __dllexport Device* create(const char* parms, size_t numThreads, size_t verbose) 
  {
    /*! stores all client connections */
    std::vector<network::socket_t> clients;
    
    /*! create copy of the server string */
    size_t len = strlen(parms);
    char* servers = new char[len+1];
    strcpy(servers,parms);
    
    /*! parse hostname:port list */
    char* host = strtok (servers," \n\r");
    while (host) {
      unsigned short port = 8484;
      char* sep = strchr(host,':');
      if (sep) {
        *sep = 0;
        port = (unsigned short) atoi(sep+1);
      }
      clients.push_back(network::connect(host,port));
      host = strtok (NULL," \n\r");
    }
    delete[] servers;
    
    return new NetworkDevice(clients);
  }

#endif

  namespace network
  {
    /*! writes a Vector3f value to the socket */
    void write(socket_t socket, Vector3f value) {
      write(socket,value.x);
      write(socket,value.y);
      write(socket,value.z);
    }
  }

  NetworkDevice::NetworkDevice () 
    : terminate(false), serverID(0), nextHandle(1) 
  {
  }

  NetworkDevice::NetworkDevice (const std::vector<network::socket_t>& servers) 
    : terminate(false), serverID(0), servers(servers), nextHandle(1) 
  {
    init();
  }

  void NetworkDevice::init ()
  {
    /* dummy 0 handle */
    counters.push_back(0);
    buffers.push_back(NULL);

    /* barrier to wait for pick results */
    pick.barrier.init(2);
    receiveBarrier.init(1+servers.size());

    /* send each server its ID and total count */
    for (size_t i=0; i<servers.size(); i++)
    {
      network::write(servers[i],(int)magick);
      network::write(servers[i],(int)EMBREE_SET_INT1);
      network::write(servers[i],(int)0);
      network::write(servers[i],std::string("serverID"));
      network::write(servers[i],(int)i);

      network::write(servers[i],(int)magick);
      network::write(servers[i],(int)EMBREE_SET_INT1);
      network::write(servers[i],(int)0);
      network::write(servers[i],std::string("serverCount"));
      network::write(servers[i],(int)servers.size());
    }
  }
  
  NetworkDevice::~NetworkDevice () 
  {
    /* enter exit mode */
    terminate = true;

    /* close connection to all servers */
    for (size_t i=0; i<servers.size(); i++)
      network::close(servers[i]);
  }

  /*******************************************************************
                     network broadcast to servers
  *******************************************************************/
  
  void NetworkDevice::broadcast(bool value) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],value);
  }
  
  void NetworkDevice::broadcast(char value) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],value);
  }
  
  void NetworkDevice::broadcast(int value) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],value);
  }
  
  void NetworkDevice::broadcast(float value) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],value);
  }
  
  void NetworkDevice::broadcast(const std::string& value) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],value);
  }

  void NetworkDevice::broadcast(const Vector3f& value) {
    broadcast(value.x);
    broadcast(value.y);
    broadcast(value.z);
  }     

  void NetworkDevice::broadcast(const void* data, size_t bytes) {
    for (size_t i=0; i<servers.size(); i++) 
      network::write(servers[i],data,bytes);
  }

  void NetworkDevice::flush() {
    for (size_t i=0; i<servers.size(); i++) 
      network::flush(servers[i]);
  }

  /*******************************************************************
                     handle ID allocations
  *******************************************************************/

  int NetworkDevice::allocHandle() 
  {
    if (pool.empty()) {
      pool.push_back(nextHandle++);
      counters.push_back(0);
      buffers.push_back(NULL);
    }
    int id = pool.back();
    counters[id] = 1;
    pool.pop_back();
    return id;
  }
  
  void NetworkDevice::incRef(int id) {
    Lock<MutexSys> lock(handleMutex);
    counters[id]++;
  }
  
  void NetworkDevice::decRef(int id) 
  {
    Lock<MutexSys> lock(handleMutex);
    if (--counters[id] == 0) {
      pool.push_back((int)id);
      buffers[id] = null;
    }
  }
  
  /*******************************************************************
                 multi-threaded receive from servers
  *******************************************************************/

  void NetworkDevice::_receive(void* data)
  {
    try { 
      ((NetworkDevice*)data)->receive();
    }
    catch (const network::Disconnect&) {
    }
    catch (const std::exception& e) {
      if (!((NetworkDevice*)data)->terminate) 
        std::cout << "Error: " << e.what() << std::endl;
    }
  }

  void NetworkDevice::receiveCommand(int id) 
  {
    /*! test if network packet is valid */
    int _magick = network::read_int(servers[id]);
    if (_magick != magick) 
      throw std::runtime_error("receiveCommand: received invalid command block: "+std::stringOf((void*)_magick)+" != "+std::stringOf((void*)magick));

    /*! decode commands */
    int cmd = network::read_int(servers[id]);
    switch (cmd) 
    {
    /*! receiving pick result */
    case EMBREE_PICK_RESULT: {
      pick.hit   = network::read_bool(servers[id]);
      pick.pos.x = network::read_float(servers[id]);
      pick.pos.y = network::read_float(servers[id]);
      pick.pos.z = network::read_float(servers[id]);
      break;
    }

    /*! receiving frames */
    case EMBREE_FRAME_DATA_NATIVE:
    case EMBREE_FRAME_DATA_RGB8: 
    case EMBREE_FRAME_DATA_RGBE8: 
    case EMBREE_FRAME_DATA_RGB_FLOAT:
    case EMBREE_FRAME_DATA_DXT1: {
      
      /*! read framebuffer to store into */
      int fbID  = network::read_int(servers[id]);
      int bufID = network::read_int(servers[id]);
      /*int width = */ network::read_int(servers[id]);
      /*int height = */ network::read_int(servers[id]);
      Ref<FrameBuffer>& buffer = buffers[fbID]->buffer(bufID);
      const std::string format = buffers[fbID]->getFormat();

      /*! read tile location */
      size_t width  = buffer->getWidth();
      size_t height = buffer->getHeight();
      const float rcp255 = 1.0f/255.0f;

      /*! read framebuffer */
      switch (cmd) 
      {
      case EMBREE_FRAME_DATA_NATIVE: {
        
        /*! calculate size of one row in bytes */
        int rowBytes = 0;
        if      (format == "RGB_FLOAT32") rowBytes = 3*width*sizeof(float);
        else if (format == "RGBA8"      ) rowBytes = 4*width*sizeof(char);
        else if (format == "RGB8"       ) rowBytes = (3*width*sizeof(char)+3)/4*4;
        else throw std::runtime_error("unknown framebuffer format");

        /*! store data directly into framebuffer */
        char* dst = (char*) buffer->getData();
        for (size_t y=0; y<height; y++, dst+=rowBytes) {
          if (((y>>2)+id) % servers.size()) continue;
          network::read(servers[id],dst,rowBytes);
        }
        break;
      }

      case EMBREE_FRAME_DATA_RGB8:
        for (size_t y=0; y<height; y++) {
          if (((y>>2)+id) % servers.size()) continue;
          for (size_t x=0; x<width; x++) { 
	    float r = float((unsigned char)network::read_char(servers[id]))*rcp255;
	    float g = float((unsigned char)network::read_char(servers[id]))*rcp255;
	    float b = float((unsigned char)network::read_char(servers[id]))*rcp255;
            buffer->set(x,y,Color(r,g,b));
          }
        }
        break;

      case EMBREE_FRAME_DATA_RGBE8:
        for (size_t y=0; y<height; y++) {
          if (((y>>2)+id) % servers.size()) continue;
          for (size_t x=0; x<width; x++) { 
            Vector3f c = decodeRGBE8(network::read_int(servers[id]));
            buffer->set(x,y,Color(c.x,c.y,c.z));
          }
        }
        break;

      case EMBREE_FRAME_DATA_RGB_FLOAT:
        for (size_t y=0; y<height; y++) {
          if (((y>>2)+id) % servers.size()) continue;
          for (size_t x=0; x<width; x++) { 
            float r = network::read_float(servers[id]);
            float g = network::read_float(servers[id]);
            float b = network::read_float(servers[id]);
            buffer->set(x,y,Color(r,g,b));
          }
        }
        break;

      default: 
        throw std::runtime_error("unknown frame encoding");
      }
      break;
    }

    default:
      throw std::runtime_error("unknown command received: "+std::stringOf(cmd));
    }
  }

  void NetworkDevice::receive()  
  {
    int id = serverID++;
    pick.barrier.wait(); // to assign IDs in correct order

    while (true) 
      receiveCommand(id);
  }

  /*******************************************************************
                      creation of objects
  *******************************************************************/

  Device::RTCamera NetworkDevice::rtNewCamera(const char* type) 
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_CAMERA);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTCamera)(long)id;
  }

  Device::RTData NetworkDevice::rtNewData(const char* type, size_t bytes, const void* data) 
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_DATA);
    broadcast((int)id);
    broadcast((int)bytes);
    broadcast(data,bytes);
    flush();

    if (!strcasecmp(type,"immutable_managed")) 
      alignedFree((void*)data);

    return (Device::RTData)(long)id;
  }

  Device::RTData NetworkDevice::rtNewDataFromFile(const char* type, const char* fileName, size_t offset, size_t bytes)
  {
    if (strcasecmp(type,"immutable"))
      throw std::runtime_error("unknown data type: "+(std::string)type);
    
    /*! load data remote */
    if (strncmp(fileName,"server:",7) == 0) 
    {
      int id = allocHandle();
      broadcast((int)magick);
      broadcast((int)EMBREE_NEW_DATA_FROM_FILE);
      broadcast((int)id);
      broadcast((std::string)type);
      broadcast((std::string)(fileName+7));
      broadcast((int)offset);
      broadcast((int)bytes);
      flush();

      return (Device::RTData)(long)id;
    }

    /*! load data locally */
    else 
    {
      /*! read data from file */
      FILE* file = fopen(fileName,"rb");
      if (!file) throw std::runtime_error("cannot open file "+(std::string)fileName);
      fseek(file,(long)offset,SEEK_SET);

      char* data = (char*) alignedMalloc(bytes);
      if (bytes != fread(data,1,sizeof(bytes),file))
        throw std::runtime_error("error filling data buffer from file");
      fclose(file);

      return rtNewData("immutable_managed", bytes, data);
    }
  }

  Device::RTImage NetworkDevice::rtNewImage(const char* type, size_t width, size_t height, const void* data, const bool copy)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_IMAGE);
    broadcast((int)id);
    broadcast((std::string)type);
    broadcast((int)width);
    broadcast((int)height);

    /*! send image data */
    size_t bytes = 0;
    if      (!strcasecmp(type,"RGB8"        )) bytes = width*height*3*sizeof(char);
    else if (!strcasecmp(type,"RGBA8"       )) bytes = width*height*4*sizeof(char);
    else if (!strcasecmp(type,"RGB_FLOAT32" )) bytes = width*height*3*sizeof(float);
    else if (!strcasecmp(type,"RGBA_FLOAT32")) bytes = width*height*4*sizeof(float);
    else throw std::runtime_error("unknown image type: "+std::string(type)); 
    broadcast(data,bytes);
    flush();

    if (!copy) free((void*)data);

    return (Device::RTImage)(long)id;
  }

  Device::RTImage NetworkDevice::rtNewImageFromFile(const char* fileName)
  {
    /*! load image locally */
    if (strncmp(fileName,"server:",7)) 
    {
#if defined(__MIC__)      
      throw std::runtime_error("rtNewImageFromFile not supported on MIC");
#else
      Ref<Image> image = loadImage(fileName);
      if (!image) throw std::runtime_error("cannot load image: "+std::string(fileName));
      else if (Ref<Image3c> cimg = image.dynamicCast<Image3c>())
        return rtNewImage("RGB8",cimg->width,cimg->height,cimg->ptr(),true);
      else if (Ref<Image4c> cimg = image.dynamicCast<Image4c>())
        return rtNewImage("RGBA8",cimg->width,cimg->height,cimg->ptr(),true);
      else if (Ref<Image3f> fimg = image.dynamicCast<Image3f>())
        return rtNewImage("RGB_FLOAT32",fimg->width,fimg->height,fimg->ptr(),true);
      else if (Ref<Image4f> fimg = image.dynamicCast<Image4f>())
        return rtNewImage("RGBA_FLOAT32",fimg->width,fimg->height,fimg->ptr(),true);
      else
        throw std::runtime_error("unknown image type");
#endif
      return (Device::RTImage) (long)0;
    }

    /*! load image remote */
    else
    {
      int id = allocHandle();
      broadcast((int)magick);
      broadcast((int)EMBREE_NEW_IMAGE_FROM_FILE);
      broadcast((int)id);
      broadcast((std::string)fileName);
      flush();

      return (Device::RTImage)(long)id;
    }
  }

  Device::RTTexture NetworkDevice::rtNewTexture(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_TEXTURE);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTTexture)(long)id;
  }

  Device::RTMaterial NetworkDevice::rtNewMaterial(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_MATERIAL);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTMaterial)(long)id;
  }

  Device::RTShape NetworkDevice::rtNewShape(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_SHAPE);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTShape)(long)id;
  }

  Device::RTLight NetworkDevice::rtNewLight(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_LIGHT);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTLight)(long)id;
  }

  Device::RTPrimitive NetworkDevice::rtNewShapePrimitive(Device::RTShape shape, 
                                                         Device::RTMaterial material, 
                                                         const float* transform)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_SHAPE_PRIMITIVE);
    broadcast((int)id);
    broadcast((int)(size_t)shape);
    broadcast((int)(size_t)material);

    if (transform) {
      for (size_t i=0; i<12; i++) {
        broadcast(transform[i]);
      }
    } else {
      broadcast(Vector3f(1.0f,0.0f,0.0f));
      broadcast(Vector3f(0.0f,1.0f,0.0f));
      broadcast(Vector3f(0.0f,0.0f,1.0f));
      broadcast(Vector3f(0.0f,0.0f,0.0f));
    }
    flush();
    
    return (Device::RTPrimitive)(long)id;
  }

  Device::RTPrimitive NetworkDevice::rtNewLightPrimitive(Device::RTLight light, 
                                                         Device::RTMaterial material, 
                                                         const float* transform)
  {

    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_LIGHT_PRIMITIVE);
    broadcast((int)id);
    broadcast((int)(size_t)light);
    broadcast((int)(size_t)material);

    if (transform) {
      for (size_t i=0; i<12; i++) 
        broadcast(transform[i]);
    } else {
      broadcast(Vector3f(1.0f,0.0f,0.0f));
      broadcast(Vector3f(0.0f,1.0f,0.0f));
      broadcast(Vector3f(0.0f,0.0f,1.0f));
      broadcast(Vector3f(0.0f,0.0f,0.0f));
    }
    flush();
    
    return (Device::RTPrimitive)(long)id;
  }

  Device::RTPrimitive NetworkDevice::rtTransformPrimitive(Device::RTPrimitive primitive, const float* transform) 
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_TRANSFORM_PRIMITIVE);
    broadcast((int)id);
    broadcast((int)(size_t)primitive);

    if (transform) {
      for (size_t i=0; i<12; i++) 
        broadcast(transform[i]);
    } else {
      broadcast(Vector3f(1.0f,0.0f,0.0f));
      broadcast(Vector3f(0.0f,1.0f,0.0f));
      broadcast(Vector3f(0.0f,0.0f,1.0f));
      broadcast(Vector3f(0.0f,0.0f,0.0f));
    }
    flush();
    
    return (Device::RTPrimitive)(long)id;
  }

  Device::RTScene NetworkDevice::rtNewScene(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_SCENE);
    broadcast((int)id);
    broadcast((std::string)type);
    flush();
    return (Device::RTScene)(long)id;
  }

  void NetworkDevice::rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_SCENE_PRIMITIVE);
    broadcast((int)(size_t)scene);
    broadcast((int)slot);
    broadcast((int)(size_t)prim);
    flush();
  }

  Device::RTToneMapper NetworkDevice::rtNewToneMapper(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_TONEMAPPER);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTToneMapper)(long)id;
  }

  Device::RTRenderer NetworkDevice::rtNewRenderer(const char* type)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_RENDERER);
    broadcast((int)id);
    broadcast((std::string)type);
    return (Device::RTRenderer)(long)id;
  }

  Device::RTFrameBuffer NetworkDevice::rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t numBuffers, void** ptrs)
  {
    int id = allocHandle();
    broadcast((int)magick);
    broadcast((int)EMBREE_NEW_FRAMEBUFFER);
    broadcast((int)id);
    broadcast((std::string)type);
    broadcast((int)width);
    broadcast((int)height);
    broadcast((int)numBuffers);
    flush();

    if      (!strcasecmp(type,"RGB_FLOAT32")) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGBFloat32::create);
    else if (!strcasecmp(type,"RGBA8"      )) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGBA8     ::create);
    else if (!strcasecmp(type,"RGB8"       )) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGB8      ::create);
    else throw std::runtime_error("unknown framebuffer type: "+std::string(type));

    return (Device::RTFrameBuffer)(long)id;
  }

  void* NetworkDevice::rtMapFrameBuffer(Device::RTFrameBuffer frameBuffer, int bufID) 
  {
    Ref<SwapChain>& swapChain = buffers[(size_t)frameBuffer];
    if (bufID < 0) bufID = swapChain->id();
    swapChain->buffer(bufID)->wait();
    return swapChain->buffer(bufID)->getData();
  }

  void NetworkDevice::rtUnmapFrameBuffer(Device::RTFrameBuffer frameBuffer, int bufID)
  {
  }

  void NetworkDevice::rtSwapBuffers(Device::RTFrameBuffer frameBuffer) 
  {
    /* broadcast swap buffer call to rendering servers */
    size_t fbID = (size_t)frameBuffer;
    broadcast((int)magick);
    broadcast((int)EMBREE_SWAP_BUFFERS);
    broadcast((int)fbID);
    flush();

    /* wait for next buffer to be available */
    for (size_t i=0; i<servers.size(); i++)
      receiveCommand(i);

    Ref<SwapChain>& swapchain = buffers[fbID];
    swapchain->swapBuffers();
  }
  
  void NetworkDevice::rtIncRef(Device::RTHandle handle)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_INCREF);
    broadcast((int)(size_t)handle);
    incRef((int)(size_t)handle);
  }

  void NetworkDevice::rtDecRef(Device::RTHandle handle)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_DECREF);
    broadcast((int)(size_t)handle);
    decRef((int)(size_t)handle);
  }

  /*******************************************************************
                    setting of parameters
  *******************************************************************/

  void NetworkDevice::rtSetBool1(Device::RTHandle handle, const char* property, bool x)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_BOOL1);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((bool)x);
  }
  
  void NetworkDevice::rtSetBool2(Device::RTHandle handle, const char* property, bool x, bool y)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_BOOL2);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((bool)x);
    broadcast((bool)y);
 }

  void NetworkDevice::rtSetBool3(Device::RTHandle handle, const char* property, bool x, bool y, bool z)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_BOOL3);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((bool)x);
    broadcast((bool)y);
    broadcast((bool)z);
  }

  void NetworkDevice::rtSetBool4(Device::RTHandle handle, const char* property, bool x, bool y, bool z, bool w)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_BOOL4);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((bool)x);
    broadcast((bool)y);
    broadcast((bool)z);
    broadcast((bool)w);
  }

  void NetworkDevice::rtSetInt1(Device::RTHandle handle, const char* property, int x)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_INT1);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)x);
  }

  void NetworkDevice::rtSetInt2(Device::RTHandle handle, const char* property, int x, int y)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_INT2);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)x);
    broadcast((int)y);
  }

  void NetworkDevice::rtSetInt3(Device::RTHandle handle, const char* property, int x, int y, int z)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_INT3);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)x);
    broadcast((int)y);
    broadcast((int)z);
  }

  void NetworkDevice::rtSetInt4(Device::RTHandle handle, const char* property, int x, int y, int z, int w)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_INT4);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)x);
    broadcast((int)y);
    broadcast((int)z);
    broadcast((int)w);
  }

  void NetworkDevice::rtSetFloat1(Device::RTHandle handle, const char* property, float x)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_FLOAT1);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((float)x);
  }

  void NetworkDevice::rtSetFloat2(Device::RTHandle handle, const char* property, float x, float y)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_FLOAT2);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((float)x);
    broadcast((float)y);
  }

  void NetworkDevice::rtSetFloat3(Device::RTHandle handle, const char* property, float x, float y, float z)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_FLOAT3);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((float)x);
    broadcast((float)y);
    broadcast((float)z);
  }

  void NetworkDevice::rtSetFloat4(Device::RTHandle handle, const char* property, float x, float y, float z, float w)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_FLOAT4);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((float)x);
    broadcast((float)y);
    broadcast((float)z);
    broadcast((float)w);
  }

  void NetworkDevice::rtSetArray(Device::RTHandle handle, const char* property, const char* type, Device::RTData data, size_t size, size_t stride, size_t ofs)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_ARRAY);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((std::string)type);
    broadcast((int)(size_t)data);
    broadcast((int)size);
    broadcast((int)stride);
    broadcast((int)ofs);
  }

  void NetworkDevice::rtSetString(Device::RTHandle handle, const char* property, const char* str)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_STRING);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((std::string)str);
  }

  void NetworkDevice::rtSetImage(Device::RTHandle handle, const char* property, Device::RTImage img)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_IMAGE);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)(size_t)img);
  }

  void NetworkDevice::rtSetTexture(Device::RTHandle handle, const char* property, Device::RTTexture tex)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_TEXTURE);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((int)(size_t)tex);
  }

  void NetworkDevice::rtSetTransform(Device::RTHandle handle, const char* property, const float* transform)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_SET_TRANSFORM);
    broadcast((int)(size_t)handle);
    broadcast((std::string)property);
    broadcast((float)transform[ 0]);
    broadcast((float)transform[ 1]);
    broadcast((float)transform[ 2]);
    broadcast((float)transform[ 3]);
    broadcast((float)transform[ 4]);
    broadcast((float)transform[ 5]);
    broadcast((float)transform[ 6]);
    broadcast((float)transform[ 7]);
    broadcast((float)transform[ 8]);
    broadcast((float)transform[ 9]);
    broadcast((float)transform[10]);
    broadcast((float)transform[11]);
  }

  void NetworkDevice::rtClear(Device::RTHandle handle)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_CLEAR);
    broadcast((int)(size_t)handle);
  }

  void NetworkDevice::rtCommit(Device::RTHandle handle)
  {
    broadcast((int)magick);
    broadcast((int)EMBREE_COMMIT);
    broadcast((int)(size_t)handle);
    flush();
  }

  /*******************************************************************
                          render calls
  *******************************************************************/

  void NetworkDevice::rtRenderFrame(Device::RTRenderer renderer, Device::RTCamera camera, Device::RTScene scene, 
                                    Device::RTToneMapper toneMapper, Device::RTFrameBuffer frameBuffer, int accumulate)
  {
    /*! send render command */
    broadcast((int)magick);
    broadcast((int)EMBREE_RENDER_FRAME);
    broadcast((int)(size_t)renderer);
    broadcast((int)(size_t)camera);
    broadcast((int)(size_t)scene);
    broadcast((int)(size_t)toneMapper);
    broadcast((int)(size_t)frameBuffer);
    broadcast((int)accumulate);
    flush();
  }

  bool NetworkDevice::rtPick(Device::RTCamera camera, float x, float y, Device::RTScene scene, float& px, float& py, float& pz) 
  {
    /*! send pick command */
    network::write(servers[0],(int)magick);
    network::write(servers[0],(int)EMBREE_PICK);
    network::write(servers[0],(int)(size_t)camera);
    network::write(servers[0],(float)x);
    network::write(servers[0],(float)y);
    network::write(servers[0],(int)(size_t)scene);
    network::flush(servers[0]);
    
    /*! wait for pick result */
    for (size_t i=0; i<servers.size(); i++)
      receiveCommand(i);

    /*! return pick data */
    px = pick.pos.x;
    py = pick.pos.y;
    pz = pick.pos.z;
    return pick.hit;
  }
}
