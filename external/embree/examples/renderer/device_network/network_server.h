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

#ifndef __EMBREE_SERVER_H__
#define __EMBREE_SERVER_H__

#include "default.h"
#include <queue>
#include <vector>
#include "sys/network.h"
#include "sys/sync/condition.h"
#include "sys/sync/mutex.h"
#include "sys/thread.h"
#include "device/device.h"

namespace embree 
{
  class NetworkServer 
  {
  public:
    
    /*! constructor */
    NetworkServer(network::socket_t socket, Device* device, int encoding, bool verbose);
    
    /*! destructor */
    ~NetworkServer();
    
  private:
    
    /*! SwapChain ============================================ */
    class SwapChain : public RefCount
    {
    public:
      
      friend class NetworkServer;
      
      /*! constructor */
      SwapChain (NetworkServer* server, std::string type, Device::RTFrameBuffer frameBuffer, 
                 int frameBufferID, size_t width, size_t height, size_t numBuffers, int serverID, int serverCount) 
        : server(server), type(type), frameBuffer(frameBuffer), frameBufferID(frameBufferID), width(width), height(height), 
        numBuffers(numBuffers), writeID(0), readID((size_t)-1), serverID(serverID), serverCount(serverCount)
      {
        server->device->rtIncRef(frameBuffer);
        encoded = (unsigned char *) malloc(width * height * 4);  
      }
      
      ~SwapChain () 
      {
        server->device->rtDecRef(frameBuffer);
        free(encoded); encoded = NULL;
      }
      
      /*! switch to next write buffer */
      void nextWriteBuffer();
      
      /*! send data to network client */
      void send();
      
      /*! switch to next read buffer */
      void nextReadBuffer();
      
    private:
      NetworkServer* server;
      std::string type;
      Device::RTFrameBuffer frameBuffer;
      int frameBufferID;
      size_t width,height;
      size_t numBuffers;
      unsigned char* encoded;  //!< helper to encode framebuffer
      
    private:
      size_t writeID;  //!< next buffer to render into
      size_t readID;   //!< next buffer to read from

    private:
      int serverID;
      int serverCount;
    };
    
    /*! render device and handle management =========================== */
    Device *device;
    std::vector<Device::RTHandle> handles;
    std::vector<Ref<SwapChain> > swapChains;
    std::vector<int> numRefs;
      
    /*! update ID to handle mapping */
    void set(int id, Device::RTHandle handle, SwapChain* swapChain = NULL);
    
    /*! return the handle associated with this ID */
    template<typename T> T get(size_t id);
    
    /*! client communication state ==================================== */
    network::socket_t socket;
    int serverID;
    int serverCount;
    
    /*! receive data from the client */
    void receive();
    
    /*! thread level function to receive data from the client */
    static void _receive(void *data);
    
    int encoding;
    bool verbose;
  };
}

#endif // __EMBREE_SERVER_H__

