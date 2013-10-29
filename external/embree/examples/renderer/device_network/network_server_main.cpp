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

#include "network_common.h"
#include "network_server.h"
#include "sys/sysinfo.h"

namespace embree 
{
  /*! encoding mode for outgoing frame buffer data */
  static int g_encoding = EMBREE_FRAME_DATA_NATIVE;
  
  /*! server port number */
  static unsigned short g_port = 8484;
  
  /*! number of threads to use during rendering */
  static int g_threads = 0;
  
  /*! verbose output mode */
  static bool g_verbose = false;

  /*! device type */
  std::string g_device_type = "default";

  /*! allow multiple subsequent connections */
  bool g_multiple_connections = true;
  
  /*! parse command line */
  void parseCommandLine(int argc, char **argv) 
  {
    /*! loop over tokens */
    for (int i = 1 ; i < argc ; i++) 
    {
      /* set device */
      if (!strcmp(argv[i],"-device")) {
        if (++i >= argc) throw std::runtime_error("no device specified");
        g_device_type = argv[i];
      }

      /*! frame buffer encoding mode */
      else if (!strcmp(argv[i], "-encode")) 
      {
        if (++i >= argc) throw std::runtime_error("no encoding specified");
        else if (!strcmp(argv[i], "native")) g_encoding = EMBREE_FRAME_DATA_NATIVE;
        else if (!strcmp(argv[i], "jpeg"  )) g_encoding = EMBREE_FRAME_DATA_JPEG;
        else if (!strcmp(argv[i], "rgb8"  )) g_encoding = EMBREE_FRAME_DATA_RGB8;
        else if (!strcmp(argv[i], "rgbe8" )) g_encoding = EMBREE_FRAME_DATA_RGBE8;
        else throw std::runtime_error("unknown encoding mode: " + (std::string) argv[i]);
      }
      
      /*! port number */
      else if (!strcmp(argv[i], "-port")) 
      {
        if (++i >= argc) throw std::runtime_error("no port number specified");
        g_port = (unsigned short) atoi(argv[i]);
        if (g_port == 0) throw std::runtime_error("invalid port specified: " + (std::string) argv[i]);
      }
      
      /*! number of threads to use during rendering */
      else if (!strcmp(argv[i], "-threads")) 
      {
        if (++i >= argc) throw std::runtime_error("no thread count specified");
        g_threads = atoi(argv[i]);
      }

       /*! terminate after first connection */
      else if (!strcmp(argv[i], "-single-connection"))
        g_multiple_connections = false;
      
      /*! verbose mode */
      else if (!strcmp(argv[i], "-verbose")) {
        g_verbose = true;
      }

      /*! invalid argument */
      else {
        std::cout << "usage: embree_network_server [-encode rgb8|rgbe8|jpeg] [-port number] [-threads number] [-verbose]" << std::endl;
        throw std::runtime_error("invalid command line argument: " + (std::string) argv[i]);
      }
    }
  }
  
  void main(int argc, char **argv) 
  {
    /*! parse command line options */
    parseCommandLine(argc, argv);
    
    /*! bind to a port */
    network::socket_t socket = network::bind(g_port);
    
    /*! handle incoming connections */
    do {
      std::cout << std::endl << "listening for connections on port " << g_port << " ... " << std::flush;
      network::socket_t client = network::listen(socket);
      std::cout << "  [CONNECTED]" << std::endl;
      new NetworkServer(client,Device::rtCreateDevice(g_device_type.c_str(), g_threads),g_encoding,g_verbose);
    } 
    while (g_multiple_connections);
  }
}

/******************************************************************************/
/*                               Main Function                                */
/******************************************************************************/

int main(int argc, char **argv) 
{
  try {
    embree::main(argc, argv);
    return(0);
  } catch (const std::exception &e) {
    std::cout << "main(): Error: " << e.what() << std::endl;
    return(1);
  } catch (...) {
    std::cout << "main(): Error: unknown exception caught." << std::endl;
    return(1);
  }
}

