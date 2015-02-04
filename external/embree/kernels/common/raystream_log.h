// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "../kernels/common/default.h"
#include "embree2/rtcore_ray.h"
#include <iostream>
#include <fstream>
//#include <pthread.h>

namespace embree
{

#if defined(__MIC__)
 #define DEFAULT_PATH_BINARY_FILES      "/home/micuser/"
#else
 #define DEFAULT_PATH_BINARY_FILES      "./"
#endif

#define DEFAULT_FILENAME_GEOMETRY      "geometry.bin"

#define DEFAULT_FILENAME_RAY16         "ray16.bin"
#define DEFAULT_FILENAME_RAY16_VERIFY  "ray16_verify.bin"

#define DEFAULT_FILENAME_RAY8          "ray8.bin"
#define DEFAULT_FILENAME_RAY8_VERIFY   "ray8_verify.bin"

#define DEFAULT_FILENAME_RAY4          "ray4.bin"
#define DEFAULT_FILENAME_RAY4_VERIFY   "ray4_verify.bin"

#define DEFAULT_FILENAME_RAY1          "ray1.bin"
#define DEFAULT_FILENAME_RAY1_VERIFY   "ray1_verify.bin"

  class RayStreamLogger
  {
  public:


    class DataStream {
    private:
      bool initialized;

      std::string filename;
      std::ofstream data;

      void open() 
      {
        data.open(filename.c_str(),std::ios::out | std::ios::binary);
        data.seekp(0, std::ios::beg);
        
        if (!data)
          {
            DBG_PRINT(filename);
            FATAL("could not open data stream");
          }
      }
      
    public:

    DataStream(std::string name) : initialized(false) 
        {
          filename = name;
        }

      ~DataStream() {
        if (initialized)
          {
            data.close();
          }
      }
    
      void write(void *ptr, const size_t size)
      {
        if (unlikely(!initialized))
          {
            open();
            initialized = true;
          }

        data.write((char*)ptr,size);
        data << std::flush;
      }

    };

  private:

    MutexSys mutex;

    DataStream *ray16;
    DataStream *ray16_verify;
    DataStream *ray8;
    DataStream *ray8_verify;
    DataStream *ray4;
    DataStream *ray4_verify;
    DataStream *ray1;
    DataStream *ray1_verify;


  public:

    enum { 
      RAY_INTERSECT = 0,
      RAY_OCCLUDED  = 1
    };

    RayStreamLogger();
    ~RayStreamLogger();

    struct __aligned(64) LogRay1  {
      unsigned int type;
      unsigned int dummy[3]; /* 16 bytes alignment */
      RTCRay ray;

      LogRay1() {
	memset(this,0,sizeof(LogRay1));
      }
    };

    struct __aligned(64) LogRay4  {
      unsigned int type;
      unsigned int m_valid;
      unsigned int numRays;
      unsigned int dummy[1]; /* 16 bytes alignment */
      RTCRay4 ray4;

      LogRay4() {
	memset(this,0,sizeof(LogRay4));
      }
    };

    struct __aligned(64) LogRay8  {
      unsigned int type;
      unsigned int m_valid;
      unsigned int numRays;
      unsigned int dummy[5]; /* 32 bytes alignment */
      RTCRay8 ray8;

      LogRay8() {
	memset(this,0,sizeof(LogRay8));
      }
    };

    struct __aligned(64) LogRay16  {
      unsigned int type;
      unsigned int m_valid;
      unsigned int numRays;
      unsigned int dummy[13]; /* 64 bytes alignment */
      RTCRay16 ray16;

      LogRay16() {
	memset(this,0,sizeof(LogRay16));
      }
    };

      
  static RayStreamLogger rayStreamLogger;

  void logRay16Intersect(const void* valid, void* scene, RTCRay16& start, RTCRay16& end);
  void logRay16Occluded (const void* valid, void* scene, RTCRay16& start, RTCRay16& end);

  void logRay8Intersect(const void* valid, void* scene, RTCRay8& start, RTCRay8& end);
  void logRay8Occluded (const void* valid, void* scene, RTCRay8& start, RTCRay8& end);

  void logRay4Intersect(const void* valid, void* scene, RTCRay4& start, RTCRay4& end);
  void logRay4Occluded (const void* valid, void* scene, RTCRay4& start, RTCRay4& end);

  void logRay1Intersect(void* scene, RTCRay& start, RTCRay& end);
  void logRay1Occluded (void* scene, RTCRay& start, RTCRay& end);

  void dumpGeometry(void* scene);
  };
};
