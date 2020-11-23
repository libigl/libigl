// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Vladimir Fonov <vladimir.fonov@gmail.com>
//               2013 Alec Jacobson <alecjacobson@gmail.com>
//               2014 Christian Schüller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef IGL_EMBREE_EMBREE_DEVICE_H
#define IGL_EMBREE_EMBREE_DEVICE_H
#include <embree3/rtcore.h>
#include <iostream>

namespace igl
{
  namespace embree
  {
     // keep track of embree device
     class EmbreeDevice
     {
       private:
        static RTCDevice g_device;
        static int       g_cntr;
       public:
        static RTCDevice get_device(void)
        {
            if(!g_device) 
            {
                g_device = rtcNewDevice (NULL);
                if(rtcGetDeviceError (g_device) != RTC_ERROR_NONE)
                    std::cerr << "Embree: An error occurred while initializing embree core!" << std::endl;
            #ifdef IGL_VERBOSE
                else
                    std::cerr << "Embree: core initialized." << std::endl;
            #endif
            }
            ++g_cntr;
            return g_device;
        }
        static void release_device(void)
        {
            if(!--g_cntr) {
                rtcReleaseDevice (g_device);
                g_device = nullptr;                
            #ifdef IGL_VERBOSE
                    std::cerr << "Embree: core released." << std::endl;
            #endif
            }
        }
     };
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "EmbreeDevice.cpp"
#endif

#endif // IGL_EMBREE_EMBREE_DEVICE_H