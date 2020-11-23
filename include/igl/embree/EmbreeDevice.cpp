#include "EmbreeDevice.h"
// Implementation

namespace igl {
    namespace embree {
     RTCDevice EmbreeDevice::g_device = nullptr;
     int EmbreeDevice::g_cntr = 0;
    }
}
