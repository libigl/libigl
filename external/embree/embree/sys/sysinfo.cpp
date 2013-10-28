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

#include "sysinfo.h"
#include "intrinsics.h"

////////////////////////////////////////////////////////////////////////////////
/// All Platforms
////////////////////////////////////////////////////////////////////////////////

namespace embree
{
  std::string getPlatformName() 
  {
#if defined(__LINUX__) && !defined(__X86_64__)
    return "Linux (32bit)";
#elif defined(__LINUX__) && defined(__X86_64__)
    return "Linux (64bit)";
#elif defined(__FREEBSD__) && !defined(__X86_64__)
    return "FreeBSD (32bit)";
#elif defined(__FREEBSD__) && defined(__X86_64__)
    return "FreeBSD (64bit)";
#elif defined(__CYGWIN__) && !defined(__X86_64__)
    return "Cygwin (32bit)";
#elif defined(__CYGWIN__) && defined(__X86_64__)
    return "Cygwin (64bit)";
#elif defined(__WIN32__) && !defined(__X86_64__)
    return "Windows (32bit)";
#elif defined(__WIN32__) && defined(__X86_64__)
    return "Windows (64bit)";
#elif defined(__MACOSX__) && !defined(__X86_64__)
    return "MacOS (32bit)";
#elif defined(__MACOSX__) && defined(__X86_64__)
    return "MacOS (64bit)";
#elif defined(__UNIX__) && !defined(__X86_64__)
    return "Unix (32bit)";
#elif defined(__UNIX__) && defined(__X86_64__)
    return "Unix (64bit)";
#else
    return "Unknown";
#endif
  }

  std::string getCPUVendor()
  {
    int cpuinfo[4]; 
    __cpuid (cpuinfo, 0); 
    int name[4];
    name[0] = cpuinfo[1];
    name[1] = cpuinfo[3];
    name[2] = cpuinfo[2];
    name[3] = 0;
    return (char*)name;
  }

  CPUModel getCPUModel() 
  {
    int out[4];
    __cpuid(out, 0);
    if (out[0] < 1) return CPU_UNKNOWN;
    __cpuid(out, 1);
    int family = ((out[0] >> 8) & 0x0F) + ((out[0] >> 20) & 0xFF);
    int model  = ((out[0] >> 4) & 0x0F) | ((out[0] >> 12) & 0xF0);
    if (family !=   6) return CPU_UNKNOWN;           // earlier than P6
    if (model == 0x0E) return CPU_CORE1;             // Core 1
    if (model == 0x0F) return CPU_CORE2;             // Core 2, 65 nm
    if (model == 0x16) return CPU_CORE2;             // Core 2, 65 nm Celeron
    if (model == 0x17) return CPU_CORE2;             // Core 2, 45 nm
    if (model == 0x1A) return CPU_CORE_NEHALEM;      // Core i7, Nehalem
    if (model == 0x1E) return CPU_CORE_NEHALEM;      // Core i7
    if (model == 0x1F) return CPU_CORE_NEHALEM;      // Core i7
    if (model == 0x2C) return CPU_CORE_NEHALEM;      // Core i7, Xeon
    if (model == 0x2E) return CPU_CORE_NEHALEM;      // Core i7, Xeon
    if (model == 0x2A) return CPU_CORE_SANDYBRIDGE;  // Core i7, SandyBridge
    return CPU_UNKNOWN;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Windows Platform
////////////////////////////////////////////////////////////////////////////////

#ifdef __WIN32__

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace embree
{
  std::string getExecutableFileName() {
    char filename[1024];
    if (!GetModuleFileName(NULL, filename, sizeof(filename))) return std::string();
    return std::string(filename);
  }

  size_t getNumberOfLogicalThreads() 
  {
#if (_WIN32_WINNT >= 0x0601)
    int groups = GetActiveProcessorGroupCount();
    int totalProcessors = 0;
    for (int i = 0; i < groups; i++) 
      totalProcessors += GetActiveProcessorCount(i);
    return totalProcessors;
#else
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#endif
  }

  int getTerminalWidth() 
  {
    HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
    if (handle == INVALID_HANDLE_VALUE) return 80;
    CONSOLE_SCREEN_BUFFER_INFO info = { 0 };
    GetConsoleScreenBufferInfo(handle, &info);
    return info.dwSize.X;
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Linux Platform
////////////////////////////////////////////////////////////////////////////////

#ifdef __LINUX__

#include <stdio.h>
#include <unistd.h>

namespace embree
{
  std::string getExecutableFileName() 
  {
    char pid[32]; sprintf(pid, "/proc/%d/exe", getpid());
    char buf[1024];
    int bytes = readlink(pid, buf, sizeof(buf)-1);
    if (bytes != -1) buf[bytes] = '\0';
    return std::string(buf);
  }
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// MacOS Platform
////////////////////////////////////////////////////////////////////////////////

#ifdef __MACOSX__

#include <mach-o/dyld.h>

namespace embree
{
  std::string getExecutableFileName()
  {
    char buf[1024];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) != 0) return std::string();
    return std::string(buf);
  }
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// Unix Platform
////////////////////////////////////////////////////////////////////////////////

#if defined(__UNIX__)

#include <unistd.h>
#include <sys/ioctl.h>

#if defined(__USE_NUMA__)
#include <numa.h>
#endif

namespace embree
{
  size_t getNumberOfLogicalThreads() {
    static int nThreads = -1;
    if (nThreads == -1)
      nThreads = sysconf(_SC_NPROCESSORS_CONF);
    return nThreads;
  }

  int getTerminalWidth() 
  {
    struct winsize info;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &info) < 0) return 80;
    return info.ws_col;
  }
}
#endif

