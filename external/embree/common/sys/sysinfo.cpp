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

#include "sysinfo.h"
#include "intrinsics.h"
#include "stl/string.h"

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
    return "Mac OS X (32bit)";
#elif defined(__MACOSX__) && defined(__X86_64__)
    return "Mac OS X (64bit)";
#elif defined(__UNIX__) && !defined(__X86_64__)
    return "Unix (32bit)";
#elif defined(__UNIX__) && defined(__X86_64__)
    return "Unix (64bit)";
#else
    return "Unknown";
#endif
  }

  std::string getCompilerName()
  {
#if defined(__INTEL_COMPILER)
    int icc_mayor = __INTEL_COMPILER / 100 % 100;
    int icc_minor = __INTEL_COMPILER % 100;
    std::string version = "Intel Compiler ";
    version += std::stringOf(icc_mayor);
    version += "." + std::stringOf(icc_minor);
#if defined(__INTEL_COMPILER_UPDATE)
    version += "." + std::stringOf(__INTEL_COMPILER_UPDATE);
#endif
    return version;
#elif defined(__clang__)
    return "CLANG " __clang_version__;
#elif defined (__GNUC__)
    return "GCC " __VERSION__;
#elif defined(_MSC_VER)
    std::string version = std::stringOf(_MSC_FULL_VER);
    version.insert(4,".");
    version.insert(9,".");
    version.insert(2,".");
    return "Visual C++ Compiler " + version;
#else
    return "Unknown Compiler";
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

  /* cpuid[eax=0].ecx */
  static const int CPU_FEATURE_BIT_SSE3   = 1 << 0;
  static const int CPU_FEATURE_BIT_SSSE3  = 1 << 9;
  static const int CPU_FEATURE_BIT_FMA3   = 1 << 12;
  static const int CPU_FEATURE_BIT_SSE4_1 = 1 << 19;
  static const int CPU_FEATURE_BIT_SSE4_2 = 1 << 20;
  static const int CPU_FEATURE_BIT_MOVBE  = 1 << 22;
  static const int CPU_FEATURE_BIT_POPCNT = 1 << 23;
  static const int CPU_FEATURE_BIT_XSAVE  = 1 << 26;
  static const int CPU_FEATURE_BIT_OXSAVE = 1 << 27;
  static const int CPU_FEATURE_BIT_AVX    = 1 << 28;
  static const int CPU_FEATURE_BIT_F16C   = 1 << 29;
  static const int CPU_FEATURE_BIT_RDRAND = 1 << 30;

  /* cpuid[eax=1].edx */
  static const int CPU_FEATURE_BIT_SSE  = 1 << 25;
  static const int CPU_FEATURE_BIT_SSE2 = 1 << 26;

  /* cpuid[eax=0x80000001].ecx */
  static const int CPU_FEATURE_BIT_LZCNT = 1 << 5;

  /* cpuid[eax=7,ecx=0].ebx */
  static const int CPU_FEATURE_BIT_BMI1  = 1 << 3;
  static const int CPU_FEATURE_BIT_AVX2  = 1 << 5;
  static const int CPU_FEATURE_BIT_BMI2  = 1 << 8;

  bool check_xcr0_ymm() 
  {
    int xcr0;
#if defined (__WIN32__)

#if defined(__INTEL_COMPILER) 
    xcr0 = (int)_xgetbv(0);
#elif (defined(_MSC_VER) && (_MSC_FULL_VER >= 160040219)) // min VS2010 SP1 compiler is required
    xcr0 = (int)_xgetbv(0); 
#else
#pragma message ("WARNING: AVX not supported by your compiler.")
    xcr0 = 0;
#endif

#else

#if defined(__INTEL_COMPILER) 
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#elif ((__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)) && (!defined(__MACOSX__) || defined(__TARGET_AVX__) || defined(__TARGET_AVX2__))
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#elif ((__clang_major__ > 3) || (__clang_major__ == 3 && __clang_minor__ >= 1))
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#else
#pragma message ("WARNING: AVX not supported by your compiler.")
    xcr0 = 0;
#endif

#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
  }

  int cpu_features = 0;

  int getCPUFeatures()
  {
    if (cpu_features) return cpu_features;
    
    int info[4]; 
    __cpuid(info, 0x00000000);
    unsigned nIds = info[0];  
    __cpuid(info, 0x80000000);
    unsigned nExIds = info[0];

    int info1[4] = { 0,0,0,0 };
    int info7[4] = { 0,0,0,0 };
    int infoe1[4] = { 0,0,0,0 };
    if (nIds >= 1) __cpuid (info1,0x00000001);
#if _WIN32
#if _MSC_VER && (_MSC_FULL_VER < 160040219)
#else
	 if (nIds >= 7) __cpuidex(info7,0x00000007,0);
#endif
#else
    if (nIds >= 7) __cpuid_count(info7,0x00000007,0);
#endif
    if (nExIds >= 0x80000001) __cpuid(infoe1,0x80000001);

    bool has_ymm = false;
    if (info1[2] & CPU_FEATURE_BIT_OXSAVE)
      has_ymm = check_xcr0_ymm();
    
    if (info1[3] & CPU_FEATURE_BIT_SSE) cpu_features |= CPU_FEATURE_SSE;
    if (info1[3] & CPU_FEATURE_BIT_SSE2) cpu_features |= CPU_FEATURE_SSE2;
    if (info1[2] & CPU_FEATURE_BIT_SSE3) cpu_features |= CPU_FEATURE_SSE3;
    if (info1[2] & CPU_FEATURE_BIT_SSSE3) cpu_features |= CPU_FEATURE_SSSE3;
    if (info1[2] & CPU_FEATURE_BIT_SSE4_1) cpu_features |= CPU_FEATURE_SSE41;
    if (info1[2] & CPU_FEATURE_BIT_SSE4_2) cpu_features |= CPU_FEATURE_SSE42;
    if (info1[2] & CPU_FEATURE_BIT_POPCNT) cpu_features |= CPU_FEATURE_POPCNT;
    if (has_ymm && info1[2] & CPU_FEATURE_BIT_AVX) cpu_features |= CPU_FEATURE_AVX;
    if (info1[2] & CPU_FEATURE_BIT_F16C) cpu_features |= CPU_FEATURE_F16C;
    if (info1[2] & CPU_FEATURE_BIT_RDRAND) cpu_features |= CPU_FEATURE_RDRAND;
    if (has_ymm && info7[1] & CPU_FEATURE_BIT_AVX2) cpu_features |= CPU_FEATURE_AVX2;
    if (has_ymm && info1[2] & CPU_FEATURE_BIT_FMA3) cpu_features |= CPU_FEATURE_FMA3;
    if (infoe1[2] & CPU_FEATURE_BIT_LZCNT) cpu_features |= CPU_FEATURE_LZCNT;
    if (info7[1] & CPU_FEATURE_BIT_BMI1) cpu_features |= CPU_FEATURE_BMI1;
    if (info7[1] & CPU_FEATURE_BIT_BMI2) cpu_features |= CPU_FEATURE_BMI2;

#if defined(__MIC__)
    cpu_features |= CPU_FEATURE_KNC;
#endif
    return cpu_features;
  }
  
  std::string stringOfCPUFeatures(int features)
  {
    std::string str;
    if (features & CPU_FEATURE_SSE   ) str += "SSE ";
    if (features & CPU_FEATURE_SSE2  ) str += "SSE2 ";
    if (features & CPU_FEATURE_SSE3  ) str += "SSE3 ";
    if (features & CPU_FEATURE_SSSE3 ) str += "SSSE3 ";
    if (features & CPU_FEATURE_SSE41 ) str += "SSE41 ";
    if (features & CPU_FEATURE_SSE42 ) str += "SSE42 ";
    if (features & CPU_FEATURE_POPCNT) str += "POPCNT ";
    if (features & CPU_FEATURE_AVX   ) str += "AVX ";
    if (features & CPU_FEATURE_F16C  ) str += "F16C ";
    if (features & CPU_FEATURE_RDRAND) str += "RDRAND ";
    if (features & CPU_FEATURE_AVX2  ) str += "AVX2 ";
    if (features & CPU_FEATURE_FMA3  ) str += "FMA3 ";
    if (features & CPU_FEATURE_LZCNT ) str += "LZCNT ";
    if (features & CPU_FEATURE_BMI1  ) str += "BMI1 ";
    if (features & CPU_FEATURE_BMI2  ) str += "BMI2 ";
    if (features & CPU_FEATURE_KNC   ) str += "KNC ";
    return str;
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

  size_t getNumberOfCores() {
    static int nCores = -1;
    if (nCores == -1) nCores = getNumberOfLogicalThreads(); // FIXME: detect if hyperthreading is enabled
    if (nCores ==  0) nCores = 1;
    return nCores;
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
/// Mac OS X Platform
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
    if (nThreads == -1) nThreads = sysconf(_SC_NPROCESSORS_CONF);
    return nThreads;
  }

  size_t getNumberOfCores() {
    static int nCores = -1;
    if (nCores == -1) nCores = sysconf(_SC_NPROCESSORS_CONF)/2; // FIXME: detect if hyperthreading is enabled
    if (nCores ==  0) nCores = 1;
    return nCores;
  }

  int getTerminalWidth() 
  {
    struct winsize info;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &info) < 0) return 80;
    return info.ws_col;
  }
}
#endif

