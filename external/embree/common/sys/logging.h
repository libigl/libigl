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

#include "platform.h"

#include <sstream>
#include <cstdio>
#include <string>
#include <ctime>
#include <stdio.h>

//#define EMBREE_LOG_SHOW_FILE               // file.cpp(line#):
//#define EMBREE_LOG_SHOW_FILE_FULL_PATH     // full path for the file name
//#define EMBREE_LOG_SHOW_TIME               // [hh:mm:ss]
//#define EMBREE_LOG_SHOW_TID                // [t####]
//#define EMBREE_LOG_SHOW_PID                // [P####]
//#define EMBREE_LOG_SHOW_MSGTYPE            // DEBUG, INFO, etc.

#define EMB_LOG_MSG_CRITICAL     0
#define EMB_LOG_MSG_ERROR        1
#define EMB_LOG_MSG_WARNING      2
#define EMB_LOG_MSG_INFO         3
#define EMB_LOG_MSG_DEBUG        4

#define EMB_LOG_MSGTYPE_TO_STRING(type) (       \
  (type == EMB_LOG_MSG_INFO)     ? "INFO"     : \
  (type == EMB_LOG_MSG_WARNING)  ? "WARNING"  : \
  (type == EMB_LOG_MSG_ERROR)    ? "ERROR"    : \
  (type == EMB_LOG_MSG_CRITICAL) ? "CRITICAL" : \
  (type == EMB_LOG_MSG_DEBUG)    ? "DEBUG"    : \
  "UNKNOWN")

#ifndef EMBREE_LOG_LEVEL
#if defined(_DEBUG) || defined(DEBUG)
  #define EMBREE_LOG_LEVEL EMB_LOG_MSG_DEBUG
#else
  #define EMBREE_LOG_LEVEL EMB_LOG_MSG_INFO
#endif
#endif

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
  #pragma warning(push)
  // for warning C4127: conditional expression is constant
  #pragma warning(disable:4127)
  // for warning C4996: 'sprintf': This function or variable may be unsafe.
  #pragma warning(disable:4996)
  #define NOMINMAX
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  #undef near
  #undef far
  #include <process.h>
  #define EMB_LOG_TID (GetCurrentThreadId())
  #define EMB_LOG_PID _getpid()
  #define EMBREE_LOG_SYSTEM_PATH_SEPARATOR "\\"
#else
  #include <sys/types.h>
  #include <sys/time.h>
  #include <unistd.h>
  #include <sys/syscall.h>
  #include <bits/syscall.h>
  #define EMB_LOG_TID syscall(__NR_gettid)
  #define EMB_LOG_PID getpid()
  #define EMBREE_LOG_SYSTEM_PATH_SEPARATOR "/"
#endif

namespace embree {

template <typename T>
class Log
{
public:
  Log() {};
  virtual ~Log();
  std::ostringstream& put(int32 level,
                          const char* file,
                          int32 lineNumber);
protected:
  std::ostringstream os;
private:
  Log(const Log&);
  Log& operator =(const Log&);

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
  std::string nowTime()
  {
    const int32 MAX_LEN = 200;
    char buffer[MAX_LEN];
    if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0,
                       "HH:mm:ss", buffer, MAX_LEN) == 0) {
      return "Error in nowTime()";
    }

    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount()) % 1000);
    return result;
  }
#else
  std::string nowTime()
  {
    char buffer[16];
    timeval tv;
    gettimeofday(&tv, NULL);
    std::strftime(buffer, sizeof(buffer), "%H:%M:%S", std::localtime(&tv.tv_sec));
    char result[16] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
  }
#endif //WIN32
};

template <typename T>
Log<T>::~Log()
{
  os << std::endl;
  T::output(os.str());
}

template <typename T>
std::ostringstream& Log<T>::put(int32 type, 
                                const char* file, 
                                int32 lineNumber)
{
#ifdef EMBREE_LOG_SHOW_FILE
  // Visual Studio debug output window will understand "filename(line#):" format
  // to allow double click and jump to line.
  #ifndef EMBREE_LOG_SHOW_FILE_FULL_PATH
    std::string filePath(file);
    size_t pos = filePath.find_last_of(EMBREE_LOG_SYSTEM_PATH_SEPARATOR);
    if (pos != std::string::npos) {
      filePath = filePath.substr(pos+1);
    }
    os << filePath;
  #else
    os << file;
  # endif
  os << "(" << lineNumber << "): ";
#endif
#ifdef EMBREE_LOG_SHOW_TIME
  os << "[" << nowTime() << "] ";
#endif 
#ifdef EMBREE_LOG_SHOW_PID || EMBREE_LOG_SHOW_TID
  os << "[";
    #ifdef EMBREE_LOG_SHOW_PID
      os << "P" << EMB_LOG_PID;
    #endif
    #ifdef EMBREE_LOG_SHOW_TID
      os << "|t" << EMB_LOG_TID;
    #endif
  os << "] ";
#endif
#ifdef EMBREE_LOG_SHOW_MSGTYPE
  os << EMB_LOG_MSGTYPE_TO_STRING(type) << ": ";
#endif

  return os;
}

class LogOutput
{
public:
  /**
   * To log output to a file, do something like this before any logging code:
   * FILE* logFile = fopen("log.txt", "a");
   * embree::LogOutput::stream() = logFile;
   */
  static FILE*& stream();
  static void output(const std::string& msg);
};

inline FILE*& LogOutput::stream()
{
  static FILE* pStream = stderr;
  return pStream;
}

inline void LogOutput::output(const std::string& msg)
{   
  FILE* pStream = stream();
  if (!pStream) {
    return;
  }
  // output in the Visual Studio debugger.
#if (defined(_DEBUG) || defined(DEBUG)) && defined(_MSC_VER)
  OutputDebugStringA(msg.c_str());
#endif
  fprintf(pStream, "%s", msg.c_str());
  fflush(pStream);
}

typedef Log<LogOutput> DefaultLogger;

#if (EMB_LOG_MSG_DEBUG <= EMBREE_LOG_LEVEL)
  #define EMB_LOG_DEBUG(message) embree::DefaultLogger().put(EMB_LOG_MSG_DEBUG, __FILE__, __LINE__) << message
#else
  #define EMB_LOG_DEBUG(message) {}
#endif

#if (EMB_LOG_MSG_INFO <= EMBREE_LOG_LEVEL)
  #define EMB_LOG_INFO(message) embree::DefaultLogger().put(EMB_LOG_MSG_INFO, __FILE__, __LINE__) << message
#else
  #define EMB_LOG_INFO(message) {}
#endif

#if (EMB_LOG_MSG_WARNING <= EMBREE_LOG_LEVEL)
  #define EMB_LOG_WARNING(message) embree::DefaultLogger().put(EMB_LOG_MSG_WARNING, __FILE__, __LINE__) << message
#else
  #define EMB_LOG_WARNING(message) {}
#endif

#if (EMB_LOG_MSG_ERROR <= EMBREE_LOG_LEVEL)
  #define EMB_LOG_ERROR(message) embree::DefaultLogger().put(EMB_LOG_MSG_ERROR, __FILE__, __LINE__) << message
#else
  #define EMB_LOG_ERROR(message) {}
#endif

#if (EMB_LOG_MSG_CRITICAL <= EMBREE_LOG_LEVEL)
  #define EMB_LOG_CRITICAL(message) embree::DefaultLogger().put(EMB_LOG_MSG_CRITICAL, __FILE__, __LINE__) << message
#else
  #define EMB_LOG_CRITICAL(message) {}
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
  #pragma warning(pop)
#endif

// fully qualified function name
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
  #define EMB_LOG_FUNC_NAME_LONG __PRETTY_FUNCTION__ << " "
#elif defined (_MSC_VER)
  #define EMB_LOG_FUNC_NAME_LONG __FUNCSIG__ << " "
#else
  #define EMB_LOG_FUNC_NAME_LONG __func__ << "() "
#endif

// short function name
#if defined (_MSC_VER)
  #define EMB_LOG_FUNC_NAME __FUNCTION__ << "() "
#else
  #define EMB_LOG_FUNC_NAME __func__ << "() "
#endif

} // embree
