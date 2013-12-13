#include "path_to_executable.h"
#ifdef __APPLE__
#  include <mach-o/dyld.h>
#endif
std::string igl::path_to_executable()
{
  // http://pastebin.com/ffzzxPzi
  using namespace std;
  std::string path;
  char buffer[1024];
  uint32_t size = sizeof(buffer);
#if defined (WIN32) || defined (WIN64)
  GetModuleFileName(buffer, &size);
  path = buffer;
#elif defined (__APPLE__)
  if(_NSGetExecutablePath(buffer, &size) == 0)
  {
    path = buffer;
  }
#elif defined(UNIX)
  if (readlink("/proc/self/exe", buffer, sizeof(buffer)) == -1)
  {
    path = buffer;
  }
#elif defined(__FreeBSD__)
  int mib[4];
  mib[0] = CTL_KERN;
  mib[1] = KERN_PROC;
  mib[2] = KERN_PROC_PATHNAME;
  mib[3] = -1;
  sysctl(mib, 4, buffer, sizeof(buffer), NULL, 0);
  path = buffer;
#elif defined(SUNOS)
  path = getexecname();
#endif
  return path;
}

