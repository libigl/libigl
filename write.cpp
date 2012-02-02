#include "write.h"

#include "writeOBJ.h"
#include "writeOFF.h"

template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::write(
                      const std::string str,
                      const Eigen::PlainObjectBase<DerivedV>& V,
                      const Eigen::PlainObjectBase<DerivedF>& F)
{
  const char* p;
  for (p = str.c_str(); *p != '\0'; p++)
    ;
  while (*p != '.')
    p--;
  
  if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
    return igl::writeOBJ(str,V,F);
  
  if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
    return igl::writeOFF(str,V,F);
}
