#include "read.h"

#include "readOBJ.h"
#include "readOFF.h"
IGL_INLINE bool igl::read(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    const char* p;
    for (p = str.c_str(); *p != '\0'; p++)
        ;
    while (*p != '.')
        p--;
    
    if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
    {
        return igl::readOBJ(str,V,F);
    }else if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
    {
        return igl::readOFF(str,V,F);
    }else
    {
      fprintf(stderr,"read() does not recognize extension: %s\n",p);
      return false;
    }
}
