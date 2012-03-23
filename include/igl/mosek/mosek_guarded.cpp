#include "mosek_guarded.h"
#include <iostream>

IGL_INLINE MSKrescodee igl::mosek_guarded(const MSKrescodee r)
{
  using namespace std;
  if(r != MSK_RES_OK)
  {
    /* In case of an error print error code and description. */      
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    MSK_getcodedesc(r,symname,desc);
    cerr<<"MOSEK ERROR ("<<r<<"): "<<symname<<" - '"<<desc<<"'"<<endl;
  }
  return r;
}

