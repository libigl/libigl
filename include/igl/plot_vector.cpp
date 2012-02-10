#include "plot_vector.h"
#include <iostream>
#include <vector>


template <typename T>
IGL_INLINE void igl::plot_vector( std::vector<T>& v)
{
  using namespace std;
  for (int i=0; i<v.size(); ++i)
    std::cerr << v[i] << " ";
  std::cerr << std::endl;
}

template <typename T>
IGL_INLINE void igl::plot_vector( std::vector< std::vector<T> >& v)
{
  using namespace std;
  for (int i=0; i<v.size(); ++i)
  {
    std::cerr << i << ": ";
    for (int j=0; j<v[i].size(); ++j)
      std::cerr << v[i][j] << " ";
    std::cerr << std::endl;
  }
}


template <typename T>
IGL_INLINE void igl::plot_vector( std::vector< std::vector< std::vector<T> > >& v)
{
  using namespace std;
  for (int m=0; m<v.size(); ++m)
  {
    std::cerr << "Matrix " << m << std::endl;

    for (int i=0; i<v[m].size(); ++i)
    {
      std::cerr << i << ": ";
      for (int j=0; j<v[m][i].size(); ++j)
        std::cerr << v[m][i][j] << " ";
      std::cerr << std::endl;
    }
    
    std::cerr << "---- end " << m << std::endl;

  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
