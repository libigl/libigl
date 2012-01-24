#ifndef IGL_PLOT_VECTOR_H
#define IGL_PLOT_VECTOR_H

#include <iostream>
namespace igl 
{
  // Not clear what these are supposed to be doing. Currently they print
  // vectors to standard error...
  template <typename T>
  inline void plot_vector( std::vector<T>& v);
  template <typename T>
  inline void plot_vector( std::vector< std::vector<T> >& v);
  template <typename T>
  inline void plot_vector( std::vector< std::vector< std::vector<T> > >& v);
}

// Implementation

template <typename T>
inline void igl::plot_vector( std::vector<T>& v)
{
  for (int i=0; i<v.size(); ++i)
    std::cerr << v[i] << " ";
  std::cerr << std::endl;
}

template <typename T>
inline void igl::plot_vector( std::vector< std::vector<T> >& v)
{
  for (int i=0; i<v.size(); ++i)
  {
    std::cerr << i << ": ";
    for (int j=0; j<v[i].size(); ++j)
      std::cerr << v[i][j] << " ";
    std::cerr << std::endl;
  }
}


template <typename T>
inline void igl::plot_vector( std::vector< std::vector< std::vector<T> > >& v)
{
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

#endif
