#ifndef IGL_PLOT_VECTOR_H
#define IGL_PLOT_VECTOR_H

namespace igl 
{
  template <typename T>
  inline void plot_vector( std::vector<T>& v)
  {
    for (int i=0; i<v.size(); ++i)
      std::cerr << v[i] << " ";
    std::cerr << std::endl;
  }

  template <typename T>
  inline void plot_vector( std::vector< std::vector<T> >& v)
  {
    for (int i=0; i<v.size(); ++i)
    {
      for (int j=0; j<v[i].size(); ++j)
        std::cerr << v[i][j] << " ";
      std::cerr << std::endl;
    }
  }

}

#endif
