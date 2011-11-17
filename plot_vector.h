#ifndef IGL_PLOT_VECTOR_H
#define IGL_PLOT_VECTOR_H

namespace igl 
{
  template <typename T>
  inline void plot_vector( vector<T>& v)
  {
    for (int i=0; i<v.size(); ++i)
      cerr << v[i] << " ";
    cerr << endl;
  }

  template <typename T>
  inline void plot_vector( vector< vector<T> >& v)
  {
    for (int i=0; i<v.size(); ++i)
    {
      for (int j=0; j<v[i].size(); ++j)
        cerr << v[i][j] << " ";
      cerr << endl;
    }
  }

}

#endif
