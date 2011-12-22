#ifndef IGL_ADJACENCY_MATRIX_H
#define IGL_ADJACENCY_MATRIX_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <plot_vector.h>
namespace igl 
{
  // Constructs the graph adjacency list of a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F       #F by dim list of mesh faces (must be triangles)
  //   sorted  flag that indicates if the list should be sorted counter-clockwise
  // Outputs: 
  //   A  vector<vector<T> > containing at row i the adjacent vertices of vertex i
  //
  // Example:
  //   // Mesh in (V,F)
  //   vector<vector<double> > A;
  //   adjacency_list(F,A);
  //
  // See also: edges, cotmatrix, diag
  template <typename T, typename M>
  inline void adjacency_list(
                             const M & F, 
                             std::vector<std::vector<T> >& A,
                             bool sorted = false
                             );
}

// Implementation
#include "verbose.h"

template <typename T, typename M>
inline void igl::adjacency_list(
                                const M & F, 
                                std::vector<std::vector<T> >& A,
                                bool sorted 
                                )
{
  A.clear(); 
  A.resize(F.maxCoeff()+1);
  
  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,(j+1)%F.cols());
      A.at(s).push_back(d);
      A.at(d).push_back(s);
    }
  }
  
  // Remove duplicates
  for(int i=0; i<A.size();++i)
  {
    std::sort(A[i].begin(), A[i].end());
    A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
  }
  
  // If needed, sort every VV
  if (sorted)
  {
	  // Loop over faces
	  
	  // for every vertex v store a set of ordered edges not incident to v that belongs to triangle incident on v.
    std::vector<std::vector<std::vector<int> > > SR; 
	  SR.resize(A.size());
	  
	  for(int i = 0;i<F.rows();i++)
	  {
      // Loop over this face
      for(int j = 0;j<F.cols();j++)
      {
        // Get indices of edge: s --> d
        int s = F(i,j);
        int d = F(i,(j+1)%F.cols());
        // Get index of opposing vertex v
        int v = F(i,(j+2)%F.cols());
        
        std::vector<int> e(2);
        e[0] = d;
        e[1] = v;
        SR[s].push_back(e);
      }
	  }
	  
    for(int v=0; v<SR.size();++v)
    {
      std::vector<int>& vv = A.at(v);
      std::vector<std::vector<int> >& sr = SR[v];
      
      std::vector<std::vector<int> > pn = sr;
      
      // Compute previous/next for every element in sr
      for(int i=0;i<sr.size();++i)
      {
        int a = sr[i][0];
        int b = sr[i][1];
        
        // search for previous
        int p = -1;
        for(int j=0;j<sr.size();++j)
          if(sr[j][1] == a)
            p = j;
        pn[i][0] = p;
        
        // search for next
        int n = -1;
        for(int j=0;j<sr.size();++j)
          if(sr[j][0] == b)
            n = j;
        pn[i][1] = n;
        
      }
      
      // assume manifoldness (look for beginning of a single chain)
      int c = 0;
      for(int j=0; j<=sr.size();++j)
        if (pn[c][0] != -1)
          c = pn[c][0];
      
      if (pn[c][0] == -1) // border case
      {
        // finally produce the new vv relation
        for(int j=0; j<sr.size();++j)
        {
          vv[j] = sr[c][0];
          if (pn[c][1] != -1)
            c = pn[c][1];
        }
        vv.back() = sr[c][1];
      }
      else
      {
        // finally produce the new vv relation
        for(int j=0; j<sr.size();++j)
        {
          vv[j] = sr[c][0];
          
          c = pn[c][1];
        }
      }
    }
  }
}

#endif
