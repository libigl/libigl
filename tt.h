//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef TT_H
#define TT_H

#include <Eigen/Core>
#include <string>
#include "ismanifold.h"

#include <vector>

namespace igl 
{

  // Preprocessing
  template<typename T> 
  inline void tt_preprocess(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT)
  {
    for(int f=0;f<F.rows();++f)
      for (int i=0;i<3;++i)
      {
        // v1 v2 f ei 
        int v1 = F(f,i);
        int v2 = F(f,(i+1)%3);
        if (v1 > v2) std::swap(v1,v2);
        std::vector<int> r(4);
        r[0] = v1; r[1] = v2;
        r[2] = f;  r[3] = i;
        TTT.push_back(r);
      }
    std::sort(TTT.begin(),TTT.end());
  }
  
  // Extract the face adjacencies
  inline void tt_extractTT(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT, Eigen::MatrixXi& TT)
  {
    TT = Eigen::MatrixXi::Constant((int)(F.rows()),3,-1);
    
    for(int i=1;i<(int)TTT.size();++i)
    {
      std::vector<int>& r1 = TTT[i-1];
      std::vector<int>& r2 = TTT[i];
      if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
      {
        TT(r1[2],r1[3]) = r2[2];
        TT(r2[2],r2[3]) = r1[2];
      }
    }
  }
  
  // Extract the face adjacencies indices (needed for fast traversal)
  inline void tt_extractTTi(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& TTT, Eigen::MatrixXi& TTi)
  {
    TTi = Eigen::MatrixXi::Constant((int)(F.rows()),3,-1);
    
    for(int i=1;i<(int)TTT.size();++i)
    {
      std::vector<int>& r1 = TTT[i-1];
      std::vector<int>& r2 = TTT[i];
      if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
      {
        TTi(r1[2],r1[3]) = r2[3];
        TTi(r2[2],r2[3]) = r1[3];
      }
    }
  }

  // Compute triangle-triangle adjacency
  template<typename T> 
  inline void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& TT)
  {
    assert(igl::isManifold(V,F));
    std::vector<std::vector<int> > TTT;
    
    tt_preprocess<T>(V,F,TTT);
    tt_extractTT(F,TTT,TT);
  }

  // Compute triangle-triangle adjacency with indices
  template<typename T> 
  inline void tt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& TT, Eigen::MatrixXi& TTi)
  {
    assert(igl::isManifold(V,F));
    std::vector<std::vector<int> > TTT;
    
    tt_preprocess<T>(V,F,TTT);
    tt_extractTT(F,TTT,TT);
    tt_extractTTi(F,TTT,TTi);
  }

}

#endif
