//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef ISMANIFOLD_H
#define ISMANIFOLD_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
    // check if the mesh is edge-manifold
    bool isManifold(Eigen::MatrixXd& V, Eigen::MatrixXd& F)
    {
        vector<vector<int> > TTT;
        for(int f=0;f<F.rows();++f)
            for (int i=0;i<3;++i)
            {
                // v1 v2 f ei 
                int v1 = F(f,i);
                int v2 = F(f,(i+1)%3);
                if (v1 > v2) std::swap(v1,v2);
                vector<int> r(4);
                r[0] = v1; r[1] = v2;
                r[2] = f;  r[3] = i;
                TTT.push_back(r);
            }
        std::sort(TTT.begin(),TTT.end());
        TT = MatrixXi::Constant((int)(F.rows()),3,-1);
        
        for(int i=2;i<TTT.size();++i)
        {
            vector<int>& r1 = TTT[i-2];
            vector<int>& r2 = TTT[i-1];
            vector<int>& r3 = TTT[i];
            if ( (r1[0] == r2[0] && r2[0] == r3[0]) 
                && 
                (r1[1] == r2[1] && r2[1] == r3[1]) )
            {
                return false;
            }
        }
        return true;
    }
}

#endif