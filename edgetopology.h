//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef EDGETOPOLOGY_H
#define EDGETOPOLOGY_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
    // Initialize Edges and their topological relations
    
    // Output:
    // E  : #Ex2, Stores the edge description as pair of indices to vertices
    // FE : #Fx3, Stores the Triangle-Edge relation
    // EF : #Ex2: Stores the Edge-Triangle relation (unsorted)

    void edgetopology(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
                      Eigen::MatrixXi& E, Eigen::MatrixXi& FE, Eigen::MatrixXi& EF)
    {
        assert(isManifold());
        vector<vector<int> > ETT;
        for(int f=0;f<F.rows();++f)
            for (int i=0;i<3;++i)
            {
                // v1 v2 f vi 
                int v1 = F(f,i);
                int v2 = F(f,(i+1)%3);
                if (v1 > v2) std::swap(v1,v2);
                vector<int> r(4);
                r[0] = v1; r[1] = v2;
                r[2] = f;  r[3] = i;
                ETT.push_back(r);
            }
        std::sort(ETT.begin(),ETT.end());
        
        // count the number of edges (assume manifoldness)
        int En = 1; // the last is always counted
        for(int i=0;i<ETT.size()-1;++i)
            if (!((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1])))
                ++En;
        
        E  = MatrixXi::Constant((int)(En),2,-1);
        FE = MatrixXi::Constant((int)(F.rows()),3,-1);
        EF = MatrixXi::Constant((int)(En),2,-1);
        En = 0;
        
        for(int i=0;i<ETT.size();++i)
        {
            
            for(int j=0;j<ETT[i].size();++j)
            {
                cerr << ETT[i][j] << "\t\t";
            }
            cerr << endl;
        }
        
        
        for(int i=0;i<ETT.size();++i)
        {
            if (i == ETT.size()-1 ||
                !((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1]))
                )
            {
                // Border edge
                vector<int>& r1 = ETT[i];
                E(En,0)         = r1[0];
                E(En,1)         = r1[1];
                EF(En,0)        = r1[2];
                FE(r1[2],r1[3]) = En;
            } 
            else
            {
                vector<int>& r1 = ETT[i];
                vector<int>& r2 = ETT[i+1];
                E(En,0)         = r1[0];
                E(En,1)         = r1[1];
                EF(En,0)        = r1[2];
                EF(En,1)        = r2[2];
                FE(r1[2],r1[3]) = En;
                FE(r2[2],r2[3]) = En;
                ++i; // skip the next one
            }
            ++En;
        }        
    }
}

#endif