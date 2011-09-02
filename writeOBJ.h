//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef WRITEOBJ_H
#define WRITEOBJ_H

#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>

namespace igl 
{
    // Write a mesh in an ascii obj file
    void writeOBJ(std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        std::ofstream s(str.c_str());
        for(int i=0;i<V.rows();++i)
            s << "v " << V(i,0) << " " << V(i,1) << " " << V(i,2) << std::endl;
        
        for(int i=0;i<F.rows();++i)
            s << "f " << F(i,0)+1 << " " << F(i,1)+1 << " " << F(i,2)+1 << std::endl;
        
        s.close();
    }
}

#endif