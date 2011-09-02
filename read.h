//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef READ_H
#define READ_H

#include <Eigen/Core>
#include <string>

#include <readOBJ.h>
#include <readOFF.h>

namespace igl 
{
    // read mesh from an ascii file with automatic detection of file format. supported: obj, off)
    void read(std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        const char* p;
        for (p = str.c_str(); *p != '\0'; p++)
            ;
        while (*p != '.')
            p--;
        
        if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
        {
            igl::readOBJ(str,V,F);
            return;
        }
        
        if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
        {
            igl::readOFF(str,V,F);
            return;
        }
    }
}

#endif