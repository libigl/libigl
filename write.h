//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef WRITE_H
#define WRITE_H

#include <Eigen/Core>
#include <string>

#include <writeOBJ.h>
#include <writeOFF.h>

namespace igl 
{
    // write mesh to an ascii file with automatic detection of file format. supported: obj, off)
    inline void write(std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        const char* p;
        for (p = str.c_str(); *p != '\0'; p++)
            ;
        while (*p != '.')
            p--;
        
        if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
            return igl::writeOBJ(str,V,F);
        
        if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
            return igl::writeOFF(str,V,F);
    }
}

#endif
