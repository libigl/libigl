//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#ifndef IGL_READOBJ_H
#define IGL_READOBJ_H

#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace igl 
{
  //! Read a mesh from an ascii obj file
  // Inputs:
  //   str  path to .obj file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  //
  // KNOWN BUG: This only knows how to read *triangle* meshes. It will probably
  // crash or give garbage on anything else.
  //
  // KNOWN BUG: This only knows how to face lines without normal or texture
  // indices. It will probably crash or give garbage on anything else.
  bool readOBJ(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}

// Implementation
bool igl::readOBJ(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::ifstream s(str.c_str());
    if (s.is_open() == false)
    {
      fprintf (stderr, "readOBJ(): could not open file %s", str.c_str());
      return false;
    }
    std::vector<Eigen::Vector3d> Vtemp;
    std::vector<Eigen::Vector3i> Ftemp;
    char buf[1000];
    while(!s.eof())
    {
        s.getline(buf, 1000);
        if (buf[0] == 'v') // vertex coordinates found
        {
            char v;
            double v1,v2,v3;
            sscanf(buf, "%c %lf %lf %lf",&v,&v1,&v2,&v3);
            Vtemp.push_back(Eigen::Vector3d(v1,v2,v3));
        }
        else if (buf[0] == 'f') // face description found
        {
            char v;
            int v1,v2,v3;
            sscanf(buf, "%c %d %d %d",&v,&v1,&v2,&v3);
            Ftemp.push_back(Eigen::Vector3i(v1-1,v2-1,v3-1));
        }
    }
    s.close();
    
    V = Eigen::MatrixXd(Vtemp.size(),3);
    for(int i=0;i<V.rows();++i)
        V.row(i) = Vtemp[i];
    
    F = Eigen::MatrixXi(Ftemp.size(),3);
    for(int i=0;i<F.rows();++i)
        F.row(i) = Ftemp[i];

    return true;
}

#endif
