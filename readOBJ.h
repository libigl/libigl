//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef READOBJ_H
#define READOBJ_H

#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace igl 
{
    //! Read a mesh from an ascii obj file
    void readOBJ(std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        std::ifstream s(str.c_str());
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
    }
}

#endif