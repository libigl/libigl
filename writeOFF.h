//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef WRITEOFF_H
#define WRITEOFF_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
    // write mesh to an ascii off file
    void writeOFF (std::string fname, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        FILE *fp = fopen (fname.c_str(), "w");
        
        if (!fp)
            fprintf (stderr, "writeOFF(): could not open file %s", fname.c_str());
        
        fprintf (fp, "OFF\n%d %d 0\n",  (int) V.rows(), (int) F.rows());
        
        for (unsigned i = 0; i < V.rows(); i++)
            fprintf (fp, "%f %f %f\n", V(i,0), V(i,1), V(i,2));
        
        for (unsigned i = 0; i < F.rows(); i++)
            fprintf (fp, "3 %d %d %d\n", F(i,0), F(i,1), F(i,2));
        
        fclose (fp);
    }
}

#endif