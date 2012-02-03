//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.
//
//
//  Example that shows the integration with matlab
//

// IMPORTANT DO NOT REMOVE OR MOVE
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET 

#include <iostream>
#include <string>
#include <read.h>

#include <matlabinterface.h>

using namespace std;

int main (int argc, const char * argv[])
{
  // This is broken
    //// read the header of matlabinterface.h for compilation instructions
    //
    //Eigen::MatrixXd V,V2;
    //Eigen::MatrixXi F,F2;
    //
    //// Read mesh from file
    //igl::read("bunny.off",V,F);
    //
    //// Send mesh to matlab
    //igl::mlsetmatrix("V",V);
    //igl::mlsetmatrix("F",F);

    //// Plot the mesh from matlab
    //igl::mleval("trimesh(F,V(:,1),V(:,2),V(:,3))");

    //// Receive mesh from matlab
    //igl::mlgetmatrix("V",V2);
    //igl::mlgetmatrix("F",F2);

    //// Plot the received mesh
    //cerr << "V " << endl << V2  << endl;
    //cerr << "F " << endl << F2  << endl;
    //
    //// It is also possible to send scalars
    //igl::mlsetscalar("s", 3);
    //cerr << "s = " << igl::mlgetscalar("s") << endl;

    //// If the program closes the matlab session is killed too..
    //getchar();
    
    return 0;
}

