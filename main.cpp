//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// IMPORTANT DO NOT REMOVE OR MOVE
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET 

#include <iostream>
#include <string>
#include <read.h>
#include <write.h>
#include <cotmatrix.h>

using namespace std;

int main (int argc, const char * argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read(string(argv[1]),V,F);

    std::cout << "Mesh loaded!\n";
    cout << "Vertex Array:" << endl;
    cout << V << endl;
    cout << "-------------" << endl;
    cout << "Face Array:" << endl;
    cout << F << endl;
    cout << "-------------" << endl;

    cout << "CotMatrix:" << endl;
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    cout << L << endl;
    cout << "-------------" << endl;
    
    igl::write(string(argv[2]),V,F);
    
    return 0;
}

