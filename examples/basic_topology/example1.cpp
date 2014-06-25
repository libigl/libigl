//
//  IGL Lib - Simple C++ mesh library
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// IMPORTANT DO NOT REMOVE OR MOVE
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
#include <string>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>

using namespace std;

int main (int argc, const char * argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("../shared/TinyTorus.obj",V,F);

    std::cout << "Mesh loaded!\n";
    cout << "Vertex Array:" << endl;
    cout << V << endl;
    cout << "-------------" << endl;
    cout << "Face Array:" << endl;
    cout << F << endl;
    cout << "-------------" << endl;

    igl::write_triangle_mesh("bunny_out.off",V,F);

    // Face Topology
    cout << "TT Topology:" << endl;
    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(V,F,TT);
    cout << TT << endl;
    cout << "-------------" << endl;

    // Edge Topology
    cout << "Edge Topology:" << endl;
    Eigen::MatrixXi EV;
    Eigen::MatrixXi FE;
    Eigen::MatrixXi EF;

    igl::edge_topology(V,F,EV,FE, EF);
    cout << EV << endl << FE << endl << EF << endl;
    cout << "-------------" << endl;


    return 0;
}
