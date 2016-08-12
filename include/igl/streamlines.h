// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_STREAMLINES_H
#define IGL_STREAMLINES_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl
{
    struct StreamlineData
    {
        Eigen::MatrixXi TT;         //  #F by #3 adjacent matrix
        Eigen::MatrixXi E;          //  #E by #3
        Eigen::MatrixXi F2E;        //  #Fx3, Stores the Triangle-Edge relation
        Eigen::MatrixXi E2F;        //  #Ex2, Stores the Edge-Triangle relation
        Eigen::MatrixXd field;      //  #F by 3N list of the 3D coordinates of the per-face vectors
                                    //      (N degrees stacked horizontally for each triangle)
        Eigen::MatrixXi match_ab;   //  #E by N matrix, describing for each edge the matching a->b, where a
                                    //      and b are the faces adjacent to the edge (i.e. vector #i of
                                    //      the vector set in a is matched to vector #mab[i] in b)
        Eigen::MatrixXi match_ba;   //  #E by N matrix, describing the inverse relation to match_ab
        int nsample;                //  #S, number of sample points
        int degree;                 //  #N, degrees of the vector field
    };

    struct StreamlineState
    {
        Eigen::MatrixXd start_point;        //  #N*S by 3 starting points of segment (stacked vertically for each degree)
        Eigen::MatrixXd end_point;          //  #N*S by 3 endpoints points of segment (stacked vertically for each degree)
        Eigen::MatrixXi current_face;       //  #S by N face indices (stacked horizontally for each degree)
        Eigen::MatrixXi current_direction;  //  #S by N field direction indices (stacked horizontally for each degree)

    };


    // Given a mesh and a field the function computes the /data/ necessary for tracing the field'
    // streamlines, and creates the initial /state/ for the tracing.
    // Inputs:
    //   V             #V by 3 list of mesh vertex coordinates
    //   F             #F by 3 list of mesh faces
    //   temp_field    #F by 3n list of the 3D coordinates of the per-face vectors
    //                    (n-degrees stacked horizontally for each triangle)
    //   treat_as_symmetric
    //              if true, adds n symmetry directions to the field (N = 2n). Else N = n
    //   percentage    [0-1] percentage of faces sampled
    // Outputs:
    //   data          struct containing topology information of the mesh and field
    //   state         struct containing the state of the tracing
    IGL_INLINE void streamlines_init(
            const Eigen::MatrixXd V,
            const Eigen::MatrixXi F,
            const Eigen::MatrixXd &temp_field,
            const bool treat_as_symmetric,
            StreamlineData &data,
            StreamlineState &state,
            double percentage = 0.3

    );

    // The function computes the next state for each point in the sample
    //   V             #V by 3 list of mesh vertex coordinates
    //   F             #F by 3 list of mesh faces
    //   data          struct containing topology information
    //   state         struct containing the state of the tracing
    IGL_INLINE void streamlines_next(
            const Eigen::MatrixXd V,
            const Eigen::MatrixXi F,
            const StreamlineData & data,
            StreamlineState & state

    );
}


#ifndef IGL_STATIC_LIBRARY
#  include "streamlines.cpp"
#endif

#endif
