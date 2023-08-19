// high level interface for MshLoader.h/.cpp

// Copyright (C) 2020 Vladimir Fonov <vladimir.fonov@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distribute
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READ_MSH_H
#define IGL_READ_MSH_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>
#include <vector>


namespace igl 
{
    /// read triangle surface mesh and tetrahedral volume mesh from .msh file
    ///
    /// @param[in] msh - file name
    /// @param[out] X  eigen double matrix of vertex positions  #X by 3
    /// @param[out] Tri  #Tri eigen integer matrix of triangular faces indices into vertex positions
    /// @param[out] Tet  #Tet eigen integer matrix of tetrahedral indices into vertex positions
    /// @param[out] TriTag #Tri eigen integer vector of tags associated with surface faces
    /// @param[out] TetTag #Tet eigen integer vector of tags associated with volume elements
    /// @param[out] XFields #XFields list of strings with field names associated with nodes
    /// @param[out] XF      #XFields list of eigen double matrices, fields associated with nodes 
    /// @param[out] EFields #EFields list of strings with field names associated with elements
    /// @param[out] TriF    #EFields list of eigen double matrices, fields associated with surface elements
    /// @param[out] TetF    #EFields list of eigen double matrices, fields associated with volume elements
    /// @return true on success
    /// \bug only version 2.2 of .msh file is supported (gmsh 3.X)
    /// \bug only triangle surface elements and tetrahedral volumetric elements are supported
    /// \bug only 3D information is supported
    /// \bug only the 1st tag per element is returned (physical) 
    /// \bug same element fields are expected to be associated with surface elements and volumetric elements
    IGL_INLINE bool readMSH(
      const std::string &msh,
      Eigen::MatrixXd &X,
      Eigen::MatrixXi &Tri,
      Eigen::MatrixXi &Tet,
      Eigen::VectorXi &TriTag,
      Eigen::VectorXi &TetTag,
      std::vector<std::string> &XFields,
      std::vector<Eigen::MatrixXd> &XF,
      std::vector<std::string> &EFields,
      std::vector<Eigen::MatrixXd> &TriF,
      std::vector<Eigen::MatrixXd> &TetF
      );
    /// \overload
    IGL_INLINE bool readMSH(const std::string &msh,
                Eigen::MatrixXd &X,
                Eigen::MatrixXi &Tri,
                Eigen::MatrixXi &Tet,
                Eigen::VectorXi &TriTag,
                Eigen::VectorXi &TetTag
                );
    /// \overload
    IGL_INLINE bool readMSH(const std::string &msh,
                Eigen::MatrixXd &X,
                Eigen::MatrixXi &Tri,
                Eigen::VectorXi &TriTag
                );
    /// \overload
    IGL_INLINE bool readMSH(const std::string &msh,
                Eigen::MatrixXd &X,
                Eigen::MatrixXi &Tri
                );

}


#ifndef IGL_STATIC_LIBRARY
#  include "readMSH.cpp"
#endif

#endif //IGL_READ_MSH_H
