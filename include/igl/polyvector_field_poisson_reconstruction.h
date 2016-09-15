// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_POLYVECTOR_FIELD_POISSON_RECONSTRUCTION
#define IGL_POLYVECTOR_FIELD_POISSON_RECONSTRUCTION
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl {
  // Poisson integration on a combed n-polyvector field, defined on a cut mesh.
  // The function finds n new integrable vector fields that are as close as possible to the
  // N-vector fields of the original combed polyvector field. Integrable means
  // that the vector fields are gradients of scalar functions. The deviation from the input
  // field measures how much the input field deviates from integrability.
  // Inputs:
  //   Vcut             #V by 3 list of the vertex positions
  //   Fcut             #F by 3 list of the faces (must be triangles)
  //   sol3D_combed     #F by 3n list of the 3D coordinates of the per-face vectors of the combed vector set field
  //                    (stacked horizontally for each triangle). Vector #1 in one face will match vector #1 in
  //                    the adjacent face.
  // Outputs:
  //   scalars          #V by n list of the per-vertex scalar functions of which the input field
  //                    is approximately the gradients
  //   sol3D_recon      #F by 3n list of the 3D coordinates of the per-face vectors of the reconstructed
  //                    vector set fields (stacked horizontally for each triangle). The fields are the
  //                    gradients of the scalar functions sF.
  //   max_error        #V by 1 list of the maximal (across the n vector fields) reconstruction error.
  //
  template <typename DerivedV, typename DerivedF, typename DerivedSF, typename DerivedS, typename DerivedE>
  IGL_INLINE double polyvector_field_poisson_reconstruction(
                                                            const Eigen::PlainObjectBase<DerivedV> &Vcut,
                                                            const Eigen::PlainObjectBase<DerivedF> &Fcut,
                                                            const Eigen::PlainObjectBase<DerivedS> &sol3D_combed,
                                                            Eigen::PlainObjectBase<DerivedSF> &scalars,
                                                            Eigen::PlainObjectBase<DerivedS> &sol3D_recon,
                                                            Eigen::PlainObjectBase<DerivedE> &max_error );
  //Wrappers with less output
  template <typename DerivedV, typename DerivedF, typename DerivedSF, typename DerivedS>
  IGL_INLINE void polyvector_field_poisson_reconstruction(
                                                            const Eigen::PlainObjectBase<DerivedV> &Vcut,
                                                            const Eigen::PlainObjectBase<DerivedF> &Fcut,
                                                            const Eigen::PlainObjectBase<DerivedS> &sol3D_combed,
                                                            Eigen::PlainObjectBase<DerivedSF> &scalars);


};


#ifndef IGL_STATIC_LIBRARY
#include "polyvector_field_poisson_reconstruction.cpp"
#endif


#endif /* defined(IGL_POLYVECTOR_FIELD_POISSON_RECONSTRUCTION) */
