# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import numpy as np
import scipy.sparse as sparse
import pyigl as igl

def p2e(m):
    if isinstance(m, np.ndarray):
        if not (m.flags['C_CONTIGUOUS'] or m.flags['F_CONTIGUOUS']):
            raise TypeError('p2e support either c-order or f-order')
        if m.dtype.type in [np.int32, np.int64]:
            return igl.eigen.MatrixXi(m.astype(np.int32))
        elif m.dtype.type in [np.float64, np.float32]:
            return igl.eigen.MatrixXd(m.astype(np.float64))
        elif m.dtype.type == np.bool:
            return igl.eigen.MatrixXb(m)
        raise TypeError("p2e only support dtype float64/32, int64/32 and bool")
    if sparse.issparse(m):
        # convert in a dense matrix with triples
        coo = m.tocoo()
        triplets = np.vstack((coo.row, coo.col, coo.data)).T

        triples_eigen_wrapper = igl.eigen.MatrixXd(triplets)

        if m.dtype.type == np.int32:
            t = igl.eigen.SparseMatrixi()
            t.fromcoo(triples_eigen_wrapper)
            return t
        elif m.dtype.type == np.float64:
            t = igl.eigen.SparseMatrixd()
            t.fromCOO(triples_eigen_wrapper)
            return t


    raise TypeError("p2e only support numpy.array or scipy.sparse")


def e2p(m):
    if isinstance(m, igl.eigen.MatrixXd):
        return np.array(m, dtype='float64', order='C')
    elif isinstance(m, igl.eigen.MatrixXi):
        return np.array(m, dtype='int32', order='C')
    elif isinstance(m, igl.eigen.MatrixXb):
        return np.array(m, dtype='bool', order='C')
    elif isinstance(m, igl.eigen.SparseMatrixd):
        coo = np.array(m.toCOO())
        I = coo[:, 0]
        J = coo[:, 1]
        V = coo[:, 2]
        return sparse.coo_matrix((V,(I,J)), shape=(m.rows(),m.cols()), dtype='float64')
    elif isinstance(m, igl.eigen.SparseMatrixi):
        coo = np.array(m.toCOO())
        I = coo[:, 0]
        J = coo[:, 1]
        V = coo[:, 2]
        return sparse.coo_matrix((V,(I,J)), shape=(m.rows(),m.cols()), dtype='int32')

def printMatrixSizes(x,xn):
    print(xn + " (" + str(x.rows()) + "," + str(x.cols()) + ")")
