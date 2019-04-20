#!/usr/bin/env python
#
# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
from iglhelpers import *

############# Dense Matrix Types #############

# Create a numpy dense array
# 2 types are supported by the wrappers: float64 and int64
dense_matrix = np.array([(1, 2, 3), (4, 5, 6), (7, 8, 9)], dtype='float64')

# libigl wrappers uses Eigen as a matrix type, you can easily convert between numpy and Eigen using
# the helper function p2e. This operation duplicates the data.
dense_matrix_eigen = p2e(dense_matrix)

# The Eigen wrappers allows you to do operations directly on this matrix,
# without having to convert back to numpy
dense_matrix_eigen_2 = dense_matrix_eigen * dense_matrix_eigen

# You can also inspect the data without converting it ...
print("Eigen Matrix: \n", dense_matrix_eigen_2, "\n", sep='')

# and access single elements
print("Eigen Matrix(0,0): ", dense_matrix_eigen_2[0, 0], "\n")

# To convert it back to a numpy array, use the helper function e2p
dense_matrix_2 = e2p(dense_matrix_eigen_2)
print("Numpy Array: \n", dense_matrix_2, "\n", sep='')

############# Sparse Matrix Types #############

# Sparse matrices are handled in a very similar way
# 2 types are supported by the wrappers: float64 and int64
sparse_matrix = sparse.rand(10, 10, 0.1)

# To convert to the eigen forma use p2e
sparse_matrix_eigen = p2e(sparse_matrix)

# They can directly be used plotted or used in computations
print("Sparse matrix Eigen: ", sparse_matrix_eigen, sep='')

# And converted back with e2p
sparse_matrix_2 = e2p(sparse_matrix_eigen)
print("Sparse matrix Numpy: ", sparse_matrix_2.todense(), sep='')
