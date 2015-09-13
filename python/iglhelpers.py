import numpy as np
import scipy.sparse as sparse
import igl

def p2e(m):
    if isinstance(m, np.ndarray):
        if m.dtype.type == np.int32:
            return igl.eigen.MatrixXi(m)
        elif m.dtype.type == np.float64:
            return igl.eigen.MatrixXd(m)
        raise TypeError("p2e only support dtype float64 or int32")
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
        return np.array(m, dtype='float64')
    elif isinstance(m, igl.eigen.MatrixXi):
        return np.array(m, dtype='int32')
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
