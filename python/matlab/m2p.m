% Converts a Matlab matrix to a python-wrapped Eigen Matrix
function [ P ] = m2p( M )
    if (isa(M, 'double'))
        % Convert the matrix to a python 1D array
        a = py.array.array('d',reshape(M,1,numel(M)));
        % Then convert it to a eigen type
        t = py.igl.eigen.MatrixXd(a.tolist());
        % Finally reshape it back
        P = t.MapMatrix(uint16(size(M,1)),uint16(size(M,2)));
    elseif (isa(M, 'integer'))
        % Convert the matrix to a python 1D array
        a = py.array.array('i',reshape(M,1,numel(M)));
        % Then convert it to a eigen type
        t = py.igl.eigen.MatrixXi(a.tolist());
        % Finally reshape it back
        P = t.MapMatrix(uint16(size(M,1)),uint16(size(M,2)));
    else
        error('Unsupported numerical type.');
    end
