% Converts a python-wrapped Eigen Matrix to a Matlab matrix
function [ M ] = p2m( P )
    if py.repr(py.type(P)) == '<class ''igl.eigen.MatrixXd''>'
        % Convert it to a python array first
        t = py.array.array('d',P);
        % Reshape it
        M = reshape(double(t),P.rows(),P.cols());
    elseif py.repr(py.type(P)) == '<class ''igl.eigen.MatrixXi''>'
        % Convert it to a python array first
        t = py.array.array('i',P);
        % Reshape it
        M = reshape(int32(t),P.rows(),P.cols());
    else
        error('Unsupported numerical type.');
    end
end

