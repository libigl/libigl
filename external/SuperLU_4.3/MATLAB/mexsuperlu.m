function [L,U,prow,pcol] = mexsuperlu(A,psparse);
% MEXSUPERLU : Supernodal LU factorization
%
%  MEXSUPERLU is the mex-file version of the supernodal factorization.
%  The user will normally call SUPERLU, which calls MEXSUPERLU.
%  See SUPERLU for a description of parameters.

% Above is the text for HELP MEXSUPERLU; the following will be executed 
% only if the mex-file appropriate for the machine can't be found.

disp('Warning:  Executable mexsuperlu.mex not found for this architecture.');
disp('The supernodal LU package seems not to be installed for Matlab.');
n = max(size(A));
L = sparse(n,n);
U = sparse(n,n);
prow = 1:n;
pcol = 1:n;
