function x = lusolve(A,b,Pcol)
% LUSOLVE : Solve linear systems by supernodal LU factorization.
% 
%  x = lusolve(A, b) returns the solution to the linear system A*x = b,
%      using a supernodal LU factorization that is faster than Matlab's 
%      builtin LU.  This m-file just calls a mex routine to do the work.
%
%  By default, A is preordered by column minimum degree before factorization.
%  Optionally, the user can supply a desired column ordering:
%
%  x = lusolve(A, b, pcol) uses pcol as a column permutation.  
%      It still returns x = A\b, but it factors A(:,pcol) (if pcol is a 
%      permutation vector) or A*Pcol (if Pcol is a permutation matrix).
%       
%  x = lusolve(A, b, 0) suppresses the default minimum degree ordering;
%      that is, it forces the identity permutation on columns.
%
%  See also SUPERLU.
%
% John Gilbert, 6 April 1995
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


[m,n] = size(A);
if m ~= n 
    error('matrix must be square'); 
end;
[mb,nb] = size(b);
if mb ~= n 
    error('right-hand side must have same row dimension as matrix');
end;
if n == 0
    x = [];
    return;
end;

% As necessary, compute the column permutation, and
% convert it from a permutation vector to a permutation matrix
% to fit the internal data structures of mexlusolve.
if nargin < 3 
    Pcol = colamd(A);
end;
if isempty(Pcol) | Pcol == 0
    Pcol = speye(n); 
end;
if min(size(Pcol)) == 1
    Pcol = sparse(1:n,Pcol,1,n,n);
end;

% Make sure the matrices are sparse and the vector is full.
if ~issparse(A)
    A = sparse(A);
end;
if issparse(b)
    b = full(b);
end;
if ~issparse(Pcol)
    Pcol = sparse(Pcol);
end;

x = mexlusolve(A,b,Pcol);
