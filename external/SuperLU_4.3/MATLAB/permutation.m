% PERMUTATION : Helpful notes on permutation matrices and vectors.
%
%
% There are two ways to represent permutations in Matlab.  
% 
% A *permutation matrix* is a square matrix P that is zero except for exactly
% one 1 in each row and each column.  A *permutation vector* is a 1 by n
% or n by 1 vector p that contains each of the integers 1:n exactly once.
% 
% 
% Let A be any n by n matrix, P be an n by n permutation matrix, and p be
% a permutation vector of length n.  Then:
% 
%   P*A is a permutation of the rows of A.
%   A*P' is the same permutation of the columns of A.  Remember the transpose!
% 
%   A(p,:) is a permutation of the rows of A.
%   A(:,p) is the same permutation of the columns of A.
% 
% 
% Matrix P and column vector p represent the same permutation if and only if
% any of the following three equivalent conditions holds:
% 
%   p = P*(1:n)'     or     P = sparse(1:n,p,1,n,n)     or     P = I(p,:),
% 
% where I = speye(n,n) is the identity matrix of the right size.
% 
% 
% The inverse of a permutation matrix is its transpose:
% 
%   Pinv = inv(P) = P'
% 
% The inverse of a permutation vector can be computed by a one-liner:
% 
%   pinv(p) = 1:n;
% 
% (This works only if pinv doesn't already exist, or exists with the correct 
% dimensions.  To be safe, say:  pinv=zeros(1,n); pinv(p)=1:n; )
% 
% 
% Products of permutations go one way for matrices and the other way for
% vectors.  If P1, P2, and P3 are permutation matrices and p1, p2, and p3
% are the corresponding permutation vectors, then
% 
%   P1 = P2 * P3     if and only if     p1 = p3(p2).
% 
% 
% Permutation vectors are a little more compact to store and efficient 
% to apply than sparse permutation matrices (and much better than full
% permutation matrices).  A sparse permutation matrix takes about twice
% the memory of a permutation vector.  The Matlab function LU returns
% a permutation matrix (sparse if the input was sparse); most other 
% built-in Matlab functions return permutation vectors, including 
% SYMRCM, SYMMMD, COLAMD, DMPERM, and ETREE.

% John Gilbert, 6 April 1995
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.

help permutation
