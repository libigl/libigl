function [L,U,prow,pcol] = superlu(A,psparse)
% SUPERLU : Supernodal LU factorization
% 
%  Executive summary:
%
%  [L,U,p] = superlu(A)          is like [L,U,P] = lu(A), but faster.
%  [L,U,prow,pcol] = superlu(A)  preorders the columns of A by min degree,
%                                    yielding A(prow,pcol) = L*U.
%
%  Details and options:
%
%  With one input and two or three outputs, SUPERLU has the same effect as LU,
%  except that the pivoting permutation is returned as a vector, not a matrix:
%
%  [L,U,p] = superlu(A) returns unit lower triangular L, upper triangular U,
%            and permutation vector p with A(p,:) = L*U.
%  [L,U] = superlu(A) returns permuted triangular L and upper triangular U
%            with A = L*U.
%
%  With a second input, the columns of A are permuted before factoring:
%
%  [L,U,prow] = superlu(A,psparse) returns triangular L and U and permutation 
%            prow with A(prow,psparse) = L*U.
%  [L,U] = superlu(A,psparse) returns permuted triangular L and triangular U 
%            with A(:,psparse) = L*U.
%  Here psparse will normally be a user-supplied permutation matrix or vector
%  to be applied to the columns of A for sparsity.  COLAMD is one way to get
%  such a permutation; see below to make SUPERLU compute it automatically.
%  (If psparse is a permutation matrix, the matrix factored is A*psparse'.)
%
%  With a fourth output, a column permutation is computed and applied:
%
%  [L,U,prow,pcol] = superlu(A,psparse)  returns triangular L and U and
%            permutations prow and pcol with A(prow,pcol) = L*U.
%            Here psparse is a user-supplied column permutation for sparsity,
%            and the matrix factored is A(:,psparse) (or A*psparse' if the
%            input is a permutation matrix).  Output pcol is a permutation
%            that first performs psparse, then postorders the etree of the 
%            column intersection graph of A.  The postorder does not affect 
%            sparsity, but makes supernodes in L consecutive.
%  [L,U,prow,pcol] = superlu(A,0) is the same as ... = superlu(A,I); it does
%            not permute for sparsity but it does postorder the etree.
%  [L,U,prow,pcol] = superlu(A) is the same as ... = superlu(A,colamd(A));
%            it uses column minimum degree to permute columns for sparsity,
%            then postorders the etree and factors.
%
%  This m-file calls the mex-file MEXSUPERLU to do the work.
%
% See also COLAMD, LUSOLVE.
%
% John Gilbert, 6 April 1995.
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

% Note on permutation matrices and vectors:
%
% Names beginning with p are permutation vectors; 
% names beginning with P are the corresponding permutation matrices.
%
% Thus  A(pfoo,pbar) is the same as Pfoo*A*Pbar'.
%
% We don't actually form any permutation matrices except Psparse,
% but matrix notation is easier to untangle in the comments.


[m,n] = size(A);
if m ~= n 
    error('matrix must be square.'); 
end;
if n == 0
    L = []; U = []; prow = []; pcol = [];
    return;
end;

% If necessary, compute the column sparsity permutation.
if nargin < 2 
    if nargout >= 4
        psparse = colamd(A);
    else
        psparse = 1:n;
    end;
end;
if max(size(psparse)) <= 1
    psparse = 1:n;
end;

% Compute the permutation-matrix version of psparse,
% which fits the internal data structures of mexsuperlu.
if min(size(psparse)) == 1
    Psparse = sparse(1:n,psparse,1,n,n);
else
    Psparse = psparse;
    psparse = Psparse*[1:n]';
end;

% Make sure the matrices are sparse.
if ~issparse(A)
    A = sparse(A);
end;
if ~issparse(Psparse)
    Psparse = sparse(Psparse);
end;


% The output permutations from the mex-file are dense permutation vectors.
[L,U,prowInv,pcolInv] = mexsuperlu(A,Psparse);
prow = zeros(1,n);
prow(prowInv) = 1:n;
pcol = zeros(1,n);
pcol(pcolInv) = 1:n;

% -- Added 12/14/2011 --
% The row indices of L & U matrices are not sorted from SuperLU, but 
% Matlab requires the matrices to be sorted.
% Now, we do double-transpose to achieve sorting. This acts like bucket sort,
% should be faster than find / sparse combination, which probably uses
% quicksort.
% [i,j,s] = find(L); L = sparse(i,j,s,m,n);

L1 = L'; L = L1';
U1 = U'; U = U1';
% ----

% We now have
%
%    Prow*A*Psparse'*Post' = L*U   (1)
%    Pcol' = Psparse'*Post'         
%
% (though we actually have the vectors, not the matrices, and
% we haven't computed Post explicitly from Pcol and Psparse).
% Now we figure out what the user really wanted returned,
% and rewrite (1) accordingly.

if nargout == 4  
    % Return both row and column permutations.  
    % This is what we've already got.
    % (1) becomes  Prow*A*Pcol' = L*U. 
elseif nargout == 3
    % Return row permutation only.  Fold the postorder perm 
    % but not the sparsity perm into it, and into L and U.
    % This preserves triangularity of L and U.
    % (1) becomes (Post'*Prow) * A * Psparse' = (Post'*L*Post) * (Post'*U*Post).
    postInv = pcolInv(psparse);
    prow = prow(postInv);
    L = L(postInv,postInv);
    U = U(postInv,postInv);
else % nargout <= 2
    % Return no permutations.  Fold the postorder perm but
    % not the sparsity perm into L and U.  Fold the pivoting
    % perm into L, which destroys its triangularity.
    % (1) becomes A * Psparse' = (Prow'*L*Post) * (Post'*U*Post).
    postInv = pcolInv(psparse);
    L = L(prowInv,postInv);
    U = U(postInv,postInv);
end;
