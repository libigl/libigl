function [r,s] = spart2(A)
% SPART2 : supernode partition
%
% [r,x] = spart2(A) partitions the columns of A according to supernode 
% definition 2:
% A supernode in A = L+U is a sequence of adjacent columns in which
% the diagonal block of L is full,
% and below the diagonal all the columns (of L) have the same row structure.
% Output: row and column partitions r and s suitable for SPYPART(A,r,s)
%
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


[nr,nc] = size(A);
A = spones(A);
A = tril(A) | speye(nr,nc);

A1 = tril([zeros(nr,1) A]);
A2 = [A ones(nr,1)];

signature = sum(xor(A1,A2));
r = find(signature);
r = r';
if nargout > 1,
    s = r;
    r = [1 nr+1];
end;
