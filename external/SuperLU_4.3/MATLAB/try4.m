function info = try4(A,pin);
% TRY4 : test SUPERLU with 4 outputs
%
% info = try4(A,pin);
% normally info is the residual norm; 
% but info is at least 10^6 if the factors are not triangular,
% or the permutations aren't permutations.
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


n = max(size(A));
info = 0;

if nargin == 1
    [l,u,pr,pc] = superlu(A);
elseif nargin == 2
    [l,u,pr,pc] = superlu(A,pin);
else 
    error('requires 1 or 2 inputs');
end;

% DOES NOT WORK IN MATLAB 5, NOT YET KNOW WHY -- XSL 5-25-98
% figure(2); [r,s] = spart2(l); spypart(l,r,s); title('L');
figure(2); spy(l); title('L');
figure(1); spy(u); title('U');
drawnow;

format compact
disp('SUPERLU with 4 outputs:');
if any(any(triu(l,1)))
    disp('L is *NOT* lower triangular.');
    info = info + 10^6;
else
    disp('L is lower triangular.');
end;
if nnz(l) == nnz(l+l)
    disp('L has no explicit zeros.');
else
    disp('L contains explicit zeros.');
    info = info+10^6;
end;
if any(any(tril(u,-1)))
    disp('U is *NOT* upper triangular.');
    info = info + 10^6;
else
    disp('U is upper triangular.');
end;
if nnz(u) == nnz(u+u)
    disp('U has no explicit zeros.');
else
    disp('U contains explicit zeros.');
    info = info+10^6;
end;
if pr == [1:n]
    disp('PROW is the identity permutation.');
elseif isperm(pr)
    disp('PROW is a non-identity permutation.');
else 
    disp('PROW is *NOT* a permutation.');
    info = info + 10^6;
end;
if pc == [1:n]
    disp('PCOL is the identity permutation.');
elseif isperm(pc)
    disp('PCOL is a non-identity permutation.');
else 
    disp('PCOL is *NOT* a permutation.');
    info = info + 10^6;
end;
if pr == pc
    disp('PROW and PCOL are the same.')
else
    disp('PROW and PCOL are different.')
end;
rnorm = norm(A(pr,pc) - l*u,inf);
fprintf(1,'||A(PROW,PCOL) -  L*U|| = %d\n', rnorm);
info = info + rnorm;
disp(' ');
