function info = try3(A,pin);
% TRY3 : test SUPERLU with 3 outputs
%
% info = try3(A,pin);
% normally info is the residual norm; 
% but info is at least 10^6 if the factors are not triangular,
% or the permutation isn't a permutation.
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.


n = max(size(A));
info = 0;

if nargin == 1
    [l,u,pr] = superlu(A);
elseif nargin == 2
    [l,u,pr] = superlu(A,pin);
else 
    error('requires 1 or 2 inputs');
end;

figure(2); spy(l); title('L');
figure(1); spy(u); title('U');
drawnow;

format compact
disp('SUPERLU with 3 outputs:');
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

if nargin == 1
    rnorm = norm(A(pr,:) - l*u,inf);
    fprintf(1,'||A(PROW,:) -  L*U|| = %d\n', rnorm);
elseif size(pin) == [n n]
    rnorm = norm(A(pr,:)*pin' - l*u,inf);
    fprintf(1,'||A(PROW,:)*PIN'' -  L*U|| = %d\n', rnorm);
elseif max(size(pin)) == n
    rnorm = norm(A(pr,pin) - l*u,inf);
    fprintf(1,'||A(PROW,PIN) -  L*U|| = %d\n', rnorm);
    if pr == pin
        disp('PROW is the same as the input permutation.');
    else
        disp('PROW is different from the input permutation.');
    end;
else
    rnorm = norm(A(pr,:) - l*u,inf);
    fprintf(1,'||A(PROW,:) -  L*U|| = %d\n', rnorm);
end;

info = info + rnorm;
disp(' ');
