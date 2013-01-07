function failed = trylusolve(A,b);
% TRYLUSOLVE : Test the Matlab interface to lusolve.
%
% failed = trylusolve;  
%          This runs several tests on lusolve, using a matrix with the
%          structure of "smallmesh".  It returns a list of failed tests.
%
% failed = trylusolve(A); 
% failed = trylusolve(A,b); 
%          This times the solution of a system with the given matrix,
%          both with lusolve and with Matlab's builtin "\" operator.
%          Usually, "\" will be faster for triangular or symmetric postive
%          definite matrices, and lusolve will be faster otherwise.
%
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

disp(' ');
disp(['Testing LUSOLVE on ' time]);
disp(' ');

format compact
resetrandoms;
if nargin == 0
    A = hbo('smallmesh');
end;
[n,m] = size(A);
if n~=m, error('matrix must be square'), end;

failed = [];
ntest = 0;
tol = 100 * eps * n^2;

if nargin >= 1 % Just try a solve with the given input matrix.
    ntest=ntest+1;
    v = spparms;
    spparms('spumoni',1);
    fprintf('Matrix size: %d by %d with %d nonzeros.\n',n,m,nnz(A));
    disp(' ');
    if nargin == 1
        b = rand(n,1);
    end;
    tic; x=A\b; t1=toc;
    fprintf('A\\b time = %d seconds.\n',t1);
    disp(' ');
    tic; x=lusolve(A,b); t2=toc;
    fprintf('LUSOLVE time = %d seconds.\n',t2);
    disp(' ');
    ratio = t1/t2;
    fprintf('Ratio = %d\n',ratio);
    residual_norm = norm(A*x-b,inf)/norm(A,inf)
    disp(' ');
    if residual_norm > tol
         disp('*** FAILED ***'); failed = [failed ntest];
    end;
    spparms(v);
    return;
end;

% With no inputs, try several tests on the small mesh.


A = sprandn(A);
p = randperm(n);
I = speye(n);
P = I(p,:);
b = rand(n,2);

ntest=ntest+1;
fprintf('Test %d: Input perm 0.\n',ntest);
x = lusolve(A,b,0);
residual_norm = norm(A*x-b,inf)
disp(' ');
if residual_norm > tol 
    disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; 
end;

ntest=ntest+1;
fprintf('Test %d: Input perm given as vector.\n',ntest);
x = lusolve(A,b,p);
residual_norm = norm(A*x-b,inf)
disp(' ');
if residual_norm > tol
    disp('*** FAILED ***^G'), disp(' '), failed = [failed ntest];
end;

ntest=ntest+1;
fprintf('Test %d: Input perm given as matrix.\n',ntest);
x = lusolve(A,b,P);
residual_norm = norm(A*x-b,inf)
disp(' ');
if residual_norm > tol 
    disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; 
end;

ntest=ntest+1;
fprintf('Test %d: No input perm (colamd done internally).\n',ntest);
x = lusolve(A,b);
residual_norm = norm(A*x-b,inf)
disp(' ');
if residual_norm > tol 
    disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; 
end;

ntest=ntest+1;
fprintf('Test %d: Empty matrix.\n',ntest);
x = lusolve([],[]);
disp(' ');
% if max(size(x))
if length(x)
    x
    disp('*** FAILED ***^G'), disp(' ');
    failed = [failed ntest];
end;

ntest=ntest+1;
fprintf('Test %d: Timing versus \\ on the 3-element airfoil matrix.\n',ntest);
A = hbo('airfoil2');
n = max(size(A));
A = sprandn(A);
b = rand(n,1);
tic; x=A\b; t1=toc;
fprintf('A\\b time = %d seconds.\n',t1);
tic; x=lusolve(A,b); t2=toc;
fprintf('LUSOLVE time = %d seconds.\n',t2);
ratio = t1/t2;
fprintf('Ratio = %d\n',ratio);
disp(' ');
if ratio < 1, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;


nfailed = length(failed);
if nfailed
    fprintf('\n%d tests failed.\n\n',nfailed);
else
    fprintf('\nAll tests passed.\n\n',nfailed);
end;
