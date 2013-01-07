function failed = trysuperlu(A);
% TRYSUPERLU : Test the Matlab interface to superlu.
%
% failed = trysuperlu;  
%          This runs several tests on superlu, using a matrix with the
%          structure of "smallmesh".  It returns a list of failed tests.
%
% failed = trysuperlu(A);
%          This just times the factorization of A.
%
% Copyright (c) 1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

disp(' ');
disp(['Testing SUPERLU on ' time]);
disp(' ');

resetrandoms;
if nargin == 0
    A = hbo('smallmesh');
end;
[n,m] = size(A);
if n~=m, error('matrix must be square'), end;

failed = [];
ntest = 0;
tol = 100 * eps * n^2;

if nargin == 1 % Just try to factor the given input matrix.
    ntest=ntest+1;
    fprintf('Matrix size: %d by %d with %d nonzeros.\n',n,m,nnz(A));
    r = trytime(A);
    if r > tol, disp('*** FAILED ***'); failed = [failed ntest]; end;
    return;
end;

% With no input argument, try several tests on the small mesh.

PivPostA = sprandn(A);
A = (n+1)*speye(n) - spones(A); % Diagonally dominant
[t,post] = etree(A'*A);

if post == [1:n], disp('note: column etree already in postorder'), end;

NoPivPostA = A(post,:);
NoPivNoPostA = A(post,post);
PivNoPostA = sprandn(NoPivNoPostA);
p = randperm(n);
pinv = zeros(1,n);
pinv(p) = 1:n;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, no postorder, 4 outputs.\n',ntest);
r = try4(NoPivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, no postorder, 3 outputs.\n',ntest);
r = try3(NoPivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, no postorder, 2 outputs.\n',ntest);
r = try2(NoPivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, no postorder, 4 outputs.\n',ntest);
r = try4(PivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, no postorder, 3 outputs.\n',ntest);
r = try3(PivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, no postorder, 2 outputs.\n',ntest);
r = try2(PivNoPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, postorder, 4 outputs.\n',ntest);
r = try4(NoPivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, postorder, 3 outputs.\n',ntest);
r = try3(NoPivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, no pivot, postorder, 2 outputs.\n',ntest);
r = try2(NoPivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, postorder, 4 outputs.\n',ntest);
r = try4(PivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, postorder, 3 outputs.\n',ntest);
r = try3(PivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm 0, pivot, postorder, 2 outputs.\n',ntest);
r = try2(PivPostA,0);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, no postorder, 4 outputs.\n',ntest);
r = try4(NoPivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, no postorder, 3 outputs.\n',ntest);
r = try3(NoPivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, no postorder, 2 outputs.\n',ntest);
r = try2(NoPivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, no postorder, 4 outputs.\n',ntest);
r = try4(PivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, no postorder, 3 outputs.\n',ntest);
r = try3(PivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, no postorder, 2 outputs.\n',ntest);
r = try2(PivNoPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, postorder, 4 outputs.\n',ntest);
r = try4(NoPivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, postorder, 3 outputs.\n',ntest);
r = try3(NoPivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, no pivot, postorder, 2 outputs.\n',ntest);
r = try2(NoPivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, postorder, 4 outputs.\n',ntest);
r = try4(PivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, postorder, 3 outputs.\n',ntest);
r = try3(PivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given, pivot, postorder, 2 outputs.\n',ntest);
r = try2(PivPostA(:,pinv),p);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: No input perm (colamd done internally), pivot, postorder, 4 outputs.\n',ntest);
r = try4(PivPostA);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: No input perm, pivot, postorder, 3 outputs.\n',ntest);
r = try3(PivPostA);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: No input perm, pivot, postorder, 2 outputs.\n',ntest);
r = try3(PivPostA);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Input perm given as matrix, pivot, postorder, 4 outputs.\n',ntest);
P = speye(n);
P = P(p,:);
r = try4(PivPostA*P,P);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

ntest=ntest+1;
fprintf('Test %d: Empty matrix.\n',ntest);
[L,U,PROW,PCOL] = superlu([]);
disp(' ');
if max([size(L) size(U) size(PROW) size(PCOL)])
    L, U, PROW, PCOL
    disp('*** FAILED ***^G'), disp(' ');
    failed = [failed ntest];
end;

ntest=ntest+1;
fprintf('Test %d: Timing versus LU on the 3-element airfoil matrix.\n',ntest);
A = hbo('airfoil2');
n = max(size(A));
A = sprandn(A);
pmmd = colamd(A);
A = A(:,pmmd);
r = trytime(A);
if r > tol, disp('*** FAILED ***'), disp(' '), failed = [failed ntest]; end;

nfailed = length(failed);
if nfailed
    fprintf('\n%d tests failed.\n\n',nfailed);
else
    fprintf('\nAll tests passed.\n\n',nfailed);
end;
