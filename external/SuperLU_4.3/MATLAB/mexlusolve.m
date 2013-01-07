function x = mexlusolve(A,b,Pcol)
% MEXLUSOLVE : Supernodal LU factor-and-solve.
%
%  MEXLUSOLVE is the mex-file version of the supernodal solver.
%  The user will normally call LUSOLVE, which calls MEXLUSOLVE.
%  See LUSOLVE for a description of parameters.

% Above is the text for HELP MEXLUSOLVE; the following will be executed 
% only if the mex-file appropriate for the machine can't be found.

disp('Warning:  Executable mexlusolve.mex not found for this architecture.');
disp('The supernodal LU package seems not to be installed (or built for Matlab).');
x = zeros(size(A,2),1);
