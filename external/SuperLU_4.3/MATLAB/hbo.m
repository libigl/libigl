function [A,xy,B,hbtype,rhstype] = hbo(file)
%HBO	Read and process a Harwell-Boeing sparse matrix file.
%	A = HBO('matfile') gets the sparse matrix from the specified .mat file.
%
%	[A,xy,B,hbtype,rhstype] = HBO('matfile') also gets:
%       xy: node coordinates (if present)
%       B:  right-hand sides, starting guesses, or exact solutions (if present)
%       hbtype: matrix type, which will be one of the following three-letter
%	        codes: CSA, PRA, PSA, PSE, PUA, RRA, RSA, RUA, RZA, UNK
%       rhstype: a description of B, for which see the user's guide.
%
%	In addition to reading the file, HBO assembles any unassembled matrices,
%	symmetrizes any symmetric matrices, and implicitizes any explicit zeros.
%
%	HBO .mat files are created from Harwell-Boeing data files by the
%	stand-alone Fortran process, HBO2MAT.  These .mat files contain:
%
%	    For assembled, type xxA, matrices:
%	        A       - the sparse matrix.
%	        hbtitle - The first 72 characters of the first "card".
%	        hbname  - The matrix name, same as file name without .mat.
%	        hbtype  - One of those three letter codes.
%	        hbfill  - If present, the value inserted in pattern matrices.
%	        hbzero  - If present, the value inserted for explicit zeros.
%	        rhstype - If present, the right hand side type.
%               xy      - If present, a matrix whose rows are node coordinates.
%               B       - If present, a matrix whose columns are right-hand
%                         sides, etc.
%
%	    For unassembled, xxE, matrices:
%	        hbtitle - The first 72 characters of the first "card".
%	        hbname  - The matrix name, same as file name without .mat.
%	        hbtype  - Only type PSE exists in current collection.
%	        varind and elptr - The indices specifying the locations
%	            of the constituent elements.
%
%  Say "help harwell" for more information about the collection.
%  See also HBOLIST, HBOFIND.

%	Cleve Moler, The MathWorks, 4/2/94.

load(file)

if ~exist('hbtype')
   hbtype = 'UNK';
   A = A;

% PSE  - Pattern symmetric unassembled
elseif strcmp(hbtype,'PSE')
   n1 = length(elptr);
   elptr(n1) = [];
   J = varind;
   I = zeros(size(J));
   I(elptr) = ones(size(elptr));
   I = cumsum(I);
   A = sparse(I,J,1);
   A = A'*A;

% RSA  - Real symmetric
elseif strcmp(hbtype,'RSA')
   A = A + A' - diag(diag(A));

% RZA  - Real skew symmetric
elseif strcmp(hbtype,'RZA')
   A = A - A';

% RUA  - Real unsymmetric
elseif strcmp(hbtype,'RUA')
   A = A;

% RRA  - Real rectangular
elseif strcmp(hbtype,'RRA')
   A = A;

% CSA  - Complex symmetric
elseif strcmp(hbtype,'CSA')
   A = A + A' - diag(diag(A));

% PSA  - Pattern symmetric
elseif strcmp(hbtype,'PSA')
   A = A + A' - diag(diag(A));

% PUA  - Pattern unsymmetric
elseif strcmp(hbtype,'PUA')
   A = A;

% PRA  - Pattern rectangular
elseif strcmp(hbtype,'PRA')
   A = A;

else
   error(['Harwell-Boeing type ' hbtype ' unexpected.'])
end

% Remove any explict zeros

if exist('hbzero')
   k = find(A == hbzero);
   A(k) = sparse(length(k),1);
end

% Get any right-hand side vectors or coordinates.

if exist('xyz')
    xy = xyz;
end;

if ~exist('xy')
    xy = [];
end;

if ~exist('B')
    B = [];
    rhstype = 'NON';
end;

if ~exist('rhstype')
    rhstype = 'UNK';
end;
