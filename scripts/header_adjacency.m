function [A,H,f2H] = header_adjacency(files)
  % HEADER_ADJACENCY given a list of .cpp/.h files find all <igl/*.h> header
  % files and determine mutual inclusions.
  %
  % Inputs:
  %   files  list of .cpp/.h files
  % Outputs:
  %   A  adjacency matrix so that A(i,j) is the number of mutual includes for
  %     H(i) and H(j)
  %   H  list of unique headers
  %
  % Example:
  %   files = textscan([ ...
  %     ls('/some/dir/*.cpp') ...
  %     ls('/some/dir/*.h')],'%s', 'delimiter', '\n' );
  %   files = files{1};
  %   [A,H] = header_adjacency(files);
  %

  mH = containers.Map();
  f2HI = [];
  f2HJ = [];
  for f = 1:numel(files)
    if isempty(regexp(files{f},'include/igl/[^/]*$'))
      fH = regexp(fileread(files{f}),'<igl/([^.]*.h)>','tokens');
    else
      fH = regexp(fileread(files{f}),'"([^.]*.h)"','tokens');
    end
    fH = [fH{:}]';
    for h = 1:numel(fH)
      h_str = fH{h};
      if ~mH.isKey(h_str)
        mH(h_str) = mH.size(1) + 1;
      end
      f2HI = [f2HI(:); f];
      f2HJ = [f2HJ(:); mH(h_str)];
    end
  end
  f2H = sparse(f2HI,f2HJ,1,numel(files),mH.size(1));
  A = f2H'*f2H;
  % diagonal is not interesting
  A = A - diag(diag(A));
  keys = mH.keys;
  values = mH.values;
  values = [values{:}];
  [~,I] = sort(values);
  H = keys(I);

end
