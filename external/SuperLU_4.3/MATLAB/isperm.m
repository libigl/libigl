function result = isperm(p)
% ISPERM        Is the argument a permutation?

result = 0;
if min(size(p)) > 1, return, end;
ds = diff(sort(p));
if any(ds ~= 1), return, end;
if min(p) ~= 1, return, end;
result = 1;
