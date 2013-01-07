function resetrandoms
% RESETRANDOMS : Initialize random number generators for repeatability.
%
% John Gilbert, 1994.
% Copyright (c) 1990-1995 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.

rand('seed',0);
randn('seed',0);

