function y = time()
%TIME   Timestamp.
%       S = TIME returns a string containing the date and time.
 
%       Copyright (c) 1984-93 by The MathWorks, Inc.
 
t = clock;
base = t(1) - rem(t(1),100);
months = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';
          'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
y = [int2str(t(3)),'-',months(t(2),:),'-',int2str(t(1)-base), ...
     '  ',int2str(t(4)),':', int2str(t(5)),':', int2str(t(6))];
