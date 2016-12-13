#!/usr/bin/env ruby

filename = nil
lineno = nil
ARGF.each do |line|
  if ARGF.filename != filename
    filename = ARGF.filename
    lineno = ARGF.lineno-1
  end
  if line =~ /^ *([A-z0-9]*) *=/
    left = $1
    if line =~ /=.*#{left}.*/
      printf "%s:%d  %s",ARGF.filename, ARGF.lineno-lineno, line
    end
  end
end
