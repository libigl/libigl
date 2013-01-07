#!/bin/csh

set ofile = stest.out			# output file
if ( -e $ofile ) then
    rm -f $ofile
endif
echo "Single-precision testing output" > $ofile

set MATRICES     = (LAPACK g20.rua)
set NVAL         = (9 19)
set NRHS         = (5)
set LWORK        = (0 10000000)

#
# Loop through all matrices ...
#
foreach m ($MATRICES)

  #--------------------------------------------
  # Test matrix types generated in LAPACK-style
  #--------------------------------------------
  if  ($m == 'LAPACK') then
      echo '== LAPACK test matrices' >> $ofile
      foreach n ($NVAL)
        foreach s ($NRHS)
          foreach l ($LWORK)
	    echo '' >> $ofile
            echo 'n='$n 'nrhs='$s 'lwork='$l >> $ofile
            ./stest -t "LA" -l $l -n $n -s $s >> $ofile
          end
        end
      end
  #--------------------------------------------
  # Test a specified sparse matrix
  #--------------------------------------------
  else
    echo '' >> $ofile
    echo '== sparse matrix:' $m >> $ofile
    foreach s ($NRHS)
        foreach l ($LWORK)
	    echo '' >> $ofile
            echo 'nrhs='$s 'lwork='$l >> $ofile
            ./stest -t "SP" -s $s -l $l < ../EXAMPLE/$m >> $ofile
        end
    end
  endif

end


