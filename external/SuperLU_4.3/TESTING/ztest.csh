#!/bin/csh

set ofile = ztest.out			# output file
if ( -e $ofile ) then
    rm -f $ofile
endif
echo "Double-precision complex testing output" > $ofile

set MATRICES     = (LAPACK cg20.cua)
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
            ./ztest -t "LA" -l $l -n $n -s $s >> $ofile
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
            ./ztest -t "SP" -s $s -l $l  < ../EXAMPLE/$m >> $ofile
        end
    end
  endif

end

