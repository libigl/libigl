      program f77_main
      integer maxn, maxnz
      parameter ( maxn = 10000, maxnz = 100000 )
      integer rowind(maxnz), colptr(maxn)
      real*8  values(maxnz), b(maxn)
      integer n, nnz, nrhs, ldb, info, iopt
      integer*8 factors
*
      call hbcode1(n, n, nnz, values, rowind, colptr)
*
      nrhs = 1
      ldb = n
      do i = 1, n
         b(i) = 1
      enddo
*
* First, factorize the matrix. The factors are stored in *factors* handle.
      iopt = 1
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      if (info .eq. 0) then
         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'INFO from factorization = ', info
      endif
*
* Second, solve the system using the existing factors.
      iopt = 2
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      if (info .eq. 0) then
         write (*,*) 'Solve succeeded'
         write (*,*) (b(i), i=1, 10)
      else
         write(*,*) 'INFO from triangular solve = ', info
      endif

* Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_fortran_dgssv( iopt, n, nnz, nrhs, values, rowind, colptr, 
     $                      b, ldb, factors, info )
*
      stop
      end


