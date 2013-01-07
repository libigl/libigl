      subroutine hbcode1(nrow, ncol, nnzero, values, rowind, colptr)

C     ================================================================
C     ... SAMPLE CODE FOR READING A SPARSE MATRIX IN STANDARD FORMAT
C     ================================================================

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 ,
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     1               NROW  , NCOL  , NNZERO, NELTVL

      INTEGER        COLPTR (*), ROWIND (*)

      REAL*8         VALUES (*)

C    ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------

      READ ( *, 1000 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO, NELTVL,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )

C     -------------------------
C     ... READ MATRIX STRUCTURE
C     -------------------------

      READ ( *, PTRFMT ) ( COLPTR (I), I = 1, NCOL+1 )

      READ ( *, INDFMT ) ( ROWIND (I), I = 1, NNZERO )

      IF  ( VALCRD .GT. 0 )  THEN

C         ----------------------
C         ... READ MATRIX VALUES
C         ----------------------

          READ ( *, VALFMT ) ( VALUES (I), I = 1, NNZERO )

      ENDIF

      return
      end
