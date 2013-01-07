CCCCC COPYRIGHT (c) 1999  Council for the Central Laboratory of the
CCCCC Research Councils.    All rights reserved.
CCCCC PACKAGE MC64A/AD
CCCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and Jacko Koster (jak@ii.uib.no)
CCCCC LAST UPDATE 20/09/99
CCCCC
C *** Conditions on external use ***
C
C The user shall acknowledge the contribution of this
C package in any publication of material dependent upon the use of
C the package. The user shall use reasonable endeavours to notify
C the authors of the package of this publication.
C
C The user can modify this code but, at no time
C shall the right or title to all or any part of this package pass
C to the user. The user shall make available free of charge
C to the authors for any purpose all information relating to any
C alteration or addition made to this package for the purposes of
C extending the capabilities or enhancing the performance of this
C package.
C
C The user shall not pass this code directly to a third party without the
C express prior consent of the authors.  Users wanting to licence their
C own copy of these routines should send email to hsl@aeat.co.uk
C
C None of the comments from the Copyright notice up to and including this
C one shall be removed or altered in any way.

C**********************************************************************
      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC64A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER ICNTL(10)
C
C  Local variables
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the defaut value, the incoming data is checked for
C     out-of-range indices and duplicates.  Setting ICNTL(4) to any
C     other will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) non-zero, if the data is
C     known to avoid these problems.
C
C    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to
C     zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
C  Purpose
C  =======
C
C This subroutine attempts to find a column permutation for an NxN 
C sparse matrix A = {a_ij} that makes the permuted matrix have N 
C entries on its diagonal.
C If the matrix is structurally nonsingular, the subroutine optionally
C returns a column permutation that maximizes the smallest element 
C on the diagonal, maximizes the sum of the diagonal entries, or 
C maximizes the product of the diagonal entries of the permuted matrix.
C For the latter option, the subroutine also finds scaling factors 
C that may be used to scale the matrix so that the nonzero diagonal 
C entries of the permuted matrix are one in absolute value and all the 
C off-diagonal entries are less than or equal to one in absolute value.
C The natural logarithms of the scaling factors u(i), i=1..N, for the 
C rows and v(j), j=1..N, for the columns are returned so that the 
C scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
C
C  Parameters
C  ==========
C
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   1 Compute a column permutation of the matrix so that the
C     permuted matrix has as many entries on its diagonal as possible. 
C     The values on the diagonal are of arbitrary size. HSL subroutine 
C     MC21A/AD is used for this. See [1].
C   2 Compute a column permutation of the matrix so that the smallest 
C     value on the diagonal of the permuted matrix is maximized.
C     See [3].
C   3 Compute a column permutation of the matrix so that the smallest
C     value on the diagonal of the permuted matrix is maximized.
C     The algorithm differs from the one used for JOB = 2 and may
C     have quite a different performance. See [2].
C   4 Compute a column permutation of the matrix so that the sum
C     of the diagonal entries of the permuted matrix is maximized.
C     See [3].
C   5 Compute a column permutation of the matrix so that the product
C     of the diagonal entries of the permuted matrix is maximized
C     and vectors to scale the matrix so that the nonzero diagonal 
C     entries of the permuted matrix are one in absolute value and 
C     all the off-diagonal entries are less than or equal to one in 
C     absolute value. See [3].
C  Restriction: 1 <= JOB <= 5.
C 
C N is an INTEGER variable which must be set by the user to the
C   order of the matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C
C NE is an INTEGER variable which must be set by the user to the
C   number of entries in the matrix. It is not altered by the 
C   subroutine.
C   Restriction: NE >= 1.
C
C IP is an INTEGER array of length N+1.
C   IP(J), J=1..N, must be set by the user to the position in array IRN 
C   of the first row index of an entry in column J. IP(N+1) must be set
C   to NE+1. It is not altered by the subroutine.
C
C IRN is an INTEGER array of length NE. 
C   IRN(K), K=1..NE, must be set by the user to hold the row indices of
C   the entries of the matrix. Those belonging to column J must be 
C   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
C   of the row indices within each column is unimportant. Repeated 
C   entries are not allowed. The array IRN is not altered by the 
C   subroutine.
C
C A is a REAL (DOUBLE PRECISION in the D-version) array of length NE. 
C   The user must set A(K), K=1..NE, to the numerical value of the 
C   entry that corresponds to IRN(K). 
C   It is not used by the subroutine when JOB = 1. 
C   It is not altered by the subroutine.
C
C NUM is an INTEGER variable that need not be set by the user.
C   On successful exit, NUM will be the number of entries on the 
C   diagonal of the permuted matrix.
C   If NUM < N, the matrix is structurally singular.
C
C CPERM is an INTEGER array of length N that need not be set by the 
C   user. On successful exit, CPERM contains the column permutation.
C   Column CPERM(J) of the original matrix is column J in the permuted 
C   matrix, J=1..N.
C
C LIW is an INTEGER variable that must be set by the user to
C   the dimension of array IW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LIW >= 5N
C     JOB = 2 :  LIW >= 4N
C     JOB = 3 :  LIW >= 10N + NE
C     JOB = 4 :  LIW >= 5N
C     JOB = 5 :  LIW >= 5N
C 
C IW is an INTEGER array of length LIW that is used for workspace.
C 
C LDW is an INTEGER variable that must be set by the user to the
C   dimension of array DW. It is not altered by the subroutine.
C   Restriction:
C     JOB = 1 :  LDW is not used
C     JOB = 2 :  LDW >= N
C     JOB = 3 :  LDW >= NE
C     JOB = 4 :  LDW >= 2N + NE
C     JOB = 5 :  LDW >= 3N + NE
C 
C DW is a REAL (DOUBLE PRECISION in the D-version) array of length LDW 
C   that is used for workspace. If JOB = 5, on return,
C   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N.
C 
C ICNTL is an INTEGER array of length 10. Its components control the 
C   output of MC64A/AD and must be set by the user before calling
C   MC64A/AD. They are not altered by the subroutine.
C
C   ICNTL(1) must be set to specify the output stream for
C   error messages. If ICNTL(1) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(2) must be set by the user to specify the output stream for
C   warning messages. If ICNTL(2) < 0, messages are suppressed.
C   The default value set by MC46I/ID is 6.
C
C   ICNTL(3) must be set by the user to specify the output stream for
C   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
C   The default value set by MC46I/ID is -1.
C
C   ICNTL(4) must be set by the user to a value other than 0 to avoid
C   checking of the input data.
C   The default value set by MC46I/ID is 0.
C    
C INFO is an INTEGER array of length 10 which need not be set by the 
C   user. INFO(1) is set non-negative to indicate success. A negative 
C   value is returned if an error occurred, a positive value if a 
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the 
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : successful entry (for structurally singular matrix).
C   +2 : the returned scaling factors are large and may cause
C        overflow when used to scale the matrix. 
C        (For JOB = 5 entry only.)
C   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : NE < 1. Value of NE held in INFO(2).
C   -4 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -5 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -6 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -7 : repeated entries are found. INFO(2) contains the index of a 
C        column in which such entries are found.
C  INFO(3) to INFO(10) are not currently used and are set to zero by 
C        the routine.
C
C References:
C  [1]  I. S. Duff, (1981),
C       "Algorithm 575. Permutations for a zero-free diagonal",
C       ACM Trans. Math. Software 7(3), 387-390.
C  [2]  I. S. Duff and J. Koster, (1998),
C       "The design and use of algorithms for permuting large
C       entries to the diagonal of sparse matrices",
C       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3]  I. S. Duff and J. Koster, (1999),
C       "On algorithms for permuting large entries to the diagonal 
C       of sparse matrices",
C       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England.

C Local variables and parameters
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
C External routines and functions
c     EXTERNAL FD05AD
c     DOUBLE PRECISION FD05AD
      EXTERNAL MC21AD,MC64BD,MC64RD,MC64SD,MC64WD, DLAMCH
      DOUBLE PRECISION DLAMCH
C Intrinsic functions
      INTRINSIC ABS,LOG

C Set RINF to largest positive real number (infinity)
c XSL    RINF = FD05AD(5)
      RINF = DLAMCH('Overflow')

C Check value of JOB
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Check value of NE
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
C Check LIW
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
C If JOB = 1, do not check
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
C Check row indices. Use IW(1:N) as workspace
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I 
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF

C Print diagnostics on input
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF

C Set components of INFO to zero
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE

C Compute maximum matching with MC21A/AD
      IF (JOB.EQ.1) THEN
C Put length of column J in IW(J)
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
C IW(N+1:5N) is workspace
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.2) THEN
C IW(1:5N), DW(1:N) are workspaces
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF

C Compute bottleneck matching
      IF (JOB.EQ.3) THEN
C Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
C Sort entries in each column by decreasing value. 
        CALL MC64RD(N,NE,IP,IW,DW)
C IW(NE+1:NE+10N) is workspace
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
C B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF

      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
C Check size of scaling factors
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE 
C       GO TO 90
      ENDIF

   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
C Matrix is structurally singular, return with warning
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
C Scaling factors are large, return with warning
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF

C Print diagnostics on output
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D) 
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL (DOUBLE PRECISION in the D-version) array of length
C   NE. A(K), K=1..NE, must be set to the value of the entry
C   that corresponds to IRN(K). It is not altered.
C IPERM is an INTEGER array of length N. On exit, it contains the 
C    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
C NUM is INTEGER variable. On exit, it contains the cardinality of the
C    matching stored in IPERM.
C IW is an INTEGER work array of length 4N.
C DW is a REAL (DOUBLE PRECISION in D-version) work array of length N.

C Local variables
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
C Local parameters 
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
C Intrinsic functions
      INTRINSIC ABS,MIN
C External subroutines and/or functions
c      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD, DLAMCH
c      DOUBLE PRECISION FD05AD, DLAMCH
      EXTERNAL MC64DD,MC64ED,MC64FD, DLAMCH
      DOUBLE PRECISION DLAMCH

C Set RINF to largest positive real number
c XSL  RINF = FD05AD(5)
      RINF = DLAMCH('Overflow')

C Initialization
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
C Scan columns of matrix;
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I 
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
C Update BV with smallest of all the largest maximum absolute values 
C of the rows.
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Rescan unassigned columns; improve initial assignment
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000

C Prepare for main loop
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's 
C algorithm for solving the single source shortest path problem 

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is cost of shortest path to any unassigned row
C ISP is matrix position of unassigned row element in shortest path
C JSP is column index of unassigned row element in shortest path
        CSP = MINONE
C Build shortest path tree starting from unassigned column JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
C Add row I to Q2
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
C Add row I to Q, and push it 
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE

        DO 150 JDUM = 1,NUM
C If Q2 is empty, extract new rows from Q
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
C End of dummy loop; this point is never reached
          ENDIF
C Move row Q0 
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
C Scan column that matches with row Q0
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
C Update D(I) 
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unassigned; update shortest path info
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
C Delete row I from Q (if necessary); add row I to Q2
                IF (DI.NE.MINONE)
     *            CALL MC64FD(L(I),QLEN,N,Q,D,L,1)
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
C Add row I to Q (if necessary); push row I up Q
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
C Update tree
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
       
C If CSP = MINONE, no augmenting path is found
  160   IF (CSP.EQ.MINONE) GO TO 190
C Update bottleneck value
        BV = MIN(BV,CSP)
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
C End of dummy loop; this point is never reached
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE 
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop

C BV is bottleneck value of final matching
      IF (NUM.EQ.N) GO TO 1000

C Matrix is structurally singular, complete IPERM.
C JPERM, PR are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = I
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables N,Q,D,L are described in MC64B/BD
C IF IWAY is equal to 1, then
C node I is pushed from its current position upwards
C IF IWAY is not equal to 1, then
C node I is pushed from its current position downwards

C Local variables and parameters
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI

      DI = D(I)
      POS = L(I)
C POS is index of current position of I in the tree
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END
 
C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64W/WD (IWAY = 2)
C The root node is deleted from the binary heap.

C Local variables and parameters
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI

C Move last element to begin of Q
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
C Exchange old last element with larger priority child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
C Exchange old last element with smaller child
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   20 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)

C Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or
C     MC64WD (IWAY = 2).
C Move last element in the heap 

      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
 
C Quick return, if possible
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF

C Move last element from queue Q to position POS0
C POS is current position of node I in the tree
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
   10   CONTINUE
C End of dummy loop; this point is never reached
   20   Q(POS) = I
        L(I) = POS
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
C End of dummy loop; this point is never reached
      ELSE
        DO 32 IDUM = 1,N
          IF (POS.LE.1) GO TO 34
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
   32   CONTINUE
C End of dummy loop; this point is never reached
   34   Q(POS) = I
        L(I) = POS
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
C End of dummy loop; this point is never reached
      ENDIF
C End of dummy if; this point is never reached
   40 Q(POS) = I
      L(I) = POS

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)

C This subroutine sorts the entries in each column of the 
C sparse matrix (defined by N,NE,IP,IRN,A) by decreasing
C numerical value.

C Local constants
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
C Local variables
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
C Local arrays
      INTEGER TODO(TDLEN)
      
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)

C Sort array roughly with partial quicksort
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
C KEY is the smallest of two values present in interval [FIRST,LAST)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
C Only one value found in interval, so it is already sorted
        TD = TD - 2
        GO TO 425

C Reorder interval [FIRST,LAST) such that entries before MID are gt KEY
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
C Both subintervals [FIRST,MID), [MID,LAST) are nonempty
C Stack the longest of the two subintervals first
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
C          TODO(TD-1) = FIRST
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2

  425   CONTINUE
        IF (TD.EQ.0) GO TO 400 
C There is still work to be done
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
C Next interval is already short enough for straightforward insertion
        TD = TD - 2
        GO TO 425

C Complete sorting with straightforward insertion
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200 
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE

  100 CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)

C N, NE, IP, IRN, are described in MC64A/AD.
C A is a REAL (DOUBLE PRECISION in the D-version) array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that 
C   corresponds to IRN(k). The entries in each column must be 
C   non-negative and ordered by decreasing value. 
C IPERM is an INTEGER array of length N. On exit, it contains the 
C   bottleneck matching: IPERM(I) - 0 or row I is matched to column 
C   IPERM(I).
C NUMX is an INTEGER variable. On exit, it contains the cardinality 
C   of the matching stored in IPERM.
C IW is an INTEGER work array of length 10N.

C FC is an integer array of length N that contains the list of 
C   unmatched columns. 
C LEN(J), LENL(J), LENH(J) are integer arrays of length N that point 
C   to entries in matrix column J.
C   In the matrix defined by the column parts IP(J)+LENL(J) we know 
C   a matching does not exist; in the matrix defined by the column 
C   parts IP(J)+LENH(J) we know one exists.
C   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
C   that is tested for a maximum matching.
C W is an integer array of length N and contains the indices of the 
C   columns for which LENL ne LENH.
C WLEN is number of indices stored in array W.
C IW is integer work array of length N.
C IW4 is integer work array of length 4N used by MC64U/UD.

      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
c      EXTERNAL FD05AD,MC64QD,MC64UD
c      DOUBLE PRECISION FD05AD
      EXTERNAL MC64QD,MC64UD, DLAMCH
      DOUBLE PRECISION DLAMCH

C BMIN and BMAX are such that a maximum matching exists for the input
C   matrix in which all entries smaller than BMIN are dropped. 
C   For BMAX, a maximum matching does not exist.
C BVAL is a value between BMIN and BMAX.
C CNT is the number of calls made to MC64U/UD so far.
C NUM is the cardinality of last matching found.

C Set RINF to largest positive real number
c XSL      RINF = FD05AD(5)
      RINF = DLAMCH('Overflow')

C Compute a first maximum matching from scratch on whole matrix.
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
C The first call to MC64U/UD
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
C Matrix is structurally singular
        BMAX = RINF
      ELSE
C Matrix is structurally nonsingular, NUM=NUMX=N;
C Set BMAX just above the smallest of all the maximum absolute 
C values of the columns
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF

C Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
C Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
C Set LEN(J), LENH(J) just after last entry in column J. 
C Set LENL(J) just after last entry in column J with value ge BMAX.
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
C Column J is empty or all entries are ge BMAX
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
C Add J to W if LENL(J) ne LENH(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE

C Main loop
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
C We have a maximum matching in IW; store IW in IPERM
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
C Keep going round this loop until matching IW is no longer maximum.
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99 
C Set LEN such that all matrix entries with value lt BVAL are 
C discarded. Store old LEN in LENH. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60 
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
C Remove entry from matching
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
C IP(J)+LEN(J)-1 is last entry in column ge BVAL
              LEN(J) = II - IP(J) + 1
C If LENH(J) = LENL(J), remove J from W
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   81     MOD = 1
        ELSE
C We do not have a maximum matching in IW. 
          BMAX = BVAL
C BMIN is the bottleneck value of a maximum matching;
C for BMAX the matching is not maximum, so BMAX>BMIN
C          IF (BMAX .EQ. BMIN) GO TO 99
C Find splitting value BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99 
C Set LEN such that all matrix entries with value ge BVAL are
C inside matrix. Store old LEN in LENL. Do this for all columns W(K).
C Each step, either K is incremented or WLEN is decremented.
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
C End of dummy loop; this point is never reached
C Set mode for next call to MC64U/UD
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))

C IW contains maximum matching of length NUM
   90 CONTINUE 
C End of dummy loop; this point is never reached

C BMIN is bottleneck value of final matching
   99 IF (NUMX.EQ.N) GO TO 1000
C The matrix is structurally singular, complete IPERM
C W, IW are work arrays
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) = J
  320 CONTINUE

 1000 RETURN
      END

C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL

C This routine searches for at most XX different numerical values 
C in the columns W(1:WLEN). XX>=2. 
C Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1 
C until XX values are found or all columns have been considered.
C On output, NVAL is the number of different values that is found 
C and SPLIT(1:NVAL) contains the values in decreasing order.
C If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2). 
C
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA

C Scan columns in W(1:WLEN). For each encountered value, if value not
C already present in SPLIT(1:NVAL), insert value such that SPLIT 
C remains sorted by decreasing value. 
C The sorting is done by straightforward insertion; therefore the use
C of this routine should be avoided for large XX (XX < 20).
      NVAL = 0 
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
C Check presence of HA in SPLIT
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
C The insertion 
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
C Exit loop if XX values are found
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
C Determine VAL
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)

      RETURN
      END
      
C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)

C PR(J) is the previous column to J in the depth first search.
C   Array PR is used as workspace in the sorting algorithm.
C Elements (I,IPERM(I)) I=1,..,N are entries at the end of the
C   algorithm unless N assignments have not been made in which case
C   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
C CV(I) is the most recent loop number (ID+JORD) at which row I
C   was visited.
C ARP(J) is the number of entries in column J which have been scanned 
C   when looking for a cheap assignment.
C OUT(J) is one less than the number of entries in column J which have 
C   not been scanned during one pass through the main loop.
C NUMX is maximum possible size of matching.

      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1

      IF (ID.EQ.1) THEN
C The first call to MC64U/UD.
C Initialize CV and ARP; parameters MOD, NUMX are not accessed
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
C Not the first call to MC64U/UD.
C Re-initialize ARP if entries were deleted since last call to MC64U/UD
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

C NUM0 is size of input matching
C NUM1 is maximum possible size of matching
C NUM2 is maximum allowed number of unassigned rows/columns
C NUM is size of current matching

C Quick return if possible
C      IF (NUM.EQ.N) GO TO 199
C NFC is number of rows/columns that could not be assigned
      NFC = 0
C Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD,
C so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1)*N 

C Main loop. Each pass round this loop either results in a new
C assignment or gives a column with no assignment

      DO 100 JORD = NUM0+1,N

C Each pass uses unique number ID1
        ID1 = ID0 + JORD
C J is unmatched column
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
C Look for a cheap assignment
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C No cheap assignment in row
          ARP(J) = LENC(J)
C Begin looking for assignment chain starting with row J
   30     OUT(J) = LENC(J) - 1
C Inner loop.  Extends chain by one or backtracks
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
C Forward scan
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
C Column J has not yet been accessed during this pass
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
C Backtracking step.
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
C No augmenting path exists for column J.
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
C A matching of maximum size NUM1 is not possible
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
C End of dummy loop; this point is never reached
   70   CONTINUE
C End of dummy loop; this point is never reached

C New assignment is made.
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
C End of dummy loop; this point is never reached

   95   IF (NUM.EQ.NUM1) THEN
C A matching of maximum size NUM1 is found
          LAST = JORD
          GO TO 101
        ENDIF
C
  100 CONTINUE

C All unassigned columns have been considered
      LAST = N

C Now, a transversal is computed or is not possible.
C Complete FC before returning.
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE

C  199 RETURN
      RETURN
      END

C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D) 
      IMPLICIT NONE
C
C *** Copyright (c) 1999  Council for the Central Laboratory of the
C     Research Councils                                             ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC64 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***
C *** Any problems?   Contact ...
C     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   ***
C
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)

C N, NE, IP, IRN are described in MC64A/AD.
C A is a REAL (DOUBLE PRECISION in the D-version) array of length NE.
C   A(K), K=1..NE, must be set to the value of the entry that
C   corresponds to IRN(K). It is not altered.
C   All values A(K) must be non-negative.
C IPERM is an INTEGER array of length N. On exit, it contains the 
C   weighted matching: IPERM(I) = 0 or row I is matched to column 
C   IPERM(I). 
C NUM is an INTEGER variable. On exit, it contains the cardinality of 
C   the matching stored in IPERM.
C IW is an INTEGER work array of length 5N.
C DW is a REAL (DOUBLE PRECISION in the D-version) array of length 2N.
C   On exit, U = D(1:N) contains the dual row variable and 
C   V = D(N+1:2N) contains the dual column variable. If the matrix 
C   is structurally nonsingular (NUM = N), the following holds:
C      U(I)+V(J) <= A(I,J)  if IPERM(I) |= J
C      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
C      U(I) = 0  if IPERM(I) = 0
C      V(J) = 0  if there is no I for which IPERM(I) = J

C Local variables
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
C Local parameters
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
C External subroutines and/or functions
c      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD
c      DOUBLE PRECISION FD05AD
      EXTERNAL MC64DD,MC64ED,MC64FD, DLAMCH
      DOUBLE PRECISION DLAMCH


C Set RINF to largest positive real number
c XSL      RINF = FD05AD(5)
      RINF = DLAMCH('Overflow')

C Initialization
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
C Initialize U(I) 
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
C Row I is not empty
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
C Assignment of column J to row I
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
C Scan unassigned columns; improve assignment
      DO 95 J = 1,N
C JPERM(J) ne 0 iff column J is already assigned
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
C Continue only if column J is not empty
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60 
          JJ = IPERM(I)
C Scan remaining part of assigned column JJ 
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
  
C Prepare for main loop
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE

C Main loop ... each pass round this loop is similar to Dijkstra's 
C algorithm for solving the single source shortest path problem 

      DO 100 JORD = 1,N

        IF (JPERM(JORD).NE.0) GO TO 100
C JORD is next unmatched column
C DMIN is the length of shortest path in the tree
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
C CSP is the cost of the shortest augmenting path to unassigned row
C IRN(ISP). The corresponding column index is JSP.
        CSP = RINF
C Build shortest path tree starting from unassigned column (root) JORD
        J = JORD
        PR(J) = -1

C Scan column J
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
C Initialize heap Q and Q2 with rows held in Q(1:QLEN)
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
C Update tree
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE

        DO 150 JDUM = 1,NUM

C If Q2 is empty, extract rows from Q 
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
C Q0 is row whose distance D(Q0) to the root is smallest
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
C Exit loop if path to Q0 is longer than the shortest augmenting path 
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1

C Scan column that matches with row Q0
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
C DNEW is new cost
            DNEW = VJ + A(K)-U(I)
C Do not update D(I) if DNEW ge cost of shortest path
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
C Row I is unmatched; update shortest path info
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
C Row I is matched; do not update D(I) if DNEW is larger
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                IF (L(I).NE.0) 
     *            CALL MC64FD(L(I),QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE   
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
C Update tree
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
       
C If CSP = RINF, no augmenting path is found
  160   IF (CSP.EQ.RINF) GO TO 190
C Find augmenting path by tracing backward in PR; update IPERM,JPERM
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J) 
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
C End of dummy loop; this point is never reached

C Update U for rows in Q(UP:N)
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE 
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE 
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE

  100 CONTINUE
C End of main loop


C Set dual column variable in D(1:N)
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE

      IF (NUM.EQ.N) GO TO 1100

C The matrix is structurally singular, complete IPERM.
C JPERM, OUT are work arrays
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = J
  320 CONTINUE
 1100 RETURN
      END

