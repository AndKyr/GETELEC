*DECK DUSRMT
      SUBROUTINE DUSRMT (I, J, AIJ, INDCAT, PRGOPT, DATTRV, IFLAG)
C***BEGIN PROLOGUE  DUSRMT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DSPLP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (USRMAT-S, DUSRMT-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   The user may supply this code
C
C***SEE ALSO  DSPLP
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DUSRMT
      DOUBLE PRECISION PRGOPT(*),DATTRV(*),AIJ
      INTEGER IFLAG(*)
C
C***FIRST EXECUTABLE STATEMENT  DUSRMT
      IF(IFLAG(1).EQ.1) THEN
C
C     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4,
C     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL.
C     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN
C     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA.
           IF(DATTRV(1).EQ.0.D0) THEN
           I = 0
           J = 0
           IFLAG(1) = 3
           ELSE
           IFLAG(2)=-DATTRV(1)
           IFLAG(3)= DATTRV(2)
           IFLAG(4)= 3
           ENDIF
C
           RETURN
      ELSE
           J=IFLAG(2)
           I=IFLAG(3)
           L=IFLAG(4)
           IF(I.EQ.0) THEN
C
C     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED.
                IFLAG(1)=3
                RETURN
           ELSE IF(I.LT.0) THEN
C
C     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN.
                J=-I
                I=DATTRV(L)
                L=L+1
           ENDIF
C
           AIJ=DATTRV(L)
C
C     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY.
           IFLAG(2)=J
           IFLAG(3)=DATTRV(L+1)
           IFLAG(4)=L+2
C
C     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE
C     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED.
           INDCAT=0
           RETURN
      ENDIF
      END
