CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE SELECTION                             CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SELECTION(JP1,JP2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/FITNESS/FITNESS(NPOPMAX),PSEL(NPOPMAX),JBEST
C
      INTEGER JP1,JP2
C
C**************************************************************************
C
C     SELECT FIRST PARENT
C
      CALL RANDOM(1,R)
      DO J=NPOP,1,-1
       IF (R.LT.PSEL(J)) JP1=J
      END DO
C
C     SELECT SECOND PARENT
C
      JP2=JP1
      DO WHILE (JP2.EQ.JP1)
       CALL RANDOM(1,R)
       DO J=NPOP,1,-1
        IF (R.LT.PSEL(J)) JP2=J
       END DO
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
