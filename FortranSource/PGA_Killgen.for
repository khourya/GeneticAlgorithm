CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE KILL                                  CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE KILLGEN()
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
      COMMON/PARAMET/PARAM(NPARMAX,NPOPMAX)
      COMMON/NEWGENE/CHILD(NPARMAX,NPOPMAX)
C
      INTEGER IPAR(NPARMAX*NBITMAX)
C
      REAL*8 B2D
C
C**************************************************************************
C
      DO J=1,NPOP-1
       III=0
       DO I=1,NPAR
        DO II=1,NBIT(I)
         III=III+1
         CALL RANDOM(1,R)
         IF (R.LT.0.5) THEN
          IPAR(III)=0
         ELSE
          IPAR(III)=1
         END IF
        END DO
       END DO
       DO I=1,NPAR
        CHILD(I,J)=B2D(IPAR,I)
       END DO    
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
