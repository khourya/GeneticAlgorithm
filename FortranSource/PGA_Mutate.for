CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE MUTATE                                CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MUTATE(ICHI)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
C
      INTEGER ICHI(NPARMAX*NBITMAX)
C
C**************************************************************************
C
C     JUMP MUTATION
C
      CALL RANDOM(1,R)
      IF (R.LT.PJMU) THEN
       III=0
       DO I=1,NPAR
        DO II=1,NBIT(I)
         III=III+1
         IF (ICHI(III).EQ.0) THEN
          ICHI(III)=1
         ELSE
          ICHI(III)=0
         END IF
        END DO
       END DO
      END IF   
C
C**************************************************************************
C
      RETURN
C
      END
