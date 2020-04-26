CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE REPRODUCT                             CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REPRODUCT(JP1,JP2,IC)
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
      INTEGER JP1,JP2
C
      REAL*8 PAR(NPARMAX)
      INTEGER IPAR1(NPARMAX*NBITMAX),IPAR2(NPARMAX*NBITMAX)
      INTEGER ICHI(NPARMAX*NBITMAX),IC
C
      REAL*8 B2D
C
C**************************************************************************
C
C     ENCODE FIRST PARENT
C
      DO I=1,NPAR
       PAR(I)=PARAM(I,JP1)
      END DO
      CALL D2B(PAR,IPAR1)
C
C     ENCODE SECOND PARENT
C
      DO I=1,NPAR
       PAR(I)=PARAM(I,JP2)
      END DO
      CALL D2B(PAR,IPAR2)
C
C     GENERATE CHILD BY UNIFORM CROSSOVER
C
      III=0
      DO I=1,NPAR
       DO II=1,NBIT(I)
        III=III+1
        CALL RANDOM(1,R)
        IF (R.LT.0.5) THEN
         ICHI(III)=IPAR1(III)
        ELSE
         ICHI(III)=IPAR2(III)
        END IF
       END DO
      END DO
C
C     PERFORM JUMP MUTATION IN NEWBORN CHILD
C
      CALL MUTATE(ICHI)      
C
C     DECODE CHILD
C
      DO I=1,NPAR
       PAR(I)=B2D(ICHI,I)
      END DO
C
C     PERFORM CREEP MUTATION IN NEWBORN CHILD
C
      CALL CREEP(PAR)
C
C     ASSIGN CHILD
C
      DO I=1,NPAR
       CHILD(I,IC)=PAR(I)
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
