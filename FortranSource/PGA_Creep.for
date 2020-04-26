CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE CREEP                                 CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CREEP(PAR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
C
      REAL*8 PAR(NPARMAX)
C
C**************************************************************************
C
C     CREEP MUTATION
C
      DO I=1,NPAR
       CALL RANDOM(1,R)
       IF (R.LT.PCMU) THEN
        CALL RANDOM(1,R)
        IF (R.LT.0.5) THEN
         PAR(I)=PAR(I)-PARRES(I)
        ELSE
         PAR(I)=PAR(I)+PARRES(I)
        END IF
        IF (PAR(I).GT.PARMAX(I)) PAR(I)=PARMAX(I)
        IF (PAR(I).LT.PARMIN(I)) PAR(I)=PARMIN(I)
       END IF
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
