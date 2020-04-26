CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE LOAD                                  CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LOAD(IDMASTER,ILOAD,FRA,NPRO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
C
      REAL*8 FRA(MAXPROC)
      INTEGER IDMASTER,NPRO
      INTEGER ILOAD(MAXPROC,NPOPMAX)
C
      INTEGER NUSED,NPOSI,NFRAC
      REAL*8 FRACTION
C
C**************************************************************************
C
      IF (NPRO.GT.1) THEN
C
       IF (NPOP.GT.1) THEN
C
        NUSED=0
        NPOSI=NUSED+1
C
C       ASSIGN AN INTEGER PROPORTION OF THE POPULATION TO MASTER PROCESSES
C
        N=IDMASTER+1
C
        NFRAC=INT(FRA(N)*DBLE(NPOP))
        FRACTION=(FRA(N)*DBLE(NPOP)-DBLE(NFRAC))
        NUSED=NUSED+NFRAC
        IF (FRACTION.GE.0.95) THEN
         NUSED=NUSED+1
         FRACTION=FRACTION-1.
        END IF
        IF (NUSED.GT.NPOP) NUSED=NPOP
        DO J=NPOSI,NUSED
         ILOAD(N,J)=1
        END DO
        NPOSI=NUSED+1
C
C       ASSIGN AN INTEGER PROPORTION OF THE POPULATION TO PROCESSES
C
        DO N=1,NPRO
         IF (N.NE.(IDMASTER+1)) THEN
          NFRAC=INT(FRA(N)*DBLE(NPOP))
          FRACTION=FRACTION+(FRA(N)*DBLE(NPOP)-DBLE(NFRAC))
          NUSED=NUSED+NFRAC
          IF (FRACTION.GE.0.5) THEN
           NUSED=NUSED+1
           FRACTION=FRACTION-1.
          END IF
          IF (NUSED.GT.NPOP) NUSED=NPOP
          DO J=NPOSI,NUSED
           ILOAD(N,J)=1
          END DO
          NPOSI=NUSED+1
         END IF
        END DO
C
C       ASSIGNS THE REMAINING (IF ANY) COLUMNS TO LAST OR FIRST PROCESS
C
        IF ((IDMASTER+1).NE.NPRO) THEN
         DO J=NPOSI,NPOP
          ILOAD(NPRO,J)=1
         END DO
        ELSE
         DO J=NPOSI,NPOP
          ILOAD(1,J)=1
         END DO
        END IF
C
       ELSE
C
        ILOAD(1,1)=1
        DO N=2,NPRO
         ILOAD(N,1)=0
        END DO
C
       END IF
C
       ELSE
C
       DO J=1,NPOP
        ILOAD(1,J)=1
       END DO
C
      END IF
C
C**************************************************************************
C
      RETURN
C
      END
