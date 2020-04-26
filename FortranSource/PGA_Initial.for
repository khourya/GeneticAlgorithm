CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE INITIAL                               CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INITIAL(ID,IDMASTER,IGEN,UNDAT,UNOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
      COMMON/PARAMET/PARAM(NPARMAX,NPOPMAX)
C
      INTEGER ID,IDMASTER,IERR
C
      REAL*8 PAR(NPARMAX)
      INTEGER IPAR(NPARMAX*NBITMAX)
      INTEGER IGEN
      INTEGER IO1,IO2,IO3
      INTEGER UNDAT,UNOUT
C
      REAL*8 B2D
C
C**************************************************************************
C
      IF (ID.EQ.IDMASTER) THEN
       OPEN (UNOUT,FILE='Model/PGA.out')
       OPEN (UNDAT,FILE='Model/PGA.dat',STATUS='OLD',IOSTAT=IO1)
C
       IF (IO1.EQ.0) THEN
C
        READ (UNDAT,*) IGEN
        DO J=1,NPOP
         READ (UNDAT,*,IOSTAT=IO2) JJ,(PARAM(I,J),I=1,NPAR)
        END DO
        CLOSE (UNDAT)
C
        DO J=1,NPOP
         DO I=1,NPAR
          IF (PARAM(I,J).GT.PARMAX(I)) PARAM(I,J)=PARMAX(I)
          IF (PARAM(I,J).LT.PARMIN(I)) PARAM(I,J)=PARMIN(I)
         END DO
        END DO
C
        DO I=1,IGEN
         READ (UNOUT,*,IOSTAT=IO3) II,FIT,(PAR(J),J=1,NPAR)
        END DO
C
       ELSE
C
        IGEN=0
        DO J=1,NPOP
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
          PARAM(I,J)=B2D(IPAR,I)
         END DO      
        END DO
C
       END IF
C
      END IF
C
C     BROADCAST PARAMETERS OVER CLUSTER
C
      DO J=1,NPOP
       DO I=1,NPAR
        PAR(I)=PARAM(I,J)
       END DO
       CALL MPI_BCAST(PAR,NPAR,MPI_DOUBLE_PRECISION,
     &                IDMASTER,MPI_COMM_WORLD,IERR)
       DO I=1,NPAR
        PARAM(I,J)=PAR(I)
       END DO
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
