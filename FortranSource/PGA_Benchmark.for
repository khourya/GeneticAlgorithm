CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                    SUBROUTINE BENCHMARK                             CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BENCHMARK(TIM,FRA,ID,IDMASTER,NPRO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
      REAL*8 STARTTIME,ENDTIME
      REAL*8 TIMESUMI
      REAL*8 TIM(MAXPROC),FRA(MAXPROC)
      INTEGER ID,IDMASTER,NPRO,IERR
C
      REAL*8 E,EPS
      INTEGER NS
C
C**************************************************************************
C
      STARTTIME=MPI_WTIME()
C
C     BENCHMARKING LOOP
C
      NS=2000
      E=0.
      DO I=1,NS
       DO J=1,NS
        K=NS*(I-1)+J
        EPS=(1./DBLE(K))**2.
        E=E+EPS
       END DO
      END DO
C
      ENDTIME=MPI_WTIME()
C
C     BENCHMARKING TIMINGS
C
      TIM(ID+1)=ENDTIME-STARTTIME
C
C     BROADCAST ALL TIMINGS
C
      DO N=1,NPRO
       CALL MPI_BCAST(TIM(N),1,MPI_DOUBLE_PRECISION,N-1,
     &                MPI_COMM_WORLD,IERR)
      END DO
C
C     COMPUTE POWER PROPORTION
C
      IF (ID.EQ.IDMASTER) THEN
       TIMESUMI=0.
       DO N=1,NPRO
        TIMESUMI=TIMESUMI+(1./TIM(N))
       END DO
       DO N=1,NPRO
        FRA(N)=(1./TIM(N))/TIMESUMI
       END DO
      END IF
C
C**************************************************************************
C
      RETURN
C
      END
