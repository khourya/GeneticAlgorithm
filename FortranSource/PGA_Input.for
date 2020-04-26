CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                         SUBROUTINE INPUT                            CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUT(ID,IDMASTER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
C
      INTEGER ID,IDMASTER,IERR
C
C**************************************************************************
C
      IF (ID.EQ.IDMASTER) THEN
C
       OPEN (14,FILE='model/PGA.inp',STATUS='OLD')
C
C      MAXIMUM NUMBER OF GENERATIONS
C
       READ (14,*) NGEN
C
C      POPULATION SIZE
C
       READ (14,*) NPOP
C
C      PROBABILITY OF JUMP MUTATION
C     
       READ (14,*) PJMU
C
C      PROBABILITY OF CREEP MUTATION
C
       READ (14,*) PCMU
C
       CLOSE (14)
C
      END IF
C
C     BROADCAST GENETIC ALGORITHM DATA OVER CLUSTER
C
      CALL MPI_BCAST(NGEN,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(NPOP,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(PJMU,1,MPI_DOUBLE_PRECISION,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(PCMU,1,MPI_DOUBLE_PRECISION,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
C**************************************************************************
C
      RETURN
C
      END
