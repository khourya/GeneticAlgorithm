CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                     SUBROUTINE FUNCTIONEVALUATE                     CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCTIONEVALUATE(ID,IDMASTER,NPRO,ILOAD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
C     USE THE FOLLOWING SPACE TO INCLUDE THE FILES FOR THE FUNCTION EVALUATION
C
      INCLUDE '../Include/Parameters.for'
      INCLUDE '../Include/Information.for'
      INCLUDE '../Include/Geometry.for'
      INCLUDE '../Include/Materials.for'
      INCLUDE '../Include/Field-S.for'
      INCLUDE '../Include/Field-M.for'
      INCLUDE '../Include/Field-E.for'
      INCLUDE '../Include/Field-P.for'
      INCLUDE '../Include/Topology.for'
      INCLUDE '../Include/Interpolation-C.for'
      INCLUDE '../Include/Interpolation-U.for'
      INCLUDE '../Include/Triangulation.for'
C
C     USE THE FOLLOWING SPACE TO INCLUDE THE COMMON BLOCKS FOR FUNCTION VARIABLES
C
      PARAMETER (NTMAX=20,NMMAX=200)
	COMMON/MEASURE/NMN,NMT,MT(NTMAX),MR(NMMAX),MN(NMMAX)
      COMMON/MEASURT/TM(NMMAX,NTMAX)
      COMMON/HOLELOC/XCH,YCH,RXH,RYH
C
C     STANDARD COMMON BLOCKS AND VARIABLES
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
      COMMON/PARAMET/PARAM(NPARMAX,NPOPMAX)
      COMMON/FITNESS/FITNESS(NPOPMAX),PSEL(NPOPMAX),JBEST
C
      INTEGER ID,IDMASTER,NPRO,IERR
      INTEGER ILOAD(MAXPROC,NPOPMAX)
C
      REAL*8 FIT
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER ISENDER,IND
C
      INTEGER IDUMMY
C
C**************************************************************************
C
      DO J=1,NPOP
       IF (ILOAD(ID+1,J).EQ.1) THEN
C
        FITNESS(J)=0.D+000
C
C       TRANSLATE PARAMETERS INTO SOURCE CHARACTERISTICS
C
        XCH=PARAM(1,J)
        YCH=PARAM(2,J)
        RXH=PARAM(3,J)
        RYH=PARAM(4,J)
C
C       INITIALIZE TIME STEPPING AND MEASUREMENT TIME
C
        ITP=0
        CTIME=0.D+000
        NT=1
C
C       INITIALIZE MODEL SETUP
C 
        CALL INITIALIZE
C
C       LOOP OVER MAXIMUM ITERATIONS
C
        DO WHILE (ITP.LT.MAXITER)
C
C        INCREASE ITERATION AND ELAPSED TIME
C
         ITP=ITP+1
         CTIME=CTIME+DT
C
C        SOLVE LEVEL-SET FIELD
C
         IF (MSL.EQ.1) CALL SOLVELEVELSET
C
C        SOLVE MOMENTUM FIELD
C
         IF (MSM.EQ.1) CALL SOLVEMOMENTUM
C
C        SOLVE ENERGY FIELD
C
         IF (MSE.EQ.1) CALL SOLVEENERGY
C
C        SOLVE PORE PRESSURE FIELD
C
         IF (MSP.EQ.1) CALL SOLVEPOREPRESSURE
C
C        SOLVE STRUCTURAL FIELD
C
         IF (MSS.EQ.1) CALL SOLVESTRUCTURAL
C
C        EVALUATE THE FITNESS OF THE Jth INDIVIDUAL GIVEN ITS PARAMETERS
C     
         IF (ITP.EQ.MT(NT)) THEN
C
C         CALCULATE RMS
C
          DO I=1,NMN
           FITNESS(J)=FITNESS(J)+(TM(I,NT)-TCC(MN(I),MR(I)))**2.D+000
CCC
CCC           WRITE (*,*) ITP,TM(I,NT),TCC(MN(I),MR(I))
CCC
          END DO
C
          NT=NT+1
C
         END IF
C
        END DO
C
C       INVERT RMS TO CALCULATE FITNESS
C
        IF (FITNESS(J).GT.EPS) THEN 
         FITNESS(J)=DSQRT(DBLE(NMN*NMT)/FITNESS(J))
        ELSE
         FITNESS(J)=1.D+020
        END IF     
CCCCC
        WRITE (*,*) "INDIVIDUAL: ",J
        WRITE (*,*) XCH,YCH
        WRITE (*,*) RXH,RYH
        WRITE (*,*) "FITNESS: ",FITNESS(J)
CCCCC
       END IF
C
      END DO
C
C**************************************************************************
C
C     SEND FITNESS TO MASTER COMPUTER
C
      CALL MPI_BCAST(IDUMMY,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
      DO J=1,NPOP
       IF (ID.NE.IDMASTER) THEN
        IF (ILOAD(ID+1,J).EQ.1) THEN
         FIT=FITNESS(J)
         CALL MPI_SEND(FIT,1,MPI_DOUBLE_PRECISION,IDMASTER,J,
     @                 MPI_COMM_WORLD,IERR)    
        END IF
       ELSE
        IF (ILOAD(ID+1,J).EQ.0) THEN
         CALL MPI_RECV(FIT,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     @                 MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
         ISENDER=STATUS(MPI_SOURCE)
         IND=STATUS(MPI_TAG)
         FITNESS(IND)=FIT 
        END IF
       END IF
      END DO
C
      CALL MPI_BCAST(IDUMMY,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
C     ACUMMULATE PROBABILITY OF SELECTION
C
      IF (ID.EQ.IDMASTER) THEN
       FITMIN=FITNESS(1)
       DO J=2,NPOP
        IF (FITNESS(J).LT.FITMIN) FITMIN=FITNESS(J)
       END DO
       IF (FITMIN.GE.0.) FITMIN=0.
       FIT=0.
       DO J=1,NPOP
        FIT=FIT+(FITNESS(J)-FITMIN)
       END DO
       DO J=1,NPOP
        PSEL(J)=(FITNESS(J)-FITMIN)/FIT
       END DO
       DO J=2,NPOP
        PSEL(J)=PSEL(J)+PSEL(J-1)
       END DO
      END IF
C
C     SELECT BEST FITNESS
C
      IF (ID.EQ.IDMASTER) THEN
       FIT=FITNESS(1)
       JBEST=1
       DO J=2,NPOP
        IF (FITNESS(J).GT.FIT) THEN
         FIT=FITNESS(J)
         JBEST=J
        END IF
       END DO
      END IF
C
C     BROADCAST BEST INDIVIDUAL
C
      CALL MPI_BCAST(JBEST,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
C**************************************************************************
C
C     ASSIGN THE BEST PARAMETERS TO FUNCTION EVALUATION VARIABLE
C
      XCH=PARAM(1,JBEST)
      YCH=PARAM(2,JBEST)
      RXH=PARAM(3,JBEST)
      RYH=PARAM(4,JBEST)
C
C**************************************************************************
C
      RETURN
C
      END
