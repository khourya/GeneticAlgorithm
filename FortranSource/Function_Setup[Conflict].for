CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                     SUBROUTINE FUNCTIONINPUT                        CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCTIONSETUP(ID,IDMASTER,NPRO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
C     USE THE FOLLOWING SPACE TO INCLUDE THE FILES FOR THE FUNCTION INPUT
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
C
      INTEGER ID,IDMASTER,NPRO,IERR
C
      CHARACTER*120 TITLE
      INTEGER IOS
C
C**************************************************************************
C
C     INPUT AND PROCESS THE NECESSARY DATA FOR FUNCTION EVALUATION
C     INCLUDE THE INPUT DATA AND PROBLEM SETUP ROUTINES
C     USE FUNCTION INPUT DATA TO STABLISH NUMBER OF PARAMETERS TO OPTIMIZE
C     AND THE MINIMUM AND MAXIMUM ALLOWABLE VALUES
C
      IF (ID.EQ.IDMASTER) THEN
C
       NPAR=4
C
C      MEASUREMENT DATA
C
       OPEN (20,FILE='Model/ALMA_a_measure.txt',STATUS='OLD',IOSTAT=IOS)
C
       IF (IOS.EQ.0) THEN       
        READ (12,*) TITLE
C
C       NUMBER OF MEASURMENT NODES AND TIMES
C
	  READ (20,*) NMN,NMT
C
C       XMIN,YMIN,XMAX,YMAX
C
        READ (20,*) XMIN,YMIN,XMAX,YMAX
C
C       X-HOLE LOCATION
C
        PARMIN(1)=XMIN
        PARMAX(1)=XMAX
        NBIT(1)=8
C
C       Y-HOLE LOCATION
C
        PARMIN(2)=YMIN
        PARMAX(2)=YMAX
        NBIT(2)=8
C
C       HOLE X-RADIUS 
C
        PARMIN(3)=(XMAX-XMIN)/1.D+002
        PARMAX(3)=(XMAX-XMIN)/1.D+001
        NBIT(3)=8
C
C       HOLE Y-RADIUS 
C
        PARMIN(4)=(YMAX-YMIN)/1.D+002
        PARMAX(4)=(YMAX-YMIN)/1.D+001
        NBIT(4)=8
C
C       MEASUREMENT REGION AND BOUNDARY NODE
C
        DO I=1,NMN
	   READ (20,*) MR(I),MN(I)
	  END DO
C
        DO NT=1,NMT
C
C        TIME-STEP
C
         READ (20,*) MT(NT)
C
C        MEASURED TEMPERATURE
C
         DO I=1,NMN
	    READ (20,*) TM(I,NT)
	   END DO
C
        END DO
C
       ELSE
C
        WRITE (*,*) '...MISSING INFORMATION FILE.'
        STOP
C
       END IF
C
       CLOSE (20)
C
C      INPUT PARAMETERS, GEOMETRY, CONNECTIVITY AND INTERPOLATION
C
       WRITE (*,*) 'READING PROBLEM DATA............................'
       CALL MESHLESS_INPUT
       WRITE (*,*) '............................................DONE'
C
C      VALIDATING THERMOPHYSICAL QUANTITIES
C
       WRITE (*,*) 'VALIDATING THERMOPHYISICAL QUANTITIES...........'
       CALL VALIDATE
       WRITE (*,*) '............................................DONE'
C
C      INPUT LEVEL-SET
C
       IF (MSL.EQ.1) THEN
        WRITE (*,*) 'READING LEVEL-SET DATA..........................'
        CALL INPUTLEVELSET
        WRITE (*,*) '............................................DONE'
       END IF
C
C      INPUT MOMENTUM
C
       IF (MSM.EQ.1) THEN
        WRITE (*,*) 'READING MOMENTUM DATA...........................'
        CALL INPUTMOMENTUM
        WRITE (*,*) '............................................DONE'
       END IF
C
C      INPUT ENERGY
C
       IF (MSE.EQ.1) THEN
        WRITE (*,*) 'READING ENERGY DATA.............................'
        CALL INPUTENERGY
        WRITE (*,*) '............................................DONE'
       END IF
C
C      INPUT PORE PRESSURE
C
       IF (MSP.EQ.1) THEN
        WRITE (*,*) 'READING PORE PRESSURE DATA......................'
        CALL INPUTPOREPRESSURE
        WRITE (*,*) '............................................DONE'
       END IF
C
C      INPUT STRUCTURAL
C
       IF (MSS.EQ.1) THEN
        WRITE (*,*) 'READING STRUCTURAL DATA.........................'
        CALL INPUTSTRUCTURAL
        WRITE (*,*) '............................................DONE'
       END IF
C
      END IF
C
C     BROADCAST MESHLESS DATA OVER CLUSTER
C
C      CALL MPI_BCAST(XXX,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C      ...
C      ...
C      ...
C      ...
C
C
C**************************************************************************
C
C     COMPUTE PARAMETER RESOLUTION AND BROADCAST OVER CLUSTER
C
      IF (ID.EQ.IDMASTER) THEN
       DO I=1,NPAR
        PARRES(I)=(PARMAX(I)-PARMIN(I))/(2.**NBIT(I)-1.)
       END DO
      END IF
C
      CALL MPI_BCAST(NPAR,1,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(NBIT,NPAR,MPI_INTEGER,IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(PARMIN,NPAR,MPI_DOUBLE_PRECISION,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(PARMAX,NPAR,MPI_DOUBLE_PRECISION,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(PARRES,NPAR,MPI_DOUBLE_PRECISION,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
C**************************************************************************
C
      RETURN
C
      END
