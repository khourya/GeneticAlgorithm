CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                     SUBROUTINE FUNCTIONOUTPUT                       CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FUNCTIONOUTPUT()
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
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
C
C**************************************************************************
C
C     OUTPUT THE FUNCTION VARIABLES
C
C
C     INITIALIZE TIME STEPPING 
C
      ITP=0
      CTIME=0.D+000
C
C     INITIALIZE MODEL SETUP
C
      CALL INITIALIZE
C
C     LOOP OVER MAXIMUM ITERATIONS
C
      DO WHILE (ITP.LT.MAXITER)
C
C      INCREASE ITERATION AND ELAPSED TIME
C
       ITP=ITP+1
       CTIME=CTIME+DT
C
C      SOLVE LEVEL-SET FIELD
C
       IF (MSL.EQ.1) CALL SOLVELEVELSET
C
C      SOLVE MOMENTUM FIELD
C
       IF (MSM.EQ.1) CALL SOLVEMOMENTUM
C
C      SOLVE ENERGY FIELD
C
       IF (MSE.EQ.1) CALL SOLVEENERGY
C
C      SOLVE PORE PRESSURE FIELD
C
       IF (MSP.EQ.1) CALL SOLVEPOREPRESSURE
C
C      SOLVE STRUCTURAL FIELD
C
       IF (MSS.EQ.1) CALL SOLVESTRUCTURAL
C
      END DO
C
C     OUTPUT FIELD
C     
      CALL MESHLESS_OUTPUT
C
C     WRITE HOLE DESCRIPTION TO FILE
C
      OPEN (16,FILE='Model\ALMA_a_Hole.out',STATUS='UNKNOWN')
	WRITE (16,'(4(2X,E12.6))') XCH,YCH,RXH,RYH
	CLOSE (16)
C
C
C**************************************************************************
C
      RETURN
C
      END
