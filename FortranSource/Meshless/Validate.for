CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE VALIDATE                               CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   VALIDATE TRANSIENT PROBLEM STABILITY USING FOURIER #'S                CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VALIDATE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Topology.for'
C
C----------------------------------------------------------------------
C
C     VALIDATING DIFFUSION PARAMETERS
C
      DO K=1,NR
C
       RRMIN=1.D+020
       DO I=NB(K)+1,NC(K)
        RR=RXA(I,K)*RXA(I,K)+RYA(I,K)*RYA(I,K)
        IF (RR.LT.RRMIN) RRMIN=RR
       END DO
C
C      CHECKING MOMENTUM FOURIER NUMBER (DIFF*DT/DX^2 < 0.25)
C
       IF ((MSM.EQ.1).AND.(KR(K).EQ.1)) THEN
        DIFF=VIE(K)/DEE(K)
        DTMAX=0.25D+000*RRMIN/DIFF
        IF (DT.GT.DTMAX) THEN
         WRITE (*,*)
     &   '...WARNING!!! MOMENTUM STABILITY LIMIT EXCEEDED'
         WRITE (*,*)
     &   '...DECREASE THE TIME STEP TO AT MOST ',DTMAX
         READ (*,*)
         STOP
        END IF
       END IF
C
C      CHECKING ENERGY FOURIER NUMBER (DIFF*DT/DX^2 < 0.5)
C
       IF (MSE.EQ.1) THEN
        DIFF=TCE(K)/(DEE(K)*SHE(K))
        DTMAX=0.25D+000*RRMIN/DIFF
        IF (DT.GT.DTMAX) THEN
         WRITE (*,*)
     &   '...WARNING!!! ENERGY STABILITY LIMIT EXCEEDED'
         WRITE (*,*)
     &   '...DECREASE THE TIME STEP TO AT MOST ',DTMAX
         READ (*,*)
         STOP
        END IF
       END IF
C
C      CHECKING PORE PRESSURE FOURIER NUMBER (DIFF*DT/DX^2 < 0.5)
C
       IF ((MSP.EQ.1).AND.(MSL.EQ.0).AND.(KR(K).EQ.2)) THEN
        DIFF=PER(K)/(VIF(K)*POR(K)*COF(K))
        DTMAX=0.25D+000*RRMIN/DIFF
        IF (DT.GT.DTMAX) THEN
         WRITE (*,*)
     &   '...WARNING!!! PORE PRESSURE STABILITY LIMIT EXCEEDED'
         WRITE (*,*)
     &   '...DECREASE THE TIME STEP TO AT MOST ',DTMAX
         READ (*,*)
         STOP
        END IF
       END IF
C
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END