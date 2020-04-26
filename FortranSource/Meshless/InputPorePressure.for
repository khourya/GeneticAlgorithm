CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUTPOREPRESSURE                      CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS PORE PRESSURE CONDITIONS AND SOLUTION                           CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUTPOREPRESSURE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
      INTEGER IOS
C
C----------------------------------------------------------------------
C
C     READ POROUS FLOW BOUNDARY CONDITIONS
C
      OPEN (11,FILE='Model/ALMA_d_por.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING POROUS FLOW BOUNDARY CONDITION FILE.'
       STOP
      END IF
C
      DO K=1,NR
       IF (KR(K).EQ.2) THEN
C
C       READ REGION NUMBER
C
        READ (11,*) KK
C
C       READ INITIAL VALUES
C
        READ (11,*) PINI(K)
C
C       READ BODY FORCES
C
        READ (11,*) QBG(K)
C
C       BOUNDARY CONDITIONS
C       KP=0:IMPERMEABLE, KP=1:PRESSURE, KP=2:VELOCITY, KP=3:RATIO
C
        DO I=1,NB(K)
         READ (11,*) II,KP(I,K),PB,AB,BB
C        IMPERMEABLE
         IF (KP(I,K).EQ.0) THEN
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=DEF(K)*(GX*XN(I,K)+GY*YN(I,K))
         END IF
C        PRESSURE
         IF (KP(I,K).EQ.1) THEN
          GP(I,1,K)=1.D+000
          GP(I,2,K)=0.D+000
          GP(I,3,K)=PB
          PCC(I,K)=PB
         END IF
C        DARCY VELOCITY
         IF (KP(I,K).EQ.2) THEN
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=0.D+000
          IF (PER(K).GT.EPS) GP(I,3,K)=DEF(K)*(GX*XN(I,K)+GY*YN(I,K))
     &                                -VIF(K)*AB/PER(K) 
         END IF
C        RATIO
         IF (KP(I,K).EQ.3) THEN
          GP(I,1,K)=VIF(K)*BB
          GP(I,2,K)=PER(K)
          GP(I,3,K)=VIF(K)*BB*PB
     &             +DEF(K)*PER(K)*(GX*XN(I,K)+GY*YN(I,K)) 
         END IF
C
        END DO
C
       ELSE
C
C       READ DUMMY VALUES
C
        READ (11,*) KDUM
        READ (11,*) PDUM
        READ (11,*) QDUM
        DO I=1,NB(K)
         READ (11,*) IDUM,KDUM,PDUM,ADUM,BDUM
        END DO
C
       END IF
      END DO
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     POROUS FLOW OUTPUT
C
      OPEN (13,FILE='Model/ALMA_o_por.bin',STATUS='OLD',IOSTAT=IOS
     &        ,FORM='UNFORMATTED')
C
      IF (IOS.EQ.0) THEN
C
       DO K=1,NR
        IF (KR(K).EQ.2) THEN
         DO I=1,NC(K)
          READ (13) PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K)
         END DO
        ELSE
         DO I=1,NC(K)
          READ (13) PDUM,DXPDUM,DYPDUM,D2PDUM
         END DO
        END IF  
       END DO
C
      ELSE
C
       DO K=1,NR
        IF (KR(K).EQ.2) THEN
         DO I=1,NC(K)
          PCC(I,K)=PINI(K)
          DXP(I,K)=0.D+000
          DYP(I,K)=0.D+000
          D2P(I,K)=0.D+000
         END DO
        END IF
       END DO
C
      END IF
C
      CLOSE (13)
C
C----------------------------------------------------------------------
C
C     INTERFACE PORE PRESSURE BOUNDARY CONDITIONS
C
      DO K=1,NR
       IF (KR(K).EQ.2) THEN
        PROP=PER(K)/(POR(K)*VIF(K))
C
C       INTERFACE BOUNDARY CONDITIONS
C       KP<-1000:INTERFACE PRESSURE, KP<-2000:INTERFACE DARCY VELOCITY
C
        DO I=1,NB(K)
C        TYPE 1 INTERFACE: IMPOSE PRESSURE
         IF ((KP(I,K).LT.(-1000)).AND.(KP(I,K).GT.(-2000))) THEN
          GP(I,1,K)=1.D+000
          GP(I,2,K)=0.D+000
          GP(I,3,K)=PCC(I,K)
         END IF
C        TYPE 2 INTERFACE: IMPOSE DARCY VELOCITY
         IF ((KP(I,K).LT.(-2000)).AND.(KP(I,K).GT.(-3000))) THEN
          VX=DXP(I,K)-GX*DEF(K)
          VY=DYP(I,K)-GX*DEF(K)
          GP(I,1,K)=0.D+000
          GP(I,2,K)=-PROP
          GP(I,3,K)=-PROP*(VX*XN(I,K)+VY*YN(I,K))
         END IF
C
        END DO
C
       END IF
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END