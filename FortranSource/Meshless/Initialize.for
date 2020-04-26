CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INITIALIZE                             CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   INITIALIZES MODEL SETUP                                               CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INITIALIZE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-S.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-E.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
C     LEVEL-SET
C
      IF (MSL.EQ.0) THEN
       DO K=1,NR
        DO I=1,NC(K)
         SCC(I,K)=0.D+000
         DXS(I,K)=0.D+000
         DYS(I,K)=0.D+000
         D2S(I,K)=0.D+000
        END DO
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     MOMENTUM
C
      IF (MSM.EQ.0) THEN
       DO K=1,NR
        IF (KR(K).EQ.1) THEN
         DO I=1,NC(K)
          UCC(I,K)=0.D+000
          DXU(I,K)=0.D+000
          DYU(I,K)=0.D+000
          D2U(I,K)=0.D+000
          VCC(I,K)=0.D+000
          DXV(I,K)=0.D+000
          DYV(I,K)=0.D+000
          D2V(I,K)=0.D+000
          PCC(I,K)=0.D+000
          DXP(I,K)=0.D+000
          DYP(I,K)=0.D+000
          D2P(I,K)=0.D+000
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     ENERGY
C
      IF (MSE.EQ.0) THEN
       DO K=1,NR
        DO I=1,NC(K)
         TCC(I,K)=0.D+000
         DXT(I,K)=0.D+000
         DYT(I,K)=0.D+000
         D2T(I,K)=0.D+000
        END DO
       END DO
      ELSE
       DO K=1,NR
        DO I=1,NC(K)
         TCC(I,K)=TINI(K)
        END DO
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     PORE PRESSURE
C
      IF (MSP.EQ.0) THEN
       DO K=1,NR
        IF (KR(K).EQ.2) THEN
         DO I=1,NC(K)
          PCC(I,K)=0.D+000
          DXP(I,K)=0.D+000
          DYP(I,K)=0.D+000
          D2P(I,K)=0.D+000
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     STRUCTURAL
C
      IF (MSS.EQ.0) THEN
       DO K=1,NR
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         DO I=1,NC(K)
          UCC(I,K)=0.D+000
          DXU(I,K)=0.D+000
          DYU(I,K)=0.D+000
          D2U(I,K)=0.D+000
          VCC(I,K)=0.D+000
          DXV(I,K)=0.D+000
          DYV(I,K)=0.D+000
          D2V(I,K)=0.D+000
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     DILATATION
C
      DO K=1,NR
       DO I=1,NC(K)
        DCC(I,K)=DXU(I,K)+DYV(I,K)
        DTD(I,K)=0.D+000
       END DO
      END DO
C
C----------------------------------------------------------------------
C
C     POROUS FLOW VELOCITY
C
      DO K=1,NR
       IF (KR(K).EQ.2) THEN
        PROP=PER(K)/(POR(K)*VIF(K))
        DO I=1,NC(K)
         IF (MSL.EQ.0) THEN
          FRONT=PROP
         ELSE
          FRONT=SCC(I,K)*PROP+(1.D+000-SCC(I,K))*1.D-003*PROP
         END IF
         VFX(I,K)=-FRONT*(DXP(I,K)-GX*DEF(K))
         VFY(I,K)=-FRONT*(DYP(I,K)-GY*DEF(K))
        END DO
       ELSE
        DO I=1,NC(K)
         VFX(I,K)=0.D+000
         VFY(I,K)=0.D+000
        END DO
       END IF
      END DO
C
C----------------------------------------------------------------------
C
C     ADJUST INTERFACE PORE PRESSURE BOUNDARY CONDITIONS
C
      IF ((MSP.EQ.1).AND.(MSL.EQ.1)) THEN
       DO K=1,NR
        IF (KR(K).EQ.2) THEN
         DO I=1,NB(K)
          IF ((KP(I,K).LT.(-2000)).AND.(KP(I,K).GT.(-3000))) THEN
           GP(I,3,K)=VFX(I,K)*XN(I,K)+VFY(I,K)*YN(I,K)
          END IF
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     ADJUST INTERFACE TRACTION BOUNDARY CONDITIONS
C
      IF (MSS.EQ.1) THEN
       DO K=1,NR
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         ELS=1.D+000-2.D+000*POE(K)
         TEX=2.D+000*VIE(K)*BEE(K)*(1.D+000+POE(K))/ELS
         DO I=1,NB(K)
          IF ((KM(I,K).LT.(-2000)).AND.(KM(I,K).GT.(-3000))) THEN
           GU(I,3,K)=GU(I,3,K)
     &              -(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*XN(I,K)
           GV(I,3,K)=GV(I,3,K)
     &              -(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*YN(I,K)
          END IF  
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     ADJUST STRUCTURAL BOUNDARY CONDITIONS WITH PORE PRESSURE
C
      IF ((MSS.EQ.1).AND.(MSP.EQ.1)) THEN
       DO K=1,NR
        IF (KR(K).EQ.2) THEN
         DO I=1,NB(K)
          IF ((KM(I,K).EQ.0).AND.(KP(I,K).EQ.1)) THEN
           GU(I,1,K)=0.D+000
           GU(I,2,K)=1.D+000
           GU(I,3,K)=GU(I,3,K)+GP(I,3,K)*XN(I,K)
           GV(I,1,K)=0.D+000
           GV(I,2,K)=1.D+000
           GV(I,3,K)=GV(I,3,K)+GP(I,3,K)*YN(I,K)
          END IF  
         END DO
        END IF
       END DO
      END IF
C
C----------------------------------------------------------------------
C
C     INITIALIZING VELOCITY AND PRESSURE FIELD
C
      IF ((ITP.EQ.0).AND.(MSM.EQ.1)) THEN
       CALL SOLVEPOTENTIALFLOW
      END IF
C
C----------------------------------------------------------------------
C
C     INITIALIZING PORE PRESSURE FIELD
C
      IF ((ITP.EQ.0).AND.(MSP.EQ.1).AND.(MSL.EQ.1)) THEN
       CALL SOLVEPOREPRESSURE
      END IF
C
C----------------------------------------------------------------------
C
C     INITIALIZING STRUCTURAL FIELD
C
      IF ((ITP.EQ.0).AND.(MSS.EQ.1)) THEN
       CALL SOLVESTRUCTURAL
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END