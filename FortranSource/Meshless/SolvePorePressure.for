CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVEPOREPRESSURE                      CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES PORE PRESSURE TRANSPORT EQUATION                               CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVEPOREPRESSURE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-S.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-P.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
C
C----------------------------------------------------------------------
C
      REAL*8 DIFF(NCMAX)
      REAL*8 RHS(NCMAX)
      REAL*8 QFP(NCMAX)
C
      REAL*8 AM(NCMAX,NFMAX),BM(NCMAX),XM(NCMAX)
C
C----------------------------------------------------------------------
C
C     INITIALIZE RESIDUAL
C
      RESP=0.D+000
C
C----------------------------------------------------------------------
C
C     SOLVE PRESSURE AS TRANSIENT PROBLEM IF NOT SOLVING FOR LEVEL-SET (MSL=0) 
C
      IF (MSL.EQ.0) THEN
      DO K=1,NR
      IF (KR(K).EQ.2) THEN
C
C      DIFFUSIVITY
C
       DO I=1,NC(K)
        DIFF(I)=PER(K)/(VIF(K)*POR(K)*COF(K))
       END DO
C
C      BOUNDARY CONDITIONS
C
       DO I=1,NB(K)
        RHS(I)=GP(I,3,K)
       END DO
C
C      BODY FORCE
C
       DO I=NB(K)+1,NC(K)
C        RHS(I)=(QBG(K)-DTD(I,K))/(POR(K)*COF(K))
        RHS(I)=QBG(K)/(POR(K)*COF(K))
       END DO
C
C      CONVECTIVE FLUX
C
       DO I=1,NC(K)
        QFP(I)=0.D+000
       END DO
C
C      TRANSPORT PORE PRESSURE
C
       CALL TRANSPORT(PCC,D2P,DXP,DYP,QFP,DIFF,RHS,GP,DT,RESPP,K)
C
C      POROUS FLOW VELOCITY
C
       PROP=PER(K)/(POR(K)*VIF(K))
       DO I=1,NC(K)
        VFX(I,K)=-PROP*(DXP(I,K)-GX*DEF(K))
        VFY(I,K)=-PROP*(DYP(I,K)-GY*DEF(K))
       END DO
C
C      ACCUMULATE PRESSURE RESIDUAL
C
       RESP=RESP+RESPP
C
      END IF
      END DO
      END IF
C
C----------------------------------------------------------------------
C
C     SOLVE PRESSURE AS STEADY PROBLEM IF SOLVING FOR LEVEL-SET (MSL=1) 
C
      IF (MSL.EQ.1) THEN
      DO K=1,NR
      IF (KR(K).EQ.2) THEN
C
C      FORM COEFFICIENT MATRIX
C
       DO I=1,NB(K)
        DO II=1,NCONN(I,K)
         AM(I,II)=GP(I,2,K)*(FXC(I,II,K)*XN(I,K)+FYC(I,II,K)*YN(I,K))
        END DO
        AM(I,1)=AM(I,1)+GP(I,1,K)
       END DO
       DO I=NB(K)+1,NC(K)
        DO II=1,NCONN(I,K)
         AM(I,II)=FXX(I,II,K)+FYY(I,II,K)
        END DO
       END DO
C
C      FORM RHS VECTOR
C
       DO I=1,NB(K)
        BM(I)=GP(I,3,K)
       END DO
C
       DO I=NB(K)+1,NC(K)
C        BM(I)=VIF(K)*(DTD(I,K)-QBG(K))/PER(K)
        BM(I)=-VIF(K)*QBG(K)/PER(K)
       END DO
C
C      SOLUTION VECTOR  
C
       DO I=1,NC(K)
        XM(I)=PCC(I,K)
       END DO
C
C      SOLVE FOR PORE PRESSURE
C
       MXP=10
       DO IIS=1,2
        CALL GMRES(AM,BM,XM,RESPP,MXP,K)
       END DO
C
       IF (ITP.EQ.0) THEN
        IIS=0
        MXP=100
        RESPP=1.D+000
        WRITE  (*,*) 'PRESSURE ITERATION  | RESIDUAL'
        DO WHILE ((IIS.LT.25).AND.(RESPP.GT.1.D-004))
         CALL GMRES(AM,BM,XM,RESPP,MXP,K)
         IIS=IIS+1
         WRITE  (*,'(12X,I6,5X,E12.6)') IIS,RESPP
        END DO
       END IF
C
C      APPROXIMATE DERIVATIVES
C
       DO I=1,NC(K)
        PCC(I,K)=XM(I)
        D2P(I,K)=0.D+000
        DXP(I,K)=0.D+000
        DYP(I,K)=0.D+000
        DO II=1,NCONN(I,K)
         D2P(I,K)=D2P(I,K)
     &           +(FXX(I,II,K)+FYY(I,II,K))*XM(ICONN(I,II,K))
         DXP(I,K)=DXP(I,K)+FXC(I,II,K)*XM(ICONN(I,II,K))
         DYP(I,K)=DYP(I,K)+FYC(I,II,K)*XM(ICONN(I,II,K))
        END DO
       END DO
C
C      POROUS FLOW VELOCITY
C
       PROP=PER(K)/(POR(K)*VIF(K))
       DO I=1,NC(K)
        FRONT=SCC(I,K)*PROP+(1.D+000-SCC(I,K))*1.D-003*PROP
        VFX(I,K)=-FRONT*(DXP(I,K)-GX*DEF(K))
        VFY(I,K)=-FRONT*(DYP(I,K)-GY*DEF(K))
       END DO
C
C      ACCUMULATE PRESSURE RESIDUAL
C
       RESP=RESP+RESPP
C
      END IF
      END DO
      END IF
C
C----------------------------------------------------------------------
C
C      AVERAGE INTERFACE PORE PRESSURE
C
      IF (NR.GT.1) THEN
       CALL INTERFACEPOREPRESSURE
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END