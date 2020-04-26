CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVEPOTENTIALFLOW                     CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES POTENTIAL FLOW AS INITIAL SOLUTION                             CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVEPOTENTIALFLOW
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-P.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
C
C----------------------------------------------------------------------
C
      REAL*8 DIFF(NCMAX)
C
      REAL*8 QFU(NCMAX)
      REAL*8 QFV(NCMAX)
C
      REAL*8 AM(NCMAX,NFMAX),BM(NCMAX),XM(NCMAX)
C
C----------------------------------------------------------------------
C
C     LOOP THROUGH SUB-REGIONS
C
      DO K=1,NR
      IF (KR(K).EQ.1) THEN
C
C      INITIALIZE POTENTIAL
C
       DO I=1,NC(K)
        PCC(I,K)=0.D+000
       END DO
C
C      DIFFUSIVITY
C
       DO I=1,NC(K)
        DIFF(I)=VIE(K)/DEE(K)
       END DO
C
C      INPUT AREA FOR PROPORTIONALITY CONSTANT
C
       AINP=0.D+000
       IINP=0
       DO I=1,NB(K)
        IF (KM(I,K).EQ.1) THEN
         IINP=1
         AINP=AINP+AR(I,K)
        END IF
       END DO
       PROP=-AINP*AINP/(8.D+000*VIE(K)) 
C
C----------------------------------------------------------------------
C
C      SOLVE POTENTIAL FLOW IF INLET PRESENT
C
       IF (IINP.EQ.1) THEN
C
C       FORM COEFFICIENT MATRIX
C
        DO I=1,NB(K)
         DO II=1,NCONN(I,K)
          AM(I,II)=GH(I,2,K)*(FXC(I,II,K)*XN(I,K)+FYC(I,II,K)*YN(I,K))
         END DO
         AM(I,1)=AM(I,1)+GH(I,1,K)
        END DO
        DO I=NB(K)+1,NC(K)
         DO II=1,NCONN(I,K)
          AM(I,II)=FXX(I,II,K)+FYY(I,II,K)
         END DO
        END DO
C
C       FORM RHS VECTOR
C
        DO I=1,NB(K)
         IF ((KM(I,K).EQ.2).OR.(KM(I,K).EQ.3)) THEN
          BM(I)=GP(I,3,K)
         ELSE
          BM(I)=(GU(I,3,K)*XN(I,K)+GV(I,3,K)*YN(I,K))/PROP
         END IF
        END DO
        DO I=NB(K)+1,NC(K)
         BM(I)=0.D+000
        END DO
C
C       SOLUTION VECTOR  
C
        DO I=1,NC(K)
         XM(I)=PCC(I,K)
        END DO
C
C       SOLVE FOR POTENTIAL
C
        MXP=100
        IIS=0
        RESP1=1.D+000
        WRITE  (*,*) 'POTENTIAL ITERATION | RESIDUAL'
        DO WHILE ((IIS.LT.40).AND.(RESP1.GT.1.D-012))
         CALL GMRES(AM,BM,XM,RESP1,MXP,K)
         IIS=IIS+1
         WRITE  (*,'(12X,I6,5X,E12.6)') IIS,RESP1
        END DO
C
C       APPROXIMATE DERIVATIVES
C
        DO I=1,NC(K)
         PCC(I,K)=XM(I)
         DXP(I,K)=0.D+000
         DYP(I,K)=0.D+000
         DO II=1,NCONN(I,K)
          DXP(I,K)=DXP(I,K)+FXC(I,II,K)*XM(ICONN(I,II,K))
          DYP(I,K)=DYP(I,K)+FYC(I,II,K)*XM(ICONN(I,II,K))
         END DO
        END DO
C
C----------------------------------------------------------------------
C
C       COMPUTE VELOCITY FROM POTENTIAL
C
        DO I=1,NC(K)
         UCC(I,K)=PROP*DXP(I,K)
         VCC(I,K)=PROP*DYP(I,K)
        END DO
C
C       VELOCITY DERIVATIVES
C
        DO I=1,NC(K)
         D2U(I,K)=0.D+000
         DXU(I,K)=0.D+000
         DYU(I,K)=0.D+000
         D2V(I,K)=0.D+000
         DXV(I,K)=0.D+000
         DYV(I,K)=0.D+000
         DO II=1,NCONN(I,K)
          D2U(I,K)=D2U(I,K)
     &            +(FXX(I,II,K)+FYY(I,II,K))*UCC(ICONN(I,II,K),K)
          DXU(I,K)=DXU(I,K)+FXC(I,II,K)*UCC(ICONN(I,II,K),K)
          DYU(I,K)=DYU(I,K)+FYC(I,II,K)*UCC(ICONN(I,II,K),K)
          D2V(I,K)=D2V(I,K)
     &            +(FXX(I,II,K)+FYY(I,II,K))*VCC(ICONN(I,II,K),K)
          DXV(I,K)=DXV(I,K)+FXC(I,II,K)*VCC(ICONN(I,II,K),K)
          DYV(I,K)=DYV(I,K)+FYC(I,II,K)*VCC(ICONN(I,II,K),K)
         END DO
        END DO
C
C       COMPUTE CONVECTIVE FLUXES
C 
        CALL UPWIND(UCC,VCC,UCC,QFU,DIFF,K)
        CALL UPWIND(UCC,VCC,VCC,QFV,DIFF,K)
C
C----------------------------------------------------------------------
C
C       SOLVE PRESSURE FIELD        
C
        IPRESSURE=0
C
C       FORM COEFFICIENT MATRIX
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
C       FORM RHS VECTOR
C
        DO I=1,NB(K)
         IF (KM(I,K).EQ.2) THEN
          IPRESSURE=1
          BM(I)=GP(I,3,K)
         ELSE
          BM(I)=-DEE(K)*(QFU(I)*XN(I,K)+QFV(I)*YN(I,K))
         END IF
        END DO
        DO I=NB(K)+1,NC(K)
         DQF=0.D+000
         DO II=1,NCONN(I,K)
          DQF=DQF+FXC(I,II,K)*QFU(ICONN(I,II,K))
     &           +FYC(I,II,K)*QFV(ICONN(I,II,K)) 
         END DO
         BM(I)=-DEE(K)*DQF
        END DO
C
C       SOLUTION VECTOR  
C
        DO I=1,NC(K)
         XM(I)=PCC(I,K)
        END DO
C
C       SOLVE FOR POTENTIAL
C
        MXP=100
        IIS=0
        RESP1=1.D+000
        WRITE  (*,*) 'PRESSURE ITERATION  | RESIDUAL'
        DO WHILE ((IIS.LT.25).AND.(RESP1.GT.1.D-010))
         CALL GMRES(AM,BM,XM,RESP1,MXP,K)
         IIS=IIS+1
         WRITE  (*,'(12X,I6,5X,E12.6)') IIS,RESP1
        END DO
C
C       APPROXIMATE DERIVATIVES
C
        DO I=1,NC(K)
         PCC(I,K)=XM(I)
         D2P(I,K)=0.D+000
         DXP(I,K)=0.D+000
         DYP(I,K)=0.D+000
         DO II=1,NCONN(I,K)
          D2P(I,K)=D2P(I,K)
     &            +(FXX(I,II,K)+FYY(I,II,K))*XM(ICONN(I,II,K))
          DXP(I,K)=DXP(I,K)+FXC(I,II,K)*XM(ICONN(I,II,K))
          DYP(I,K)=DYP(I,K)+FYC(I,II,K)*XM(ICONN(I,II,K))
         END DO
        END DO
C
C       FIX PRESSURE IF NOT IMPOSED
C
        IF (IPRESSURE.EQ.0) THEN
         PREF=PCC(1,K)
         DO I=1,NC(K)
          PCC(I,K)=PCC(I,K)-PREF
         END DO
        END IF
C
       END IF
C
C----------------------------------------------------------------------
C
C      RE-CALCULATE DILATATION AND TIME DERIVATIVE
C
       DO I=1,NC(K)
        DCC(I,K)=DXU(I,K)+DYV(I,K)
        DTD(I,K)=0.D+000
       END DO
C
C----------------------------------------------------------------------
C
      END IF
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END