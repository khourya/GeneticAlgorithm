CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVEMOMENTUM                          CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES MOMENTUM TRANSPORT EQUATIONS                                   CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVEMOMENTUM
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
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
C
C----------------------------------------------------------------------
C
      REAL*8 PMIN,PMAX,RESP
C
      REAL*8 DENS(NCMAX)
      REAL*8 DIFF(NCMAX)
C
      REAL*8 RHX(NCMAX),RHY(NCMAX)
      REAL*8 QFU(NCMAX),QFV(NCMAX)
C
      REAL*8 BCP(NBMAX)
C
      REAL*8 AM(NCMAX,NFMAX),BM(NCMAX),XM(NCMAX)
C
      REAL*8 UCO(NCMAX),VCO(NCMAX),PCO(NCMAX),DCO(NCMAX)
C
C----------------------------------------------------------------------
C
C     INITIALIZE RESIDUALS
C
      RESC=0.D+000
      RESX=0.D+000
      RESY=0.D+000
C
C----------------------------------------------------------------------
C
C     LOOP OVER SUBREGIONS FOR FIELD SOLUTION
C
      DO K=1,NR
      IF (KR(K).EQ.1) THEN
C
C----------------------------------------------------------------------
C
C      STORE CURRENT FIELD
C
       DO I=1,NC(K)
        UCO(I)=UCC(I,K)
        VCO(I)=VCC(I,K)
        PCO(I)=PCC(I,K)
        DCO(I)=DXU(I,K)+DYV(I,K)
       END DO
C
C----------------------------------------------------------------------
C
C      DENSITY
C
       DO I=1,NC(K)
        DENS(I)=DEE(K)
       END DO
C
C      DIFFUSIVITY
C
       DO I=1,NC(K)
        DIFF(I)=VIE(K)/DENS(I)
       END DO
C
C      COMPUTE CONVECTIVE FLUXES
C
       CALL UPWIND(UCC,VCC,UCC,QFU,DIFF,K)
       CALL UPWIND(UCC,VCC,VCC,QFV,DIFF,K)
C
C      PRESSURE BOUNDARY CONDITIONS
C
       DO I=1,NB(K)
        BCP(I)=DIFF(I)*(D2U(I,K)*XN(I,K)+D2V(I,K)*YN(I,K))
       END DO
C
C----------------------------------------------------------------------
C
C      SOLVE VELOCITY FIELD
C
       DO I=1,NB(K)
C
C       NO-SLIP WALLS OR INLETS
C
        IF ((KM(I,K).LE.0).OR.(KM(I,K).EQ.1)) THEN
         RHX(I)=GU(I,3,K)
         RHY(I)=GV(I,3,K)
        END IF
C
C       OUTLETS: PRESSURE AND OUTFLOW
C
        IF ((KM(I,K).EQ.2).OR.(KM(I,K).EQ.3)) THEN
         RHX(I)=0.D+000
         RHY(I)=0.D+000
        END IF
C
C       SLIP (SHEAR) WALLS
C
        IF (KM(I,K).EQ.4) THEN
         RHX(I)=GU(I,3,K)/VIE(K)
         RHY(I)=GV(I,3,K)/VIE(K)
        END IF
C
       END DO
C
C      BODY FORCES
C      
       DO I=NB(K)+1,NC(K)
        RHX(I)=(FBX(K)-DXP(I,K))/DENS(I)-GX*BEE(K)*(TCC(I,K)-TRE(K))
        RHY(I)=(FBY(K)-DYP(I,K))/DENS(I)-GY*BEE(K)*(TCC(I,K)-TRE(K))
       END DO
C
C      TRANSPORT
C
       CALL TRANSPORT(UCC,D2U,DXU,DYU,QFU,DIFF,RHX,GU,DT,RESX1,K)
       CALL TRANSPORT(VCC,D2V,DXV,DYV,QFV,DIFF,RHY,GV,DT,RESY1,K)
C
C      ACCUMULATE RESIDUALS
C
       RESX=RESX+RESX1
       RESY=RESY+RESY1
C
C----------------------------------------------------------------------
C
C      SOLVE FOR HELMHOLTZ POTENTIAL: HOMOGENEOUS BOUNDARY CONDITIONS
C
C      FORM COEFFICIENT MATRIX
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
C      FORM RHS VECTOR
C
       DO I=1,NB(K)
        BM(I)=0.D+000
       END DO
       DO I=NB(K)+1,NC(K)
        BM(I)=DENS(I)*(DXU(I,K)+DYV(I,K))/DT
       END DO
C
C      INITIAL CONDITION FOR HELMHOLTZ POTENTIAL
C
       DO I=1,NC(K)
        XM(I)=0.D+000
       END DO
C
C      SOLVE FOR HELMHOLTZ POTENTIAL
C
       MXP=10
       IF (ITP.EQ.1) MXP=100
       DO IIS=1,2
        CALL GMRES(AM,BM,XM,RESP1,MXP,K)
       END DO
C
C      APPROXIMATE DERIVATIVES
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
C      CORRECT VELOCITY
C
       DO I=1,NC(K)
        UCC(I,K)=UCC(I,K)-DT*DXP(I,K)/DENS(I)
        VCC(I,K)=VCC(I,K)-DT*DYP(I,K)/DENS(I)
       END DO
C
C      CALCULATE FLOW DEFECT TO ENFORCE CONTINUITY
C
       QINP=0.D+000
       QOUT=0.D+000
       AOUT=0.D+000
       DO I=1,NB(K)
        VN=UCC(I,K)*XN(I,K)+VCC(I,K)*YN(I,K)
        IF (KM(I,K).EQ.1) THEN
         QINP=QINP+VN*AR(I,K)
        END IF
        IF ((KM(I,K).EQ.2).OR.(KM(I,K).EQ.3)) THEN
         QOUT=QOUT+VN*AR(I,K)
         AOUT=AOUT+AR(I,K)
        END IF
       END DO
       DEFECT=QOUT+QINP
C
C      CORRECT OUTLET VELOCITIES
C
       DO I=1,NB(K)
        IF ((KM(I,K).EQ.2).OR.(KM(I,K).EQ.3)) THEN
         UCC(I,K)=UCC(I,K)-XN(I,K)*DEFECT/AOUT
         VCC(I,K)=VCC(I,K)-YN(I,K)*DEFECT/AOUT
        END IF
       END DO
C
C      VELOCITY DERIVATIVES
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
     &           +(FXX(I,II,K)+FYY(I,II,K))*UCC(ICONN(I,II,K),K)
         DXU(I,K)=DXU(I,K)+FXC(I,II,K)*UCC(ICONN(I,II,K),K)
         DYU(I,K)=DYU(I,K)+FYC(I,II,K)*UCC(ICONN(I,II,K),K)
         D2V(I,K)=D2V(I,K)
     &           +(FXX(I,II,K)+FYY(I,II,K))*VCC(ICONN(I,II,K),K)
         DXV(I,K)=DXV(I,K)+FXC(I,II,K)*VCC(ICONN(I,II,K),K)
         DYV(I,K)=DYV(I,K)+FYC(I,II,K)*VCC(ICONN(I,II,K),K)
        END DO
       END DO
C
C----------------------------------------------------------------------
C
C      COMPUTE CONTINUITY RESIDUAL 
C
       IF (DABS(QINP).GT.EPS) THEN
        RESC1=DABS(DEFECT/QINP)
       ELSE
        RESC1=EPS
       END IF
C
C      ACCUMULATE CONTINUITY RESIDUAL
C
       RESC=RESC+RESC1
C
C----------------------------------------------------------------------
C
C      CORRECT PRESSURE
C
       DO I=1,NC(K)
        PCC(I,K)=PCO(I)+PCC(I,K)
       END DO
C
C----------------------------------------------------------------------
C
C      SOLVE PRESSURE FIELD        
C
       IPRESSURE=0
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
        IF (KM(I,K).EQ.2) THEN
         IPRESSURE=1
         BM(I)=GP(I,3,K)
        ELSE
         BCP(I)=BCP(I)
     &         -(QFU(I)*XN(I,K)+QFV(I)*YN(I,K))
     &         +(FBX(K)*XN(I,K)+FBY(K)*YN(I,K))/DENS(I)
     &         -BEE(K)*(TCC(I,K)-TRE(K))*(GX*XN(I,K)+GY*YN(I,K))
     &         -((UCC(I,K)-UCO(I))*XN(I,K)+(VCC(I,K)-VCO(I))*YN(I,K))/DT
         BM(I)=DENS(I)*BCP(I)
        END IF
       END DO
C
       DO I=NB(K)+1,NC(K)
        DQF=0.D+000
        DO II=1,NCONN(I,K)
         DQF=DQF+FXC(I,II,K)*QFU(ICONN(I,II,K))
     &          +FYC(I,II,K)*QFV(ICONN(I,II,K)) 
        END DO
        BM(I)=-DENS(I)*(DQF+BEE(K)*(GX*DXT(I,K)+GY*DYT(I,K)))
       END DO
C
C      SOLUTION VECTOR  
C
       DO I=1,NC(K)
        XM(I)=PCC(I,K)
       END DO
C
C      SOLVE FOR PRESSURE
C
       MXP=10
       IF (ITP.EQ.1) MXP=100
       DO IIS=1,2
        CALL GMRES(AM,BM,XM,RESP1,MXP,K)
       END DO
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
C      FIX PRESSURE IF NOT IMPOSED
C
       IF (IPRESSURE.EQ.0) THEN
        PREF=PCC(1,K)
        DO I=1,NC(K)
         PCC(I,K)=PCC(I,K)-PREF
        END DO
       END IF
C
C----------------------------------------------------------------------
C
C      FIX BOUNDARY VELOCITY
C
       DO I=1,NB(K)
        IF ((KM(I,K).LE.0).OR.(KM(I,K).EQ.1)) THEN
         UCC(I,K)=GU(I,3,K)
         VCC(I,K)=GV(I,3,K)
        END IF
       END DO
C
C----------------------------------------------------------------------
C
C      RE-CALCULATE DILATATION AND TIME DERIVATIVE
C
       DO I=1,NC(K)
        DCC(I,K)=DXU(I,K)+DYV(I,K)
        DTD(I,K)=(DCC(I,K)-DCO(I))/DT
       END DO
C
C----------------------------------------------------------------------
C
C     END SUBREGION LOOP
C
      END IF
      END DO
C
C----------------------------------------------------------------------
C
C      AVERAGE INTERFACE MOMENTUM INTERACTION
C
      IF (NR.GT.1) THEN
       CALL INTERFACEMOMENTUM
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END