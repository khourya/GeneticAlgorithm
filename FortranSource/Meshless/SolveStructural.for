CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVESTRUCTURAL                        CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES STRUCTURAL PROBLEM                                             CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVESTRUCTURAL
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
      REAL*8 DCO(NCMAX)
C
      REAL*8 RESS(NDIM)
C
      REAL*8 AM(NDIM*NCMAX,NDIM*NFMAX),BM(NDIM*NCMAX),XM(NDIM*NCMAX)
C
C----------------------------------------------------------------------
C
C     INITIALIZE RESIDUAL
C
      RESU=0.D+000
      RESV=0.D+000
C
C----------------------------------------------------------------------
C
C     LOOP THROUGH SUB-REGIONS
C
      DO K=1,NR
      IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
C
C      RESTORE FIELD POINT POSITIONS
C
       DO I=1,NC(K)
        XC(I,K)=XC(I,K)-UCC(I,K)
        YC(I,K)=YC(I,K)-VCC(I,K)
       END DO
C
C      OLD DILATATION
C
       DO I=1,NC(K)
        DCO(I)=DXU(I,K)+DYV(I,K)
       END DO
C
C      COEFFICIENTS
C
       ELS=1.D+000-2.D+000*POE(K)
       TEX=2.D+000*VIE(K)*BEE(K)*(1.D+000+POE(K))/ELS
       RAS=(1.D+000-POE(K))/ELS
C
C----------------------------------------------------------------------
C
C      BOUNDARY RIGHT-HAND SIDE TERMS
C
       DO I=1,NB(K)
        IDIM=NDIM*(I-1)
C       DEFORMATIONS
        IF (KM(I,K).EQ.1) THEN
         BM(IDIM+1)=GU(I,3,K)
         BM(IDIM+2)=GV(I,3,K)
        END IF
C       TRACTIONS
        IF ((KM(I,K).EQ.0).OR.(KM(I,K).EQ.2)) THEN
         BM(IDIM+1)=(GU(I,3,K)
     &             +(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*XN(I,K))/VIE(K)
         BM(IDIM+2)=(GV(I,3,K)
     &             +(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*YN(I,K))/VIE(K)
        END IF
C       SLIDING WALLS
        IF (KM(I,K).EQ.3) THEN
         BM(IDIM+1)=GU(I,3,K)/VIE(K)
         BM(IDIM+2)=GV(I,3,K)/VIE(K)
        END IF
C       TYPE 1 INTERFACE: DEFORMATIONS
        IF ((KM(I,K).LT.(-1000)).AND.(KM(I,K).GT.(-2000))) THEN
         BM(IDIM+1)=GU(I,3,K)
         BM(IDIM+2)=GV(I,3,K)
        END IF
C       TYPE 2 INTERFACE: TRACTIONS
        IF ((KM(I,K).LT.(-2000)).AND.(KM(I,K).GT.(-3000))) THEN
         BM(IDIM+1)=(GU(I,3,K)
     &             +(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*XN(I,K))/VIE(K)
         BM(IDIM+2)=(GV(I,3,K)
     &             +(PCC(I,K)+TEX*(TCC(I,K)-TRE(K)))*YN(I,K))/VIE(K)
        END IF
C
       END DO
C
C      INTERNAL RIGHT-HAND SIDE TERMS
C
       DO I=NB(K)+1,NC(K)
        IDIM=NDIM*(I-1)
        BM(IDIM+1)=(DXP(I,K)+TEX*DXT(I,K)-FBX(K))/VIE(K)
        BM(IDIM+2)=(DYP(I,K)+TEX*DYT(I,K)-FBY(K))/VIE(K)
       END DO
C
C----------------------------------------------------------------------
C
C      BOUNDARY COEFFICIENT MATRIX
C
       DO I=1,NB(K)
        IDIM=NDIM*(I-1)
C       DEFORMATIONS
        IF (KM(I,K).EQ.1) THEN
         AM(IDIM+1,1)=1.D+000
         AM(IDIM+1,2)=0.D+000
         AM(IDIM+2,1)=0.D+000
         AM(IDIM+2,2)=1.D+000
         DO II=2,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=0.D+000
          AM(IDIM+1,JDIM+2)=0.D+000
          AM(IDIM+2,JDIM+1)=0.D+000
          AM(IDIM+2,JDIM+2)=0.D+000
         END DO
        END IF
C       TRACTIONS
        IF ((KM(I,K).EQ.0).OR.(KM(I,K).EQ.2)) THEN
         DO II=1,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=2.D+000*XN(I,K)*RAS*FXC(I,II,K)
     &                     +YN(I,K)*FYC(I,II,K)
          AM(IDIM+1,JDIM+2)=2.D+000*POE(K)*XN(I,K)*FYC(I,II,K)/ELS
     &                     +YN(I,K)*FXC(I,II,K) 
          AM(IDIM+2,JDIM+1)=2.D+000*POE(K)*YN(I,K)*FXC(I,II,K)/ELS
     &                     +XN(I,K)*FYC(I,II,K)
          AM(IDIM+2,JDIM+2)=2.D+000*YN(I,K)*RAS*FYC(I,II,K)
     &                     +XN(I,K)*FXC(I,II,K)
         END DO
        END IF
C       SLIDING WALLS
        IF (KM(I,K).EQ.3) THEN
         DO II=1,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=FXC(I,II,K)*XN(I,K)+FYC(I,II,K)*YN(I,K)
          AM(IDIM+1,JDIM+2)=0.D+000
          AM(IDIM+2,JDIM+1)=0.D+000
          AM(IDIM+2,JDIM+2)=FXC(I,II,K)*XN(I,K)+FYC(I,II,K)*YN(I,K)
         END DO
        END IF
C       TYPE 1 INTERFACE: DEFORMATIONS
        IF ((KM(I,K).LT.(-1000)).AND.(KM(I,K).GT.(-2000))) THEN
         AM(IDIM+1,1)=1.D+000
         AM(IDIM+1,2)=0.D+000
         AM(IDIM+2,1)=0.D+000
         AM(IDIM+2,2)=1.D+000
         DO II=2,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=0.D+000
          AM(IDIM+1,JDIM+2)=0.D+000
          AM(IDIM+2,JDIM+1)=0.D+000
          AM(IDIM+2,JDIM+2)=0.D+000
         END DO
        END IF
C       TYPE 2 INTERFACE: TRACTIONS
        IF ((KM(I,K).LT.(-2000)).AND.(KM(I,K).GT.(-3000))) THEN
         DO II=1,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=2.D+000*XN(I,K)*RAS*FXC(I,II,K)
     &                     +YN(I,K)*FYC(I,II,K)
          AM(IDIM+1,JDIM+2)=2.D+000*POE(K)*XN(I,K)*FYC(I,II,K)/ELS
     &                     +YN(I,K)*FXC(I,II,K) 
          AM(IDIM+2,JDIM+1)=2.D+000*POE(K)*YN(I,K)*FXC(I,II,K)/ELS
     &                     +XN(I,K)*FYC(I,II,K)
          AM(IDIM+2,JDIM+2)=2.D+000*YN(I,K)*RAS*FYC(I,II,K)
     &                     +XN(I,K)*FXC(I,II,K)
         END DO
        END IF
C
       END DO
C
C      INTERNAL COEFFICIENT MATRIX
C
       DO I=NB(K)+1,NC(K)
        IDIM=NDIM*(I-1)
         DO II=1,NCONN(I,K)
          JDIM=NDIM*(II-1)
          AM(IDIM+1,JDIM+1)=FXX(I,II,K)+FYY(I,II,K)+FXX(I,II,K)/ELS
          AM(IDIM+1,JDIM+2)=FXY(I,II,K)/ELS
          AM(IDIM+2,JDIM+1)=FXY(I,II,K)/ELS
          AM(IDIM+2,JDIM+2)=FXX(I,II,K)+FYY(I,II,K)+FYY(I,II,K)/ELS
         END DO
       END DO
C
C----------------------------------------------------------------------
C
C      SOLUTION VECTOR
C
       DO I=1,NC(K)
        IDIM=NDIM*(I-1)
        XM(IDIM+1)=UCC(I,K)
        XM(IDIM+2)=VCC(I,K)
       END DO
C
C----------------------------------------------------------------------
C
C      SOLVE DEFORMATION FIELD
C
       MXP=10
       DO IIS=1,4
        CALL GMRESMULTI(AM,BM,XM,RESS,MXP,K)
       END DO
C
       IF (ITP.EQ.0) THEN
        IIS=0
        MXP=200
        RESS(1)=1.D+000
        RESS(2)=1.D+000
        WRITE  (*,*) 'STRUCTURE ITERATION | RESIDUAL-X | RESIDUAL-Y |'
        DO WHILE ((IIS.LT.25).AND.
     &            ((RESS(1).GT.1.D-003).OR.(RESS(2).GT.1.D-003)))
         CALL GMRESMULTI(AM,BM,XM,RESS,MXP,K)
         IIS=IIS+1
         WRITE  (*,'(12X,I6,5X,E10.5,3X,E10.5)') IIS,RESS(1),RESS(2)
        END DO
       END IF
C
C----------------------------------------------------------------------
C
C      FIELD DERIVATIVES
C
       DO I=1,NC(K)
        IDIM=NDIM*(I-1)
        UCC(I,K)=XM(IDIM+1)
        VCC(I,K)=XM(IDIM+2)
        D2U(I,K)=0.D+000
        DXU(I,K)=0.D+000
        DYU(I,K)=0.D+000
        D2V(I,K)=0.D+000
        DXV(I,K)=0.D+000
        DYV(I,K)=0.D+000
        DO II=1,NCONN(I,K)
         IJDIM=NDIM*(ICONN(I,II,K)-1)
         D2U(I,K)=D2U(I,K)+(FXX(I,II,K)+FYY(I,II,K))*XM(IJDIM+1)
         DXU(I,K)=DXU(I,K)+FXC(I,II,K)*XM(IJDIM+1)
         DYU(I,K)=DYU(I,K)+FYC(I,II,K)*XM(IJDIM+1)
         D2V(I,K)=D2V(I,K)+(FXX(I,II,K)+FYY(I,II,K))*XM(IJDIM+2)
         DXV(I,K)=DXV(I,K)+FXC(I,II,K)*XM(IJDIM+2)
         DYV(I,K)=DYV(I,K)+FYC(I,II,K)*XM(IJDIM+2)
        END DO
       END DO
C
C----------------------------------------------------------------------
C
C      ACCUMULATE RESIDUAL
C
       RESU=RESU+RESS(1)
       RESV=RESV+RESS(2)
C
C----------------------------------------------------------------------
C
C      UPDATE FIELD POINT POSITIONS
C
       DO I=1,NC(K)
        XC(I,K)=XC(I,K)+UCC(I,K)
        YC(I,K)=YC(I,K)+VCC(I,K)
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
      END IF
      END DO
C
C----------------------------------------------------------------------
C
C      AVERAGE INTERFACE STRUCTURAL INTERACTION
C
      IF (NR.GT.1) THEN
       CALL INTERFACESTRUCTURAL
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END