CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC   SUBROUTINE TRANSPORT                                                  CCC
CCC                                                                         CCC
CCC   ADVANCES EQUATION: dH/dt+V.Grad(H)=Diff.L(H)+BOD                      CCC
CCC                                                                         CCC
CCC   EXPLICITLY IN TIME                                                    CCC
CCC                                                                         CCC
CCC   WHERE: V.Grad(H)=QFH                                                  CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE TRANSPORT(HCC,D2H,DXH,DYH,QFH,DIFF,RHS,GH,DTH,RESH,K)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
C
C----------------------------------------------------------------------
C
      REAL*8 HCC(NCMAX,NRMAX),D2H(NCMAX,NRMAX)
      REAL*8 DXH(NCMAX,NRMAX),DYH(NCMAX,NRMAX)
      REAL*8 QFH(NCMAX)
C
      REAL*8 DIFF(NCMAX)
      REAL*8 RHS(NCMAX)
      REAL*8 GH(NBMAX,3,NRMAX)
      REAL*8 DTH
      INTEGER K
C     
      REAL*8 HCO(NCMAX)
      REAL*8 HMIN,HMAX,RESH
C
C----------------------------------------------------------------------
C
C     PREVIOUS FIELD
C
      DO I=1,NC(K)
       HCO(I)=HCC(I,K)
      END DO
C
C----------------------------------------------------------------------
C
C     SUB-LEVEL ITERATIONS
C
      DO IS=1,ISUB
C
C      APPLY BOUNDARY CONDITIONS AT BOUNDARY NODES
C
       DO I=1,NB(K)
        GN1=GH(I,1,K)
        GN2=GH(I,2,K)
        GN3=RHS(I)
        DN1=FXC(I,1,K)*XN(I,K)+FYC(I,1,K)*YN(I,K)
        DNH=DXH(I,K)*XN(I,K)+DYH(I,K)*YN(I,K)
        HCC(I,K)=(GN3+GN2*DN1*HCC(I,K)-GN2*DNH)/(GN1+GN2*DN1)
       END DO
C
C      APPLY GOVERNING EQUATION AT INTERNAL NODES
C
       DO I=NB(K)+1,NC(K)
        HCC(I,K)=HCO(I)+DTH*(DIFF(I)*D2H(I,K)+RHS(I)-QFH(I))
       END DO
C
C      APPROXIMATE DERIVATIVES
C
       DO I=1,NC(K)
        D2H(I,K)=0.D+000
        DXH(I,K)=0.D+000
        DYH(I,K)=0.D+000
        DO II=1,NCONN(I,K)
         D2H(I,K)=D2H(I,K)
     &           +(FXX(I,II,K)+FYY(I,II,K))*HCC(ICONN(I,II,K),K)
         DXH(I,K)=DXH(I,K)+FXC(I,II,K)*HCC(ICONN(I,II,K),K)
         DYH(I,K)=DYH(I,K)+FYC(I,II,K)*HCC(ICONN(I,II,K),K)
        END DO
       END DO
C
      END DO
C
C----------------------------------------------------------------------
C
C     COMPUTE RESIDUAL
C
cc      RESH=0.D+000
cc      HMIN= 1.D+020
cc      HMAX=-1.D+020
cc      DO I=NB(K)+1,NC(K)
cc       RESH=RESH+(HCC(I,K)-HCO(I))*(HCC(I,K)-HCO(I))
cc       IF (HCC(I,K).LT.HMIN) HMIN=HCC(I,K)
cc       IF (HCC(I,K).GT.HMAX) HMAX=HCC(I,K)
cc      END DO
C
C     NON-DIMENSIONALIZE RESIDUAL
C
cc     IF ((HMAX-HMIN).GT.EPS) THEN
cc       RESH=DSQRT(RESH/DBLE(NI(K)))/(HMAX-HMIN)
cc      ELSE
       RESH=EPS
cc      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END