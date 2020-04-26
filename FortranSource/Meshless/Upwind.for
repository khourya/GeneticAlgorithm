CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE UPWIND                                 CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   PERFORMS 1ST OR 2ND ORDER UPWINDING FOR LAGRANGE DERIVATIVES          CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPWIND(UCC,VCC,HCC,QFH,DIFF,K)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
      INCLUDE '../../Include/Interpolation-U.for'
C
C----------------------------------------------------------------------
C
      REAL*8 UCC(NCMAX,NRMAX),VCC(NCMAX,NRMAX)
      REAL*8 HCC(NCMAX,NRMAX),QFH(NCMAX)
      REAL*8 DIFF(NCMAX)
      INTEGER K
C
      REAL*8 FXU(NFMAX),FYU(NFMAX)
C
      ICOU=1
      COUMAX=0.D+000
C
C----------------------------------------------------------------------
C
C     LOOP OVER DATA CENTERS
C
      DO I=1,NC(K)
C
       DO II=1,NCONN(I,K)
        FXU(II)=FXC(I,II,K)
        FYU(II)=FYC(I,II,K)
       END DO
C
C----------------------------------------------------------------------
C
C      LOCAL COURANT NUMBER
C
       COU=DABS(UCC(I,K))*DT/RXA(I,K)+DABS(VCC(I,K))*DT/RYA(I,K)
       IF (COU.GT.COUMAX) THEN
        ICOU=I
        COUMAX=COU
       END IF
C
C----------------------------------------------------------------------
C
C      LOCAL X AND Y REYNOLD NUMBERS
C
       REYX=UCC(I,K)*RXA(I,K)/DIFF(I)
       REYY=VCC(I,K)*RYA(I,K)/DIFF(I)
C
C----------------------------------------------------------------------
C
C      DETERMINE UPWINDING ORDER BLENDING
C
       IF ((DABS(REYX).GT.1.D+000).OR.(DABS(REYY).GT.1.D+000)) THEN
C
        THETA=1.0+000*THU
C
        NMAX=2
        NMIN=2
        PHIMAX=HCC(ICONN(I,2,K),K)
        PHIMIN=HCC(ICONN(I,2,K),K)
        DO II=3,NCONN(I,K)
         IF(HCC(ICONN(I,II,K),K).GT.PHIMAX) PHIMAX=HCC(ICONN(I,II,K),K)
         IF(HCC(ICONN(I,II,K),K).LT.PHIMIN) PHIMIN=HCC(ICONN(I,II,K),K)
        END DO
        DEL=(PHIMAX-PHIMIN)/1.D+002
        IF (DEL.LT.1.D-010) DEL=1.D-010
        IF ((HCC(I,K)-PHIMAX).LT.-DEL) THETA=0.D+000*THU
        IF ((HCC(I,K)-PHIMIN).GT. DEL) THETA=0.D+000*THU
C
C       X-UPWINDING
C
        IF (REYX.GT.1.D+000) THEN
         DO II=1,NCONN(I,K)
          FXU(II)=(1.D+000-THETA)*FXW(I,II,K)+THETA*SXW(I,II,K)
         END DO
        END IF
C
        IF (REYX.LT.-1.D+000) THEN
         DO II=1,NCONN(I,K)
          FXU(II)=(1.D+000-THETA)*FXE(I,II,K)+THETA*SXE(I,II,K)
         END DO
        END IF
C
C       Y-UPWINDING
C
        IF (REYY.GT.1.D+000) THEN
         DO II=1,NCONN(I,K)
          FYU(II)=(1.D+000-THETA)*FYS(I,II,K)+THETA*SYS(I,II,K)
         END DO
        END IF
C
        IF (REYY.LT.-1.D+000) THEN
         DO II=1,NCONN(I,K)
          FYU(II)=(1.D+000-THETA)*FYN(I,II,K)+THETA*SYN(I,II,K)
         END DO
        END IF
C
       END IF
C
C----------------------------------------------------------------------
C
C      CONVECTIVE FLUXES
C
       QHX=0.D+000
       QHY=0.D+000
       DO II=1,NCONN(I,K)
        QHX=QHX+FXU(II)*HCC(ICONN(I,II,K),K)
        QHY=QHY+FYU(II)*HCC(ICONN(I,II,K),K)
       END DO
       QFH(I)=UCC(I,K)*QHX+VCC(I,K)*QHY
C
C----------------------------------------------------------------------
C
      END DO
C
C----------------------------------------------------------------------
C
C     COURANT NUMBER WARNING
C
      IF (COUMAX.GT.1.D+000) THEN
       I=ICOU
       DTMAX=1.D+000/(DABS(UCC(I,K))/RXA(I,K)+DABS(VCC(I,K))/RYA(I,K))
       WRITE (*,*) '...WARNING!!! COURANT LIMIT EXCEEDED'
       WRITE (*,*) '...DECREASE THE TIME STEP TO AT MOST ',DTMAX
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END