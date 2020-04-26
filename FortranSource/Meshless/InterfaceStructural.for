CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INTERFACESTRUCTURAL                    CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   INTERPOLATES DEFORMATION AND TRACTION VALUES ACCROSS INTERFACES       CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERFACESTRUCTURAL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-E.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
C     LOOP OVER SUBREGIONS TO AVERAGE DATA
C
      DO KK=1,NR
      IF ((KR(KK).EQ.0).OR.(KR(KK).EQ.2)) THEN
C
       DO II=1,NB(KK)
       IF (KM(II,KK).LT.0) THEN
C
C       SHIFT BOUNDARY CONDITION FLAG
C
        IF ((KM(II,KK).LT.(-1000)).AND.(KM(II,KK).GT.(-2000))) THEN
         KMIIKK=KM(II,KK)-1000
        END IF
        IF ((KM(II,KK).LT.(-2000)).AND.(KM(II,KK).GT.(-3000))) THEN
         KMIIKK=KM(II,KK)+1000
        END IF
C
C       FIND CLOSEST POINT IN OPOSITE REGION
C
        RMIN=1.D+040
        IJ=II
        KJ=KK
        DO K=1,NR
         IF (K.NE.KK) THEN
          DO I=1,NB(K)
           IF (KM(I,K).EQ.KMIIKK) THEN
            XX=XC(I,K)-XC(II,KK)
            YY=YC(I,K)-YC(II,KK)
            RR=XX*XX+YY*YY
            IF (RR.LT.RMIN) THEN
             IJ=I
             KJ=K
             RMIN=RR
            END IF
           END IF
          END DO
         END IF
        END DO
C
C       DETERMINE VALUES AT CLOSEST POINT
C
        IF ((KR(KJ).EQ.0).OR.(KR(KJ).EQ.2)) THEN
         UCCIJKJ=UCC(IJ,KJ)
         VCCIJKJ=VCC(IJ,KJ)
         ELS=1.D+000-2.D+000*POE(KJ)
         TEX=2.D+000*VIE(KJ)*BEE(KJ)*(1.D+000+POE(KJ))/ELS
         SXX=2.D+000*VIE(KJ)*DXU(IJ,KJ)
     &      +2.D+000*VIE(KJ)*POE(KJ)*DCC(IJ,KJ)/ELS
     &      -PCC(IJ,KJ)-TEX*(TCC(IJ,KJ)-TRE(KJ))
         SYY=2.D+000*VIE(KJ)*DYV(IJ,KJ)
     &      +2.D+000*VIE(KJ)*POE(KJ)*DCC(IJ,KJ)/ELS
     &      -PCC(IJ,KJ)-TEX*(TCC(IJ,KJ)-TRE(KJ))
         SXY=VIE(KJ)*(DYU(IJ,KJ)+DXV(IJ,KJ))
         TRXIJKJ=SXX*XN(IJ,KJ)+SXY*YN(IJ,KJ)
         TRYIJKJ=SXY*XN(IJ,KJ)+SYY*YN(IJ,KJ)
        END IF
C
        IF ((KR(KJ).EQ.1).OR.(KR(KJ).EQ.3)) THEN
         UCCIJKJ=UCC(II,KK)
         VCCIJKJ=VCC(II,KK)
         SXX=2.D+000*VIE(KJ)*DXU(IJ,KJ)-PCC(IJ,KJ)
         SYY=2.D+000*VIE(KJ)*DYV(IJ,KJ)-PCC(IJ,KJ)
         SXY=VIE(KJ)*(DYU(IJ,KJ)+DXV(IJ,KJ))
         TRXIJKJ=SXX*XN(IJ,KJ)+SXY*YN(IJ,KJ)
         TRYIJKJ=SXY*XN(IJ,KJ)+SYY*YN(IJ,KJ)
        END IF
C
C       AVERAGE BOUNDARY VALUES 
C
        IF ((KM(II,KK).LT.(-1000)).AND.(KM(II,KK).GT.(-2000))) THEN
         UCC(II,KK)=(1.D+000-THI)*UCC(II,KK)+THI*UCCIJKJ
         VCC(II,KK)=(1.D+000-THI)*VCC(II,KK)+THI*VCCIJKJ
         GU(II,3,KK)=UCC(II,KK)
         GV(II,3,KK)=VCC(II,KK)
        END IF
C
        IF ((KM(II,KK).LT.(-2000)).AND.(KM(II,KK).GT.(-3000))) THEN
         ELS=1.D+000-2.D+000*POE(KK)
         TEX=2.D+000*VIE(KK)*BEE(KK)*(1.D+000+POE(KK))/ELS
         SXX=2.D+000*VIE(KK)*DXU(II,KK)
     &      +2.D+000*VIE(KK)*POE(KK)*DCC(II,KK)/ELS
     &      -PCC(II,KK)-TEX*(TCC(II,KK)-TRE(KK))
         SYY=2.D+000*VIE(KK)*DYV(II,KK)
     &      +2.D+000*VIE(KK)*POE(KK)*DCC(II,KK)/ELS
     &      -PCC(II,KK)-TEX*(TCC(II,KK)-TRE(KK))
         SXY=VIE(KK)*(DYU(II,KK)+DXV(II,KK))
         TRXIIKK=SXX*XN(II,KK)+SXY*YN(II,KK)
         TRYIIKK=SXY*XN(II,KK)+SYY*YN(II,KK)
         GU(II,3,KK)=(1.D+000-THI)*TRXIIKK+THI*(-TRXIJKJ)
         GV(II,3,KK)=(1.D+000-THI)*TRYIIKK+THI*(-TRYIJKJ)
        END IF
C
       END IF
       END DO
C
      END IF
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END