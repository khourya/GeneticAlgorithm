CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INTERFACEPOREPRESSURE                  CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   INTERPOLATES PORE PRESSURE AND DARCY VELOCITY ACCROSS INTERFACES      CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERFACEPOREPRESSURE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
C     LOOP OVER SUBREGIONS TO AVERAGE DATA
C
      DO KK=1,NR
      IF (KR(KK).EQ.2) THEN
C
       DO II=1,NB(KK)
       IF (KP(II,KK).LT.0) THEN
C
C       FIND CLOSEST POINT IN OPOSITE REGION
C
        RMIN=1.D+040
        IJ=II
        KJ=KK
        DO K=1,NR
         IF (K.NE.KK) THEN
          DO I=1,NB(K)
           XX=XC(I,K)-XC(II,KK)
           YY=YC(I,K)-YC(II,KK)
           RR=XX*XX+YY*YY
           IF (RR.LT.RMIN) THEN
            IJ=I
            KJ=K
            RMIN=RR
           END IF
          END DO
         END IF
        END DO
C
C       DETERMINE VALUES AT CLOSEST POINT
C
        IF (KR(KJ).EQ.0) THEN
         PCCIJKJ=PCC(II,KK)
         VNCIJKJ=0.D+000
        END IF
C
        IF ((KR(KJ).EQ.1).OR.(KR(KJ).EQ.3)) THEN
         PCCIJKJ=PCC(IJ,KJ)
         VNCIJKJ=UCC(IJ,KJ)*XN(IJ,KJ)+VCC(IJ,KJ)*YN(IJ,KJ)
        END IF
C
        IF (KR(KJ).EQ.2) THEN
         PCCIJKJ=PCC(IJ,KJ)
         VNCIJKJ=VFX(IJ,KJ)*XN(IJ,KJ)+VFY(IJ,KJ)*YN(IJ,KJ)
        END IF
C
C       AVERAGE BOUNDARY VALUES 
C
        IF ((KP(II,KK).LT.(-1000)).AND.(KP(II,KK).GT.(-2000))) THEN
         PCC(II,KK)=(1.D+000-THI)*PCC(II,KK)+THI*PCCIJKJ
         GP(II,3,KK)=PCC(II,KK)
        END IF
C
        IF ((KP(II,KK).LT.(-2000)).AND.(KP(II,KK).GT.(-3000))) THEN
         VNCIIKK=VFX(II,KK)*XN(II,KK)+VFY(II,KK)*YN(II,KK)
         GP(II,3,KK)=(1.D+000-THI)*VNCIIKK+THI*(-VNCIJKJ)
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