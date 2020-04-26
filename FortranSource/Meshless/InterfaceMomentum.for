CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INTERFACEMOMENTUM                      CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   INTERPOLATES DEFORMATION AND PRESSURE VALUES ACCROSS INTERFACES       CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERFACEMOMENTUM
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
      IF ((KR(KK).EQ.1).OR.(KR(KK).EQ.3)) THEN
C
       DO II=1,NB(KK)
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
C       MATCHING INTERFACE
C
        IF ((KM(IJ,KJ).LT.0).OR.(KP(IJ,KJ).LT.0)) THEN
C
C        DEFORMATIONS
C
         XC(II,KK)=XC(IJ,KJ)
         YC(II,KK)=YC(IJ,KJ)
C
C        INLETS
C
         IF (KM(II,KK).EQ.1) THEN
          IF (KR(KJ).EQ.2) THEN
           UCCIJKJ=VFX(IJ,KJ)
           VCCIJKJ=VFY(IJ,KJ)
           UCC(II,KK)=(1.D+000-THI)*UCC(II,KK)+THI*UCCIJKJ
           VCC(II,KK)=(1.D+000-THI)*VCC(II,KK)+THI*VCCIJKJ
           GU(II,3,KK)=UCC(II,KK)
           GV(II,3,KK)=VCC(II,KK)
          END IF
         END IF
C
C        OUTLETS
C
         IF ((KM(II,KK).EQ.2).OR.(KM(II,KK).EQ.3)) THEN
          IF (KR(KJ).EQ.2) THEN
           PCCIJKJ=PCC(IJ,KJ)
           PCC(II,KK)=(1.D+000-THI)*PCC(II,KK)+THI*PCCIJKJ
           GP(II,3,KK)=PCC(II,KK)
          END IF
         END IF
C
        END IF
C
C----------------------------------------------------------------------
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