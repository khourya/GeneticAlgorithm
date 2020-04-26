CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INTERFACEENERGY                        CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   INTERPOLATES TEMPERATURE AND FLUX VALUES ACCROSS INTERFACES           CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERFACEENERGY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-E.for'
C
C----------------------------------------------------------------------
C
C     LOOP OVER SUBREGIONS TO AVERAGE DATA
C
      DO KK=1,NR
C
       DO II=1,NB(KK)
       IF (KE(II,KK).LT.0) THEN
C
C       SHIFT BOUNDARY CONDITION FLAG
C
        IF ((KE(II,KK).LT.(-1000)).AND.(KE(II,KK).GT.(-2000))) THEN
         KEIIKK=KE(II,KK)-1000
        END IF
        IF ((KE(II,KK).LT.(-2000)).AND.(KE(II,KK).GT.(-3000))) THEN
         KEIIKK=KE(II,KK)+1000
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
           IF (KE(I,K).EQ.KEIIKK) THEN
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
        TCCIJKJ=TCC(IJ,KJ)
        QX=DXT(IJ,KJ)
        QY=DYT(IJ,KJ)
        QNCIJKJ=-TCE(KJ)*(QX*XN(IJ,KJ)+QY*YN(IJ,KJ))
C
C       AVERAGE BOUNDARY VALUES 
C
        IF ((KE(II,KK).LT.(-1000)).AND.(KE(II,KK).GT.(-2000))) THEN
         TCC(II,KK)=(1.D+000-THI)*TCC(II,KK)+THI*TCCIJKJ
         GT(II,3,KK)=TCC(II,KK)
        END IF
C
        IF ((KE(II,KK).LT.(-2000)).AND.(KE(II,KK).GT.(-3000))) THEN
         QX=DXT(II,KK)
         QY=DYT(II,KK)
         QNCIIKK=-TCE(KK)*(QX*XN(II,KK)+QY*YN(II,KK))
         GT(II,3,KK)=(1.D+000-THI)*QNCIIKK+THI*(-QNCIJKJ)
        END IF
C
       END IF
       END DO
C
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END