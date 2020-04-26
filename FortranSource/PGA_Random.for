CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                         SUBROUTINE RANDOM                           CCC
CCC                                                                     CCC
CCC   PURPOSE:                                                          CCC
CCC       GENERATES A RANDOM NUMBER BETWEEN 0 AND 1                     CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RANDOM(IDUM,RAND)
C
C     NEGATIVE IDUM REINITIALIZES THE SEQUENCE
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      SAVE
C
      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
C
      DIMENSION MA(55)
C
      DATA IFF /0/
C
C**************************************************************************
C
      IF ((IDUM.LT.0).OR.(IFF.EQ.0)) THEN
       IFF=1
       MJ=MSEED-DBLE(IABS(IDUM))
       MJ=DMOD(MJ,MBIG)
       MA(55)=MJ
       MK=1
       DO I=1,54
        II=MOD(21*I,55)
        MA(II)=MK
        MK=MJ-MK
        IF (MK.LT.MZ) MK=MK+MBIG
        MJ=MA(II)
       END DO
       DO K=1,4
        DO I=1,55
         MA(I)=MA(I)-MA(1+MOD(I+30,55))
         IF (MA(I).LT.MZ) MA(I)=MA(I)+MBIG
        END DO
       END DO
       INEXT=0
       INEXTP=31
      END IF
C
      INEXT=INEXT+1
      IF (INEXT.EQ.56) INEXT=1
      INEXTP=INEXTP+1
      IF (INEXTP.EQ.56) INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAND=MJ*FAC
C
C**************************************************************************
C
      RETURN
C
      END