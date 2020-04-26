CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                         SUBROUTINE D2B                              CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE D2B(PAR,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
C
      REAL*8 PAR(NPARMAX)
      INTEGER IPAR(NPARMAX*NBITMAX)
      INTEGER IVAL
C
C**************************************************************************
C
C     CODE THE PARAMETER ARRAY INTO A BINARY STRING
C
      III=0
      DO I=1,NPAR
       IVAL=INT((2.**NBIT(I)-1.)*(PAR(I)-PARMIN(I))
     @                          /(PARMAX(I)-PARMIN(I)))
       DO II=NBIT(I),1,-1
        IPAR(III+II)=MOD(IVAL,2)
        IVAL=IVAL/2
       END DO
       III=III+NBIT(I)
      END DO
C
C**************************************************************************
C
      RETURN
C
      END
