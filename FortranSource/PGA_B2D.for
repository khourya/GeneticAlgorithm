CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     CCC
CCC                         SUBROUTINE B2D                              CCC
CCC                                                                     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION B2D(IPAR,I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
C
      INTEGER IPAR(NPARMAX*NBITMAX),I
C
C**************************************************************************
C
C     DECODE THE Ith PARAMETER OF THE BINARY STRING IPAR
C
      III=0
      DO J=1,I
       III=III+NBIT(J)
      END DO
      ISUM=0
      DO II=NBIT(I),1,-1
       ISUM=ISUM+IPAR(III)*(2**(NBIT(I)-II))
       III=III-1
      END DO
C
      B2D=PARMIN(I)+(PARMAX(I)-PARMIN(I))*DBLE(ISUM)/DBLE(2**NBIT(I)-1)
C
C**************************************************************************
C
      RETURN
C
      END
