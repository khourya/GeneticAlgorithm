CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUTLEVELSET                          CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS LEVEL-SET CONDITIONS AND SOLUTION                               CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUTLEVELSET
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Field-S.for'
C
C----------------------------------------------------------------------
C
      INTEGER IOS
C
C----------------------------------------------------------------------
C
C     READ LEVEL-SET BOUNDARY CONDITIONS
C
      OPEN (11,FILE='Model/ALMA_d_lev.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING LEVEL-SET BOUNDARY CONDITION FILE.'
       STOP
      END IF
C
      DO K=1,NR
C
C      READ REGION NUMBER
C
       READ (11,*) KK
C
C      READ INITIAL VALUES
C
       READ (11,*) SINI(K)
C
C      BOUNDARY CONDITIONS
C      KS=0:IMPERMEABLE, KS=1:LEVEL (0-1), KS=2:FLUX, KS=3:RATIO
C
       DO I=1,NB(K)
        READ (11,*) II,KS(I,K),SB,CB,DB
C       IMPERMEABLE
        IF (KS(I,K).LE.0) THEN
         GS(I,1,K)=0.D+000
         GS(I,2,K)=1.D+000
         GS(I,3,K)=0.D+000
        END IF
C       LEVEL
        IF (KS(I,K).EQ.1) THEN
         GS(I,1,K)=1.D+000
         GS(I,2,K)=0.D+000
         GS(I,3,K)=SB
        END IF
C       FLUX
        IF (KS(I,K).EQ.2) THEN
         GS(I,1,K)=0.D+000
         GS(I,2,K)=1.D+000
         GS(I,3,K)=CB
        END IF
C       RATIO
        IF (KS(I,K).EQ.3) THEN
         GS(I,1,K)=DB
         GS(I,2,K)=1.D+000
         GS(I,3,K)=DB*SB 
        END IF
       END DO
C
      END DO
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     LEVEL-SET OUTPUT
C
      OPEN (13,FILE='Model/ALMA_o_lev.bin',STATUS='OLD',IOSTAT=IOS
     &        ,FORM='UNFORMATTED')
C
      IF (IOS.EQ.0) THEN
C
       DO K=1,NR
        DO I=1,NC(K)
         READ (13) SCC(I,K),DXS(I,K),DYS(I,K),D2S(I,K)
        END DO
       END DO
C
      ELSE
C
       DO K=1,NR
        DO I=1,NC(K)
         SCC(I,K)=SINI(K)
         DXS(I,K)=0.D+000
         DYS(I,K)=0.D+000
         D2S(I,K)=0.D+000
        END DO
       END DO
C
      END IF
C
      CLOSE (13)
C
C----------------------------------------------------------------------
C
      RETURN
      END