CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUTMOMENTUM                          CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS MOMENTUM CONDITIONS AND SOLUTION                                CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUTMOMENTUM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
      INTEGER IOS
C
C----------------------------------------------------------------------
C
C     READ MOMENTUM BOUNDARY CONDITIONS
C
      OPEN (11,FILE='Model/ALMA_d_mom.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING MOMENTUM BOUNDARY CONDITION FILE.'
       STOP
      END IF
C
      DO K=1,NR
       IF (KR(K).EQ.1) THEN
C
C       READ REGION NUMBER
C
        READ (11,*) KK
C
C       READ INITIAL VALUES
C
        READ (11,*) UINI(K),VINI(K)
C
C       READ BODY FORCES
C
        READ (11,*) FBX(K),FBY(K)
C
C       BOUNDARY CONDITIONS: 
C       KM=0: NO-SLIP WALL
C       KM=1: FLOW INLET
C       KM=2: PRESSURE OUTLET
C       KM=3: OUTFLOW
C       KM=4: SLIP WALL (PB = SHEAR)
C
        DO I=1,NB(K)
         READ (11,*) II,KM(I,K),UB,VB,PB
C        SOLID NO-SLIP WALL
         IF (KM(I,K).LE.0) THEN
          GU(I,1,K)=1.D+000
          GU(I,2,K)=0.D+000
          GU(I,3,K)=UB
          GV(I,1,K)=1.D+000
          GV(I,2,K)=0.D+000
          GV(I,3,K)=VB
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=0.D+000
          GH(I,1,K)=0.D+000
          GH(I,2,K)=1.D+000
          GH(I,3,K)=0.D+000
         END IF
C        INLET 
         IF (KM(I,K).EQ.1) THEN
          GU(I,1,K)=1.D+000
          GU(I,2,K)=0.D+000
          GU(I,3,K)=UB
          GV(I,1,K)=1.D+000
          GV(I,2,K)=0.D+000
          GV(I,3,K)=VB
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=0.D+000
          GH(I,1,K)=0.D+000
          GH(I,2,K)=1.D+000
          GH(I,3,K)=0.D+000
         END IF
C        PRESSURE OUTLET
         IF (KM(I,K).EQ.2) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)=0.D+000
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=0.D+000
          GP(I,1,K)=1.D+000
          GP(I,2,K)=0.D+000
          GP(I,3,K)=PB
          GH(I,1,K)=1.D+000
          GH(I,2,K)=0.D+000
          GH(I,3,K)=0.D+000
         END IF
C        OUTFLOW
         IF (KM(I,K).EQ.3) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)=0.D+000
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=0.D+000
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=0.D+000
          GH(I,1,K)=1.D+000
          GH(I,2,K)=0.D+000
          GH(I,3,K)=0.D+000
         END IF
C        SOLID SLIP WALL
         IF (KM(I,K).EQ.4) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)= PB*YN(I,K)
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=-PB*XN(I,K)
          GP(I,1,K)=0.D+000
          GP(I,2,K)=1.D+000
          GP(I,3,K)=0.D+000
          GH(I,1,K)=0.D+000
          GH(I,2,K)=1.D+000
          GH(I,3,K)=0.D+000
         END IF
C
        END DO
C
       ELSE
C
C       READ DUMMY VALUES
C
        READ (11,*) KDUM
        READ (11,*) UDUM,VDUM
        READ (11,*) XDUM,YDUM
        DO I=1,NB(K)
         READ (11,*) IDUM,KDUM,UDUM,VDUM,RDUM
        END DO
C
       END IF
      END DO
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     MOMENTUM OUTPUT
C
      OPEN (13,FILE='Model/ALMA_o_mom.bin',STATUS='OLD',IOSTAT=IOS
     &        ,FORM='UNFORMATTED')
C
      IF (IOS.EQ.0) THEN
C
       DO K=1,NR
        IF (KR(K).EQ.1) THEN
         DO I=1,NC(K)
          READ (13) UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K)
     &             ,PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K)
         END DO
        ELSE
         DO I=1,NC(K)
          READ (13) UDUM,DXUDUM,DYUDUM,D2UDUM
     &             ,VDUM,DXVDUM,DYVDUM,D2VDUM
     &             ,PDUM,DXPDUM,DYPDUM,D2PDUM
         END DO
        END IF
       END DO
C
      ELSE
C
       DO K=1,NR
        IF (KR(K).EQ.1) THEN
         DO I=1,NC(K)
          UCC(I,K)=UINI(K)
          DXU(I,K)=0.D+000
          DYU(I,K)=0.D+000
          D2U(I,K)=0.D+000
          VCC(I,K)=VINI(K)
          DXV(I,K)=0.D+000
          DYV(I,K)=0.D+000
          D2V(I,K)=0.D+000
          PCC(I,K)=0.D+000
          DXP(I,K)=0.D+000
          DYP(I,K)=0.D+000
          D2P(I,K)=0.D+000
         END DO
        END IF
       END DO
C
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END