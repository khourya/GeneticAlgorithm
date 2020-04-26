CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUTSTRUCTURAL                        CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS STRUCTURAL CONDITIONS AND SOLUTION                              CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUTSTRUCTURAL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-M.for'
C
C----------------------------------------------------------------------
C
      INTEGER IOS
C
C----------------------------------------------------------------------
C
C     READ STRUCTURAL BOUNDARY CONDITIONS
C
      OPEN (11,FILE='Model/ALMA_d_str.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING STRUCTURAL BOUNDARY CONDITION FILE.'
       STOP
      END IF
C
      DO K=1,NR
       IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
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
C       KM=0:PRESSURE (RB = NORMAL TRACTION)
C       KM=1:DEFORMATION
C       KM=2:TRACTION
C       KM=3:SLIDING (RB = SHEAR TRACTION)
C
        DO I=1,NB(K)
         READ (11,*) II,KM(I,K),UB,VB,PB
C        PRESSURE 
         IF (KM(I,K).EQ.0) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)=-PB*XN(I,K)
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=-PB*YN(I,K)
         END IF
C        DEFORMATION 
         IF (KM(I,K).EQ.1) THEN
          GU(I,1,K)=1.D+000
          GU(I,2,K)=0.D+000
          GU(I,3,K)=UB
          GV(I,1,K)=1.D+000
          GV(I,2,K)=0.D+000
          GV(I,3,K)=VB
         END IF
C        TRACTION
         IF (KM(I,K).EQ.2) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)=UB
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=VB
         END IF
C        SLIDING WALL (RB = SHEAR TRACTION)
         IF (KM(I,K).EQ.3) THEN
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)= PB*YN(I,K)
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=-PB*XN(I,K)
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
C     STRUCTURAL OUTPUT
C
      OPEN (13,FILE='Model/ALMA_o_str.bin',STATUS='OLD',IOSTAT=IOS
     &        ,FORM='UNFORMATTED')
C
      IF (IOS.EQ.0) THEN
C
       DO K=1,NR
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         DO I=1,NC(K)
          READ (13) XC(I,K),UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,YC(I,K),VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K)
         END DO
        ELSE
         DO I=1,NC(K)
          READ (13) XDUM,UDUM,DXUDUM,DYUDUM,D2UDUM
     &             ,YDUM,VDUM,DXVDUM,DYVDUM,D2VDUM
         END DO
        END IF        
       END DO
C
      ELSE
C
       DO K=1,NR
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         DO I=1,NC(K)
          UCC(I,K)=UINI(K)
          XC(I,K)=XC(I,K)+UCC(I,K)
          DXU(I,K)=0.D+000
          DYU(I,K)=0.D+000
          D2U(I,K)=0.D+000
          VCC(I,K)=VINI(K)
          YC(I,K)=YC(I,K)+VCC(I,K)
          DXV(I,K)=0.D+000
          DYV(I,K)=0.D+000
          D2V(I,K)=0.D+000
         END DO
        END IF
       END DO
C
      END IF
C
      CLOSE (13)
C
C----------------------------------------------------------------------
C
C     INTERFACE STRUCTURAL BOUNDARY CONDITIONS
C
      DO K=1,NR
       IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
C
C       INTERFACE BOUNDARY CONDITIONS
C       KM<-1000:DEFORMATIONS, KM<-2000:TRACTIONS
C
        DO I=1,NB(K)
C        TYPE 1 INTERFACE: IMPOSE DEFORMATIONS
         IF ((KM(I,K).LT.(-1000)).AND.(KM(I,K).GT.(-2000))) THEN
          GU(I,1,K)=1.D+000
          GU(I,2,K)=0.D+000
          GU(I,3,K)=UCC(I,K)
          GV(I,1,K)=1.D+000
          GV(I,2,K)=0.D+000
          GV(I,3,K)=VCC(I,K)
         END IF
C        TYPE 2 INTERFACE: IMPOSE TRACTIONS
         IF ((KM(I,K).LT.(-2000)).AND.(KM(I,K).GT.(-3000))) THEN
          ELS=1.D+000-2.D+000*POE(K)
          TEX=2.D+000*VIE(K)*BEE(K)*(1.D+000+POE(K))/ELS
          DIL=POE(K)*(DXU(I,K)+DYV(I,K))/ELS
          SXX=2.D+000*VIE(K)*DXU(I,K)+2.D+000*VIE(K)*DIL
          SYY=2.D+000*VIE(K)*DYV(I,K)+2.D+000*VIE(K)*DIL
          SXY=VIE(K)*(DYU(I,K)+DXV(I,K))
          GU(I,1,K)=0.D+000
          GU(I,2,K)=1.D+000
          GU(I,3,K)=SXX*XN(I,K)+SXY*YN(I,K)
          GV(I,1,K)=0.D+000
          GV(I,2,K)=1.D+000
          GV(I,3,K)=SXY*XN(I,K)+SYY*YN(I,K)
         END IF
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