CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUTENERGY                            CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS ENERGY CONDITIONS AND SOLUTION                                  CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INPUTENERGY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-E.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
C
C----------------------------------------------------------------------
C
      INTEGER IOS
	REAL*8 TCV(NCMAX)
C
C----------------------------------------------------------------------
C
C     READ ENERGY BOUNDARY CONDITIONS
C
      OPEN (11,FILE='Model/ALMA_d_ene.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING ENERGY BOUNDARY CONDITION FILE.'
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
       READ (11,*) TINI(K)
C
C      READ BODY FORCES
C
       READ (11,*) UBG(K)
C
C      BOUNDARY CONDITIONS
C      KE=0:INSULATED, KE=1:DIRICHLET, KE=2:NEUMANN, KE=3:ROBIN
C      KE<-1000:INTERFACE TEMPERATURE, KE<-2000:INTERFACE FLUX
C
       DO I=1,NB(K)
        READ (11,*) II,KE(I,K),TB,QB,HB
C       INSULATED WALL
        IF (KE(I,K).EQ.0) THEN
         GT(I,1,K)=0.D+000
         GT(I,2,K)=1.D+000
         GT(I,3,K)=0.D+000
        END IF
C       DIRICHLET (TEMP)
        IF (KE(I,K).EQ.1) THEN
         GT(I,1,K)=1.D+000
         GT(I,2,K)=0.D+000
         GT(I,3,K)=TB
        END IF
C       NEUMANN (FLUX)
        IF (KE(I,K).EQ.2) THEN
         GT(I,1,K)=0.D+000
         GT(I,2,K)=-TCE(K)
         GT(I,3,K)=QB
        END IF
C       ROBIN (CONVECTION)
        IF (KE(I,K).EQ.3) THEN
         GT(I,1,K)=HB
         GT(I,2,K)=TCE(K)
         GT(I,3,K)=HB*TB
        END IF
C
       END DO
C
      END DO
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     ENERGY OUTPUT
C
      OPEN (13,FILE='Model/ALMA_o_ene.bin',STATUS='OLD',IOSTAT=IOS
     &        ,FORM='UNFORMATTED')
C
      IF (IOS.EQ.0) THEN
C
       DO K=1,NR
        DO I=1,NC(K)
         READ (13) TCC(I,K),DXT(I,K),DYT(I,K),D2T(I,K)
        END DO
       END DO
C
      ELSE
C
       DO K=1,NR
        DO I=1,NC(K)
         TCC(I,K)=TINI(K)
         DXT(I,K)=0.D+000
         DYT(I,K)=0.D+000
         D2T(I,K)=0.D+000
        END DO
       END DO
C
      END IF
C
      CLOSE (13)
C
C----------------------------------------------------------------------
C
C     INTERFACE ENERGY BOUNDARY CONDITIONS
C
      DO K=1,NR
C
C      INTERFACE BOUNDARY CONDITIONS
C      KE<-1000:INTERFACE TEMPERATURE, KE<-2000:INTERFACE HEAT FLUX
C
       DO I=1,NB(K)
C       TYPE 1 INTERFACE: IMPOSE TEMPERATURE
        IF ((KE(I,K).LT.(-1000)).AND.(KE(I,K).GT.(-2000))) THEN
         GT(I,1,K)=1.D+000
         GT(I,2,K)=0.D+000
         GT(I,3,K)=TCC(I,K)
        END IF
C       TYPE 2 INTERFACE: IMPOSE HEAT FLUX
        IF ((KE(I,K).LT.(-2000)).AND.(KE(I,K).GT.(-3000))) THEN
         GT(I,1,K)=0.D+000
         GT(I,2,K)=-TCE(K)
         GT(I,3,K)=-TCE(K)*(DXT(I,K)*XN(I,K)+DYT(I,K)*YN(I,K))
        END IF
C
       END DO
C
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END