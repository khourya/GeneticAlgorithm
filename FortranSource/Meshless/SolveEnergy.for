CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVEENERGY                            CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES ENERGY TRANSPORT EQUATION                                      CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVEENERGY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Field-S.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-E.for'
C
C----------------------------------------------------------------------
C
C     USE THE FOLLOWING SPACE TO INCLUDE THE COMMON BLOCKS FOR FUNCTION VARIABLES
C
      PARAMETER (NTMAX=20,NMMAX=200)
	COMMON/MEASURE/NMN,NMT,MT(NTMAX),MR(NMMAX),MN(NMMAX)
      COMMON/MEASURT/TM(NMMAX,NTMAX)
      COMMON/HOLELOC/XCH,YCH,RXH,RYH
C
C----------------------------------------------------------------------
C
      REAL*8 DIFF(NCMAX)
      REAL*8 RHS(NCMAX)
      REAL*8 QFT(NCMAX)
C
C----------------------------------------------------------------------
C
C     INITIALIZE RESIDUAL
C
      RESE=0.D+000
C
C----------------------------------------------------------------------
C
C     LOOP OVER SUBREGIONS FOR FIELD SOLUTION
C
      DO K=1,NR
C
C      DIFFUSIVITY
C
       DO I=1,NC(K)
        DIFF(I)=TCE(K)/(DEE(K)*SHE(K))
CCC
CCC
CCC
        ELX=(XC(I,K)-XCH)/(1.33D+000*RXH)
        ELY=(YC(I,K)-YCH)/(1.33D+000*RYH)
        IF ((K.EQ.1).AND.((ELX*ELX+ELY*ELY).LT.1.D+000)) THEN
         DIFF(I)=30.D+000*DIFF(I)
C         WRITE (*,*) '*'
        END IF
CCC
CCC
CCC
       END DO
C
C      BOUNDARY CONDITIONS
C
       DO I=1,NB(K)
        RHS(I)=GT(I,3,K)
       END DO
C
C      BODY FORCE
C
       DO I=NB(K)+1,NC(K)
        RHS(I)=UBG(K)/(DEE(K)*SHE(K))
       END DO
C
C      CONVECTIVE ENERGY FLUX
C
       IF (KR(K).EQ.0) THEN
        DO I=1,NC(K)
         QFT(I)=0.D+000
        END DO
       END IF
C
cc       IF (KR(K).EQ.1) THEN
cc        CALL UPWIND(UCC,VCC,TCC,QFT,DIFF,K)
cc       END IF
C
cc       IF (KR(K).EQ.2) THEN
cc        CALL UPWIND(VFX,VFY,TCC,QFT,DIFF,K)
cc        DO I=1,NC(K)
cc         QFT(I)=DEF(K)*SHF(K)*QFT(I)/(DEE(K)*SHE(K))
cc        END DO
cc       END IF
C
C      TRANSPORT ENERGY
C
       CALL TRANSPORT(TCC,D2T,DXT,DYT,QFT,DIFF,RHS,GT,DT,REST,K)
C
C      ACCUMULATE ENERGY RESIDUAL
C
cc       RESE=RESE+REST
C
C      END SUBREGION LOOP
C
      END DO
C
C----------------------------------------------------------------------
C
C      AVERAGE INTERFACE ENERGY
C
      IF (NR.GT.1) THEN
       CALL INTERFACEENERGY
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END