CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE SOLVELEVELSET                          CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   SOLVES LEVEL-SET TRANSPORT EQUATION                                   CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLVELEVELSET
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Field-S.for'
      INCLUDE '../../Include/Field-M.for'
      INCLUDE '../../Include/Field-P.for'
C
C----------------------------------------------------------------------
C
      REAL*8 DIFF(NCMAX)
      REAL*8 RHS(NCMAX)
      REAL*8 QFS(NCMAX)
C
C----------------------------------------------------------------------
C
C     INITIALIZE RESIDUAL
C
      RESL=0.D+000
C
C----------------------------------------------------------------------
C
C      LOOP OVER SUBREGIONS FOR FIELD SOLUTION
C
      DO K=1,NR
C
C      DIFFUSIVITY
C
       DO I=1,NC(K)
        DIFF(I)=EPS
       END DO
C
C      BOUNDARY CONDITIONS
C
       DO I=1,NB(K)
        RHS(I)=GS(I,3,K)
       END DO
C
C      BODY FORCE
C
       DO I=NB(K)+1,NC(K)
        RHS(I)=0.D+000
       END DO
C
C      CONVECTIVE FLUX
C
       CALL UPWIND(VFX,VFY,SCC,QFS,DIFF,K)
C
C      TRANSPORT LEVEL-SET
C
       CALL TRANSPORT(SCC,D2S,DXS,DYS,QFS,DIFF,RHS,GS,DT,RESS,K)
C
C      ACCUMULATE LEVEL-SET RESIDUAL
C
       RESL=RESL+RESS
C
      END DO
C
C----------------------------------------------------------------------
C
      RETURN
      END