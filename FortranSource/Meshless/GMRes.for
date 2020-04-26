CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC   SUBROUTINE GMRES                                                      CCC
CCC                                                                         CCC
CCC   SOLVES EQUATION: [A]{X}={B}                                           CCC
CCC                                                                         CCC
CCC   USING A GMRES ALGORITHM                                               CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GMRES(AM,BM,XM,RESH,MXNT,K)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NLOOPMAX=100)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Topology.for'
C
C----------------------------------------------------------------------
C
      REAL*8 AM(NCMAX,NFMAX),BM(NCMAX),XM(NCMAX)
      REAL*8 RESH
      INTEGER MXNT,K
C
      REAL*8 R(NCMAX),V(NCMAX,NLOOPMAX+1)
      REAL*8 H(NLOOPMAX+1,NLOOPMAX+1),G(NLOOPMAX+1)
      REAL*8 COSTH(NLOOPMAX+1),SINTH(NLOOPMAX+1)
C
C----------------------------------------------------------------------
C
C     LOOP SIZE
C
      IF (MXNT.GT.NLOOPMAX) MXNT=NLOOPMAX
C
C     APPLY PRECONDITIONING TO [A] AND {B}
C
      DO I=1,NC(K)
       HT1=AM(I,1)
       DO IJ=1,NCONN(I,K)
        AM(I,IJ)=AM(I,IJ)/HT1
       END DO
       BM(I)=BM(I)/HT1
      END DO
C
C     COMPUTE INITIAL RESIDUAL VECTOR AND ITS NORM
C
      RES=0.D+000
      DO I=1,NC(K)
       R(I)=0.D+000
       DO IJ=1,NCONN(I,K)
        R(I)=R(I)+AM(I,IJ)*XM(ICONN(I,IJ,K))
       END DO
       R(I)=BM(I)-R(I)
       RES=RES+R(I)*R(I)
      END DO
      IF (RES.GT.EPS) THEN
       RES=DSQRT(RES)
      ELSE
       RES=EPS
      END IF
C
C----------------------------------------------------------------------
C
      IF (RES.GT.EPS) THEN
C
C     INITIALIZE Gj AND {Vj}
C
      JJ=1
      G(JJ)=RES
      DO I=1,NC(K)
       V(I,JJ)=R(I)/RES
      END DO
C
C     BEGIN MAIN LOOP
C
      DO JJ=1,MXNT
C
C      COMPUTE [Hij]=TRANS([A]{Vj}){Vi}
C
       DO II=1,JJ
        H(II,JJ)=0.D+000
        DO I=1,NC(K)
         HT1=0.D+000
         DO IJ=1,NCONN(I,K)
          HT1=HT1+AM(I,IJ)*V(ICONN(I,IJ,K),JJ)
         END DO
         H(II,JJ)=H(II,JJ)+HT1*V(I,II)
        END DO
       END DO
C
C      SUM {Vj+1}={Vj+1}-[Hij]{Vj}
C
       DO I=1,NC(K)
        V(I,JJ+1)=0.D+000
        DO II=1,JJ
         V(I,JJ+1)=V(I,JJ+1)-H(II,JJ)*V(I,II)
        END DO
       END DO
C
C      COMPUTE {Vj+1}={Vj+1}+[A]{Vj}
C         
       DO I=1,NC(K)
        DO IJ=1,NCONN(I,K)
         V(I,JJ+1)=V(I,JJ+1)+AM(I,IJ)*V(ICONN(I,IJ,K),JJ)
        END DO
       END DO
C
C      COMPUTE [Hj+1,j]=||{Vj+1}||
C      
       H(JJ+1,JJ)=0.D+000
       DO I=1,NC(K)
        H(JJ+1,JJ)=H(JJ+1,JJ)+V(I,JJ+1)*V(I,JJ+1)
       END DO
       H(JJ+1,JJ)=DSQRT(H(JJ+1,JJ))
C
C      COMPUTE {Vj+1}={Vj+1}/||{Vj+1}||
C      
       DO I=1,NC(K)
        V(I,JJ+1)=V(I,JJ+1)/H(JJ+1,JJ)
       END DO
C
C      APPLY ROTATIONS L=1,....,JJ-1 TO Jth COLUMN OF [H]
C
       DO II=1,JJ-1
        HT1=COSTH(II)*H(II,JJ)-SINTH(II)*H(II+1,JJ)
        HT2=SINTH(II)*H(II,JJ)+COSTH(II)*H(II+1,JJ)
        H(II,JJ)  =HT1
        H(II+1,JJ)=HT2
       END DO
C
C      CALCULATE CURRENT ROTATIONS
C
       COSTH(JJ)=   H(JJ,JJ)
     &          /DSQRT(H(JJ,JJ)*H(JJ,JJ)+H(JJ+1,JJ)*H(JJ+1,JJ))
       SINTH(JJ)=-H(JJ+1,JJ)
     &          /DSQRT(H(JJ,JJ)*H(JJ,JJ)+H(JJ+1,JJ)*H(JJ+1,JJ))
C
C      APPLY ROTATIONS TO [Hj,j] AND [Hj+1,j]
C
       HT1=COSTH(JJ)*H(JJ,JJ)-SINTH(JJ)*H(JJ+1,JJ)
       HT2=SINTH(JJ)*H(JJ,JJ)+COSTH(JJ)*H(JJ+1,JJ)
       H(JJ,JJ)  =HT1
       H(JJ+1,JJ)=HT2
C
C      APPLY ROTATIONS TO Gj AND Gj+1
C
       G(JJ+1)=0.D+000
       HT1=COSTH(JJ)*G(JJ)-SINTH(JJ)*G(JJ+1)
       HT2=SINTH(JJ)*G(JJ)+COSTH(JJ)*G(JJ+1)
       G(JJ)  =HT1
       G(JJ+1)=HT2
C
C      APPROXIMATE RESIDUAL
C
       RES=DABS(G(JJ+1))
C
      END DO
C
C     SOLVE [H]{Y}={G} BY BACK SUBSTITUTION
C
      DO II=MXNT,1,-1   
       SINTH(II)=0.D+000
       DO KK=II+1,MXNT
        SINTH(II)=SINTH(II)+H(II,KK)*COSTH(KK)
       END DO
       COSTH(II)=(G(II)-SINTH(II))/H(II,II)
      END DO
C
C     FORM NEW SOLUTION VECTOR {X}={X}+([C]^-1)[V]{Y}
C
      DO I=1,NC(K)
       R(I)=0.D+000
       DO II=1,MXNT
        R(I)=R(I)+V(I,II)*COSTH(II)
       END DO
       XM(I)=XM(I)+R(I)
      END DO
C
C     COMPUTE RESIDUAL
C
      RESH=0.D+000
      RHSH=0.D+000
      DO I=1,NC(K)
       DHSH=AM(I,1)*XM(I)
       DO II=2,NCONN(I,K)
        DHSH=DHSH+AM(I,II)*XM(ICONN(I,II,K))
       END DO
       RESH=RESH+(DHSH-BM(I))*(DHSH-BM(I))
       RHSH=RHSH+BM(I)*BM(I)
      END DO
C
C     NON-DIMENSIONALIZE RESIDUAL
C
      IF (RHSH.GT.EPS) THEN
       RESH=DSQRT(RESH/RHSH)
      ELSE
       RESH=EPS
      END IF
C
C----------------------------------------------------------------------
C
      ELSE
C
C     DROP RESIDUAL
C
       RESH=EPS
C
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END