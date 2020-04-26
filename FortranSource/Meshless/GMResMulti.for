CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC   SUBROUTINE GMRESMULTI                                                 CCC
CCC                                                                         CCC
CCC   SOLVES EQUATION: [A]{X}={B} IN MULTIDIMENSIONS                        CCC
CCC                                                                         CCC
CCC   USING A GMRES ALGORITHM                                               CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GMRESMULTI(AM,BM,XM,RESH,MXNT,K)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NLOOPMAX=200)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Topology.for'
C
C----------------------------------------------------------------------
C
      REAL*8 AM(NDIM*NCMAX,NDIM*NFMAX),BM(NDIM*NCMAX),XM(NDIM*NCMAX)
      REAL*8 RESH(NDIM)
      INTEGER MXNT,K
C
      REAL*8 R(NDIM*NCMAX),V(NDIM*NCMAX,NLOOPMAX+1)
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
       DO ID=1,NDIM
        IDIM=NDIM*(I-1)+ID
        HT1=AM(IDIM,ID)
        DO IJ=1,NCONN(I,K)
         DO JD=1,NDIM
          JDIM=NDIM*(IJ-1)+JD
          AM(IDIM,JDIM)=AM(IDIM,JDIM)/HT1
         END DO
        END DO
        BM(IDIM)=BM(IDIM)/HT1
       END DO
      END DO
C
C     COMPUTE INITIAL RESIDUAL VECTOR AND ITS NORM
C
      RES=0.D+000
      DO I=1,NC(K)
       DO ID=1,NDIM
        IDIM=NDIM*(I-1)+ID
        R(IDIM)=0.D+000
        DO IJ=1,NCONN(I,K)
         DO JD=1,NDIM
          JDIM=NDIM*(IJ-1)+JD
          IJDIM=NDIM*(ICONN(I,IJ,K)-1)+JD
          R(IDIM)=R(IDIM)+AM(IDIM,JDIM)*XM(IJDIM)
         END DO
        END DO
        R(IDIM)=BM(IDIM)-R(IDIM)
        RES=RES+R(IDIM)*R(IDIM)
       END DO
      END DO
C
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
       DO ID=1,NDIM
        IDIM=NDIM*(I-1)+ID
        V(IDIM,JJ)=R(IDIM)/RES
       END DO
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
         DO ID=1,NDIM
          IDIM=NDIM*(I-1)+ID
          HT1=0.D+000
          DO IJ=1,NCONN(I,K)
           DO JD=1,NDIM
            JDIM=NDIM*(IJ-1)+JD
            IJDIM=NDIM*(ICONN(I,IJ,K)-1)+JD
            HT1=HT1+AM(IDIM,JDIM)*V(IJDIM,JJ)
           END DO
          END DO
          H(II,JJ)=H(II,JJ)+HT1*V(IDIM,II)
         END DO
        END DO
       END DO
C
C      SUM {Vj+1}={Vj+1}-[Hij]{Vj}
C
       DO I=1,NC(K)
        DO ID=1,NDIM
         IDIM=NDIM*(I-1)+ID
         V(IDIM,JJ+1)=0.D+000
         DO II=1,JJ
          V(IDIM,JJ+1)=V(IDIM,JJ+1)-H(II,JJ)*V(IDIM,II)
         END DO
        END DO
       END DO
C
C      COMPUTE {Vj+1}={Vj+1}+[A]{Vj}
C         
       DO I=1,NC(K)
        DO ID=1,NDIM
         IDIM=NDIM*(I-1)+ID
         DO IJ=1,NCONN(I,K)
          DO JD=1,NDIM
           JDIM=NDIM*(IJ-1)+JD
           IJDIM=NDIM*(ICONN(I,IJ,K)-1)+JD
           V(IDIM,JJ+1)=V(IDIM,JJ+1)+AM(IDIM,JDIM)*V(IJDIM,JJ)
          END DO
         END DO
        END DO
       END DO
C
C      COMPUTE [Hj+1,j]=||{Vj+1}||
C      
       H(JJ+1,JJ)=0.D+000
       DO I=1,NC(K)
        DO ID=1,NDIM
         IDIM=NDIM*(I-1)+ID
         H(JJ+1,JJ)=H(JJ+1,JJ)+V(IDIM,JJ+1)*V(IDIM,JJ+1)
        END DO
       END DO
       H(JJ+1,JJ)=DSQRT(H(JJ+1,JJ))
C
C      COMPUTE {Vj+1}={Vj+1}/||{Vj+1}||
C      
       DO I=1,NC(K)
        DO ID=1,NDIM
         IDIM=NDIM*(I-1)+ID
         V(IDIM,JJ+1)=V(IDIM,JJ+1)/H(JJ+1,JJ)
        END DO
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
       DO ID=1,NDIM
        IDIM=NDIM*(I-1)+ID
        R(IDIM)=0.D+000
        DO II=1,MXNT
         R(IDIM)=R(IDIM)+V(IDIM,II)*COSTH(II)
        END DO
        XM(IDIM)=XM(IDIM)+R(IDIM)
       END DO
      END DO
C
C     COMPUTE RESIDUAL
C
      DO ID=1,NDIM
       RESH(ID)=0.D+000
      END DO
      RHSH=0.D+000
C
      DO I=1,NC(K)
       DO ID=1,NDIM
        IDIM=NDIM*(I-1)+ID
        DHSH=0.D+000
        DO IJ=1,NCONN(I,K)
         DO JD=1,NDIM
          JDIM=NDIM*(IJ-1)+JD
          IJDIM=NDIM*(ICONN(I,IJ,K)-1)+JD
          DHSH=DHSH+AM(IDIM,JDIM)*XM(IJDIM)
         END DO
        END DO
        RESH(ID)=RESH(ID)+(DHSH-BM(IDIM))*(DHSH-BM(IDIM))
        RHSH=RHSH+BM(IDIM)*BM(IDIM)
       END DO
      END DO
C
C     NON-DIMENSIONALIZE RESIDUAL
C
      DO ID=1,NDIM
       IF (RHSH.GT.EPS) THEN
        RESH(ID)=DSQRT(RESH(ID)/RHSH)
       ELSE
        RESH(ID)=EPS
       END IF
      END DO
C
C----------------------------------------------------------------------
C
      ELSE
C
C     DROP RESIDUAL
C
      DO ID=1,NDIM
       RESH(ID)=EPS
      END DO
C
      END IF
C
C----------------------------------------------------------------------
C
      RETURN
      END