CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE INPUT                                  CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   READS INPUT, DATA, VECTOR, AND TRIANGULATION                          CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MESHLESS_INPUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../../Include/Parameters.for'
      INCLUDE '../../Include/Information.for'
      INCLUDE '../../Include/Geometry.for'
      INCLUDE '../../Include/Materials.for'
      INCLUDE '../../Include/Topology.for'
      INCLUDE '../../Include/Interpolation-C.for'
      INCLUDE '../../Include/Interpolation-U.for'
      INCLUDE '../../Include/Triangulation.for'
C
C----------------------------------------------------------------------
C
      CHARACTER*120 TITLE
      INTEGER IOS
C
C----------------------------------------------------------------------
C
C     READ ITERATION PARAMETERS
C
      OPEN (12,FILE='Model/ALMA_i_inf.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.EQ.0) THEN
C
C      SOLVE MOMENTUM, ENERGY, STRUCTURAL, PORE PRESSURE, LEVEL-SET
C
       READ (12,*) TITLE
       READ (12,*) MSM,MSE,MSS,MSP,MSL
C
C      TIME STEP, NUMBER OF STEPS, OUTPUT FREQUENCY, RESIDUAL FREQUENCY, SUB-LEVEL ITERATIONS
C
       READ (12,*) TITLE
       READ (12,*) DT,MAXITER,IWRITE,IRES,ISUB
C
C      RELAXATION FOR: POTENTIAL, UPWIND, INTERFACE
C
       READ (12,*) TITLE
       READ (12,*) THP,THU,THI
C
C      GRAVITY
C
       READ (12,*) TITLE
       READ (12,*) GX,GY
C
      ELSE
C
       WRITE (*,*) '...MISSING INFORMATION FILE.'
       STOP
C
      END IF
C
      CLOSE (12)
C
C----------------------------------------------------------------------
C
C     VALIDATE ITERATION PARAMETERS
C
      IF (DT.LT.EPS) THEN
       WRITE (*,*)
     & '...TIME STEP IS TOO SMALL.....................<STOP>'
       STOP
      END IF
C
      IF (IWRITE.LT.1) IWRITE=MAXITER
      IF (IRES.LT.1) IRES=1
      IF (IRES.GT.IWRITE) IRES=IWRITE
      IF (ISUB.LT.1) ISUB=1
C
      IF (THP.LT.0.0D+000) THP=0.0D+000
      IF (THP.GE.1.0D+000) THP=1.0D+000
C
      IF (THU.LT.0.0D+000) THU=0.0D+000
      IF (THU.GT.1.0D+000) THU=1.0D+000
C
      IF (THI.LT.0.0D+000) THI=0.0D+000
      IF (THI.GT.1.0D+000) THI=1.0D+000
C
C----------------------------------------------------------------------
C
C     READ GEOMETRY 
C
      OPEN (11,FILE='Model/ALMA_d_geo.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.EQ.0) THEN
C
       READ (11,*) TITLE
C
C      NUMBER OF REGIONS
C
       READ (11,*) NR
C
       DO K=1,NR
C
C       NUMBER OF BOUNDARY POINTS AND INTERNAL POINTS
C
        READ (11,*) KK,NB(K),NI(K)
C
C       BOUNDARY GEOMETRY
C
        DO I=1,NB(K)
         READ (11,*) II,XC(I,K),YC(I,K),AR(I,K),XN(I,K),YN(I,K)
        END DO
C
C       INTERNAL POINTS
C
        NC(K)=NB(K)+NI(K)
        DO I=NB(K)+1,NC(K)
         READ (11,*) II,XC(I,K),YC(I,K)
        END DO
C
       END DO
C
      ELSE
C
       WRITE (*,*) '...MISSING GEOMETRIC DATA FILE.'
       STOP
C
      END IF
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     READ MATERIALS FILE
C
      OPEN (11,FILE='Model/ALMA_i_mat.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.NE.0) THEN
C
       WRITE (*,*) '...MISSING MATERIALS FILE.'
       STOP
C
      ELSE
C
       DO K=1,NR
C
C       READ REGION NUMBER AND MATERIAL TYPE
C
        READ (11,*) KK,KR(K)
C
C       KR = 0 : SOLID
C
        IF (KR(K).EQ.0) THEN
C
         READ (11,*) TITLE
         READ (11,*) DES(K)
         READ (11,*) VIS(K),POS(K)
         READ (11,*) TCS(K),SHS(K)
         READ (11,*) BES(K),TRS(K)
C
         DEE(K)=DES(K)
         VIE(K)=VIS(K)
         POE(K)=POS(K)
         TCE(K)=TCS(K)
         SHE(K)=SHS(K)
         BEE(K)=BES(K)
         TRE(K)=TRS(K)
C
        END IF
C
C       KR = 1 : FLUID
C
        IF (KR(K).EQ.1) THEN
C
         READ (11,*) TITLE
         READ (11,*) DEF(K)
         READ (11,*) VIF(K),POF(K)
         READ (11,*) TCF(K),SHF(K)
         READ (11,*) BEF(K),TRF(K)
C
         DEE(K)=DEF(K)
         VIE(K)=VIF(K)
         POE(K)=POF(K)
         TCE(K)=TCF(K)
         SHE(K)=SHF(K)
         BEE(K)=BEF(K)
         TRE(K)=TRF(K)
C
        END IF
C
C       KR = 2 : POROUS MEDIUM  
C
        IF (KR(K).EQ.2) THEN
C
         READ (11,*) TITLE
         READ (11,*) DES(K)
         READ (11,*) VIS(K),POS(K)
         READ (11,*) TCS(K),SHS(K)
         READ (11,*) BES(K),TRS(K)
         READ (11,*) PER(K),POR(K)
C
         READ (11,*) TITLE
         READ (11,*) DEF(K)
         READ (11,*) VIF(K),POF(K)
         READ (11,*) TCF(K),SHF(K)
         READ (11,*) BEF(K),TRF(K)
         READ (11,*) COF(K)
C
         DEE(K)=DES(K)
         VIE(K)=VIS(K)
         POE(K)=POS(K)
         TCE(K)=TCS(K)
         SHE(K)=SHS(K)
         BEE(K)=BES(K)
         TRE(K)=TRS(K)
C
        END IF
C
       END DO
C
      END IF
C
      CLOSE (11)
C
C----------------------------------------------------------------------
C
C     READ INTERPOLATION VECTORS
C
      OPEN (21,FILE='Model/ALMA_d_vec.bin',FORM='UNFORMATTED'
     &        ,STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.EQ.0) THEN
C
       READ (21) NRC
       IF (NRC.EQ.NR) THEN
        DO K=1,NR
         READ (21) NCC
         IF (NCC.EQ.NC(K)) THEN
          DO I=1,NC(K)
           READ (21) NCONN(I,K)
           READ (21) RXA(I,K),RYA(I,K)
           DO II=1,NCONN(I,K)
            READ (21) ICONN(I,II,K),
     &                FXC(I,II,K),FYC(I,II,K), 
     &                FXX(I,II,K),FYY(I,II,K),FXY(I,II,K),
     &                FXE(I,II,K),FXW(I,II,K),FYN(I,II,K),FYS(I,II,K),
     &                SXE(I,II,K),SXW(I,II,K),SYN(I,II,K),SYS(I,II,K)
           END DO
          END DO
         ELSE
          IOS=1
         END IF
        END DO
       ELSE
        IOS=1
       END IF
      END IF
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING OR CORRUPTED INTERPOLATION VECTOR FILE.'
       STOP
      END IF
C
      CLOSE (21)
C
C----------------------------------------------------------------------
C
C     READ TRIANGULATION FROM FILE
C
      OPEN (35,FILE='Model/ALMA_d_tri.txt',STATUS='OLD',IOSTAT=IOS)
C
      IF (IOS.EQ.0) THEN
C
      READ (35,*) NRR
       IF (NRR.NE.NR) IOS=1
       DO K=1,NRR
        READ (35,*) NCC,NMESH(K)
        IF (NCC.NE.NC(K)) IOS=1
        DO I=1,NMESH(K)
         READ (35,*) (MESH(I,II,K),II=1,3)
        END DO
       END DO
C
      END IF
C
      IF (IOS.NE.0) THEN
       WRITE (*,*) '...MISSING OR CORRUPTED TRIANGULATION FILE.'
       DO K=1,NR
        NMESH(K)=0
       END DO
      END IF
C
      CLOSE (35)
C
C----------------------------------------------------------------------
C
C     OUTPUT INFORMATION
C
C      OPEN (12,FILE='Model/ALMA_o_inf.txt',STATUS='OLD',IOSTAT=IOS)
C
C      IF (IOS.EQ.0) THEN
C
C       READ (12,*) ITP,CTIME
C
C      ELSE
C
       ITP=0
       CTIME=0.D+000
C
C      END IF
C
C      CLOSE (12)
C
C----------------------------------------------------------------------
C
C     INCREASE MAXIMUM ITERATIONS TO ACCOUNT FOR PREVIOUS RESULTS
C
      MAXITER=MAXITER+ITP
C
C----------------------------------------------------------------------
C
      RETURN
      END