CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                         CCC
CCC                       SUBROUTINE OUTPUT                                 CCC
CCC                                                                         CCC
CCC   PURPOSE:                                                              CCC
CCC                                                                         CCC
CCC   WRITES SOLUTION, BINARY AND TECPLOT FILES                             CCC
CCC                                                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MESHLESS_OUTPUT
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
      INCLUDE '../../Include/Field-P.for'
      INCLUDE '../../Include/Triangulation.for'
C
C----------------------------------------------------------------------
C
      CHARACTER*200 TITLE
C
C----------------------------------------------------------------------
C
C     OUTPUT TO TEXT, BINARY, AND TECPLOT FILES
C
      N1=(ITP)/100000
      N2=(ITP-100000*N1)/10000
      N3=(ITP-100000*N1-10000*N2)/1000
      N4=(ITP-100000*N1-10000*N2-1000*N3)/100
      N5=(ITP-100000*N1-10000*N2-1000*N3-100*N4)/10
      N6=(ITP-100000*N1-10000*N2-1000*N3-100*N4-10*N5)
C
C----------------------------------------------------------------------
C
C     OUTPUT TIME INFO
C
C      OPEN (12,FILE='Model/ALMA_o_inf.txt')
C
C      WRITE (12,'(I7,2X,E12.6)') ITP,CTIME
C
C      CLOSE (12)
C
C----------------------------------------------------------------------
C
C     OUTPUT MOMENTUM
C
      IF (MSM.EQ.1) THEN
C
       TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                       //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                       //'_mom.bin'
C
       OPEN (13,FILE=TITLE,FORM='UNFORMATTED')
       OPEN (14,FILE='Model/ALMA_o_mom.bin',FORM='UNFORMATTED')
C
       DO K=1,NR
        DO I=1,NC(K)
         WRITE (13) UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K) 
     &             ,PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K) 
         WRITE (14) UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K) 
     &             ,PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K) 
        END DO
       END DO
C
       CLOSE (13)
       CLOSE (14)
C
      END IF
C
C----------------------------------------------------------------------
C
C     OUTPUT STRUCTURAL
C
      IF (MSS.EQ.1) THEN
C
       TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                       //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                       //'_str.bin'
C
       OPEN (13,FILE=TITLE,FORM='UNFORMATTED')
       OPEN (14,FILE='Model/ALMA_o_str.bin',FORM='UNFORMATTED')
C
       DO K=1,NR
        DO I=1,NC(K)
         WRITE (13) XC(I,K),UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,YC(I,K),VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K) 
         WRITE (14) XC(I,K),UCC(I,K),DXU(I,K),DYU(I,K),D2U(I,K)
     &             ,YC(I,K),VCC(I,K),DXV(I,K),DYV(I,K),D2V(I,K) 
        END DO
       END DO
C
       CLOSE (13)
       CLOSE (14)
C
      END IF
C
C----------------------------------------------------------------------
C
C     OUTPUT ENERGY
C
      IF (MSE.EQ.1) THEN
C
       TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                       //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                       //'_ene.bin'
C
c       OPEN (13,FILE=TITLE,FORM='UNFORMATTED')
       OPEN (14,FILE='Model/ALMA_o_ene.bin',FORM='UNFORMATTED')
C
       DO K=1,NR
        DO I=1,NC(K)
c         WRITE (13) TCC(I,K),DXT(I,K),DYT(I,K),D2T(I,K)
         WRITE (14) TCC(I,K),DXT(I,K),DYT(I,K),D2T(I,K)
        END DO
       END DO
C
c       CLOSE (13)
       CLOSE (14)
C
      END IF
C
C----------------------------------------------------------------------
C
C     OUTPUT PORE PRESSURE
C
      IF (MSP.EQ.1) THEN
C
       TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                       //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                       //'_por.bin'
C
       OPEN (13,FILE=TITLE,FORM='UNFORMATTED')
       OPEN (14,FILE='Model/ALMA_o_por.bin',FORM='UNFORMATTED')
C
       DO K=1,NR
        DO I=1,NC(K)
         WRITE (13) PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K)
         WRITE (14) PCC(I,K),DXP(I,K),DYP(I,K),D2P(I,K)
        END DO
       END DO
C
       CLOSE (13)
       CLOSE (14)
C
      END IF
C
C----------------------------------------------------------------------
C
C     OUTPUT LEVEL-SET
C
      IF (MSL.EQ.1) THEN
C
       TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                       //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                       //'_lev.bin'
C
       OPEN (13,FILE=TITLE,FORM='UNFORMATTED')
       OPEN (14,FILE='Model/ALMA_o_lev.bin',FORM='UNFORMATTED')
C
       DO K=1,NR
        DO I=1,NC(K)
         WRITE (13) SCC(I,K),DXS(I,K),DYS(I,K),D2S(I,K)
         WRITE (14) SCC(I,K),DXS(I,K),DYS(I,K),D2S(I,K)
        END DO
       END DO
C
       CLOSE (13)
       CLOSE (14)
C
      END IF
C
C----------------------------------------------------------------------
C
C     TECPLOT OUTPUT
C
      TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                      //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                      //'_tec.plt'
c      OPEN (17,FILE=TITLE)
      OPEN (18,FILE='Model/ALMA_o_tec.plt')
C
      TITLE='VARIABLES = "X","Y","Nx","Ny","L"'
     &    //',"T","Qx","Qy","|Q|"'
     &    //',"P","Vx","Vy","|V|"'
     &    //',"Ux","Uy","|U|"'
     &    //',"Sxx","Sxy","Syy","|S|"'
C
c      WRITE (17,'(A200)') TITLE
      WRITE (18,'(A200)') TITLE
C
C     LOOP OVER SUBREGIONS
C
      DO K=1,NR
c       WRITE (17,*)
       WRITE (18,*) 
C
C      TECPLOT HEADER
C
c       WRITE (17,'("ZONE T=",A,"r:",I2," t:",E12.6,A,",N=",I5,
c     &             ",E=",I5,",F=FEPOINT,ET=TRIANGLE")') 
c     &             CHAR(34),K,CTIME,CHAR(34),NC(K),NMESH(K)
       WRITE (18,'("ZONE T=",A,"r:",I2," t:",E12.6,A,",N=",I5,
     &             ",E=",I5,",F=FEPOINT,ET=TRIANGLE")') 
     &             CHAR(34),K,CTIME,CHAR(34),NC(K),NMESH(K)
C
C      COEFFICIENTS
C
       ELS=1.D+000-2.D+000*POE(K)
       TEX=2.D+000*VIE(K)*BEE(K)*(1.D+000+POE(K))/ELS
C
C      WRITE SOLUTION FIELD
C
       DO I=1,NC(K)
        XCC=XC(I,K)
        YCC=YC(I,K)
        XNN=0.D+000
        YNN=0.D+000
        IF (I.LE.NB(K)) THEN
         XNN=XN(I,K)
         YNN=YN(I,K)
        END IF
        SC=SCC(I,K)
        TC=TCC(I,K)
        QX=-TCE(K)*DXT(I,K)
        QY=-TCE(K)*DYT(I,K)
        QM=QX*QX+QY*QY
        IF (QM.GT.EPS) QM=DSQRT(QM)
        PC=PCC(I,K)
        IF (KR(K).EQ.0) THEN
         VX=0.D+000
         VY=0.D+000
         VM=0.D+000
        END IF
        IF (KR(K).EQ.1) THEN
         VX=UCC(I,K)
         VY=VCC(I,K)
         VM=VX*VX+VY*VY
         IF (VM.GT.EPS) VM=DSQRT(VM)
        END IF
        IF (KR(K).EQ.2) THEN
         VX=VFX(I,K)
         VY=VFY(I,K)
         VM=VX*VX+VY*VY
         IF (VM.GT.EPS) VM=DSQRT(VM)
        END IF
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         UX=UCC(I,K)
         UY=VCC(I,K)
         UM=UX*UX+UY*UY
         IF (UM.GT.EPS) UM=DSQRT(UM)
         SXX=2.D+000*VIE(K)*DXU(I,K)
     &      +2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SYY=2.D+000*VIE(K)*DYV(I,K)
     &      +2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SZZ=2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SXY=VIE(K)*(DYU(I,K)+DXV(I,K))
         SXZ=0.D+000
         SYZ=0.D+000
         SMM=0.5D+000*(SXX*SXX+SYY*SYY+SZZ*SZZ)+SXY*SXY+SXZ*SXZ+SYZ*SYZ
         IF (SMM.GT.EPS) SMM=DSQRT(SMM)
        END IF
        IF (KR(K).EQ.1) THEN
         UX=0.D+000
         UY=0.D+000
         UM=0.D+000
         SXX=2.D+000*VIE(K)*DXU(I,K)
         SYY=2.D+000*VIE(K)*DYV(I,K)
         SZZ=0.D+000
         SXY=VIE(K)*(DYU(I,K)+DXV(I,K))
         SXZ=0.D+000
         SYZ=0.D+000
         SMM=0.5D+000*(SXX*SXX+SYY*SYY+SZZ*SZZ)+SXY*SXY+SXZ*SXZ+SYZ*SYZ
         IF (SMM.GT.EPS) SMM=DSQRT(SMM)
        END IF
C
c        WRITE (17,103) XC(I,K),YC(I,K),XNN,YNN,SC,TC,QX,QY,QM
c     &                ,PC,VX,VY,VM,UX,UY,UM,SXX,SXY,SYY,SMM   
        WRITE (18,103) XC(I,K),YC(I,K),XNN,YNN,SC,TC,QX,QY,QM
     &                ,PC,VX,VY,VM,UX,UY,UM,SXX,SXY,SYY,SMM   
       END DO
C
C      WRITE TRIANGULATION CONNECTIVITY
C
c       WRITE (17,*)
       WRITE (18,*)
       DO I=1,NMESH(K)
c        WRITE (17,*) MESH(I,1,K),MESH(I,2,K),MESH(I,3,K)
        WRITE (18,*) MESH(I,1,K),MESH(I,2,K),MESH(I,3,K)
       END DO
C
      END DO
C
c      CLOSE (17)
      CLOSE (18)
C
  103 FORMAT(10(E12.4,1X),5(1X,4(E12.4,1X)))
C
C----------------------------------------------------------------------
C
C     PARAVIEW OUTPUT
C
      TITLE='Model/ALMA_ox_'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)
     &                      //CHAR(N4+48)//CHAR(N5+48)//CHAR(N6+48)
     &                      //'_pvw.vtk'
c      OPEN (19,FILE=TITLE)
      OPEN (20,FILE='Model/ALMA_o_pvw.vtk')
C
C     PARAVIEW HEADER
C
c      WRITE (19,'("# vtk DataFile Version 2.0")') 
      WRITE (20,'("# vtk DataFile Version 2.0")') 
c      WRITE (19,'("ALMA 2D Output. Time =",X,E12.6)') CTIME
      WRITE (20,'("ALMA 2D Output. Time =",X,E12.6)') CTIME
c      WRITE (19,'("ASCII")') 
      WRITE (20,'("ASCII")') 
c      WRITE (19,'("DATASET UNSTRUCTURED_GRID")')            
      WRITE (20,'("DATASET UNSTRUCTURED_GRID")')            
C
C     COUNT NUMBER OF POINTS AND CELLS
C
      NPK=0
      NCK=0
      DO K=1,NR
       NPK=NPK+NC(K)
       NCK=NCK+NMESH(K)
      END DO
C
C     OUTPUT POINTS
C      
c      WRITE (19,'("POINTS",X,I6,X,"float")') NPK
      WRITE (20,'("POINTS",X,I6,X,"float")') NPK
      DO K=1,NR
       DO I=1,NC(K)
        XCC=XC(I,K)
        YCC=YC(I,K)
        ZCC=0.D+000
c        WRITE (19,'(3(E12.6,X))') XCC,YCC,ZCC            
        WRITE (20,'(3(E12.6,X))') XCC,YCC,ZCC
       END DO
      END DO
C
C     OUTPUT CELLS
C      
c      WRITE (19,'("CELLS",X,I6,X,I6)') NCK,4*NCK
      WRITE (20,'("CELLS",X,I6,X,I6)') NCK,4*NCK
      NCKC=0
      DO K=1,NR
       DO I=1,NMESH(K)
        MMH1=NCKC+MESH(I,1,K)-1
        MMH2=NCKC+MESH(I,2,K)-1
        MMH3=NCKC+MESH(I,3,K)-1
c        WRITE (19,'(I1,3(X,I6))') 3,MMH1,MMH2,MMH3
        WRITE (20,'(I1,3(X,I6))') 3,MMH1,MMH2,MMH3
       END DO
       NCKC=NCKC+NC(K)
      END DO
C
C     OUTPUT CELL TYPES
C      
c      WRITE (19,'("CELL_TYPES",X,I6)') NCK
      WRITE (20,'("CELL_TYPES",X,I6)') NCK
      DO K=1,NR
       DO I=1,NMESH(K)
c        WRITE (19,'(I1)') 5            
        WRITE (20,'(I1)') 5            
       END DO
      END DO
C
C     OUTPUT FIELD HEADER
C
c      WRITE (19,'("POINT_DATA",X,I6)') NPK
      WRITE (20,'("POINT_DATA",X,I6)') NPK
C
C     OUTPUT NORMALS
C
c      WRITE (19,'("NORMALS N float")') 
      WRITE (20,'("NORMALS N float")') 
      DO K=1,NR
       DO I=1,NB(K)
        XNN=XN(I,K)
        YNN=YN(I,K)
        ZNN=0.D+000
c        WRITE (19,'(3(E12.6,X))') XNN,YNN,ZNN            
        WRITE (20,'(3(E12.6,X))') XNN,YNN,ZNN
       END DO
       DO I=NB(K)+1,NC(K)
        XNN=0.D+000
        YNN=0.D+000
        ZNN=1.D+000
c        WRITE (19,'(3(E12.6,X))') XNN,YNN,ZNN            
        WRITE (20,'(3(E12.6,X))') XNN,YNN,ZNN            
       END DO
      END DO
C
C     OUTPUT LEVEL
C
c      WRITE (19,'("SCALARS L float 1")') 
      WRITE (20,'("SCALARS L float 1")') 
c      WRITE (19,'("LOOKUP_TABLE default")') 
      WRITE (20,'("LOOKUP_TABLE default")') 
      DO K=1,NR
       DO I=1,NC(K)
c        WRITE (19,'(E12.6)') SCC(I,K)            
        WRITE (20,'(E12.6)') SCC(I,K)            
       END DO
      END DO
C
C     OUTPUT TEMPERATURE
C
c      WRITE (19,'("SCALARS T float 1")') 
      WRITE (20,'("SCALARS T float 1")') 
c      WRITE (19,'("LOOKUP_TABLE default")') 
      WRITE (20,'("LOOKUP_TABLE default")') 
      DO K=1,NR
       DO I=1,NC(K)
c        WRITE (19,'(E12.6)') TCC(I,K)            
        WRITE (20,'(E12.6)') TCC(I,K)            
       END DO
      END DO
C
C     OUTPUT HEAT FLUX
C
c      WRITE (19,'("VECTORS Q float")') 
      WRITE (20,'("VECTORS Q float")') 
      DO K=1,NR
       DO I=1,NC(K)
        QX=-TCE(K)*DXT(I,K)
        QY=-TCE(K)*DYT(I,K)
        QZ=0.D+000
c        WRITE (19,'(3(E12.6,X))') QX,QY,QZ            
        WRITE (20,'(3(E12.6,X))') QX,QY,QZ
       END DO
      END DO
C
C     OUTPUT PORE PRESSURE
C
c      WRITE (19,'("SCALARS P float 1")') 
      WRITE (20,'("SCALARS P float 1")') 
c      WRITE (19,'("LOOKUP_TABLE default")') 
      WRITE (20,'("LOOKUP_TABLE default")') 
      DO K=1,NR
       DO I=1,NC(K)
c        WRITE (19,'(E12.6)') PCC(I,K)            
        WRITE (20,'(E12.6)') PCC(I,K)            
       END DO
      END DO
C
C     OUTPUT POROUS VELOCITY
C
c      WRITE (19,'("VECTORS V float")') 
      WRITE (20,'("VECTORS V float")') 
      DO K=1,NR
       DO I=1,NC(K)
        IF (KR(K).EQ.0) THEN
         VX=0.D+000
         VY=0.D+000
         VZ=0.D+000
        END IF
        IF (KR(K).EQ.1) THEN
         VX=UCC(I,K)
         VY=VCC(I,K)
         VZ=0.D+000
        END IF
        IF (KR(K).EQ.2) THEN
         VX=VFX(I,K)
         VY=VFY(I,K)
         VZ=0.D+000
        END IF
c        WRITE (19,'(3(E12.6,X))') VX,VY,VZ           
        WRITE (20,'(3(E12.6,X))') VX,VY,VZ
       END DO
      END DO
C
C     OUTPUT DISPLACEMENTS
C
c      WRITE (19,'("VECTORS D float")') 
      WRITE (20,'("VECTORS D float")') 
      DO K=1,NR
       DO I=1,NC(K)
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         UX=UCC(I,K)
         UY=VCC(I,K)
         UZ=0.D+000
        END IF
        IF (KR(K).EQ.1) THEN
         UX=0.D+000
         UY=0.D+000
         UZ=0.D+000
        END IF
c        WRITE (19,'(3(E12.6,X))') UX,UY,UZ            
        WRITE (20,'(3(E12.6,X))') UX,UY,UZ
       END DO
      END DO
C
C     OUTPUT STRESSES
C
c      WRITE (19,'("TENSORS S float")') 
      WRITE (20,'("TENSORS S float")') 
      DO K=1,NR
       ELS=1.D+000-2.D+000*POE(K)
       TEX=2.D+000*VIE(K)*BEE(K)*(1.D+000+POE(K))/ELS
       DO I=1,NC(K)
        IF ((KR(K).EQ.0).OR.(KR(K).EQ.2)) THEN
         SXX=2.D+000*VIE(K)*DXU(I,K)
     &      +2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SYY=2.D+000*VIE(K)*DYV(I,K)
     &      +2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SZZ=2.D+000*VIE(K)*POE(K)*DCC(I,K)/ELS
     &      -PCC(I,K)-TEX*(TCC(I,K)-TRE(K))
         SXY=VIE(K)*(DYU(I,K)+DXV(I,K))
         SXZ=0.D+000
         SYZ=0.D+000
        END IF
        IF (KR(K).EQ.1) THEN
         SXX=2.D+000*VIE(K)*DXU(I,K)
         SYY=2.D+000*VIE(K)*DYV(I,K)
         SZZ=0.D+000
         SXY=VIE(K)*(DYU(I,K)+DXV(I,K))
         SXZ=0.D+000
         SYZ=0.D+000
        END IF
c        WRITE (19,'(3(E12.6,X))') SXX,SXY,SXZ            
c        WRITE (19,'(3(E12.6,X))') SXY,SYY,SYZ 
c        WRITE (19,'(3(E12.6,X))') SXZ,SYZ,SZZ            
        WRITE (20,'(3(E12.6,X))') SXX,SXY,SXZ            
        WRITE (20,'(3(E12.6,X))') SXY,SYY,SYZ 
        WRITE (20,'(3(E12.6,X))') SXZ,SYZ,SZZ            
       END DO
      END DO
C
c      CLOSE (19)
      CLOSE (20)
C
C----------------------------------------------------------------------
C
      RETURN
      END