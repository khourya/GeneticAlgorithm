CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     
CCC  PROGRAM NDE_iVoF_MMGA                         
CCC                                                                     
CCC  Version 1.0: Parallel (mpich)                              
CCC                                                                     
CCC  Non-Destructive Evaluation                       
CCC                    
CCC  Inverse Volume-of-Fluid
CCC
CCC  Meshless Method
CCC
CCC  Genetic Algorithm    
CCC                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  
CCC   Embry-Riddle Aeronautical University
CCC
CCC   MDBL: Multi-Disciplinary Bioengineering Lab
CCC
CCC   Eduardo Divo
CCC
CCC   Hussein Saad
CCC
CCC   September 23, 2014
CCC 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                     
CCC                             MAIN PROGRAM                            
CCC                                                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PGA_Parameters.h'
      INCLUDE 'mpif.h'
C
      COMMON/PGAINFO/NGEN,NPOP,NPAR,NBIT(NPARMAX)
      COMMON/PGAMUTA/PJMU,PCMU
      COMMON/PARINFO/PARMIN(NPARMAX),PARMAX(NPARMAX),PARRES(NPARMAX)
      COMMON/PARAMET/PARAM(NPARMAX,NPOPMAX)
      COMMON/NEWGENE/CHILD(NPARMAX,NPOPMAX)
      COMMON/FITNESS/FITNESS(NPOPMAX),PSEL(NPOPMAX),JBEST
C
      REAL*8 STARTTIME,ENDTIME
      REAL*8 TIMEMY,TIMETO
      REAL*8 TIM(MAXPROC),FRA(MAXPROC)
      INTEGER ID,IDMASTER,NPRO,IERR
      INTEGER ILOAD(MAXPROC,NPOPMAX)
C
      INTEGER IGEN,IG
      INTEGER JP1,JP2,IC
      INTEGER IK
      INTEGER UNDAT,UNOUT
C
C**********************************************************************
C
C
C     START OF PARALLEL BENCHMARKING
C
      CALL MPI_INIT(IERR)
C
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,ID,IERR)
C
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPRO,IERR)
C
C     SPECIFY MASTER PROCESS
C
      IDMASTER=0
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*)
       WRITE(*,*) 'PGA: MULTIVARIABLE PARALLEL GENETIC ALGORITHM'
       WRITE(*,*) '     OPTIMIZATION TOOL'
       WRITE(*,*) 
       WRITE(*,*) 'VERSION 1.1'
       WRITE(*,*) 
       WRITE(*,*) 'EMBRY-RIDDLE AERONAUTICAL UNIVERSITY'
       WRITE(*,*) 
       WRITE(*,*)
       WRITE(*,*) 'NUMBER OF PROCESSORS...................:',NPRO
       WRITE(*,*)
      END IF
C
C**********************************************************************
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
       WRITE(*,*) 
       WRITE(*,*) 'CLUSTER BENCHMARKING PROCESS STARTED..........: [OK]'
       WRITE(*,*)
       WRITE(*,*)
      END IF
C
      CALL BENCHMARK(TIM,FRA,ID,IDMASTER,NPRO)
C
      IF (ID.EQ.IDMASTER) THEN
       DO N=1,NPRO
        WRITE (*,'("   PROCESS, BENCHMARK TIME & FRACTION...: ",I3,
     &         2X,F6.2,2X,F6.4)') N,TIM(N),FRA(N)
       END DO
       WRITE(*,*) 
       WRITE(*,*) 'CLUSTER BENCHMARKING PROCESS ENDED............: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
C
C**********************************************************************
C
C
C     INPUT DATA FILE AND INITIAL PARAMETERS
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'INPUT GENETIC ALGORITHM PARAMETERS STARTED....: [OK]'
       WRITE(*,*)
      END IF 
C
      CALL INPUT(ID,IDMASTER)
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,'("   NUMBER OF GENERATIONS.............: ",I5)') NGEN
       WRITE(*,'("   POPULATION SIZE...................: ",I5)') NPOP
       WRITE(*,'("   PROBABILITY OF JUMP MUTATION......:  ",F4.2)') PJMU
       WRITE(*,'("   PROBABILITY OF CREEP MUTATION.....:  ",F4.2)') PCMU
       WRITE(*,*) 
       WRITE(*,*) 'INPUT GENETIC ALGORITHM PARAMETERS ENDED......: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
C
C**********************************************************************
C
C     LOAD BALANCING ALGORITHM
C
C
C     INITIALIZE RANDOM GENERATOR
C
      CALL RANDOM(-1000-100*ID,R)
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'LOAD BALANCING PROCESS STARTED................: [OK]'
C
       CALL LOAD(IDMASTER,ILOAD,FRA,NPRO)
C
       WRITE(*,*) 
       WRITE(*,*) '  FINAL POPULATION BALANCE OVER PROCESSORS....:'
       WRITE(*,*) 
       DO N=1,NPRO
        WRITE(*,'("  PROCESS & FRACTION: ",I3,1X,F6.4,2X,100(I1))') 
     @                           N,FRA(N),(ILOAD(N,K),K=1,NPOP)
       END DO
C
       WRITE(*,*) 
       WRITE(*,*) 'LOAD BALANCING PROCESS ENDED..................: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
C
C     BROADCASTING CLUSTER WORK LOAD
C
      CALL MPI_BCAST(ILOAD,MAXPROC*NPOPMAX,MPI_INTEGER,
     &               IDMASTER,MPI_COMM_WORLD,IERR)
C
C**********************************************************************
C
C     INPUT DATA FOR FUNCTION EVALUATION
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'FUNCTION EVALUATION DATA INPUT STARTED........: [OK]'
      END IF
C
      CALL FUNCTIONSETUP(ID,IDMASTER,NPRO)
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,'("   NUMBER OF PARAMETERS..............: ",I5)') NPAR
       WRITE(*,*) 
       WRITE(*,*) 'FUNCTION EVALUATION DATA INPUT ENDED..........: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
C
C**********************************************************************
C
C     GENERATE INITIAL POPULATION
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'INITIAL POPULATION GENERATION STARTED.........: [OK]'
      END IF
C
C     DATA FILES UNIT NUMBERS
C
      UNDAT=8
      UNOUT=9
C
      CALL INITIAL(ID,IDMASTER,IGEN,UNDAT,UNOUT)
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'INITIAL POPULATION GENERATION ENDED...........: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
ccc
C      PARAM(1,1)=0.1
C      PARAM(2,1)=0.2
C      PARAM(3,1)=0.025
C      PARAM(4,1)=0.05
ccc
C
C**********************************************************************
C
C     MAIN OPTIMIZATION LOOP
C
      IF (ID.EQ.IDMASTER) THEN
       STARTTIME=MPI_WTIME()
       WRITE(*,*) 
       WRITE(*,*) 'MAIN GENETIC OPTIMIZATION LOOP STARTED........: [OK]'
       FITNESSMAX=0.
       IK=0
      END IF
C
      CALL FUNCTIONEVALUATE(ID,IDMASTER,NPRO,ILOAD)
C
      DO IG=1,NGEN
C
       DO IC=1,NPOP-1
C
        IF (ID.EQ.IDMASTER) CALL SELECTION(JP1,JP2)
C
        IF (ID.EQ.IDMASTER) CALL REPRODUCT(JP1,JP2,IC)
C
       END DO
C
       IF (ID.EQ.IDMASTER) THEN
        IK=IK+1
        IF (FITNESS(JBEST).GT.FITNESSMAX) THEN
         WRITE(*,010) IGEN+IG,FITNESS(JBEST)
         FITNESSMAX=FITNESS(JBEST)
         IK=0
        END IF
        IF (IK.GE.50) THEN
         CALL KILLGEN()
         WRITE(*,010) IGEN+IG,FITNESS(JBEST)
         IK=0
        END IF
       END IF
C
       CALL NEWGEN(ID,IDMASTER,NPRO,IGEN,IG,ILOAD,UNDAT,UNOUT)
C
       CALL FUNCTIONEVALUATE(ID,IDMASTER,NPRO,ILOAD)
C
      END DO
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,010) IGEN+IG-1,FITNESS(JBEST)
  010  FORMAT ('  GENERATION: ',I6,'    BEST FITNESS: ',E10.4)
       WRITE(*,*) 
       WRITE(*,*) 'MAIN GENETIC OPTIMIZATION LOOP ENDED..........: [OK]'
       CLOSE (UNOUT)
       ENDTIME=MPI_WTIME()
      END IF
C
C**********************************************************************
C
C     OUTPUT FUNCTION DATA
C
      IF (ID.EQ.IDMASTER) THEN
       WRITE(*,*) 
       WRITE(*,*) 'FUNCTION DATA OUTPUT STARTED..................: [OK]'
C
       CALL FUNCTIONOUTPUT()
C
       WRITE(*,*) 
       WRITE(*,*) 'FUNCTION DATA OUTPUT ENDED....................: [OK]'
       WRITE(*,*) 
       WRITE(*,*) '****************************************************'
      END IF
C
C**********************************************************************
C
C     OUTPUT TIMES
C
      IF (ID.EQ.IDMASTER) THEN
       TOTALTIME=ENDTIME-STARTTIME
       OPEN (21,FILE='Model/PGA.time')
       WRITE (21,*) 'NUMBER OF GENERATIONS PERFORMED..: ',NGEN
       WRITE (21,*) 'TOTAL TIME ELAPSED...............: ',TOTALTIME
       TOTALTIME=TOTALTIME/DBLE(NGEN)
       WRITE (21,*) 'TIME ELAPSED PER GENERATION......: ',TOTALTIME
       CLOSE (21)
      END IF
C
C**********************************************************************
C
      CALL MPI_FINALIZE(IERR)
C
      END
