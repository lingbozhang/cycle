! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Main program for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
program DISCOVERY
    !--------------------------------------------------
    !   FEM PROGRAM (MAIN)
    !--------------------------------------------------
    USE INTERPRETER
    USE FILEM
    USE INPUTM
    USE WAREHOUSE
    USE OPERATIONM
    USE GENERATEM
    USE SOLVERM
    USE OUTPUTM
!    USE ODEM
    IMPLICIT NONE
    !
    INTEGER :: IIN
    INTEGER :: IOUT,ADBC,ADELE,ADSURF,ADK,ADM
    INTEGER :: ADDISP,ADSTRE,ADDISPD
    !
    REAL(8) :: TT(2)
    INTEGER,PARAMETER :: NumKW=7
    CHARACTER(20),DIMENSION(NumKW) :: KeyCtrl
    !   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    ! TEST
    !--------------------------------------------------
    !   VARIABLES FOR DYNAMIC SOLVER
    !--------------------------------------------------
    !
    !------------------------
    !   Variables definition
    !------------------------
    KeyCtrl(1) = '\START'
    KeyCtrl(2) = '\INPUT'
    KeyCtrl(3) = '\GENERATE'
    KeyCtrl(4) = '\STATIC_SOLVER'
    KeyCtrl(5) = '\DYNAMIC_SOLVER'
    KeyCtrl(6) = '\OUTPUT'
    KeyCtrl(7) = '\END'
    ! INPUT FILE NUMBER
    !   OPEN FILES
    WRITE(*,*) "PLEASE INPUT NUMBER OF FILES!"
    READ(*,*) NFILE
!
    DO J=1,NFILE
        write(FILE_NAME,'(I3.3)') J
!
        CALL OPENFILES()  
        !   READ CONTROL FILE ADDRESS
        IIN=InputAD()
        CALL ReadFile(IIN)  !LineNumber&WordNumber are resetted to 0
        !
        I=KeyWord(KeyCtrl,NumKW)
        DO WHILE(I<=NumKW) 
            SELECT CASE(I)
             CASE(1) 
                    IOUT=OutputAD()
                    ADBC=OutputBC()
                    ADELE=OutputEle()
                    ADSURF=OutputSurf()
                    ADK=OutputK()
                    ADM=OutputMa()
                    !
                    ADDISP =OutputDisp()
                    ADDISPD=OutputDispD()
                    ADSTRE =OutputStress()
                    WRITE(IOUT,*) KeyCtrl(1)
             CASE(2)
                    WRITE(IOUT,*) KeyCtrl(2)
                    WRITE(*,*) "*********************"
                    WRITE(*,*) "INPUT PROCESSING"
                    WRITE(*,*) "*********************"
!
                    CALL CPU_TIME(TT(1))
                    CALL INPUT(IIN,IOUT,ADBC,ADELE,ADSURF)
                    CALL COL_ALLOCATE()
                    CALL CPU_TIME(TT(2))
!
                    WRITE(*,*)  "---------------------------------------------------"
                    WRITE(*,'(A20,ES14.7)')  "INPUT TIME =" ,TT(2)-TT(1)
                    WRITE(*,*)  "---------------------------------------------------"
             CASE(3)
                    WRITE(IOUT,*) KeyCtrl(3)
                    WRITE(*,*) "*********************"
                    WRITE(*,*) "GENERATE PROCESSING"
                    WRITE(*,*) "*********************"
!
                    CALL CPU_TIME(TT(1))
                    IF (DYNAMIC.EQ.0) THEN
                        CALL GENERATE(ADK,ADM,ADELE,ADBC)
                    ELSE
                        IF (IND_INITIAL.EQ.0) THEN
                            WRITE(*,*) "ERROR, INITIAL CONDITIONS REQUIRED!"
                            STOP
                        ENDIF
                        CALL GENERATE(ADK,ADM,ADELE,ADBC)
                    ENDIF
!                    
                    CALL CPU_TIME(TT(2))
!                   CALCULATE HALF_BANDWIDTH WITH ZERO ELIMINATED
                    CALL HALF_BANDWIDTH(K_COLSOL(1),NEQ,NWK_COLSOL)
                    WRITE(*,*)  "HALF-BANDWIDTH IS ",HALF_BAND
!
                    WRITE(*,*)  "---------------------------------------------------"
                    WRITE(*,'(A20,ES14.7)')  "GENERATION TIME =" ,TT(2)-TT(1)
                    WRITE(*,*)  "---------------------------------------------------"
             CASE(4)
                    WRITE(IOUT,*) KeyCtrl(4)
             ! ---- STATIC_SOLVER
             !{{
                    WRITE(*,*) "*********************"  
                    WRITE(*,*) "STATIC SOLVER"  
                    WRITE(*,*) "*********************"  
!
                    CALL CPU_TIME(TT(1))
                    CALL COLSOL_PRE()
                    CALL COLSOL(K_COLSOL(1),DISP(1),MAXA_COLSOL(1),NEQ,NWK_COLSOL,NEQ+1,1)
!                   CALCULATE HALF_BANDWIDTH WITH ZERO ELIMINATED AFTER FACTORIZATION
!                    CALL HALF_BANDWIDTH(K_COLSOL(1),NEQ,NWK_COLSOL)
!                    WRITE(*,*)  "HALF-BANDWIDTH (AFTER FACTORIZATION) IS ",HALF_BAND
!
                    CALL COLSOL(K_COLSOL(1),DISP(1),MAXA_COLSOL(1),NEQ,NWK_COLSOL,NEQ+1,2)
                    CALL CPU_TIME(TT(2))

!
                    WRITE(*,*)  "---------------------------------------------------"
                    WRITE(*,'(A20,ES14.7)')  "SOLUTION TIME =" ,TT(2)-TT(1)
                    WRITE(*,*)  "---------------------------------------------------"
                  !  WRITE(*,*) "DISPX-1",DISP(ID(1,1))
                  !  WRITE(*,*) "DISPY-1",DISP(ID(2,1))
                  !  WRITE(*,*) "DISPZ-1",DISP(ID(3,1))

!
                    ESTRAIN=0
                    HESTRAIN=MATRIXB_VECTOR(KMATRIX,MAXA_COLSOL,MHT,DISP,NWK_COLSOL,NEQ)
                    DO K=1,NEQ
                       ESTRAIN=ESTRAIN+DISP(K)*HESTRAIN(K)
                    ENDDO
                    ESTRAIN=ESTRAIN/2.0
                    WRITE(*,*) "STRAIN ENERGY = ", ESTRAIN
             !}}}
             CASE(5)
               WRITE(*,*) "ERROR, ODEM IS NOT INCLUDED YET!"
               STOP
    !         ! ---- DYNAMIC_SOLVER
    !         !{{{
    !                IF (DYNAMIC.NE.1) THEN
    !                    WRITE(*,*) "ERROR, DYNAMIC IS .NE. 1"
    !                    STOP
    !                ENDIF
    !                CALL ODE_IN(IOUT,R(1),RODE(1),ACC(1))
    !                OUTPUT_ODE=.FALSE.
    !                HALF_ODE=.FALSE.
    !                TIME=TPARA(1)
    !                DTIME=TPARA(3)
    !                NUM_STEP=0
    !                WRITE(*,2000) NUM_STEP,DTIME,TIME
    !                IF (CONTACT.EQ.0) THEN
    !                  WRITE(ADDISPD,1000)
    !                  CALL OUTPUTD(IND=2,IDISP=ADDISPD)
    !                  DO WHILE(.NOT.OUTPUT_ODE)
    !                    IF((TPARA(2)-TIME).GT.2.0*DTIME) THEN
    !                      CONTINUE
    !                    ELSEIF (((TPARA(2)-TIME).GT.DTIME).AND.(.NOT.HALF_ODE)) THEN
    !                      DTIME=0.5*DTIME
    !                      HALF_ODE=.TRUE.
    !                    ELSE
    !                      DTIME=TPARA(2)-TIME
    !                      OUTPUT_ODE=.TRUE.
    !                    END IF
    !                    TIME=TIME+DTIME
    !                    !
    !                    CALL ODE_SOLVER(KMATRIX(1),MMATRIX(1),R(1),NRF(1),U(1),V(1),ACC(1),TIME,DTIME)
    !                    !
    !                    CALL OUTPUTD(IND=2,IDISP=ADDISPD)
    !                    NUM_STEP=NUM_STEP+1
    !                    WRITE(*,2000) NUM_STEP,DTIME,TIME
    !                  END DO
    !                ELSE
    !                  CONTINUE
    !                ENDIF
    !         !}}}
             CASE(6)
                    WRITE(IOUT,*) KeyCtrl(6)
                    WRITE(*,*) "*********************"                     
                    WRITE(*,*) "OUTPUT PROCESSING"
                    WRITE(*,*) "*********************"                     
                    CALL CPU_TIME(TT(1))
                    CALL OUTPUT(IIN,IOUT,ADDISP,ADSTRE,ADELE,ADBC) 
                    CALL CPU_TIME(TT(2))
!
                    WRITE(*,*)  "---------------------------------------------------"
                    WRITE(*,'(A20,ES14.7)')  "OUTPUT TIME =" ,TT(2)-TT(1)
                    WRITE(*,*)  "---------------------------------------------------"                 
             CASE(7)
                  CALL CLOSEFILES
                  CALL RELEASE
!
                  EXIT
             CASE DEFAULT
                WRITE(*,*) "ERROR, WRONG INPUT IN CONTROL FILE"
                STOP
             END SELECT
             I=KeyWord(KeyCtrl,NumKW)
        END DO
    ENDDO
    print *, char(7)
    STOP
!
end program DISCOVERY
