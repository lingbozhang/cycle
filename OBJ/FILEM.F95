! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    File module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
      module Warehouse_FILE
!{{{
       implicit none
       TYPE FileADType
         INTEGER :: IIN = 3
         INTEGER :: IOUT  = 7
         INTEGER :: IOUT2 = 8      ! K MATRIX
         INTEGER :: IOUT3 = 9      ! M MATRIX
         INTEGER :: IOUT4 = 10     ! C MATRIX
         INTEGER :: IOUT5 = 11     ! STATIC DISPLACEMENT DATA
         INTEGER :: IOUT8 = 15     ! DYNAMIC DISPLACEMENT DATA
         INTEGER :: IOUT6 = 12     ! STRESS DATA
         INTEGER :: IOUT7 = 13     ! TEST DATA
!
         INTEGER :: ITEM1 = 14     ! TEMPORAL ELEMENT DATA
         INTEGER :: ITEM3 = 17     ! TEMPORAL BOUNDARY DATA
         INTEGER :: ITEM2 = 16     ! TEMPORAL SURFACE DATA
       END TYPE FileADType
       TYPE(FileADType) :: FileAD
!}}}
      end module Warehouse_FILE
!
      MODULE FILEM
        IMPLICIT NONE
!--------------------------------------------------
!   INTERFACE FOR FILEM
!          OPENFILES()
!          CLOSEFILES()
!          InputAD()
!          OutputAD()
!--------------------------------------------------
        PUBLIC  :: OPENFILES
        PUBLIC  :: CLOSEFILES
        PUBLIC  :: InputAD
        PUBLIC  :: OutputAD
        PUBLIC  :: OUTPUTELE
        PUBLIC  :: OUTPUTBC
        PUBLIC  :: OUTPUTSURF
        PUBLIC  :: OUTPUTK
        PUBLIC  :: OUTPUTMA
        PUBLIC  :: OUTPUTDISP
        PUBLIC  :: OUTPUTDISPD
        PUBLIC  :: OUTPUTSTRESS
!  
        CONTAINS
! * * * * * * * * * * * * * * * * * * *
! *     SUBROUTINE OPENFILES          *
! * * * * * * * * * * * * * * * * * * *
        SUBROUTINE OPENFILES()
        USE WAREHOUSE, ONLY : FILE_NAME
!{{{
        USE Warehouse_FILE
        IMPLICIT NONE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
        LOGICAL :: EX
!--IF DISCOVERY.IN DOES NOT EXIST: STOP
!     
        WRITE(*,*) 'DISCOVERY'//FILE_NAME//'.IN'
        INQUIRE(FILE ='DISCOVERY'//FILE_NAME//'.IN',EXIST=EX)
        IF (.NOT.EX) THEN
             PRINT *, "*** STOP *** FILE DISCOVERY.IN DOES NOT EXIST !"
             STOP
        END IF
!       
!--OPEN FILES: STORE INPUT DATA
        OPEN(FileAD%IIN   ,  FILE ='DISCOVERY'//FILE_NAME//'.IN', STATUS="OLD")
!--OPEN FILES: SYSTEM DATA
        OPEN(FileAD%IOUT  ,  FILE = "DISCOVERY.OUT", STATUS = "REPLACE")
        OPEN(FileAD%IOUT2 ,  FILE = "K.OUT",     STATUS = "REPLACE")
        OPEN(FileAD%IOUT3 ,  FILE = "M.OUT",     STATUS = "REPLACE")
        OPEN(FileAD%IOUT4 ,  FILE = "C.OUT",     STATUS = "REPLACE")
        OPEN(FileAD%IOUT5 ,  FILE = "DISP.OUT",  STATUS = "REPLACE")
        OPEN(FileAD%IOUT8 ,  FILE = "DISPD.OUT", STATUS = "REPLACE")
        OPEN(FileAD%IOUT6 ,  FILE = "STRESS.OUT",STATUS = "REPLACE")
        OPEN(FileAD%IOUT7 ,  FILE = "TEST.OUT",  STATUS = "REPLACE")
!
        OPEN(FileAD%ITEM1 ,  FORM = "UNFORMATTED", STATUS = "SCRATCH")
        OPEN(FileAD%ITEM2 ,  FORM = "UNFORMATTED", STATUS = "SCRATCH")
        OPEN(FileAD%ITEM3 ,  FORM = "UNFORMATTED", STATUS = "SCRATCH")       
!}}}
        END SUBROUTINE OPENFILES

! * * * * * * * * * * * * * * * * * * *
! *     SUBROUTINE CLOSEFILES         *
! * * * * * * * * * * * * * * * * * * *
        SUBROUTINE CLOSEFILES()
!{{{
        USE Warehouse_FILE
        IMPLICIT NONE
!--CLOSE FILES: STORE INPUT DATA
        CLOSE(FileAD%IIN)
!--CLOSE FILES: SYSTEM DATA
        CLOSE(FileAD%IOUT)
        CLOSE(FileAD%IOUT2)
        CLOSE(FileAD%IOUT3)
        CLOSE(FileAD%IOUT4)
        CLOSE(FileAD%IOUT5)
        CLOSE(FileAD%IOUT6)
        CLOSE(FileAD%IOUT7)
        CLOSE(FileAD%IOUT8)
!
        CLOSE(FileAD%ITEM1)
        CLOSE(FileAD%ITEM2)
        CLOSE(FileAD%ITEM3)
!}}}
        END SUBROUTINE CLOSEFILES
!
        INTEGER FUNCTION InputAD() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IIN
!}}}
        END FUNCTION InputAD
!
        INTEGER FUNCTION OutputEle() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%ITEM1
!}}}
        END FUNCTION OutputEle
!
        INTEGER FUNCTION OutputBC() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%ITEM3
!}}}
        END FUNCTION OutputBC
!
        INTEGER FUNCTION OutputSURF() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%ITEM2
!}}}
        END FUNCTION OutputSURF
!
        INTEGER FUNCTION OutputAD() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE
        IMPLICIT NONE
!
        I=FileAD%IOUT
!}}}
        END FUNCTION OutputAD
!
        INTEGER FUNCTION OutputK() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IOUT2
!}}}
        END FUNCTION OutputK
!
        INTEGER FUNCTION OutputMa() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IOUT3
!}}}
        END FUNCTION OutputMa
!
        INTEGER FUNCTION OutputDisp() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IOUT5
!}}}
        END FUNCTION OutputDisp
!
        INTEGER FUNCTION OutputDispD() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IOUT8
!}}}
        END FUNCTION OutputDispD    
!
        INTEGER FUNCTION OutputStress() RESULT(I)
!{{{
!-----------------------------------------------------
!       OUTPUT INPUT DOCUMENT ADDRESS
!-----------------------------------------------------
        USE Warehouse_FILE                
        IMPLICIT NONE
        I=FileAD%IOUT6
!}}}
        END FUNCTION OutputStress
!-----------------------------------------------------
        END MODULE FILEM

