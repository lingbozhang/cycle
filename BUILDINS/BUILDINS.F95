! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            LOOP  V0.1                                 
! .                                                                       
! .    Build_ins for LOOP 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
      MODULE BUILDINS
!{{{
      PUBLIC :: L_OPEN
      PUBLIC :: IFILE
!
      CONTAINS
! * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINES AND FUNCTIONS FOR FILE  *
! * * * * * * * * * * * * * * * * * * * * * *
! --------------------------------------------------
!       INPUTS:
!          FILENAME; 
!          INDICATOR:    
!              1. IN   -- STATUS="IN" 
!              2. OUT  -- STATUS="OUT" 
!              3. TEMP -- STATUS="SCRATCH"
!       FUNCTION:
!          OPEN FILE;
! --------------------------------------------------
      SUBROUTINE L_OPEN(FILENAME,INDICATOR)
!{{{
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: INDICATOR   
!     CHARACTER FOR STATUS
      CHARACTER(LEN=4) :: IND="OUT "
      CHARACTER(LEN=4),PARAMETER :: IND_1_C="IN  "
      CHARACTER(LEN=4),PARAMETER :: IND_2_C="OUT "
      CHARACTER(LEN=4),PARAMETER :: IND_3_C="TEMP"
!
      CHARACTER(LEN=30) :: FILE_NAME
      CHARACTER(LEN=30) :: IND_NAME
      INTEGER :: I,J,K
      INTEGER :: FILE_NUMBER
      INTEGER :: FILE_LENGTH
      LOGICAL :: EX
!
      FILE_NAME=ADJUSTL(FILENAME)
      FILE_LENGTH=LEN_TRIM(FILE_NAME)      
!     GENERATE FILE_NUMBER
      K=IACHAR("0")
      FILE_NUMBER=0
      DO I=1,FILE_LENGTH
         J=IACHAR(FILE_NAME(I:I))-K
         FILE_NUMBER=FILE_NUMBER+36*(I-1)*J
      ENDDO
!
      IF(PRESENT(INDICATOR)) THEN
         IND_NAME=ADJUSTL(INDICATOR)
         IND=IND_NAME(1:4)
      ENDIF
!
      SELECT CASE(IND)
      CASE(IND_1_C)
        INQUIRE(FILE=FILE_NAME(1:FILE_LENGTH)//"."//"IN",EXIST=EX)
        IF(.NOT.EX) THEN
           PRINT *, "*** STOP *** FILE DOES NOT EXIST! --- "&
                    ,FILE_NAME(1:FILE_LENGTH)//"."//"IN"
           STOP
        ENDIF
        OPEN(FILE_NUMBER, FILE=FILE_NAME(1:FILE_LENGTH)//"."//"IN",STATUS="OLD")
      CASE(IND_2_C)
        OPEN(FILE_NUMBER, FILE=FILE_NAME(1:FILE_LENGTH)//"."//"OUT",STATUS="REPLACE")
      CASE(IND_3_C)
        OPEN(FILE_NUMBER, FORM="UNFORMATTED", STATUS="SCRATCH")
      CASE DEFAULT
        PRINT *,"*** STOP *** /BUILDINS/L_OPEN/SELECT CASE"
        STOP
      END SELECT
!}}}
      END SUBROUTINE L_OPEN
! --------------------------------------------------
!       INPUTS
!          FILENAME;
!       FUNCTION:
!          CONVERT FILENAME TO INTEGER
! --------------------------------------------------
      FUNCTION IFILE(FILENAME)
!{{{
      IMPLICIT NONE
      INTEGER :: IFILE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
!
      CHARACTER(LEN=30) :: FILE_NAME
      INTEGER :: I,J,K
      INTEGER :: FILE_LENGTH
!
      FILE_NAME=ADJUSTL(FILENAME)
      FILE_LENGTH=LEN_TRIM(FILE_NAME) 
!     GENERATE FILE_NUMBER
      K=IACHAR("0")
      IFILE=0
      DO I=1,FILE_LENGTH
         J=IACHAR(FILE_NAME(I:I))-K
         IFILE=IFILE+36*(I-1)*J
      ENDDO
!}}}
      END FUNCTION IFILE
!}}}
      END MODULE BUILDINS

