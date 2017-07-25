! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.1                                 
! .                                                                       
! .    Build_ins for CYCLE
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
      MODULE BUILDINS
!{{{
      PUBLIC :: OPEN_FILE
      PUBLIC :: HASH
!     PRECISION CONTROL
      INTEGER, PARAMETER :: DP     = SELECTED_REAL_KIND(15,307)
      INTEGER, PARAMETER :: SP     = SELECTED_REAL_KIND(6,37)
      INTEGER, PARAMETER :: SP_INT = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: DP_INT = SELECTED_INT_KIND(9)
!     MATH CONSTANTS
      REAL(DP),PARAMETER :: DEG_PER_RAD = 57.295779513082320876798155_DP
      REAL(DP),PARAMETER :: RAD_PER_DEG = 0.017453292519943295769237_DP
!
      REAL(DP),PARAMETER :: E_VALUE     = 2.718281828459045235360287_DP
      REAL(DP),PARAMETER :: E_RECIP     = 0.3678794411714423215955238_DP
      REAL(DP),PARAMETER :: E_SQUARED   = 7.389056098930650227230427_DP
      REAL(DP),PARAMETER :: LOG10_OF_E  = 0.4342944819032518276511289_DP
!
      REAL(DP),PARAMETER :: LN_2        = 0.693147180559945309417232121_DP
      REAL(DP),PARAMETER :: LN_10       = 2.3025850929940456840179915_DP
      REAL(DP),PARAMETER :: LOG10_OF_2  = 0.3010299956639811952137389_DP
!
      REAL(DP),PARAMETER :: PI_VALUE    = 3.141592653589793238462643_DP
      REAL(DP),PARAMETER :: PI_LN       = 1.144729885849400174143427_DP
      REAL(DP),PARAMETER :: PI_LOG10    = 0.4971498726941338543512683_DP
      REAL(DP),PARAMETER :: PI_OVER_2   = 1.570796326794896619231322_DP
      REAL(DP),PARAMETER :: PI_OVER_3   = 1.047197551196597746154214_DP
      REAL(DP),PARAMETER :: PI_OVER_4   = 0.7853981633974483096156608_DP
      REAL(DP),PARAMETER :: PI_RECIP    = 0.3183098861837906715377675_DP
      REAL(DP),PARAMETER :: PI_SQUARED  = 9.869604401089358618834491_DP
      REAL(DP),PARAMETER :: PI_SQ_ROOT  = 1.772453850905516027298167_DP
!
      REAL(DP),PARAMETER :: SQ_ROOT_OF_2= 1.4142135623730950488_DP
      REAL(DP),PARAMETER :: SQ_ROOT_OF_3= 1.7320508075688772935_DP
!--------------------------------------------------------------------------

      CONTAINS

!--------------------------------------------------
! Open file function (07/25/2017) 
! From Lingbo Zhang
! --------------------------------------------------
! --------------------------------------------------
!       INPUTS:
!          FILENAME; 
!          FLAG:    
!              1. IN   -- STATUS="IN" 
!              2. OUT  -- STATUS="OUT" 
!              3. TEMP -- STATUS="SCRATCH"
!       FUNCTION:
!          OPEN FILE;
! --------------------------------------------------
      FUNCTION OPEN_FILE(FILENAME,FLAG)
!{{{
      IMPLICIT NONE
      INTEGER(SP_INT) :: OPEN_FILE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FLAG   
!     CHARACTER FOR STATUS
      CHARACTER(LEN=4) :: STATUS_DEFAULT="OUT "
      CHARACTER(LEN=4),PARAMETER :: STATUS_1="IN  "
      CHARACTER(LEN=4),PARAMETER :: STATUS_2="OUT "
      CHARACTER(LEN=4),PARAMETER :: STATUS_3="TEMP"
!
      INTEGER(SP_INT) :: FILE_ID
      LOGICAL         :: EX
!     Pre-set value FILE_OPEN
      OPEN_FILE=0     
!
      FILE_ID=HASH(FILENAME)
      print *, FILE_ID
!
      IF(PRESENT(FLAG)) THEN
         STATUS_DEFAULT=ADJUSTL(FLAG)
         STATUS_DEFAULT=STATUS_DEFAULT(1:4)
      ENDIF
!
      SELECT CASE(STATUS_DEFAULT)
      CASE(STATUS_1)
        INQUIRE(FILE=FILENAME//"."//"IN",EXIST=EX)
        IF(.NOT.EX) THEN
           PRINT *, "*** STOP *** FILE DOES NOT EXIST! --- "&
                    ,FILENAME//"."//"IN"
           STOP
        ENDIF
        OPEN(FILE_ID, FILE=FILENAME//"."//"IN",STATUS="OLD")
      CASE(STATUS_2)
        OPEN(FILE_ID, FILE=FILENAME//"."//"OUT",STATUS="REPLACE")
      CASE(STATUS_3)
        OPEN(FILE_ID, FORM="UNFORMATTED", STATUS="SCRATCH")
      CASE DEFAULT
        PRINT *,"*** STOP *** /BUILDINS/FILE_OPEN/SELECT CASE"
        STOP
      END SELECT
!     File opens successfully and set value FILE_OPEN=1
      OPEN_FILE=1
!}}}
      END FUNCTION OPEN_FILE
!-------------------------------------------------------------------------
! Hashing function (07/25/2017) 
! FROM Lingbo Zhang
! Warnings:
!         The length of the input string should be less than or equal to 6
! --------------------------------------------------------------------------
      function HASH(str) result(HASH_VALUE)
!{{{
      implicit none
      character(len=*),intent(in) :: str
      integer(DP_INT) :: HASH_VALUE
      integer(SP_INT) :: i

      HASH_VALUE = 0_DP_INT

      do i=1,len(str)
         HASH_VALUE = (ishft(HASH_VALUE,5) + HASH_VALUE) + ichar(str(i:i))
      end do
!}}}
      end function HASH
!------------------------------------------      
!}}}
      END MODULE BUILDINS

