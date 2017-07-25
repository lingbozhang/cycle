! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            LOOP  V0.1                                 .
! .                                                                       .
! .    MAIN PROGRAM FOR TEST                                              .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
      PROGRAM MAIN
      USE BUILDINS
      IMPLICIT NONE
!
      INTEGER(SP_INT) :: A,B
      CHARACTER :: C
!
      A=OPEN_FILE("HELLO   ")
      write(*,*) A
!
      END PROGRAM MAIN
