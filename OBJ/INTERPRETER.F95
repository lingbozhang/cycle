! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Test parsing module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

      MODULE Warehouse_INT
!{{{
!------------------------------------------------------------------
!   purpose: warehouse for storing variables for interpreter
!------------------------------------------------------------------
           IMPLICIT NONE
           INTEGER :: InFile    ! INPUT FILE
           PUBLIC   :: LineDate
           TYPE LineDate
              INTEGER          :: LineNumber, WordNumber, ReadNumber
              CHARACTER(256)   :: Line_original
              CHARACTER(20)    :: Line_words(15)
           END TYPE LineDate
           TYPE(LineDate) :: Line
!}}}
      END module Warehouse_INT
!
      MODULE INTERPRETER
!----------------------------------------------------------------------
!     INTERPERTER FOR FORTRAN LANGUAGE
!     FUNCTIONS:
!        ReadFile(address)  :: Read File information
!        KeyWord(kW,nbkw)   :: read one line from INPUT file and identify the keyword
!        GetString(string)  :: Get a Linedate
!        GetInt()           :: Get a Integer number
!        GetReal()          :: Get a Real number
!        NewLine()          :: Start a new line
!----------------------------------------------------------------------

      IMPLICIT NONE
        PUBLIC   :: ReadFile
        PUBLIC   :: KeyWord
        PUBLIC   :: GetString
        PUBLIC   :: GetInt
        PUBLIC   :: GetReal
        PUBLIC   :: NewLine
        PRIVATE  :: pcomp
        PRIVATE  :: ReadLine
        PRIVATE  :: isNumber
!
        CONTAINS
!
        Integer function ReadLine()
!{{{
!------------------------------------------------------------------
!       purpose: read one line from INFILE file and parse it
!       inputs:
!           InFile
!           Line
!       outputs:
!           ReadLine    - number of word in this line
!------------------------------------------------------------------
!           Line_words() is filled with words in this line
!------------------------------------------------------------------
        USE Warehouse_INT
        Implicit none
        Integer :: i, j, Line_length
!   initialize Line_words
        Line%Line_words = ' '
!   read one line from Document.IN
        read(InFile,'(a)',err=901,end=902) Line%Line_original
        Line%LineNumber = Line%LineNumber + 1
!   Strip horizontal tab character (Ctrl-I = ASCII Character 9)
        do i = 1,256
            if(ichar(Line%Line_original(i:i)).eq.9) Line%Line_original(i:i) = ' '
        end do
!   Strip comments
        do i = 1, 256
             if(Line%Line_original(i:i).eq.'!') then
                  Line%Line_original(i:256) = ' '
                  exit
             end if
        end do
!   Strip leading blanks
        Line%Line_original = adjustl(Line%Line_original)
        Line_length = len_trim(Line%Line_original)
        if(Line_length.eq.0) then   ! blank line or comment line
            ReadLine = 0
            return
        else
            ReadLine = 1
        end if
!   Remove extra blanks before separator
        i = 1
        do while(i.lt.Line_length)
             do while(Line%Line_original(i:i).eq.' ' .and. (Line%Line_original(i+1:i+1).eq.' ' &
                        .or. Line%Line_original(i+1:i+1).eq.','))
                  Line%Line_original(i:Line_length-1) = Line%Line_original(i+1:Line_length)
                  Line%Line_original(Line_length:Line_length) = ' '
                  Line_length = Line_length - 1
             end do
             i = i + 1
        end do
!   Remove extra blanks after ',' and parse Line%Line_original to Line%Line_words
        j = 1
        i = 1
        do while(i.lt.Line_length) !i = 1, Line_length-2
             if(Line%Line_original(i:i).eq.',' .and. Line%Line_original(i+1:i+1).eq.' ') then
                  Line%Line_original(i+1:Line_length-1) = Line%Line_original(i+2:Line_length)
                  Line%Line_original(Line_length:Line_length) = ' '
                  Line_length = Line_length - 1
             end if
             if(Line%Line_original(i+1:i+1).eq.',' .or. Line%Line_original(i+1:i+1).eq.' ') then
                  Line%Line_words(ReadLine) = Line%Line_original(j:i)
                  j = i + 2
                  ReadLine = ReadLine + 1
             end if
             i = i + 1
        end do
        Line%Line_words(ReadLine) = Line%Line_original(j:Line_length)
!
        return
!   read error encountered
  901   write(*,"('error on reading, line ',i5)") Line%LineNumber
        stop 'error on reading'
!   EOF encountered
  902   print *, 'End of File!'
        stop 'encounter EOF when reading without seeing END'
!}}}
        end function ReadLine
!
        Integer function KeyWord(kw,nbkw)
!{{{
!------------------------------------------------------------------
!   purpose: read one line from INPUT file and identify the keyword
!   inputs:
!       kw      - keyword set
!       nbkw    - number of keyword
!   outputs:
!       Keyword - position of inputed keyword
!------------------------------------------------------------------
        USE Warehouse_INT
        implicit none
        character(20),intent(in) :: kw(*)
        Integer,intent(in) :: nbkw
        Integer :: i,j
!
        KeyWord = -1
!
        do while(Line%WordNumber.eq.0)
             Line%WordNumber = ReadLine()
             Line%ReadNumber = 0
        end do
        j=LEN_TRIM(Line%Line_words(Line%ReadNumber+1))
        do i = 1, nbkw
             if(pcomp(Line%Line_words(Line%ReadNumber+1),kw(i),j)) then
                  KeyWord = i
                  Line%WordNumber = Line%WordNumber - 1
                  Line%ReadNumber = Line%ReadNumber + 1
                  exit
             end if
        end do
        if(KeyWord.eq.-1) then
             print *, '*** Warning *** invalid keyword ',Line%Line_words(Line%ReadNumber+1),'on line ',Line%LineNumber
             STOP
        end if
!}}}
        end function KeyWord
!
        logical function pcomp(a,b,n)
!{{{
!c      * * F E A P * * A Finite Element Analysis Program
!c....  Copyright (c) 1984-2000: Robert L. Taylor
!c-----[--.----+----.----+----.-----------------------------------------]
!c      Purpose: Compare character strings for match
!c               Ignores upper/lower case differences.
!c      Inputs:
!c         a(*)   - Character string 1
!c         b(*)   - Character string 2
!c         n      - Number of characters to compare
!c      Outputs:
!c         pcomp  - Flag, true if a = b
!c-----[--.----+----.----+----.-----------------------------------------]
        implicit  none
        integer   n, inc, i, ia,ib,j
        character a*(*),b*(*)

        pcomp = .false.

!   Compute increment between an upper and lower case letter
        inc = ichar('A') - ichar('a')

!   Compare for match
        j=LEN_TRIM(b)
        if(n.NE.j) return
        do i = 1,n
             ia = ichar(a(i:i))
             ib = ichar(b(i:i))
!   Test all permutations of characters for match
             if(ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
        end do

        pcomp = .true.
!}}}
        end function pcomp

        subroutine GetString(mystring)
!{{{
!------------------------------------------------------------------
!   purpose: Get a string
!------------------------------------------------------------------
        USE Warehouse_INT
        Implicit none
        CHARACTER(256),Intent(INOUT) :: mystring
!
        do while(Line%WordNumber.eq.0)
             Line%WordNumber = ReadLine()
             Line%ReadNumber = 0
        end do
!
        mystring = Line%Line_original
!}}}
        end subroutine GetString
!
        integer function GetInt()
!{{{
!------------------------------------------------------------------
!   purpose: Get an integer number
!   inputs:
!   outputs:
!       GetInt  - the integer number
!------------------------------------------------------------------
        USE Warehouse_INT
        implicit none

        do while(Line%WordNumber.eq.0)
             Line%WordNumber = ReadLine()
             Line%ReadNumber = 0
        end do

        if(isNumber(Line%Line_words(Line%ReadNumber+1))) then
             read(Line%Line_words(Line%ReadNumber+1),*) GetInt ! change string to digit
             Line%WordNumber = Line%WordNumber - 1
             Line%ReadNumber = Line%ReadNumber + 1
        else
             print *, '*** Warning *** Non Numeric Field ',Line%Line_words(Line%ReadNumber+1),'on line ',Line%LineNumber
             Line%WordNumber = Line%WordNumber - 1
             Line%ReadNumber = Line%ReadNumber + 1
        end if
!}}}
        end function GetInt

        real(8) function GetReal()
!{{{
!------------------------------------------------------------------
!   purpose: Get a real number
!   inputs:
!   outputs:
!       GetReal - the real number
!------------------------------------------------------------------
        USE Warehouse_INT
        implicit none

        do while(Line%WordNumber.eq.0)
             Line%WordNumber = ReadLine()
             Line%ReadNumber = 0
        end do

        if(isNumber(Line%Line_words(Line%ReadNumber+1))) then
             read(Line%Line_words(Line%ReadNumber+1),*) GetReal
             Line%WordNumber = Line%WordNumber - 1
             Line%ReadNumber = Line%ReadNumber + 1
        else
             print *, '*** Warning *** Non Numeric Field ',  &
                Line%Line_words(Line%ReadNumber+1),'on line ',Line%LineNumber
                Line%WordNumber = Line%WordNumber - 1
                Line%ReadNumber = Line%ReadNumber + 1
        end if
!}}}
        end function GetReal
!
        logical function isNumber(num)
!{{{
!------------------------------------------------------------------
!   purpose: verify a numeric field
!------------------------------------------------------------------
        USE Warehouse_INT
        implicit none
        character(20),intent(In) :: num
        character a

        a = num(1:1)
        isNumber = .false.
        if(a.eq.'-' .or. (ichar(a).gt.47 .and. ichar(a).lt.58)) isNumber = .true.
!}}}
        end function isNumber
!
        SUBROUTINE NewLine()
!{{{
!------------------------------------------------------------------
!   purpose: start a new line
!------------------------------------------------------------------
        USE Warehouse_INT
        IMPLICIT NONE
             Line%WordNumber = ReadLine()
             Line%ReadNumber = 0
!}}}
        END SUBROUTINE NewLine
!
        SUBROUTINE ReadFile(Address)
!------------------------------------------------------------------
!   purpose: ReadFile Information
!------------------------------------------------------------------
!{{{ 
        USE Warehouse_INT
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: Address
!
        InFile = Address
        Line%LineNumber = 0
        Line%WordNumber = 0
!}}}
        END SUBROUTINE ReadFile
!        
        end module INTERPRETER
! 


