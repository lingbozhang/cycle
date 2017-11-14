! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Input module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
module INPUTM            
    USE INTERPRETER
    USE WAREHOUSE
    USE SOLVERM
    IMPLICIT NONE
    PUBLIC  :: INPUT    
    PRIVATE :: GLOBALD
    PRIVATE :: NODED
    PRIVATE :: BOUNDD
    PRIVATE :: LOADD
    PRIVATE :: INITLD
    PRIVATE :: MATERIALD
!
    PRIVATE :: BOUNDARY_IN
    PRIVATE :: LINED
    PRIVATE :: SURFACED
!   ELEMENTS
    PRIVATE :: ELEMENT_IN
    PRIVATE :: TRUSSD
    PRIVATE :: MFS_TRUSS_ASSD
    PRIVATE :: QUADSD
!
    PRIVATE :: TRUSS_MFSD
    PRIVATE :: QUADS_MFSD
    PRIVATE :: MFS_2D_ASSD
!
    PRIVATE :: OVLP_QUADSD
    PRIVATE :: OVLP_QUADS_ASSD
    PRIVATE :: COLHT_OVLP_QUADS
!
    PRIVATE :: OFE_TRUSSD
    PRIVATE :: COLHT_OFE_TRUSS
!
    PRIVATE :: OFE1_TRUSSD
!
    PRIVATE :: OFE1_QUADSD
    PRIVATE :: COLHT_OFE1_QUADS
    PRIVATE :: COLHT_OFE1_COUPLE
    PRIVATE :: COLHT_OFE1_COUPLE2
!   FINITE ELEMENTS WITH COVERS
    PRIVATE :: FEMC_QUADSD
    PRIVATE :: COLHT_FEMC_QUADS
!   PRIVATE :: OFE1_LINED
    PRIVATE :: TRUSS_OVLPD
    PRIVATE :: TRUSS_OVLPD2
    PRIVATE :: TRUSS_OVLPD3
    PRIVATE :: OVLP_ASSD
    PRIVATE :: OVLP_ASSD3
    PRIVATE :: COLHT_OVLP
    PRIVATE :: COLHT_OVLP3
    PRIVATE :: WEIND_CAL
!
    PRIVATE :: SURFACE_IN
    PRIVATE :: PTPTD
!   OVERLAPPING RECTANGULAR ELEMENTS
    PRIVATE :: OFE_RECT
    PRIVATE :: COLHT_OFE_RECT
    PRIVATE :: COLHT_OFE_CPL

!------------------------------------------
!   3D ELEMENTS
!------------------------------------------
!   HEX8
    PRIVATE :: HEX8
!
    PRIVATE :: OFE_TETRA
    PRIVATE :: COLHT_OFE_TETRA
! 
!    PRIVATE :: SURFACE_IN
!    PRIVATE :: POINTD
    CONTAINS
!
    SUBROUTINE INPUT(IIN,IOUT,ITEM3,ITEM1,ITEM2)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *     INPUTM                                        
! *      READ AND STORE INPUT DATA                    
! *      INPUT:                                       
! *            1. GLOBAL   DATA                                 
! *            2. NODE     DATA                         
! *            3. LOAD     DATA  
! *            4. INITIAL  DATA
! *            5. MATERIAL DATA
! *            6. LINEAR ELEMENT DATA
! *            7. SURFACE  DATA
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *         
!
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: IIN,IOUT,ITEM3,ITEM1,ITEM2
!        
        INTEGER :: I
        INTEGER,PARAMETER :: KWNUM=10
        CHARACTER(LEN=20),DIMENSION(KWNUM) :: KWINPUT ! KEYWORDS OF DIFFERENT INPUT PHASE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
        KWINPUT(1) = '\GLOBAL'
        KWINPUT(2) = '\NODE'
        KWINPUT(3) = '\LOAD'
        KWINPUT(4) = '\BOUND'
        KWINPUT(5) = '\INITIAL'
        KWINPUT(6) = '\MATERIAL'
        KWINPUT(7) = '\BOUNDARY'
        KWINPUT(8) = '\ELEMENT'
        KWINPUT(9) = '\SURFACE'
        KWINPUT(10) = '\END{INPUT}'     
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        I = KeyWord(KWINPUT,KWNUM)
        DO WHILE(I<=KWNUM)
           WRITE(IOUT,1000) KWINPUT(I)
           SELECT CASE(I)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE GLOBAL DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(1)
              CALL GLOBALD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE NODE   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(2)
              CALL NODED(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE LOAD   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(3)
              CALL LOADD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE LOAD   DATA               *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(4)
              CALL BOUNDD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE INITIAL DATA              *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(5)
              IND_INITIAL=1
              CALL INITLD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE MATERIAL DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(6)
              CALL MATERIALD(IIN,IOUT)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE BOUNDARY  DATA            *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(7)
              CALL BOUNDARY_IN(IIN,IOUT,ITEM3)
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE ELEMENT  DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(8)
              CALL ELEMENT_IN(IIN,IOUT,ITEM1)   
!     * * * * * * * * * * * * * * * * * * * * * *
!     *         STORE MATERIAL DATA             *
!     * * * * * * * * * * * * * * * * * * * * * *
           CASE(9)
              IND_SURFACE=1
              CALL SURFACE_IN(IIN,IOUT,ITEM2) 
           CASE(10)
              RETURN          
           CASE DEFAULT
              WRITE(*,*) "ERROR, INPUT OPERATION PHASE DOES NOT EXIST!"
              STOP
           END SELECT
           I = KeyWord(KWINPUT,KWNUM)
        END DO
!
        RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 1000   FORMAT(/,A20)
!}}}
    END SUBROUTINE INPUT
!
    SUBROUTINE GLOBALD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE GLOBALD 
! *      READ AND STORE GLOBAL DATA                            
! *      INPUT:                                       
! *    	   HED    -   The master heading information for use in labeling the output
! *        NUMNP  -   Total number of nodal points
! *				0 : program stop
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT NONE 
        INTEGER,INTENT(IN)    :: IIN,IOUT
!
        INTEGER :: I,K
        INTEGER,PARAMETER :: NumKwd = 16
        CHARACTER(LEN=20),DIMENSION(NUMKWD) :: KWGLOBAL ! KEYWORDS OF DIFFERENT  PHASE 
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
        KWGLOBAL(1) = 'PROJECT:'
        KWGLOBAL(2) = 'NUMNP:'
        KWGLOBAL(3) = 'NLOAD:'
        KWGLOBAL(4) = 'NBOUND:'
        KWGLOBAL(5) = 'NBGROUP:'
        KWGLOBAL(6) = 'NMAT:'
        KWGLOBAL(7) = 'NEGROUP:'
        KWGLOBAL(8) = 'NSPAIR:'
        KWGLOBAL(9) = 'DYNAMIC:'
        KWGLOBAL(10) = 'CONTACT:'
        KWGLOBAL(11) = 'NITP_IN:'
        KWGLOBAL(12) = 'NIDIV_IN:'
        KWGLOBAL(13) = 'NITP_BC:'
        KWGLOBAL(14) = 'NIDIV_BC:'
        KWGLOBAL(15) = 'PENALTY:'
        KWGLOBAL(16) = '\END'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
        I = KeyWord(KWGLOBAL,NumKwd)
        DO WHILE(I<=NumKwd)
          SELECT CASE(I)
          CASE(1)
             CALL GetString(HEAD)
             CALL NewLine()
             WRITE(IOUT,1000) HEAD
          CASE(2)
             K = GetInt()
             NUMNP = K
          CASE(3)
             K = GetInt()
             NLOAD = K
          CASE(4)
             K = GetInt()
             NBOUND = K
          CASE(5)
             K = GetInt()
             NBGROUP = K
          CASE(6)
             K = GetInt()
             NMAT = K
          CASE(7)
             K = GetInt()
             NEGROUP = K
          CASE(8)
             K = GetInt()
             NSPAIR = K
          CASE(9)
             K = GetInt()
             DYNAMIC = K
          CASE(10)
             K = GetInt()
             CONTACT = K
          CASE(11)
             K = GetInt()
             NITP_IN=K
          CASE(12)
             K = GetInt()
             NIDIV_IN=K
          CASE(13)
             K = GetInt()
             NITP_BC=K
          CASE(14)
             K = GetInt()
             NIDIV_BC=K
          CASE(15)
             K = GetInt()
             PENALTY=K 
          CASE(16)
             EXIT
          CASE DEFAULT
             WRITE(*,*) "ERROR, UNEXPECTED GLOBAL DATA"
             STOP
          END SELECT
          I = KeyWord(KWGLOBAL,NumKwd)
        END DO 
        IF(NUMNP.EQ.0) THEN
           WRITE(*,*) "ERROR,NUMBER OF NODES .EQ. ZERO!"
           STOP
        END IF
!
        WRITE(IOUT,2000) NUMNP,NLOAD,NBOUND,NMAT,NBGROUP,NEGROUP,NSPAIR,PENALTY,DYNAMIC,CONTACT
        WRITE(IOUT,2010) NITP_IN,NIDIV_IN,NITP_BC,NIDIV_BC
!
        RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A40)
 2000   FORMAT(/,&          
          ' C O N T R O L   I N F O R M A T I O N',//,&
          '      NUMBER OF NODAL POINTS ',10(' .'),' (NUMNP)   = ',I5,//,&
          '      NUMBER OF LOAD     ',12(' .'),    ' (NLOAD)   = ',I5,//,&
          '      NUMBER OF EBC      ',12(' .'),    ' (NBOUND)  = ',I5,//,&
          '      NUMBER OF MATERIAL SET ',10(' .'),' (NMAT)    = ',I5,//,&
          '      NUMBER OF BOUNDARY GROUP',10(' .'),'(NBGROUP) = ',I5,//,&
          '      NUMBER OF ELEMENT GROUP',10(' .'),'(NEGROUP)  = ',I5,//,&
          '      NUMBER OF SURFACE PAIR ',10(' .'),'(NSPAIR)   = ',I5,//,&
          '      PENALTY FOR DIRICHLET BCS',10(' .'),'(PENALTY)= ',I5,//,&
          '                      0 --- PENALTY MODE OFF          ',/,&
          '                      1 --- PENALTY MODE ON           ',/,&
          '      DYNAMIC SIMULATION MODE',10(' .'),'(DYNAMIC)  = ',I5,//,&
          '                      0 --- DYNAMIC MODE OFF          ',/,&
          '                      1 --- DYNAMIC MODE ON           ',/,&
          '      CONTACT SIMULATION MODE',10(' .'),'(CONTACT)  = ',I5,//,&
          '                      0 --- CONTACT MODE OFF          ',/,&
          '                      1 --- CONTACT MODE ON           ',//)
 2010   FORMAT(/,&
          ' I N T E G R A T I O N   I N F O R M A T I O N (MFS)',//,&
          '      NUMBER OF GAUSS-LEGENDRE INTEGRATION POINTS (INNER) ',5(' .'),' (NITP_IN)   = ',I5,//,&
          '      NUMBER OF INTEGRATION SECTIONS (INNER)              ',5(' .'),' (NIDIV_IN)  = ',I5,//,&
          '      NUMBER OF GAUSS-LEGENDRE INTEGRATION POINTS (BC)    ',5(' .'),' (NITP_BC)   = ',I5,//,&
          '      NUMBER OF INTEGRATION SECTIONS (BC)                 ',5(' .'),' (NIDIV_BC)  = ',I5,//)  
!}}}
    END SUBROUTINE GLOBALD
!
    SUBROUTINE NODED(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *     SUBROUTINE NODED
! *     INPUT AND STORE NODE DATA                            
! *      INPUT:     
! *        ID(3,NUMNP) : Boundary condition codes (0=free,1=deleted)
! *        X(NUMNP)    : X coordinates
! *	   Y(NUMNP)    : Y coordinates
! *	   Z(NUMNP)    : Z coordinates
! *        NPOR(NUMNP) : Number of the polynomial order
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
       INTEGER :: I, J
! Internal variables for node input data
       INTEGER :: N,KN
       INTEGER :: IDD(3)
       REAL(8) :: XX,YY,ZZ
       INTEGER :: NDOFX,NDOFY,NDOFZ
! Assistant variables
       INTEGER :: NCOUNT
       INTEGER :: KNOLD,NOLD
       REAL(8) :: XNUM, DX, DY, DZ
       INTEGER :: NUM, NUMN
       INTEGER :: K,KK
       INTEGER :: ALLocateStatus
!-------- Allocate memory
       ALLOCATE (ID(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (X(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (Y(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (Z(NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (NDOF(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! RESET NODE DATA
       DO I=1,NUMNP
          X(I)=0.0
          Y(I)=0.0
          Z(I)=0.0
          DO J=1,3
             ID(J,I)=0
             NDOF(J,I)=0
          END DO
       END DO  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2000)
       WRITE (IOUT,2010)
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       NCOUNT=0                  !RECORD NUMBER OF INPUTTED NODES     
       KNOLD=0
       NOLD=0
!
       N = 0
       READ (IIN,*) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
       IF(NDOFX.GT.NPOR_MAX) NPOR_MAX=NDOFX
       IF(NDOFY.GT.NPOR_MAX) NPOR_MAX=NDOFY
       IF(NDOFZ.GT.NPOR_MAX) NPOR_MAX=NDOFZ
!
       DO  WHILE(N.NE.0)  
          IF(N.GT.NUMNP) THEN
            WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP!"
            STOP
          END IF
          ID(:,N)=IDD
          X(N)=XX
          Y(N)=YY
          Z(N)=ZZ
          NDOF(1,N)=NDOFX
          NDOF(2,N)=NDOFY
          NDOF(3,N)=NDOFZ
          NCOUNT=NCOUNT+1
          WRITE (IOUT,2030) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
            IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                 XNUM=NUM
                 DX=(X(N)-X(NOLD))/XNUM
                 DY=(Y(N)-Y(NOLD))/XNUM
                 DZ=(Z(N)-Z(NOLD))/XNUM
                 K=NOLD
                 DO J=1,NUMN
                    KK=K
                    K=K + KNOLD
                    X(K) = X(KK) + DX
                    Y(K) = Y(KK) + DY
                    Z(K) = Z(KK) + DZ
                    NDOF(1,K)=NDOF(1,KK)
                    NDOF(2,K)=NDOF(2,KK)
                    NDOF(3,K)=NDOF(3,KK)
                    NCOUNT=NCOUNT+1
                    DO I=1,3
                       ID(I,K) = ID(I,KK)
                    END DO
                 END DO
              END IF
            ENDIF
!        
           NOLD=N
           KNOLD=KN
           READ (IIN,*) N,(IDD(I),I=1,3),XX,YY,ZZ,NDOFX,NDOFY,NDOFZ,KN
           IF(NDOFX.GT.NPOR_MAX) NPOR_MAX=NDOFX
           IF(NDOFY.GT.NPOR_MAX) NPOR_MAX=NDOFY
           IF(NDOFZ.GT.NPOR_MAX) NPOR_MAX=NDOFZ
!
       END DO
       IF(NCOUNT.NE.NUMNP) THEN
         WRITE(*,*) "ERROR, NUMBER OF NODES .LT. NUMNP!"
         STOP
       END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write complete nodal data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015)
       WRITE (IOUT,2020)
       DO N=1,NUMNP
          WRITE (IOUT,2030) N,(ID(I,N),I=1,3),X(N),Y(N),Z(N),(NDOF(I,N),I=1,3),KN
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Number unknowns
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       NEQ=0
       DO N=1,NUMNP
          DO I=1,3           
              IF (ID(I,N) .EQ. 0) THEN
                 NEQ=NEQ + 1
                 ID(I,N)=NEQ
                 IF(NDOF(I,N).EQ.0) THEN
                    WRITE(*,*) "ERROR, NDOF SHOULD .GE. 0!", N
                    STOP
                 ELSE
                    NEQ=NEQ+NDOF(I,N)-1
                 ENDIF
              ELSE 
                 ID(I,N)=0 
                 IF(NDOF(I,N).NE.0) THEN
                    WRITE(*,*) "ERROR, NDOF DOES NOT EQUAL TO 0!", N
                    STOP
                 ENDIF
              END IF             
          END DO
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write equation numbers
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2040) (N,(ID(I,N),I=1,3),N=1,NUMNP)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-------- Allocate memory for load vector 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ALLOCATE (R(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (NRF(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (RESS(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (RNAT(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NEQ
          R(I)=0.0D0
          NRF(I)=1
       END DO     
!
       DO I=1,NUMNP
          DO J=1,3
             RNAT(J,I)=0.
             RESS(J,I)=0.
          ENDDO
       ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2000  FORMAT(//,' N O D A L   P O I N T   D A T A',/)
 2010  FORMAT(' INPUT NODAL DATA',//)
 2015  FORMAT(//,' GENERATED NODAL DATA',//)
 2020  FORMAT('  NODE',10X,'BOUNDARY',25X,'NODAL POINT',35X,'NODE',/,&
       ' NUMBER     CONDITION  CODES',21X,'COORDINATES',32X,&
       'GENERATING',/,95X,'CODE',/,&
       15X,'X    Y    Z',15X,'X',12X,'Y',12X,'Z',11X,'NDOF(3,N)',8X,'KN')
 2030  FORMAT (I5,6X,3I5,6X,ES14.7,1X,ES14.7,1X,ES14.7,3X,3I5,6X,I4)
 2040  FORMAT(//,' EQUATION NUMBERS',//,'   NODE',9X,&
       'DEGREES OF FREEDOM',/,'  NUMBER',//,&
       '     N',13X,'X      Y      Z',/,(1X,I5,9X,I5,2X,I5,2X,I5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE NODED
!
    SUBROUTINE LOADD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE LOADD 
! *      READ AND STORE LOAD DATA                            
! *      INPUT: 
! *       LOAD CASE
! *         IDIRN(NLOAD): LOAD DIRECTION CODES (1:X DIRECTION,2:Y DIRECTION; 3: Z DIRECTION)
! *         FLOAD(NLOAD): MAGNITUDE OF LOAD
! *         NFUNC(NLOAD): LOAD FUNCTION CODES 
! *               1: CONSTANT
! *               2: LINEAR
! *               3: TRIGONOMETRIC
! *         R(NEQ)         : NODAL LOAD DATA
! *         NRF(NEQ)    : LOAD FUNCTION CODES                               
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
       INTEGER :: I,J
!--------------------------------------------------
!      Internal variables
!--------------------------------------------------
       INTEGER :: N
       INTEGER :: IDIRNN,NFUNCC
       REAL(8) :: FLOADD
!--------------------------------------------------
!      Assistant variables
!--------------------------------------------------
       INTEGER :: NumNEQ=0
!
       WRITE (IOUT,2015) NLOAD
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       N = 0
       I=0
       READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       DO WHILE(N.NE.0)
          IF(N.GT.NUMNP) THEN
            WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP IN LOAD_DATA!"
            STOP
          END IF
          I=I+1        
          IF(I.GT.NLOAD) THEN
            WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .GT. NLOAD!"
            STOP
          END IF
          IF(NFUNCC.EQ.0) NFUNCC=1      
          WRITE (IOUT,2030) N,IDIRNN,FLOADD,NFUNCC
!
          NumNEQ=ID(IDIRNN,N)
          RNAT(IDIRNN,N)=FLOADD
          IF (NumNEQ.NE.0) THEN
             R(NumNEQ)=FLOADD
             NRF(NumNEQ)=NFUNCC
!             IF(NDOF(IDIRNN,N).GE.2) THEN
!                DO J=1,NDOF(IDIRNN,N)-1
!                  R(NumNEQ+J)=FLOADD
!                  NRF(NumNEQ+J)=NFUNCC
!                ENDDO
!             ENDIF
          ENDIF
!
          READ (IIN,*) N,IDIRNN,FLOADD,NFUNCC
       END DO
       IF(I.NE.NLOAD) THEN
          WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .LT. NLOAD!"
          STOP
       END IF
!
       WRITE(IOUT,2040)
       DO I=1,NEQ
          WRITE(IOUT,2050) I,R(I),NRF(I)  
       END DO
!
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2015  FORMAT(/,'NUMBER OF LOADS =',I5,/)
 2020  FORMAT(' NODE         IDRIN         LOAD           NFUNC',/,&
              'NUMBER        NUMBER      MAGNITUDE        NUMBER')
 2030  FORMAT(I5,9X,I5,3X,F13.3,5X,I6)
 2040  FORMAT(/,'EQUATION        LOAD           FUNCTION',/&
              ' NUMBER       MAGNITUDE          CODES ')
 2050  FORMAT(I5,5X,ES18.6,7X,I6)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE LOADD
!
    SUBROUTINE BOUNDD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE BOUNDD 
! *      READ AND STORE DIRICHLET BOUNDARY DATA                            
! *      INPUT: 
! *     BOUND CASE
! *         IDIRN_BC(NBOUND): BOUNDARY DIRECTION CODES (1:X DIRECTION,2:Y DIRECTION; 3: Z DIRECTION)
! *         FBOUND(NBOUND): MAGNITUDE OF ESSENTIAL BOUNDARY
! *         NFUNC_BC(NBOUND): BOUNDARY FUNCTION CODES 
! *               1: CONSTANT
! *               2: LINEAR
! *               3: TRIGONOMETRIC                            
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
       INTEGER :: I,J
!--------------------------------------------------
!      Internal variables
!--------------------------------------------------
       INTEGER :: N
       INTEGER :: IDIRNN_BC,NFUNCC_BC
       REAL(8) :: FBOUNDD
!--------------------------------------------------
!      Assistant variables
!--------------------------------------------------
       INTEGER :: NumNEQ=0    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015) NBOUND
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       N = 0
       I=0
       READ (IIN,*) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
       DO WHILE(N.NE.0)
         IF(N.GT.NUMNP) THEN
           WRITE(*,*) "ERROR,NODE NUMBER .GT. NUMNP IN LOAD_DATA!"
           STOP
         END IF
         I=I+1        
         IF(I.GT.NBOUND) THEN
            WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .GT. NLOAD!"
            STOP
         END IF
         IF(NFUNCC_BC.EQ.0) NFUNCC_BC=1
!       
! write
!
!         WRITE (IOUT,2030) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
!
         NumNEQ=ID(IDIRNN_BC,N)         
         RESS(IDIRNN_BC,N)=FBOUNDD
         IF (NumNEQ.NE.0) THEN
            R(NumNEQ)=FBOUNDD
            NRF(NumNEQ)=NFUNCC_BC
            IF(NDOF(IDIRNN_BC,N).GE.2) THEN
               DO J=1,NDOF(IDIRNN_BC,N)-1
                  R(NumNEQ+J)=FBOUNDD
                  NRF(NumNEQ+J)=NFUNCC_BC
               ENDDO
            ENDIF
         END IF
!
         READ (IIN,*) N,IDIRNN_BC,FBOUNDD,NFUNCC_BC
       END DO
       IF(I.NE.NBOUND) THEN
          WRITE(*,*) "ERROR, NUMBER OF INPUTTED LOADS .LT. NLOAD!"
          STOP
       END IF
!
       WRITE(IOUT,2040)
       DO I=1,NEQ
!          WRITE(IOUT,2050) I,R(I),NRF(I)  
       END DO
!
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2015  FORMAT(/,'NUMBER OF BOUNDARY CONDITIONS =',I5,/)
 2020  FORMAT(' NODE        IDRIN_BC     BOUNDARY        NFUNC_BC',/,&
              'NUMBER        NUMBER      MAGNITUDE        NUMBER')
 2030  FORMAT(I5,9X,I5,3X,F13.3,5X,I6)
 2040  FORMAT(/,'EQUATION        LOAD           FUNCTION',/&
              ' NUMBER       MAGNITUDE          CODES ')
 2050  FORMAT(I5,5X,F13.3,7X,I6)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE BOUNDD
!
    SUBROUTINE INITLD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE INITLDT
! *      READ AND STORE INITIAL CONDITIONS                           
! *      INPUT: 
! *        UI(3,NUMNP)    : X Y Z DISPLACEMENT
! *	   VI(3,NUMNP)    : X Y Z VELOCITY
! *        U(NEQ)         : DISPLACEMENT VARIABLES
! *        V(NEQ)         : VELOCITY VARIABLES
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IIN,IOUT
!
       INTEGER :: I,K,J,KK
       REAL(8) :: UXX,UYY,UZZ,VXX,VYY,VZZ
       REAL(8) :: DUX,DUY,DUZ,DVX,DVY,DVZ
       INTEGER :: KN,KNOLD,N,NOLD
       INTEGER :: NUMN,NUM
!--------------------------------------------------
!       Assistant variables
!--------------------------------------------------
       INTEGER :: NCOUNT
       INTEGER :: ALLocateStatus
!-------- Allocate memory
       ALLOCATE (UI(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (VI(3,NUMNP),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (U(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (V(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
       ALLOCATE (ACC(NEQ),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DO I=1,NUMNP
          DO J=1,3
             UI(J,I)=0
             VI(J,I)=0
          END DO
       END DO
!
       DO I=1,NEQ
          U(I)=0
          V(I)=0
          ACC(I)=0
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2000)
       WRITE (IOUT,2010)
       WRITE (IOUT,2020)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
       N = 0
       NCOUNT=0                  !RECORD NUMBER OF INPUTTED NODES     
       KNOLD=0
       NOLD=0

       READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ,KN  
       DO WHILE (N.NE.0)
          IF(N.GT.NUMNP) THEN
             WRITE(*,*) "ERROR, NODE NUMBER .GT. NUMNP!"
             STOP
          END IF
          UI(1,N)=UXX
          UI(2,N)=UYY
          UI(3,N)=UZZ
          VI(1,N)=VXX
          VI(2,N)=VYY
          VI(3,N)=VZZ  
          NCOUNT=NCOUNT+1
!----------------write out initial conditions
          WRITE (IOUT,2030) N,UI(1,N),UI(2,N),UI(3,N),VI(1,N),VI(2,N),VI(3,N)
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                 DUX=(UI(1,N)-UI(1,NOLD))/NUM
                 DUY=(UI(2,N)-UI(2,NOLD))/NUM
                 DUZ=(UI(3,N)-UI(3,NOLD))/NUM
                 DVX=(VI(1,N)-VI(1,NOLD))/NUM
                 DVY=(VI(2,N)-VI(2,NOLD))/NUM
                 DVZ=(VI(3,N)-VI(3,NOLD))/NUM
                 K=NOLD
                 DO J=1,NUMN
                    KK=K
                    K=K + KNOLD
                    UI(1,K) = UI(1,KK) + DUX
                    UI(2,K) = UI(2,KK) + DUY
                    UI(3,K) = UI(3,KK) + DUZ                
                    VI(1,K) = VI(1,KK) + DVX
                    VI(2,K) = VI(2,KK) + DVY
                    VI(3,K) = VI(3,KK) + DVZ
                    NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
!        
          NOLD=N
          KNOLD=KN
          READ (IIN,*) N,UXX,UYY,UZZ,VXX,VYY,VZZ,KN   
       END DO
       IF(NCOUNT.NE.NUMNP) THEN
          WRITE(*,*) "ERROR,NUMBER OF INPUTTED INITIAL DATA .LT. NUMNP!"
          STOP
       END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Write complete initial data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       WRITE (IOUT,2015)
       WRITE (IOUT,2020)
       DO N=1,NUMNP
         WRITE (IOUT,2030) N,UI(1,N),UI(2,N),UI(3,N),VI(1,N),VI(2,N),VI(3,N)
       END DO       
!----------------data transform
       DO N=1,NUMNP
          DO I=1,3
            IF(ID(I,N).NE.0) THEN
              K=ID(I,N)
              U(K)=UI(I,N)
              V(K)=VI(I,N)
            END IF
          END DO
       END DO
!----------------Write new Initial conditions
       WRITE(IOUT,2040)
       DO I=1,NEQ
          WRITE(IOUT,2050) I,U(I),V(I)  
       END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 2000  FORMAT(//,' N O D A L   I N I T I A L   D A T A',/)
 2010  FORMAT(' INPUT NODAL DATA',//)
 2015  FORMAT(//,' GENERATED INITIAL DATA',//)
 2020  FORMAT(' NODE',14X,'UX',11X,'UY',11X,'UZ',11X,'VX',11X,'VY',11X,&
              'VZ',/,' NUMBER')
 2030  FORMAT (I5,6X,6F13.3)
 2040  FORMAT(/,'EQUATION     INITIAL_DISP           INITIAL_VEL',/&
               ' NUMBER                                        ')
 2050  FORMAT(I5,5X,F13.3,7X,F13.3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!}}}
    END SUBROUTINE INITLD 
!
    SUBROUTINE MATERIALD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
! *     SUBROUTINE MATERIALD
! *      READ AND STORE MATERIAL DATA                           
! *      INPUT: 
! *        MATP(7,NMAT)  :  MATERIAL PROPERTIES
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: IIN,IOUT
! Internal variables        
        INTEGER :: N=0
        REAL(8) :: E=0,PR=0,DEN=0,PAR1=0
! Assistant variables
        INTEGER :: I,J
        INTEGER :: NCOUNT
        INTEGER :: AllocateStatus
!-------- Allocate memory
        ALLOCATE (MATP(7,NMAT),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
        DO I=1,NMAT
           DO J=1,7
              MATP(J,I)=0
           END DO
        END DO    
! SET DEFAULT VALUE
        NCOUNT=0
!
        WRITE (IOUT,2000)
        WRITE (IOUT,2010)
        WRITE (IOUT,2020)
!
        READ(IIN,*) N,E,PR,DEN,PAR1
        DO WHILE(N.NE.0)
           IF(N.GT.NMAT) THEN
                WRITE(*,*) "ERROR, MATERIAL SET NUMBER .GT. NMAT"
                STOP
           END IF
           NCOUNT=NCOUNT+1
           MATP(1,N)=E
           MATP(2,N)=PR
           MATP(3,N)=DEN
           MATP(4,N)=PAR1
           WRITE(IOUT,2030) N, E, PR,DEN,PAR1
           READ(IIN,*) N,E,PR,DEN,PAR1
        END DO
        IF(NCOUNT.NE.NMAT) THEN
            WRITE(*,*) "ERROR, NUMBER OF MATERIAL SET .LT. NMAT"
            STOP
        END IF
        RETURN
!
 2000  FORMAT(//,' M A T E R I A L    D A T A',/)
 2010  FORMAT(' INPUT MATERIAL DATA',//)
 2020  FORMAT(' SET ',22X,' E ',14X,' PR ',14X,'DEN',14X,'PAR1',/,&
              'NUMBER')
 2030  FORMAT(I6,10X,ES13.3,5X,F13.3,5X,F13.3,5X,F13.3)
!}}}
    END SUBROUTINE MATERIALD
!
    SUBROUTINE BOUNDARY_IN(IIN,IOUT,ITEM3)
!{{{
!--------------------------------------------------
!   SELECT BOUNDARY TYPE
!   NPAR()
!   LINE
!       NPAR(1) = 1
!       NPAR(2) = NUMBC
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM3
!   ASSISTANT VARIABLES
    INTEGER :: I,J,L,M
    INTEGER,PARAMETER :: KwBC=6
    CHARACTER(20),DIMENSION(KwBC) :: KwBoundary 
!
!---  ASSISTANT VARIABLES    
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwBoundary(1)='\LINE'
    KwBoundary(2)='\OFE1_LINE'
    KwBoundary(3)='\QUADS_LINE'
    KwBoundary(4)='\FEMC_LINE'
    KwBoundary(5)='\OFE_RECT'
    KwBoundary(6)='\OFE_TETRA'
!
    WRITE(IOUT,2000)
!
    IF(NBGROUP<1e-7) THEN
       WRITE(*,*) "ERROR, NBGROUP EQUALS TO ZEROS"
       STOP
    ENDIF
!
    DO J=1,NBGROUP
       I = KeyWord(KwBoundary,KwBC)
       WRITE(IOUT,1000) KwBoundary(I)
       SELECT CASE(I)
       CASE(1)
          NPAR=0
          NPAR(1)=1
!
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
   ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE(2)
          NPAR=0
          NPAR(1)=2
          OFE1BCONOFF=.TRUE.
!
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
   ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE(3)
          NPAR=0
          NPAR(1)=3
          QUADSBCONOFF=.TRUE.
!
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
   ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE(4)
          NPAR=0
          NPAR(1)=4
          FEMCBCONOFF=.TRUE.
!
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
   ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE(5)
          NPAR=0
          NPAR(1)=5
!
          CALL LINED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,2),M=1,NPAR(2)),((BCIND(L,M),L=1,4),M=1,NPAR(2)),((BCPAR(M,L),M=1,2),L=1,NPAR(2))
   ! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE(6)
          NPAR=0
          NPAR(1)=6
          CALL SURFACED(IIN,IOUT)
          WRITE(ITEM3) (NPAR(L),L=1,7)
          WRITE(ITEM3) ((BCED(L,M),L=1,3),M=1,NPAR(2)),((BCIND(L,M),L=1,6),M=1,NPAR(2)),((BCPAR(M,L),M=1,3),L=1,NPAR(2))
! Deallocate variables
          DEALLOCATE (BCED)
          DEALLOCATE (BCIND)
          DEALLOCATE (BCPAR) 
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
!
    END DO
!
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000 FORMAT(/,A20)
 2000 FORMAT(//,' B O U N D A R Y   I N F O R M A T I O N',///)
!}}}
    END SUBROUTINE BOUNDARY_IN
!
    SUBROUTINE LINED(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE LINED
! *      READ AND STORE LINE DATA
! *      LINE
! *         NPAR(1) = 1
! *         NPAR(2) = NUMBC
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    INTEGER :: N,II,JJ,KN,INDI,INDJ,FUNCI,FUNCJ
    INTEGER :: OINDI,OINDJ,OFUNCI,OFUNCJ
    REAL(8) :: FLOADI,FLOADJ
    REAL(8) :: OFLOADI,OFLOADJ
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    INTEGER :: NCOUNT
    INTEGER :: DII
!
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(2) :: KWLINE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWLINE(1) = 'NUMBC:'
    KWLINE(2) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
!
    I = KeyWord(KWLINE,2)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
       CASE(1)
           J = GetInt()
           NPAR(2)=J
       CASE(2)
           EXIT
       CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN LINED"
           STOP
       END SELECT
       I = KeyWord(KWLINE,2)
    END DO
!   
    IF(NPAR(2).EQ.0) THEN
       WRITE(*,*) "ERROR, NUMBER OF LINE BOUNDARY .EQ.0"
       STOP
    END IF
!
    WRITE(IOUT,1000) NPAR(2)
    NumBc=NPAR(2)
!
!-------- Allocate memory
!
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
    DO I=1,NumBc
       DO J=1,4
          BCIND(J,I)=0
       END DO
       BCED(1,I)=0
       BCED(2,I)=0
       BCPAR(1,I)=0.0D0
       BCPAR(2,I)=0.0D0
    END DO
!
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumBc) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       BCED(1,N)=II 
       BCED(2,N)=JJ
!       
       BCPAR(1,N)=FLOADI
       BCPAR(2,N)=FLOADJ
!
       BCIND(1,N)=INDI
       BCIND(2,N)=INDJ
       BCIND(3,N)=FUNCI
       BCIND(4,N)=FUNCJ
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              DII=(II-IOLD)/NUM
             DO J=1,NUMN
                K=K + KNOLD
!
                NII=NII+DII
                NJJ=NJJ+DII
!
                BCED(1,K)=NII 
                BCED(2,K)=NJJ
!
                BCPAR(1,K)=OFLOADI
                BCPAR(2,K)=OFLOADJ
!
                BCIND(1,K)=OINDI
                BCIND(2,K)=OINDJ
                BCIND(3,K)=OFUNCI
                BCIND(4,K)=OFUNCJ
!
               WRITE(IOUT,2050) K,NII,NJJ,OINDI,OINDJ,OFLOADI,OFLOADJ,OFUNCI,OFUNCJ                 
               NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
        WRITE(IOUT,2050) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
!         
        OINDI=INDI
        OINDJ=INDJ
        OFLOADI=FLOADI
        OFLOADJ=FLOADJ
        OFUNCI=FUNCI
        OFUNCJ=FUNCJ
!
        READ(IIN,*)  N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
    END DO
!
    IF(NCOUNT.NE.NumBc) THEN
        WRITE(*,*) "ERROR, BOUNDARY NUMBER .LT. NUMBC", NCOUNT
        STOP
    END IF
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE BOUNDARY TYPE: LINE ',/,&
              '    THE NUMBER OF LINE BOUNDARY      =',I5,/)
 2000 FORMAT(/,&
              ' BOUNDARY   ELEMENT  ELEMENT       IND(2,N)       FLOAD(2,N)       FUNC(2,N)',/,&
              ' NUMBER-N      I        J       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,2I5,2X,2ES13.3,2X,2I5)
!}}}
    END SUBROUTINE LINED
!
    SUBROUTINE SURFACED(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE LINED
! *      READ AND STORE LINE DATA
! *      LINE
! *         NPAR(1) = 1
! *         NPAR(2) = NUMBC
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    INTEGER :: N,II,JJ,KK
    INTEGER :: INDI,INDJ,INDK,FUNCI,FUNCJ,FUNCK
    REAL(8) :: FLOADI,FLOADJ,FLOADK
    INTEGER :: NCOUNT
!
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(2) :: KWLINE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWLINE(1) = 'NUMBC:'
    KWLINE(2) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
!
    I = KeyWord(KWLINE,2)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
       CASE(1)
           J = GetInt()
           NPAR(2)=J
       CASE(2)
           EXIT
       CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN LINED"
           STOP
       END SELECT
       I = KeyWord(KWLINE,2)
    END DO
!   
    IF(NPAR(2).EQ.0) THEN
       WRITE(*,*) "ERROR, NUMBER OF LINE BOUNDARY .EQ.0"
       STOP
    END IF
!
    WRITE(IOUT,1000) NPAR(2)
    NumBc=NPAR(2)
!
!-------- Allocate memory
!
    ALLOCATE (BCED(3,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(6,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(3,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
    BCED=0
    BCIND=0
    BCPAR=0.0D0
!
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, KK, INDI, INDJ, INDK, FLOADI, FLOADJ, FLOADK,  FUNCI, FUNCJ, FUNCK
    DO WHILE(N.NE.0)
       IF(N.GT.NumBc) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       BCED(1,N)=II 
       BCED(2,N)=JJ
       BCED(3,N)=KK
!       
       BCPAR(1,N)=FLOADI
       BCPAR(2,N)=FLOADJ
       BCPAR(3,N)=FLOADK
!
       BCIND(1,N)=INDI
       BCIND(2,N)=INDJ
       BCIND(3,N)=INDK
!
       BCIND(4,N)=FUNCI
       BCIND(5,N)=FUNCJ
       BCIND(6,N)=FUNCK
       WRITE(IOUT,2050) N, II, JJ, KK, INDI, INDJ,INDK, FLOADI, FLOADJ,FLOADK, FUNCI, FUNCJ,FUNCK
!
       READ(IIN,*)  N, II, JJ, KK, INDI, INDJ, INDK,FLOADI, FLOADJ, FLOADK,  FUNCI, FUNCJ, FUNCK 
    END DO
!
    IF(NCOUNT.NE.NumBc) THEN
        WRITE(*,*) "ERROR, BOUNDARY NUMBER .LT. NUMBC", NCOUNT
        STOP
    END IF
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE BOUNDARY TYPE: LINE ',/,&
              '    THE NUMBER OF LINE BOUNDARY      =',I5,/)
 2000 FORMAT(/,&
              ' BOUNDARY   ELEMENT  ELEMENT  ELEMENT      IND(3,N)       FLOAD(3,N)       FUNC(3,N)',/,&
              ' NUMBER-N      I        J       K ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4x,I5,2X,3I5,2X,3ES13.3,2X,3I5)
!}}}
    END SUBROUTINE SURFACED
!
!    SUBROUTINE OFE1_LINED(IIN,IOUT)
!!{{{
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!! *      SUBROUTINE OFE1_LINED
!! *      READ AND STORE LINE DATA
!! *      LINE
!! *         NPAR(1) = 2
!! *         NPAR(2) = NUMBC
!! *      INPUT:
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!    IMPLICIT NONE
!    INTEGER,INTENT(IN) :: IIN,IOUT
!!   ASSISTANT VARIABLES
!    INTEGER :: I,J,K
!    INTEGER :: NumBc
!    INTEGER :: N,II,JJ,KN,INDI,INDJ,FUNCI,FUNCJ
!    INTEGER :: OINDI,OINDJ,OFUNCI,OFUNCJ
!    REAL(8) :: FLOADI,FLOADJ
!    REAL(8) :: OFLOADI,OFLOADJ
!    INTEGER :: KNOLD,NOLD,NUM,NUMN
!    INTEGER :: IOLD,JOLD,NII,NJJ
!    INTEGER :: NCOUNT
!    INTEGER :: DII
!!
!    INTEGER :: AllocateStatus
!!
!    CHARACTER(20),DIMENSION(2) :: KWLINE
!! - - - - - - - - - - - - - - -  
!! -     VARIABLE DEFINITION   -
!! - - - - - - - - - - - - - - -
!    KWLINE(1) = 'NUMBC:'
!    KWLINE(2) = '\END'
!!   SET DEFAULT VALUE
!    NPAR(2)=0
!!
!    I = KeyWord(KWLINE,2)
!    DO WHILE(.TRUE.)
!       SELECT CASE(I)
!       CASE(1)
!           J = GetInt()
!           NPAR(2)=J
!       CASE(2)
!           EXIT
!       CASE DEFAULT
!           WRITE(*,*) "ERROR, ERROR INPUT IN OFE1_LINED"
!           STOP
!       END SELECT
!       I = KeyWord(KWLINE,2)
!    END DO
!!   
!    IF(NPAR(2).EQ.0) THEN
!       WRITE(*,*) "ERROR, NUMBER OF LINE BOUNDARY .EQ.0"
!       STOP
!    END IF
!!
!    WRITE(IOUT,1000) NPAR(2)
!    NumBc=NPAR(2)
!!
!!-------- Allocate memory
!!
!    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
!! ----------------------------------------------------------------
!    BCED=0
!    BCIND=0
!    BCPAR=0.0D0
!!
!    WRITE(IOUT,2000)
!    NCOUNT=0
!    READ(IIN,*) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
!    DO WHILE(N.NE.0)
!       IF(N.GT.NumBc) THEN
!          WRITE(*,*) "ERROR N .GT. NumEle"
!          STOP
!       END IF
!       NCOUNT=NCOUNT+1
!!
!       BCED(1,N)=II 
!       BCED(2,N)=JJ
!!       
!       BCPAR(1,N)=FLOADI
!       BCPAR(2,N)=FLOADJ
!!
!       BCIND(1,N)=INDI
!       BCIND(2,N)=INDJ
!       BCIND(3,N)=FUNCI
!       BCIND(4,N)=FUNCJ
!!
!       IF (KNOLD.NE.0) THEN
!          NUM=(N-NOLD)/KNOLD
!          NUMN=NUM - 1
!          IF (NUMN.GE.1) THEN
!             K=NOLD
!             NII=IOLD
!             NJJ=JOLD
!             DII=(II-IOLD)/NUM
!             DO J=1,NUMN
!                K=K + KNOLD
!!
!                NII=NII+DII
!                NJJ=NJJ+DII
!!
!                BCED(1,K)=NII 
!                BCED(2,K)=NJJ
!!
!                BCPAR(1,K)=OFLOADI
!                BCPAR(2,K)=OFLOADJ
!!
!                BCIND(1,K)=OINDI
!                BCIND(2,K)=OINDJ
!                BCIND(3,K)=OFUNCI
!                BCIND(4,K)=OFUNCJ
!!
!!
!                WRITE(IOUT,2050) K,NII,NJJ,OINDI,OINDJ,OFLOADI,OFLOADJ,OFUNCI,OFUNCJ                 
!               NCOUNT=NCOUNT+1
!             END DO
!          END IF
!        ENDIF
!        WRITE(IOUT,2050) N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ
!!     
!        NOLD=N
!        KNOLD=KN
!        IOLD=II
!        JOLD=JJ
!!         
!        OINDI=INDI
!        OINDJ=INDJ
!        OFLOADI=FLOADI
!        OFLOADJ=FLOADJ
!        OFUNCI=FUNCI
!        OFUNCJ=FUNCJ
!!
!        READ(IIN,*)  N, II, JJ, INDI, INDJ, FLOADI, FLOADJ, FUNCI, FUNCJ, KN
!    END DO
!!
!    IF(NCOUNT.NE.NumBc) THEN
!        WRITE(*,*) "ERROR, BOUNDARY NUMBER .LT. NUMBC", NCOUNT
!        STOP
!    END IF
!!
!    RETURN
!!--------------------------------------------------
! 1000 Format(/,&
!              '    THE BOUNDARY TYPE: LINE ',/,&
!              '    THE NUMBER OF LINE BOUNDARY      =',I5,/)
! 2000 FORMAT(/,&
!              ' BOUNDARY   ELEMENT  ELEMENT       IND(2,N)       FLOAD(2,N)       FUNC(2,N)',/,&
!              ' NUMBER-N      I        J       ',/)
! 2050 FORMAT(I5,6X,I5,4X,I5,2X,2I5,2X,2ES13.3,2X,2I5)
!!}}}
!    END SUBROUTINE OFE1_LINED
!
!**************************************************
!
!   SUBROUTINES FOR ELEMENT INPUT
!
!**************************************************
    SUBROUTINE ELEMENT_IN(IIN,IOUT,ITEM1)
!{{{
!--------------------------------------------------
!   SELECT ELEMENT TYPE
!   NPAR()
!   TRUSS
!       NPAR(1) = 1
!       NPAR(2) = MTYPE
!       NPAR(3) = NUMELE
!   QUADS
!       NPAR(1) = 2
!       NPAR(2) = ELEMENT TYPE
!                 0: AXISYMMETRIC
!                 1: PLANE STRAIN
!                 2: PLANE STRESS
!       NPAR(3) = MTYPE
!       NPAR(4) = NUMELE
!       NPAR(5) = NINTG
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM1
!
    INTEGER :: I,J,L,M
    INTEGER,PARAMETER :: KwEle=16
    CHARACTER(20),DIMENSION(KwEle) :: KwElement
!---  ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
    INTEGER(KIND=8) :: EFFORT
!---- ALLOCATE VARIABLES
    ALLOCATE (MHT(NEQ),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
    ALLOCATE (MHT_ND(NUMNP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
    MHT=0 
    MHT_ND=0    
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwElement(1)='\TRUSS'
    KwElement(2)='\QUADS' 
    KwElement(3)='\TRUSS_MFS'
    KwElement(4)='\QUADS_MFS'
    KwElement(5)='\COUP1'
    KwElement(6)='\TRUSS_OVLP'
    KwElement(7)='\TRUSS_OVLP2'
    KwElement(8)='\TRUSS_OVLP3'
    KwElement(9)='\OVLP_QUADS'
    KwElement(10)='\OFE_TRUSS'
    KwElement(11)='\OFE1_TRUSS'
    KwElement(12)='\OFE1_QUADS'
    KwElement(13)='\FEMC_QUADS'
    KwElement(14)='\OFE_TETRA'
    KwElement(15)='\OFE_RECT'
    KwElement(16)='\HEX8'
!
    WRITE(IOUT,2000)
!
    IF(NEGROUP.EQ.0) THEN
       WRITE(*,*) "ERROR, NEGROUP EQUALS TO ZEROS"
       STOP
    ENDIF
!
    DO J=1,NEGROUP
       I = KeyWord(KwElement,KwEle)
       WRITE(IOUT,1000) KwElement(I)
       SELECT CASE(I)
       CASE(1)
          NPAR=0
          NPAR(1)=1
          CALL TRUSSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((XYZ(L,M),L=1,6),M=1,NPAR(3)),((LM(L,M),L=1,6),M=1,NPAR(3)),(ELEPAR(1,L),L=1,NPAR(3))
          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
          DEALLOCATE (ELEPAR)
       CASE(2)
          NPAR=0
          NPAR(1)=2
          CALL QUADSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((XYZ(L,M),L=1,8),M=1,NPAR(4)),((LM(L,M),L=1,8),M=1,NPAR(4)),(ELEPAR(1,L),L=1,NPAR(4))
          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
          DEALLOCATE (ELEPAR)
        CASE(3)
          NPAR=0
          NPAR(1)=3
          CALL TRUSS_MFSD(IIN,IOUT)         
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(M,L),M=1,3),L=1,NPAR(3)),((ELEPAR(M,L),M=1,2),L=1,NPAR(3)),((LM(M,L),M=1,NPAR(4)),L=1,NPAR(3))
          ! Deallocate variables
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM)
       CASE(4)
          NPAR=0
          NPAR(1)=4
          CALL QUADS_MFSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(M,L),M=1,8),L=1,NPAR(4)),((ELEPAR(M,L),M=1,6),L=1,NPAR(4)),((LM(M,L),M=1,NPAR(5)),L=1,NPAR(4))
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM)
       CASE(5)
          NPAR=0
          NPAR(1)=5
          CALL COUPE1D(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(M,L),M=1,8),L=1,NPAR(4)),((ELEPAR(M,L),M=1,4),L=1,NPAR(4))          
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
       CASE(6)
          NPAR=0
          NPAR(1)=6
          CALL TRUSS_OVLPD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,2),M=1,NPAR(4)),(ELEPAR(1,M),M=1,NPAR(4)),((LM(L,M),L=1,NPAR(5)),M=1,NPAR(4)),&
          ((WEIND(L,M),L=1,NPAR(5)),M=1,NPAR(4)),((NDINT(L,M),L=1,2),M=1,NPAR(3)),(INTPAR(1,M),M=1,NPAR(3))
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM) 
          DEALLOCATE(WEIND)          
       CASE(7)
          NPAR=0
          NPAR(1)=7
          CALL TRUSS_OVLPD2(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,2),M=1,NPAR(4)),(ELEPAR(1,M),M=1,NPAR(4)),((LM(L,M),L=1,NPAR(5)),M=1,NPAR(4)),&
          ((NDINT(L,M),L=1,2),M=1,NPAR(3)),(INTPAR(1,M),M=1,NPAR(3))
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM) 
       CASE(8)
          NPAR=0
          NPAR(1)=8
          CALL TRUSS_OVLPD3(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,2),M=1,NPAR(4)),(ELEPAR(1,M),M=1,NPAR(4)),((LM(L,M),L=1,NPAR(5)),M=1,NPAR(4)),&
          ((NDINT(L,M),L=1,4),M=1,NPAR(3)),(INTPAR(1,M),M=1,NPAR(3)),(NumLM(L),L=1,NPAR(4))
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM) 
          DEALLOCATE(NumLM)
       CASE(9)
          NPAR=0
          NPAR(1)=9
          CALL OVLP_QUADSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,6),M=1,NPAR(4)),((ELEPAR(L,M),L=1,4),M=1,NPAR(4)),((LM(L,M),L=1,NPAR(5)),M=1,NPAR(4)),&
          ((NDINT(L,M),L=1,3),M=1,NPAR(3)),(INTPAR(1,M),M=1,NPAR(3)),(NumLM(L),L=1,NPAR(4))
!
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(LM) 
          DEALLOCATE(NumLM)
       CASE(10)
!------------------------------------------
          NPAR=0
          NPAR(1)=10
          CALL OFE_TRUSSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDINT(L,M),L=1,4),M=1,NPAR(3)),((INTPAR(L,M),L=1,2),M=1,NPAR(3))
!          
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
!------------------------------------------
       CASE(11)
!------------------------------------------
          NPAR=0
          NPAR(1)=11
          CALL OFE1_TRUSSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDINT(L,M),L=1,4),M=1,NPAR(3)),((INTPAR(L,M),L=1,2),M=1,NPAR(3))
!          
          DEALLOCATE(NDINT)
          DEALLOCATE(INTPAR)
       CASE(12)
          NPAR=0
          NPAR(1)=12
          CALL OFE1_QUADSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,5),M=1,NPAR(4)),((ELEPAR(L,M),L=1,4),M=1,NPAR(4))
!
          IF(NPAR(3).GT.0) THEN
            WRITE(ITEM1) ((NDINT(L,M),L=1,3),M=1,NPAR(3))
            DEALLOCATE(NDINT)
          ENDIF
!
          IF(NPAR(8).GT.0) THEN
             WRITE(ITEM1) ((CPLINT(L,M),L=1,4),M=1,NPAR(8))
             DEALLOCATE(CPLINT)
          ENDIF
!
          IF(NPAR(10).GT.0) THEN
             WRITE(ITEM1) ((CPLINT(L,M),L=1,6),M=1,NPAR(10))
             DEALLOCATE(CPLINT)
          ENDIF
!
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
       CASE(13)
          NPAR=0
          NPAR(1)=13
          CALL FEMC_QUADSD(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,3),M=1,NPAR(4)),(ELEPAR(1,M),M=1,NPAR(4))
          WRITE(ITEM1) ((NDINT(L,M),L=1,3),M=1,NPAR(3))
!
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(NDINT)
       CASE(14)
          NPAR=0
          NPAR(1)=14
          CALL OFE_TETRA(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,7),M=1,NPAR(3)),((ELEPAR(L,M),L=1,5),M=1,NPAR(3))
          WRITE(ITEM1) ((NDINT(L,M),L=1,4),M=1,NPAR(4))
!
          DEALLOCATE(NDMFS)
          DEALLOCATE(ELEPAR)
          DEALLOCATE(NDINT)          
       CASE(15)
          NPAR=0
          NPAR(1)=15
          CALL OFE_RECT(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((NDMFS(L,M),L=1,5),M=1,NPAR(3))
          WRITE(ITEM1) ((NDINT(L,M),L=1,4),M=1,NPAR(4))
          WRITE(ITEM1) ((CPLINT(L,M),L=1,4),M=1,NPAR(6))

          DEALLOCATE(NDMFS)
          DEALLOCATE(NDINT)
          DEALLOCATE(CPLINT)
       CASE(16)
          NPAR=0
          NPAR(1)=16
          CALL HEX8(IIN,IOUT)
          WRITE(ITEM1) (NPAR(L),L=1,12)
          WRITE(ITEM1) ((XYZ(L,M),L=1,24),M=1,NPAR(3)),((LM(L,M),L=1,24),M=1,NPAR(3))          ! Deallocate variables
          DEALLOCATE (XYZ)
          DEALLOCATE (LM)
!------------------------------------------
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
    END DO  
!
    EFFORT=0
    DO I=1,NEQ
       EFFORT=MHT(I)*MHT(I)+EFFORT
    ENDDO        
!
    WRITE(*,*) "SolE EFFORT", EFFORT/2.0D0
!
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000 FORMAT(/,A20)
 2000 FORMAT(//,' E L E M E N T   I N F O R M A T I O N',///)
!}}}  
    END SUBROUTINE ELEMENT_IN
!
    SUBROUTINE HEX8(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE HEX8
! *      READ AND STORE TRUSS DATA
! *  QUADS
! *      NPAR(1) = 16
! *      NPAR(2) = MTYPE
! *      NPAR(3) = NUMELE
! *      NPAR(4) = NINTG
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS
! *        NINTG    -  NUMBER OF NUMERICAL INTEGRATION POINTS (x-direction)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle
    REAL(8) :: THIC
    INTEGER :: N,I1,I2,I3,I4,I5,I6,I7,I8
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(4) :: KWQUADS
!- - - - - - - - - - - - - - -  
!-     VARIABLE DEFINITION   -
!- - - - - - - - - - - - - - -
    KWQUADS(1) = 'MTYPE:'
    KWQUADS(2) = 'NUMELE:'
    KWQUADS(3) = 'NINTG:'
    KWQUADS(4) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
!  
    I = KeyWord(KWQUADS,4)
    DO WHILE(I<=5)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
            NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN HEX8"
           STOP
       END SELECT
       I = KeyWord(KWQUADS,4)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF HEX8 ELEMENT .EQ.0"
        STOP
    END IF
    NumEle=NPAR(3)
!
!-------- Allocate memory
!
        ALLOCATE (XYZ(24,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (LM(24,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! ----------------------------------------------------------------
        XYZ=0
        LM=0
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, I1, I2, I3, I4, I5, I6, I7, I8
        DO WHILE(N.NE.0)
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           XYZ(1,N)=X(I1)
           XYZ(9,N)=Y(I1)
           XYZ(17,N)=Z(I1)
!
           XYZ(2,N)=X(I2)
           XYZ(10,N)=Y(I2)
           XYZ(18,N)=Z(I2)
!
           XYZ(3,N)=X(I3)
           XYZ(11,N)=Y(I3)
           XYZ(19,N)=Z(I3)
!
           XYZ(4,N)=X(I4)
           XYZ(12,N)=Y(I4)
           XYZ(20,N)=Z(I4)
!
           XYZ(5,N)=X(I5)
           XYZ(13,N)=Y(I5)
           XYZ(21,N)=Z(I5)
!
           XYZ(6,N)=X(I6)
           XYZ(14,N)=Y(I6)
           XYZ(22,N)=Z(I6)
!
           XYZ(7,N)=X(I7)
           XYZ(15,N)=Y(I7)
           XYZ(23,N)=Z(I7)
!
           XYZ(8,N)=X(I8)
           XYZ(16,N)=Y(I8)
           XYZ(24,N)=Z(I8)
!
           DO I=1,3
             LM(I  ,N)=ID(I,I1)
             LM(I+3,N)=ID(I,I2)
             LM(I+6,N)=ID(I,I3)
             LM(I+9,N)=ID(I,I4)
             LM(I+12,N)=ID(I,I5)
             LM(I+15,N)=ID(I,I6)
             LM(I+18,N)=ID(I,I7)
             LM(I+21,N)=ID(I,I8)
           END DO
           CALL COLHT(24,LM(1,N))
           WRITE(IOUT,2050) N,I1,I2,I3,I4,I5,I6,I7,I8
!
           READ(IIN,*) N, I1, I2, I3, I4, I5, I6, I7, I8
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
            WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
            STOP
        END IF
    RETURN
!--------------------------------------------------
 1000 Format(/,&                 
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF HEX8 ELEMENTS       =',I5,/,&
              '    GAUSS NUMERICAL INTEGRATION ORDER IN X-DIRECTION =',I5,/)
 2000 FORMAT(/,&
              ' ELEMENT   N1   N2  N3  N4   N5   N6   N7   N8 ',/,&
              ' NUMBER-N                                      ',/)
 2050 FORMAT(I5,6X,8I5)
!}}}
    END SUBROUTINE HEX8
!
    SUBROUTINE OFE_RECT(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OFE_RECT
! *      READ AND STORE QUADS DATA
! *      QUADS
! *         NPAR(1) = 15
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE :  TOTAL NUMBER OF ELEMENTS
! *         NPAR(4) = NUMOLR :  TOTAL NUMBER OF OVERLAP REGIONS
! *         NPAR(5) = NGDOF  :  DEGREE OF PRECISION OF NUMERICAL INTEGRATION
! *         NPAR(6) = NUMCPL :  TOTAL NUMBER OF COUPLING ELEMENTS
! *         NPAR(7) = NING   :  NUMBER OF INTEGRATION POINTS IN X-DIRECTION
! *
! *         NPAR(8) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumOFE,NumOLR,NumCPL
    INTEGER :: N,II,JJ,KK,LL,MM,KN,IND
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,KOLD,LOLD,NII,NJJ,NKK
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
    CHARACTER(20),DIMENSION(8) :: KWQUADS_OFE
!
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_OFE(1) = 'MTYPE:'
    KWQUADS_OFE(2) = 'NUMELE:'
    kWQUADS_OFE(3) = 'NUMOLR:'
    KWQUADS_OFE(4) = 'NGDOF:'
    KWQUADS_OFE(5) = 'NUMCPL:'
    KWQUADS_OFE(6) = 'NING:'
    KWQUADS_OFE(7) = 'ITYPE:'
    KWQUADS_OFE(8) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=4
    NPAR(6)=0
    NPAR(7)=3
    NPAR(8)=1
!   
    I = KeyWord(KWQUADS_OFE,8)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           J = GetInt()
           NPAR(7)=J
         CASE(7)
           J = GetInt()
           NPAR(8)=J
         CASE(8)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OFE1_QUADSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_OFE,8)
    END DO
!
    WRITE(IOUT,1000) NPAR(8),NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6),NPAR(7)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAPPING ELEMENTS .EQ.0"
        STOP
    END IF
    NumOFE=NPAR(3)
    NumOLR=NPAR(4)
    NumCPL=NPAR(6)
!-------- Allocate memory
    ALLOCATE (NDMFS(5,NumOFE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    NDMFS=0
!
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, IND, II,JJ,KK,LL
    DO WHILE(N.NE.0)
       IF(N.GT.NumOFE) THEN
          WRITE(*,*) "ERROR N .GT. NUMOFE"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=IND
       NDMFS(2,N)=II
       NDMFS(3,N)=JJ        
       NDMFS(4,N)=KK     
       NDMFS(5,N)=LL
!
       WRITE(IOUT,2250) N,IND,II,JJ,KK,LL
!
       READ(IIN,*) N,IND,II,JJ,KK,LL
    END DO
!---------------------------------------------------------------
!
    IF(NCOUNT.NE.NumOFE) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumOFE!"
       STOP
    ENDIF
!
    ALLOCATE (NDINT(3,NumOLR),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, KK,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumOLR) THEN
          WRITE(*,*) "ERROR N .GT. NUMOLR"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       NDINT(3,N)=KK

       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              NKK=KOLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD
                 NKK=NKK+KNOLD

                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 NDINT(3,K)=NKK

                 WRITE(IOUT,2050) K,NII,NJJ,NKK                  
                 CALL COLHT_OFE_RECT(NII,NJJ,NKK)
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF

        WRITE(IOUT,2050) N,II,JJ,KK
        CALL COLHT_OFE_RECT(II,JJ,KK)
   
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        KOLD=KK

        READ(IIN,*) N,II,JJ,KK,KN
    END DO

    IF(NCOUNT.NE.NumOLR) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumOLR"
       STOP
    END IF
!
    ALLOCATE (CPLINT(4,NumCPL),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2010)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, KK,LL
    DO WHILE(N.NE.0)
       IF(N.GT.NumCPL) THEN
          WRITE(*,*) "ERROR N .GT. NUMCPL"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       CPLINT(1,N)=II
       CPLINT(2,N)=JJ          
       CPLINT(3,N)=KK
       CPLINT(4,N)=LL

       WRITE(IOUT,2060) N,II,JJ,KK,LL
       CALL COLHT_OFE_CPL(II,JJ,KK,LL)

       READ(IIN,*) N,II,JJ,KK,LL
    END DO

    IF(NCOUNT.NE.NumCPL) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumCPL"
       STOP
    END IF
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE:QUADS            =',I5,/,&
              '               1: PLANE STRAIN          ',/,&
              '               2: PLANE STRESS          ',/,&       
              '    THE MATERIAL SET NUMBER            =',I5,/,&
              '    THE NUMBER OF OVERLAPPING ELEMENTS =',I5,/,&
              '    THE NUMBER OF OVERLAP REGIONS      =',I5,/,&
              '    THE DEGREE PRECISION OF THE GAUSS INTEGRATION     =',I5,/,&
              '    THE NUMBER OF COUPLING ELEMENTS     =',I5,/,&
              '    NUMBER OF INTEGRATION POINTS IN X-DIRECTION     =',I5)
 2200 FORMAT(/,&
              ' ELEMENT  IND  ND1   ND2   ND3   ND4 ',/,&
              ' NUMBER-N                                ',/)
 2250 FORMAT(I5,2X,I5,4X,4I5)

 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE           ',/,&
              ' NUMBER-N      I        J       K             ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5)
 2010 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE   ELE        ',/,&
              ' NUMBER-N      I        J       K     L        ',/)
 2060 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5)
!}}}
    END SUBROUTINE OFE_RECT
!
    SUBROUTINE COLHT_OFE_RECT(II,JJ,KK)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, JJ, KK
!   IN-FUNCTION
    INTEGER,DIMENSION(12):: ND
    INTEGER :: IDD,LS,NN,NI  
    INTEGER :: I,J,K,KKK
    INTEGER :: ME
!   
    ND(1:4)=NDMFS(2:5,II)
    ND(5:8)=NDMFS(2:5,JJ)
    ND(9:12)=NDMFS(2:5,KK)
!
    NN=MINVAL(ND)
    LS=ID(1,NN)
!------------------------------------------
    DO I=1,12
       NI=ND(I)
       DO J=1,2
          IDD=ID(J,NI)
          IF(IDD.NE.0) THEN
              DO K=1,NDOF(J,NI)
                 KKK=IDD+K-1
                 ME=KKK-LS
                 IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
              ENDDO
          ENDIF
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE_RECT
!
    SUBROUTINE COLHT_OFE_CPL(II,JJ,KK,LL)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, JJ, KK,LL
!   IN-FUNCTION
    INTEGER,DIMENSION(4):: ND
    INTEGER :: IDD,LS,NN,NI
    INTEGER :: I,J,K
    INTEGER :: KKK
    INTEGER :: ME
!   
    ND(1)=NDMFS(2,II)
    ND(2)=NDMFS(2,JJ)
    ND(3)=NDMFS(2,KK)
    ND(4)=NDMFS(2,LL)
!
    NN=MINVAL(ND)
    LS=ID(1,NN)
!------------------------------------------
    DO I=1,12
       NI=ND(I)
       DO J=1,2
          IDD=ID(J,NI)
          IF(IDD.NE.0) THEN
              DO K=1,NDOF(J,NI)
                 KKK=IDD+K-1
                 ME=KKK-LS
                 IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
              ENDDO
          ENDIF
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE_CPL
!
    SUBROUTINE OFE_TETRA(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OFE_3D
! *      READ AND STORE QUADS DATA
! *      QUADS
! *         NPAR(1) = 14
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE :  TOTAL NUMBER OF ELEMENTS
! *         NPAR(4) = NUMOLR :  TOTAL NUMBER OF OVERLAP REGIONS
! *         NPAR(5) = NGDOF  :  DEGREE OF PRECISION OF NUMERICAL INTEGRATION
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    CHARACTER(20),DIMENSION(5) :: KWQUADS_OFE
    INTEGER :: NUMELE,NUMOLR
    INTEGER :: FIX_X,FIX_Y,FIX_Z,FIX_X_OLD,FIX_Y_OLD,FIX_Z_OLD
    INTEGER :: I1,I1OLD,I2,I2OLD,I3,I3OLD,I4,I4OLD
    INTEGER :: NUM,NUMN,KN,KNOLD,N,NOLD
    INTEGER :: NI1,NI2,NI3,NI4
    INTEGER :: NCOUNT
    REAL(8),DIMENSION(3) :: XYZ1,XYZ2,XYZ3,XYZ4
    REAL(8),DIMENSION(3) :: XYZ12,XYZ13,XYZ14,XYZ15,XYZ23,XYZ24,XYZ34
    REAL(8) :: L12,L13,L14,L23,L24,L34
    REAL(8) :: R1,R2,R3,R4
    INTEGER :: AllocateStatus
    INTEGER :: NOD1,NOD2,NOD3,NOD4
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_OFE(1) = 'MTYPE:'
    KWQUADS_OFE(2) = 'NUMELE:'
    kWQUADS_OFE(3) = 'NUMOLR:'
    KWQUADS_OFE(4) = 'NGDOF:'
    KWQUADS_OFE(5) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=6
!   
    I = KeyWord(KWQUADS_OFE,5)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OFE_TETRA"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_OFE,5)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAPPING ELEMENTS .EQ.0"
        STOP
    END IF
    NumELE=NPAR(3)
    NumOLR=NPAR(4)
!-------- Allocate memory
    ALLOCATE (NDMFS(7,NumELE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(5,NumELE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    NDMFS=0
    ELEPAR(1,:)=0.0    ! Store the largest tetrahedron length
    ELEPAR(2,:)=10000  ! Store the smallest tetrahedron length
    ELEPAR(3:5,:)=0.0D0
!
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, I1, FIX_X,FIX_Y,FIX_Z,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumELE) THEN
          WRITE(*,*) "ERROR N .GT. NUMELE"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=I1
       NDMFS(2,N)=FIX_X        
       NDMFS(3,N)=FIX_Y     
       NDMFS(4,N)=FIX_Z
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NI1=I1OLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NI1=NI1+KNOLD
!
                 NDMFS(1,K)=NI1
                 NDMFS(2,K)=FIX_X_OLD
                 NDMFS(3,K)=FIX_Y_OLD
                 NDMFS(4,K)=FIX_Z_OLD
                 WRITE(IOUT,2250) K,NI1,FIX_X_OLD,FIX_Y_OLD,FIX_Z_OLD               
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,I1,FIX_X,FIX_Y,FIX_Z
!     
        NOLD=N
        KNOLD=KN
        I1OLD=I1
        FIX_X_OLD=FIX_X
        FIX_Y_OLD=FIX_Y
        FIX_Z_OLD=FIX_Z
!
        READ(IIN,*) N, I1, FIX_X,FIX_Y,FIX_Z,KN
    END DO
!---------------------------------------------------------------
    IF(NCOUNT.NE.NumELE) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!
    ALLOCATE (NDINT(4,NumOLR),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!-------------------------------------------------------------
    CALL NewLine()
!-------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, I1, I2, I3, I4, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumOLR) THEN
          WRITE(*,*) "ERROR N .GT. NUMOLR"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=I1
       NDINT(2,N)=I2          
       NDINT(3,N)=I3
       NDINT(4,N)=I4

       NOD1=NDMFS(1,I1)
       NOD2=NDMFS(1,I2)
       NOD3=NDMFS(1,I3)
       NOD4=NDMFS(1,I4)

       XYZ1(1)=X(NOD1)
       XYZ1(2)=Y(NOD1)
       XYZ1(3)=Z(NOD1)

       XYZ2(1)=X(NOD2)
       XYZ2(2)=Y(NOD2)
       XYZ2(3)=Z(NOD2)

       XYZ3(1)=X(NOD3)
       XYZ3(2)=Y(NOD3)
       XYZ3(3)=Z(NOD3)

       XYZ4(1)=X(NOD4)
       XYZ4(2)=Y(NOD4)
       XYZ4(3)=Z(NOD4)

       XYZ12=XYZ1-XYZ2
       XYZ13=XYZ1-XYZ3
       XYZ14=XYZ1-XYZ4
       XYZ23=XYZ2-XYZ3
       XYZ24=XYZ2-XYZ4
       XYZ34=XYZ3-XYZ4

       L12=SQRT(DOT_PRODUCT(XYZ12,XYZ12))
       L13=SQRT(DOT_PRODUCT(XYZ13,XYZ13))
       L14=SQRT(DOT_PRODUCT(XYZ14,XYZ14))
       L23=SQRT(DOT_PRODUCT(XYZ23,XYZ23))
       L24=SQRT(DOT_PRODUCT(XYZ24,XYZ24))
       L34=SQRT(DOT_PRODUCT(XYZ34,XYZ34))

       R1=MAX(L12,L13,L14)
       R2=MAX(L12,L23,L24)
       R3=MAX(L13,L23,L34)
       R4=MAX(L14,L24,L34)
       IF(ELEPAR(1,I1).LT.R1) ELEPAR(1,I1)=R1
       IF(ELEPAR(1,I2).LT.R2) ELEPAR(1,I2)=R2
       IF(ELEPAR(1,I3).LT.R3) ELEPAR(1,I3)=R3
       IF(ELEPAR(1,I4).LT.R4) ELEPAR(1,I4)=R4

       R1=MIN(L12,L13,L14)
       R2=MIN(L12,L23,L24)
       R3=MIN(L13,L23,L34)
       R4=MIN(L14,L24,L34)
       IF(ELEPAR(2,I1).GT.R1) ELEPAR(2,I1)=R1
       IF(ELEPAR(2,I2).GT.R2) ELEPAR(2,I2)=R2
       IF(ELEPAR(2,I3).GT.R3) ELEPAR(2,I3)=R3
       IF(ELEPAR(2,I4).GT.R4) ELEPAR(2,I4)=R4

       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NI1=I1OLD
              NI2=I2OLD
              NI3=I3OLD
              NI4=I4OLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NI1=NI1+KNOLD
                 NI2=NI2+KNOLD
                 NI3=NI3+KNOLD
                 NI4=NI4+KNOLD
 
                 NDINT(1,K)=NI1
                 NDINT(2,K)=NI2
                 NDINT(3,K)=NI3
                 NDINT(4,K)=NI4

                 NOD1=NDMFS(1,NI1)
                 NOD2=NDMFS(1,NI2)
                 NOD3=NDMFS(1,NI3)
                 NOD4=NDMFS(1,NI4)

                 XYZ1(1)=X(NOD1)
                 XYZ1(2)=Y(NOD1)
                 XYZ1(3)=Z(NOD1)
       !
                 XYZ2(1)=X(NOD2)
                 XYZ2(2)=Y(NOD2)
                 XYZ2(3)=Z(NOD2)
       !
                 XYZ3(1)=X(NOD3)
                 XYZ3(2)=Y(NOD3)
                 XYZ3(3)=Z(NOD3)
       !
                 XYZ4(1)=X(NOD4)
                 XYZ4(2)=Y(NOD4)
                 XYZ4(3)=Z(NOD4)
       !
                 XYZ12=XYZ1-XYZ2
                 XYZ13=XYZ1-XYZ3
                 XYZ14=XYZ1-XYZ4
                 XYZ23=XYZ2-XYZ3
                 XYZ24=XYZ2-XYZ4
                 XYZ34=XYZ3-XYZ4
       !
                 L12=SQRT(DOT_PRODUCT(XYZ12,XYZ12))
                 L13=SQRT(DOT_PRODUCT(XYZ13,XYZ13))
                 L14=SQRT(DOT_PRODUCT(XYZ14,XYZ14))
                 L23=SQRT(DOT_PRODUCT(XYZ23,XYZ23))
                 L24=SQRT(DOT_PRODUCT(XYZ24,XYZ24))
                 L34=SQRT(DOT_PRODUCT(XYZ34,XYZ34))
       !
                 R1=MAX(L12,L13,L14)
                 R2=MAX(L12,L23,L24)
                 R3=MAX(L13,L23,L34)
                 R4=MAX(L14,L24,L34)
                 IF(ELEPAR(1,NI1).LT.R1) ELEPAR(1,NI1)=R1
                 IF(ELEPAR(1,NI2).LT.R2) ELEPAR(1,NI2)=R2
                 IF(ELEPAR(1,NI3).LT.R3) ELEPAR(1,NI3)=R3
                 IF(ELEPAR(1,NI4).LT.R4) ELEPAR(1,NI4)=R4
       !
                 R1=MIN(L12,L13,L14)
                 R2=MIN(L12,L23,L24)
                 R3=MIN(L13,L23,L34)
                 R4=MIN(L14,L24,L34)
                 IF(ELEPAR(2,NI1).GT.R1) ELEPAR(2,NI1)=R1
                 IF(ELEPAR(2,NI2).GT.R2) ELEPAR(2,NI2)=R2
                 IF(ELEPAR(2,NI3).GT.R3) ELEPAR(2,NI3)=R3
                 IF(ELEPAR(2,NI4).GT.R4) ELEPAR(2,NI4)=R4

                 WRITE(IOUT,2050) K,NI1,NI2,NI3,NI4                 
                 CALL COLHT_OFE_TETRA(NI1,NI2,NI3,NI4)
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF

        WRITE(IOUT,2050) N,I1,I2,I3,I4
        CALL COLHT_OFE_TETRA(I1,I2,I3,I4)
   
        NOLD=N
        KNOLD=KN
        I1OLD=I1
        I2OLD=I2
        I3OLD=I3
        I4OLD=I4

        READ(IIN,*) N,I1,I2,I3,I4,KN
    END DO

    IF(NCOUNT.NE.NumOLR) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumOLR"
       STOP
    END IF
!{{{
!    IF(NumCPL.GT.0) THEN
!       ALLOCATE (CPLINT(8,NumCPL),STAT = ALLocateStatus)
!       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!!---------------------------------------------------------------
!       CALL NewLine()
!!---------------------------------------------------------------
!       WRITE(IOUT,2010)
!       NCOUNT=0
!       READ(IIN,*) N, I1, I2, I3, I4, I5, I6, I7, I8, KN
!       DO WHILE(N.NE.0)
!          IF(N.GT.NumCPL) THEN
!             WRITE(*,*) "ERROR N .GT. NumCPL"
!             STOP
!          END IF
!          NCOUNT=NCOUNT+1
!!   
!          CPLINT(1,N)=I1
!          CPLINT(2,N)=I2          
!          CPLINT(3,N)=I3
!          CPLINT(4,N)=I4
!          CPLINT(5,N)=I5
!          CPLINT(6,N)=I6
!          CPLINT(7,N)=I7
!          CPLINT(8,N)=I8
!!
!          NOD1=NDMFS(1,I1)
!          NOD2=NDMFS(1,I2)
!          NOD3=NDMFS(1,I3)
!          NOD5=NDMFS(1,I5)          
!!
!          XYZ1(1)=X(NOD1)
!          XYZ1(2)=Y(NOD1)
!          XYZ1(3)=Z(NOD1)
!!
!          XYZ2(1)=X(NOD2)
!          XYZ2(2)=Y(NOD2)
!          XYZ2(3)=Z(NOD2)
!!
!          XYZ3(1)=X(NOD3)
!          XYZ3(2)=Y(NOD3)
!          XYZ3(3)=Z(NOD3)
!!
!          XYZ5(1)=X(NOD5)
!          XYZ5(2)=Y(NOD5)
!          XYZ5(3)=Z(NOD5)
!!
!          XYZ12=XYZ1-XYZ2
!          XYZ23=XYZ2-XYZ3
!          XYZ15=XYZ1-XYZ5
!!
!          L12=DOT_PRODUCT(XYZ12,XYZ12)
!          L23=DOT_PRODUCT(XYZ23,XYZ23)
!          L15=DOT_PRODUCT(XYZ15,XYZ15)
!!
!          R1=SQRT(L12+L23+L15)
!          IF(ELEPAR(1,I1).LT.R1) ELEPAR(1,I1)=R1
!          IF(ELEPAR(1,I2).LT.R1) ELEPAR(1,I2)=R1
!          IF(ELEPAR(1,I3).LT.R1) ELEPAR(1,I3)=R1
!          IF(ELEPAR(1,I4).LT.R1) ELEPAR(1,I4)=R1
!          IF(ELEPAR(1,I5).LT.R1) ELEPAR(1,I5)=R1
!          IF(ELEPAR(1,I6).LT.R1) ELEPAR(1,I6)=R1
!          IF(ELEPAR(1,I7).LT.R1) ELEPAR(1,I7)=R1
!          IF(ELEPAR(1,I8).LT.R1) ELEPAR(1,I8)=R1
!
!          R1=MIN(SQRT(L12),SQRT(L23),SQRT(L15))
!          IF(ELEPAR(2,I1).GT.R1) ELEPAR(2,I1)=R1
!          IF(ELEPAR(2,I2).GT.R1) ELEPAR(2,I2)=R1
!          IF(ELEPAR(2,I3).GT.R1) ELEPAR(2,I3)=R1
!          IF(ELEPAR(2,I4).GT.R1) ELEPAR(2,I4)=R1
!          IF(ELEPAR(2,I5).GT.R1) ELEPAR(2,I5)=R1
!          IF(ELEPAR(2,I6).GT.R1) ELEPAR(2,I6)=R1
!          IF(ELEPAR(2,I7).GT.R1) ELEPAR(2,I7)=R1
!          IF(ELEPAR(2,I8).GT.R1) ELEPAR(2,I8)=R1
!!
!          IF (KNOLD.NE.0) THEN
!             NUM=(N-NOLD)/KNOLD
!             NUMN=NUM - 1
!             IF (NUMN.GE.1) THEN
!                 K=NOLD
!                 NI1=I1OLD
!                 NI2=I2OLD
!                 NI3=I3OLD
!                 NI4=I4OLD
!                 NI5=I5OLD
!                 NI6=I6OLD
!                 NI7=I7OLD
!                 NI8=I8OLD
!                 DO J=1,NUMN
!                    K=K + KNOLD
!                    
!                    NI1=NI1+KNOLD
!                    NI2=NI2+KNOLD
!                    NI3=NI3+KNOLD
!                    NI4=NI4+KNOLD
!                    NI5=NI5+KNOLD
!                    NI6=NI6+KNOLD
!                    NI7=NI7+KNOLD
!                    NI8=NI8+KNOLD
!  
!                    CPLINT(1,K)=NI1
!                    CPLINT(2,K)=NI2          
!                    CPLINT(3,K)=NI3
!                    CPLINT(4,K)=NI4
!                    CPLINT(5,K)=NI5
!                    CPLINT(6,K)=NI6          
!                    CPLINT(7,K)=NI7
!                    CPLINT(8,K)=NI8
!
!                    NOD1=NDMFS(1,NI1)
!                    NOD2=NDMFS(1,NI2)
!                    NOD3=NDMFS(1,NI3)
!                    NOD5=NDMFS(1,NI5)          
!          !
!                    XYZ1(1)=X(NOD1)
!                    XYZ1(2)=Y(NOD1)
!                    XYZ1(3)=Z(NOD1)
!          !
!                    XYZ2(1)=X(NOD2)
!                    XYZ2(2)=Y(NOD2)
!                    XYZ2(3)=Z(NOD2)
!          !
!                    XYZ3(1)=X(NOD3)
!                    XYZ3(2)=Y(NOD3)
!                    XYZ3(3)=Z(NOD3)
!          !
!                    XYZ5(1)=X(NOD5)
!                    XYZ5(2)=Y(NOD5)
!                    XYZ5(3)=Z(NOD5)
!          !
!                    XYZ12=XYZ1-XYZ2
!                    XYZ23=XYZ2-XYZ3
!                    XYZ15=XYZ1-XYZ5
!          !
!                    L12=DOT_PRODUCT(XYZ12,XYZ12)
!                    L23=DOT_PRODUCT(XYZ23,XYZ23)
!                    L15=DOT_PRODUCT(XYZ15,XYZ15)
!          !
!                    R1=SQRT(L12+L23+L15)
!                    IF(ELEPAR(1,NI1).LT.R1) ELEPAR(1,NI1)=R1
!                    IF(ELEPAR(1,NI2).LT.R1) ELEPAR(1,NI2)=R1
!                    IF(ELEPAR(1,NI3).LT.R1) ELEPAR(1,NI3)=R1
!                    IF(ELEPAR(1,NI4).LT.R1) ELEPAR(1,NI4)=R1
!                    IF(ELEPAR(1,NI5).LT.R1) ELEPAR(1,NI5)=R1
!                    IF(ELEPAR(1,NI6).LT.R1) ELEPAR(1,NI6)=R1
!                    IF(ELEPAR(1,NI7).LT.R1) ELEPAR(1,NI7)=R1
!                    IF(ELEPAR(1,NI8).LT.R1) ELEPAR(1,NI8)=R1
!          
!                    R1=MIN(SQRT(L12),SQRT(L23),SQRT(L15))
!                    IF(ELEPAR(2,NI1).GT.R1) ELEPAR(2,NI1)=R1
!                    IF(ELEPAR(2,NI2).GT.R1) ELEPAR(2,NI2)=R1
!                    IF(ELEPAR(2,NI3).GT.R1) ELEPAR(2,NI3)=R1
!                    IF(ELEPAR(2,NI4).GT.R1) ELEPAR(2,NI4)=R1
!                    IF(ELEPAR(2,NI5).GT.R1) ELEPAR(2,NI5)=R1
!                    IF(ELEPAR(2,NI6).GT.R1) ELEPAR(2,NI6)=R1
!                    IF(ELEPAR(2,NI7).GT.R1) ELEPAR(2,NI7)=R1
!                    IF(ELEPAR(2,NI8).GT.R1) ELEPAR(2,NI8)=R1
!!
!                    WRITE(IOUT,2060) K,NI1,NI2,NI3,NI4,NI5,NI6,NI7,NI8                 
!                    CALL COLHT_OFE_3D(NI1,NI2,NI3,NI4,NI5,NI6,NI7,NI8)
!                    NCOUNT=NCOUNT+1
!                 END DO
!             END IF
!           ENDIF
!!
!           WRITE(IOUT,2060) N,I1,I2,I3,I4,I5,I6,I7,I8
!           CALL COLHT_OFE_3D(I1,I2,I3,I4,I5,I6,I7,I8)
!!     
!           NOLD=N
!           KNOLD=KN
!           I1OLD=I1
!           I2OLD=I2
!           I3OLD=I3
!           I4OLD=I4
!           I5OLD=I5
!           I6OLD=I6
!           I7OLD=I7
!           I8OLD=I8
!!
!           READ(IIN,*) N,I1,I2,I3,I4,I5,I6,I7,I8,KN
!       END DO
!!
!       IF(NCOUNT.NE.NumCPL) THEN
!          WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumCPL"
!          STOP
!       ENDIF
!    ENDIF
!}}}
    ELEPAR(2,:)=ELEPAR(1,:)/ELEPAR(2,:)
    ELEPAR(1,:)=ELEPAR(2,:)**3*ELEPAR(1,:)
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&   
              '    THE MATERIAL SET NUMBER            =',I5,/,&
              '    THE NUMBER OF OVERLAPPING ELEMENTS =',I5,/,&
              '    THE NUMBER OF OVERLAP REGIONS      =',I5,/,&
              '    DEGREE OF PRECISION OF NUMERICAL INTEGRATION  =',I5)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND(3,N)  ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,3I5)

 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE   ELE         ',/,&
              ' NUMBER-N      1        2       3     4           ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5)
! 2010 FORMAT(/,&
!              ' INTEGRATION   ELE      ELE      ELE     ELE    ELE   ELE   ELE   ELE     ',/,&
!              ' NUMBER-N      I1       I2       I3      I4     I5    I6    I7    I8      ',/)
! 2060 FORMAT(I5,6X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
!}}}
    END SUBROUTINE OFE_TETRA
!
    SUBROUTINE COLHT_OFE_TETRA(I1,I2,I3,I4)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE COLHT_OFE_3D
! *      PURPOSE:  BUILD MHT MATRIX 
! *      INPUT
! *         ELEMENT NUMBER
! *         I1,I2,I3,I4
! *      OUTPUT:
! *         UPDATE MHT MATRIX
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: I1,I2,I3,I4
!   IN-FUNCTION
    INTEGER :: N1,N2,N3,N4,NN
    INTEGER :: IDD,LS   
    INTEGER :: I,J
    INTEGER :: KKK,ME
!
    N1=NDMFS(1,I1)
    N2=NDMFS(1,I2)
    N3=NDMFS(1,I3)
    N4=NDMFS(1,I4)
!
    LS=1E9
    NN=MIN(N1,N2,N3,N4)
    DO J=1,3
       IDD=ID(J,N1)
       IF (IDD.NE.0) THEN
          IF (LS>IDD) LS=IDD
       ENDIF
       IDD=ID(J,N2)
       IF (IDD.NE.0) THEN
          IF (LS>IDD) LS=IDD
       ENDIF
       IDD=ID(J,N3)
       IF (IDD.NE.0) THEN
          IF (LS>IDD) LS=IDD
       ENDIF
       IDD=ID(J,N4)
       IF (IDD.NE.0) THEN
          IF (LS>IDD) LS=IDD
       ENDIF
    ENDDO
!
    DO J=1,3
       IDD=ID(J,N1)
       IF(IDD.GT.0) THEN
           DO I=1,NDOF(J,N1)
              KKK=IDD+I-1
              ME=KKK-LS
              IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
           ENDDO
       ENDIF
!
       IDD=ID(J,N2)
       IF(IDD.GT.0) THEN
          DO I=1,NDOF(J,N2)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,N3)
       IF(IDD.GT.0) THEN
          DO I=1,NDOF(J,N3)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,N4)
       IF(IDD.GT.0) THEN
          DO I=1,NDOF(J,N4)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
   !
    ENDDO
!------------------------------------------
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE_TETRA
!
    SUBROUTINE FEMC_QUADSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OFE1_QUADSD
! *      READ AND STORE QUADS DATA
! *      QUADS
! *         NPAR(1) = 12
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMOFE_Triangle
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NGDOF
! *         NPAR(6) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumOFE,NumEle
    INTEGER :: N,II,JJ,KK,LL,KN
    INTEGER :: IOLD,JOLD,KOLD
    INTEGER :: NII,NJJ,NKK
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IND1,IND2,IND1OLD,IND2OLD
    REAL(8) :: RAD1,R1OLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
    CHARACTER(20),DIMENSION(6) :: KWQUADS_OFE
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_OFE(1) = 'MTYPE:'
    KWQUADS_OFE(2) = 'NUMOFE:'
    kWQUADS_OFE(3) = 'NUMELE:'
    KWQUADS_OFE(4) = 'NGDOF:'
    KWQUADS_OFE(5) = 'ITYPE:'
    KWQUADS_OFE(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=4
    NPAR(6)=1
!   
    I = KeyWord(KWQUADS_OFE,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OFE1_QUADSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_OFE,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(6),NPAR(2),NPAR(4),NPAR(3),NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAPPING ELEMENTS .EQ.0"
        STOP
    END IF
    NumOFE=NPAR(3)
    NumEle=NPAR(4)
!-------- Allocate memory
    ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND1,IND2,RAD1,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NUMELE"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND1        
       NDMFS(3,N)=IND2     
       ELEPAR(1,N)=RAD1
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=IND1OLD
                 NDMFS(3,K)=IND2OLD
                 ELEPAR(1,K)=R1OLD
                 WRITE(IOUT,2250) K,NII,IND1OLD,IND2OLD,R1OLD               
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND1,IND2,RAD1
        NOLD=N
        KNOLD=KN
        IOLD=II
        IND1OLD=IND1
        IND2OLD=IND2
        R1OLD=RAD1
!
        READ(IIN,*) N,II,IND1,IND2,RAD1,KN
    END DO
!---------------------------------------------------------------
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!
    ALLOCATE (NDINT(3,NumOFE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, KK,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumOFE) THEN
          WRITE(*,*) "ERROR N .GT. NUMOFE"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       NDINT(3,N)=KK
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              NKK=KOLD
              DO J=1,NUMN
                 K=K + KNOLD
!
                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD
                 NKK=NKK+KNOLD
!
                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 NDINT(3,K)=NKK
!
                 WRITE(IOUT,2050) K,NII,NJJ,NKK                  
                 CALL COLHT_FEMC_QUADS(NII,NJJ,NKK)
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF

        WRITE(IOUT,2050) N,II,JJ,KK
        CALL COLHT_FEMC_QUADS(II,JJ,KK)
   
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        KOLD=KK
!
        READ(IIN,*) N,II,JJ,KK,KN
    END DO
!
    IF(NCOUNT.NE.NumOFE) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumOFE"
       STOP
    END IF
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE:QUADS            =',I5,/,&
              '               1: PLANE STRAIN          ',/,&
              '               2: PLANE STRESS          ',/,&       
              '    THE MATERIAL SET NUMBER            =',I5,/,&
              '    THE NUMBER OF OVERLAPPING ELEMENTS =',I5,/,&
              '    THE NUMBER OF OVERLAP REGIONS      =',I5,/,&
              '    THE DEGREE PRECISION OF THE GAUSS INTEGRATION     =',I5)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND(2,N)     RADIUS(1,N) ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,2I5,2X,ES13.3)

 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE           ',/,&
              ' NUMBER-N      I        J       K             ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5)
!}}}
    END SUBROUTINE FEMC_QUADSD
!
    SUBROUTINE COLHT_FEMC_QUADS(II,JJ,KK)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, JJ, KK
!   IN-FUNCTION
    INTEGER :: NI,NJ,NK
    INTEGER :: IDD,LS   
    INTEGER :: I,J
    INTEGER :: KKK,ME
!       
    NI=NDMFS(1,II)
    NJ=NDMFS(1,JJ)
    NK=NDMFS(1,KK)
!
    LS=100000
!
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
    ENDDO
!------------------------------------------
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           DO I=1,NDOF(J,NI)
              KKK=IDD+I-1
              ME=KKK-LS
              IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
           ENDDO
       ENDIF
!
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NJ)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NK)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_FEMC_QUADS
!
    SUBROUTINE OFE_TRUSSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OEL_TRUSSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 10
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMOLP
! *         NPAR(4) = NITG
! *         NPAR(5) = NWEI
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NUMOLP
    INTEGER :: N,ELEI,ELEJ,INDI,INDJ,KN
    INTEGER :: NOLD,EIOLD,EJOLD,INDIOLD,INDJOLD,KNOLD
    REAL(8) :: RADI,RADJ
    REAL(8) :: RIOLD,RJOLD
    INTEGER :: NUMN,NUM
    INTEGER :: NEI,NEJ
    INTEGER :: AllocateStatus
    INTEGER :: NCOUNT
    CHARACTER(20),DIMENSION(6) :: OFE_TRUSS_KW
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    OFE_TRUSS_KW(1) = 'MTYPE:'
    OFE_TRUSS_KW(2) = 'NUMOLP:'
    OFE_TRUSS_KW(3) = 'NITG:'
    OFE_TRUSS_KW(4) = 'NWEI:'
    OFE_TRUSS_KW(5) = '\END'
    OFE_TRUSS_KW(6) = 'NDIV:'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=NITP_IN
    NPAR(5)=1
    NPAR(6)=1
!   
    I = KeyWord(OFE_TRUSS_KW,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE(6)
           J = GetInt()
           NPAR(6)=J
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OEL_TRUSSD"
           STOP
       END SELECT
       I = KeyWord(OFE_TRUSS_KW,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(6),NPAR(5)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAP REGION .EQ.0"
        STOP
    END IF
    NUMOLP=NPAR(3)
!-------- Allocate memory
    ALLOCATE (NDINT(4,NUMOLP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(2,NUMOLP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!-------- INITIALIZATION ------------------
    NDINT=0
    INTPAR=0.0D0
!
    WRITE(IOUT,2000)
!---------------------------------------------------------------
    NCOUNT=0
    READ(IIN,*) N, ELEI, ELEJ, INDI, INDJ, RADI, RADJ, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NUMOLP) THEN
          WRITE(*,*) "ERROR N .GT. NUMOLP"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDINT(1,N)=ELEI
       NDINT(2,N)=ELEJ
       NDINT(3,N)=INDI
       NDINT(4,N)=INDJ
!
       INTPAR(1,N)=RADI
       INTPAR(2,N)=RADJ
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NEI=EIOLD
              NEJ=EJOLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NEI=NEI+KNOLD
                 NEJ=NEJ+KNOLD
!
                 NDINT(1,K)=NEI
                 NDINT(2,K)=NEJ
                 NDINT(3,K)=INDIOLD
                 NDINT(4,K)=INDJOLD
                 
                 INTPAR(1,K)=RIOLD
                 INTPAR(2,K)=RJOLD 
!
                 WRITE(IOUT,2010) K,NEI,NEJ,INDIOLD,INDJOLD,RIOLD,RJOLD              
                 NCOUNT=NCOUNT+1
                 CALL COLHT_OFE_TRUSS(NEI,NEJ)
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2010) N,ELEI,ELEJ,INDI,INDJ,RADI,RADJ
        CALL COLHT_OFE_TRUSS(ELEI,ELEJ)
!     
        NOLD=N
        KNOLD=KN
        EIOLD=ELEI
        EJOLD=ELEJ
        INDIOLD=INDI
        INDJOLD=INDJ
        RIOLD=RADI
        RJOLD=RADJ
!
        READ(IIN,*) N, ELEI, ELEJ, INDI, INDJ, RADI, RADJ, KN
    END DO
!
    IF(NCOUNT.NE.NUMOLP) THEN
       WRITE(*,*) "ERROR, OVERLAP REGION NUMBER .LT. NUMOLP!"
       STOP
    ENDIF
!--------------------------------------------------
    RETURN
!--------------------------------------------------
 1000 Format(/,&     
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF INTEGRATION REGION     =',I5,/,&
              '    THE NUMBER OF GAUSS INTEGRATION POINTS   =',I5,/,&
              '    THE NUMBER OF SUBINTERVALS               =',I5,/,&
              '    THE WEIGHT FUNCTION NUMBER       =',I5,/,&
              '          1: 1-6S^2+8S^3-3S^4         ',/)
 2000 FORMAT(/,&
              'OVERLAP    ELE       ELE     IND     IND     RAD     RAD ',/,&
              'NUMBER       I         J       I       J       I       J ',/)
 2010 FORMAT(I5,1X,I5,1X,I5,1X,I5,1X,I5,2X,ES13.3,1X,ES13.3)
!
!}}}
    END SUBROUTINE OFE_TRUSSD
!
    SUBROUTINE OFE1_TRUSSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OEL_TRUSSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 10
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMOLP
! *         NPAR(4) = NITG
! *         NPAR(5) = NWEI
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NUMOLP
    INTEGER :: N,ELEI,ELEJ,INDI,INDJ,KN
    INTEGER :: NOLD,EIOLD,EJOLD,INDIOLD,INDJOLD,KNOLD
    REAL(8) :: RADI,RADJ
    REAL(8) :: RIOLD,RJOLD
    INTEGER :: NUMN,NUM
    INTEGER :: NEI,NEJ
    INTEGER :: AllocateStatus
    INTEGER :: NCOUNT
    CHARACTER(20),DIMENSION(6) :: OFE_TRUSS_KW
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    OFE_TRUSS_KW(1) = 'MTYPE:'
    OFE_TRUSS_KW(2) = 'NUMOLP:'
    OFE_TRUSS_KW(3) = 'NITG:'
    OFE_TRUSS_KW(4) = 'NWEI:'
    OFE_TRUSS_KW(5) = '\END'
    OFE_TRUSS_KW(6) = 'NDIV:'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=NITP_IN
    NPAR(5)=1
    NPAR(6)=1
!   
    I = KeyWord(OFE_TRUSS_KW,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE(6)
           J = GetInt()
           NPAR(6)=J
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OEL_TRUSSD"
           STOP
       END SELECT
       I = KeyWord(OFE_TRUSS_KW,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(6),NPAR(5)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAP REGION .EQ.0"
        STOP
    END IF
    NUMOLP=NPAR(3)
!-------- Allocate memory
    ALLOCATE (NDINT(4,NUMOLP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(2,NUMOLP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!-------- INITIALIZATION ------------------
    NDINT=0
    INTPAR=0.0D0
!
    WRITE(IOUT,2000)
!---------------------------------------------------------------
    NCOUNT=0
    READ(IIN,*) N, ELEI, ELEJ, INDI, INDJ, RADI, RADJ, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NUMOLP) THEN
          WRITE(*,*) "ERROR N .GT. NUMOLP"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDINT(1,N)=ELEI
       NDINT(2,N)=ELEJ
       NDINT(3,N)=INDI
       NDINT(4,N)=INDJ
!
       INTPAR(1,N)=RADI
       INTPAR(2,N)=RADJ
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NEI=EIOLD
              NEJ=EJOLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NEI=NEI+KNOLD
                 NEJ=NEJ+KNOLD
!
                 NDINT(1,K)=NEI
                 NDINT(2,K)=NEJ
                 NDINT(3,K)=INDIOLD
                 NDINT(4,K)=INDJOLD
                 
                 INTPAR(1,K)=RIOLD
                 INTPAR(2,K)=RJOLD 
!
                 WRITE(IOUT,2010) K,NEI,NEJ,INDIOLD,INDJOLD,RIOLD,RJOLD              
                 NCOUNT=NCOUNT+1
                 CALL COLHT_OFE_TRUSS(NEI,NEJ)
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2010) N,ELEI,ELEJ,INDI,INDJ,RADI,RADJ
        CALL COLHT_OFE_TRUSS(ELEI,ELEJ)
!     
        NOLD=N
        KNOLD=KN
        EIOLD=ELEI
        EJOLD=ELEJ
        INDIOLD=INDI
        INDJOLD=INDJ
        RIOLD=RADI
        RJOLD=RADJ
!
        READ(IIN,*) N, ELEI, ELEJ, INDI, INDJ, RADI, RADJ, KN
    END DO
!
    IF(NCOUNT.NE.NUMOLP) THEN
       WRITE(*,*) "ERROR, OVERLAP REGION NUMBER .LT. NUMOLP!"
       STOP
    ENDIF
!--------------------------------------------------
    RETURN
!--------------------------------------------------
 1000 Format(/,&     
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF INTEGRATION REGION     =',I5,/,&
              '    THE NUMBER OF GAUSS INTEGRATION POINTS   =',I5,/,&
              '    THE NUMBER OF SUBINTERVALS               =',I5,/,&
              '    THE WEIGHT FUNCTION NUMBER       =',I5,/,&
              '          1: 1-6S^2+8S^3-3S^4         ',/)
 2000 FORMAT(/,&
              'OVERLAP    ELE       ELE     IND     IND     RAD     RAD ',/,&
              'NUMBER       I         J       I       J       I       J ',/)
 2010 FORMAT(I5,1X,I5,1X,I5,1X,I5,1X,I5,2X,ES13.3,1X,ES13.3)
!
!}}}
    END SUBROUTINE OFE1_TRUSSD
!
    SUBROUTINE COLHT_OFE_TRUSS(NI,NJ)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NI,NJ
!   IN-FUNCTION
    INTEGER :: IDD,LS
    INTEGER :: I
    INTEGER :: KK,ME
!
    LS=100000
!
    IDD=ID(1,NI)
    IF(IDD.NE.0) THEN
        IF(IDD-LS.LT.0) LS=IDD
    ENDIF
!
    IDD=ID(1,NJ)
    IF(IDD.NE.0) THEN
        IF(IDD-LS.LT.0) LS=IDD
    ENDIF
!
    IDD=ID(1,NI)
    IF(IDD.NE.0) THEN
        DO I=1,NDOF(1,NI)
           KK=IDD+I-1
           ME=KK-LS
           IF(ME.GT.MHT(KK)) MHT(KK) = ME
        ENDDO
    ENDIF
!
    IDD=ID(1,NJ)
    IF(IDD.NE.0) THEN
       DO I=1,NDOF(1,NJ)
          KK=IDD+I-1
          ME=KK-LS
          IF(ME.GT.MHT(KK)) MHT(KK) = ME
       ENDDO
    ENDIF
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE_TRUSS
!
    SUBROUTINE OVLP_QUADSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_QUADS
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 9
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMINT
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NEleAss
! *         NPAR(6) = NITG
! *         NPAR(7) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: IDD
    INTEGER :: NumInt,NumEle,NEleAss
    INTEGER :: N,M,II,JJ,KK,KN,ME
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,KOLD,NII,NJJ,NKK
    REAL(8) :: THIC,THICOLD
    INTEGER :: IND1,IND2,IND1OLD,IND2OLD
    REAL(8) :: RAD1,R1OLD,RAD2,R2OLD
    INTEGER :: INDSP,INDSPOLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
    INTEGER :: NDOFF
!
    CHARACTER(20),DIMENSION(7) :: KWTRUSS_OVLP
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_OVLP(1) = 'MTYPE:'
    KWTRUSS_OVLP(2) = 'NUMINT:'
    kWTRUSS_OVLP(3) = 'NUMELE:'
    KWTRUSS_OVLP(4) = 'ELEASS:'
    KWTRUSS_OVLP(5) = 'NGDOF:'
    KWTRUSS_OVLP(6) = 'ITYPE:'
    KWTRUSS_OVLP(7) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=0
    NPAR(6)=NITP_IN
    NPAR(7)=1
!   
    I = KeyWord(KWTRUSS_OVLP,7)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           J = GetInt()
           NPAR(7)=J
         CASE(7)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_OVLP,7)
    END DO
!
    WRITE(IOUT,1000) NPAR(7),NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF INTEGRATION REGION .EQ.0"
        STOP
    END IF
    NumInt=NPAR(3)
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF ELEMENT NUMBER .EQ.0"
        STOP
    END IF
    NumEle=NPAR(4)
    NEleAss=NPAR(5)
!-------- Allocate memory
    ALLOCATE (NDINT(3,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(1,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (NDMFS(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"         
    ALLOCATE (NumLM(NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!---------------------------------------------------------------
    DO I=1,NumInt
       DO J=1,3
          NDINT(J,I)=0
       ENDDO
       INTPAR(1,I)=0.0D0
    END DO
!
    DO I=1,NumEle
       NumLM(I)=0
       DO J=1,NEleAss
          LM(J,I)=0
       END DO
!
       DO J=1,4
          ELEPAR(J,I)=0.0D0
       ENDDO
!
       DO J=1,6
          NDMFS(J,I)=0
       ENDDO
    END DO
!---------------------------------------------------------------
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND1,IND2,RAD1,RAD2,INDSP,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND1        
       NDMFS(3,N)=IND2
       NDMFS(4,N)=INDSP       
       ELEPAR(1,N)=RAD1
       ELEPAR(2,N)=RAD2
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=IND1OLD
                 NDMFS(3,K)=IND2OLD
                 NDMFS(4,K)=INDSPOLD
                 ELEPAR(1,K)=R1OLD
                 ELEPAR(2,K)=R2OLD
                 WRITE(IOUT,2250) K,NII,IND1OLD,IND2OLD,R1OLD,R2OLD,INDSPOLD               
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND1,IND2,RAD1,RAD2,INDSP
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        IND1OLD=IND1
        IND2OLD=IND2
        R1OLD=RAD1
        R2OLD=RAD2
        INDSPOLD=INDSP
!
        READ(IIN,*) N,II,IND1,IND2,RAD1,RAD2,INDSP,KN
    END DO
!
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, KK, THIC, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumInt) THEN
          WRITE(*,*) "ERROR N .GT. NumInt"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       NDINT(3,N)=KK
       INTPAR(1,N)=THIC
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              NKK=KOLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD
                 NKK=NKK+KNOLD

                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 NDINT(3,K)=NKK
                 INTPAR(1,K)=THICOLD
                 WRITE(IOUT,2050) K,NII,NJJ,NKK,THICOLD                   
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2050) N,II,JJ,KK,THIC
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        KOLD=KK
        THICOLD=THIC
!
        READ(IIN,*) N,II,JJ,KK,THIC,KN
    END DO
!
    IF(NCOUNT.NE.NumInt) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumInt"
       STOP
    END IF
!--------------------------------------------------
    CALL OVLP_QUADS_ASSD(NumInt,NEleAss,NumEle)
!
!   GENERATE MHT
!
    DO I=1,NumInt
       II=NDINT(1,I)
       JJ=NDINT(2,I)
       KK=NDINT(3,I)
       CALL COLHT_OVLP_QUADS(NEleAss,II,JJ,KK,I)
    ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE:QUADS          =',I5,/,&
              '               1: PLANE STRAIN        ',/,&
              '               2: PLANE STRESS        ',/,&       
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF QUADS ELEMENTS     =',I5,/,&
              '    THE NUMBER OF QUADS NODES        =',I5,/,&
              '    THE MAX ASSOCIATION NUMBER       =',I5,/,&
              '    THE DEGREE PRECISION OF THE GAUSS INTEGRATION     =',I5)
 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE    THIC   ',/,&
              ' NUMBER-N      I        J       K             ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5,2X,ES13.3)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND(2,N)     RADIUS(2,N)   INDSP ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,2I5,2X,2ES13.3,2X,I5)
!}}}
    END SUBROUTINE OVLP_QUADSD
!
    SUBROUTINE OFE1_QUADSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE OFE1_QUADSD
! *      READ AND STORE QUADS DATA
! *      QUADS
! *         NPAR(1) = 12
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMOFE_Triangle
! *         NPAR(8) = NUMOFE_Rectangle
! *         NPAR(10)= NUMOFE_Triangle2
! *         NPAR(9) = NINTG
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NGDOF
! *         NPAR(6) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *         NPAR(7) = NWEI 
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumOFE,NumEle,NumCpl,NumCpl2
    INTEGER :: N,II,JJ,KK,LL,KN
    INTEGER :: I1,I2,I3,I4,I5,I6
    INTEGER :: NI1,NI2,NI3,NI4,NI5,NI6,I1OLD,I2OLD,I3OLD,I4OLD,I5OLD,I6OLD
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,KOLD,LOLD,NII,NJJ,NKK,NLL
    INTEGER :: IND1,IND2,IND1OLD,IND2OLD
    REAL(8) :: RAD1,R1OLD,RAD2,R2OLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
    CHARACTER(20),DIMENSION(10) :: KWQUADS_OFE
    LOGICAL :: INTIND=.false.
!
    INTEGER :: NOD1,NOD2,NOD3
    REAL(8),DIMENSION(2) :: XY1,XY2,XY3,XYC1,XYC2,XYC3
    REAL(8) :: L1,L2,L3,R1,R2,R3
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_OFE(1) = 'MTYPE:'
    KWQUADS_OFE(2) = 'NUMOFE:'
    kWQUADS_OFE(3) = 'NUMELE:'
    KWQUADS_OFE(4) = 'NGDOF:'
    KWQUADS_OFE(5) = 'ITYPE:'
    KWQUADS_OFE(6) = 'NWEI:'
    KWQUADS_OFE(7) = '\END'
    KWQUADS_OFE(8) = 'NUMCPL:'
    KWQUADS_OFE(9) = 'NINTG:'
    KWQUADS_OFE(10)= 'NUMCPL2:'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=6
    NPAR(6)=1
    NPAR(7)=1
    NPAR(8)=0
    NPAR(9)=3
    NPAR(10)=0
!   
    I = KeyWord(KWQUADS_OFE,10)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           J = GetInt()
           NPAR(7)=J
         CASE(7)
           EXIT
         CASE(8)
           J = GetInt()
           NPAR(8)=J
         CASE(9)
           J = GetInt()
           NPAR(9)=J
         CASE(10)
           J = GetInt()
           NPAR(10)=J
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN OFE1_QUADSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_OFE,10)
    END DO
!
    WRITE(IOUT,1000) NPAR(6),NPAR(2),NPAR(4),NPAR(3),NPAR(8),NPAR(10),NPAR(9),NPAR(7),NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF OVERLAPPING ELEMENTS .EQ.0"
        STOP
    END IF
    NumOFE=NPAR(3)
    NumEle=NPAR(4)
    NumCpl=NPAR(8)
    NumCpl2=NPAR(10)
!-------- Allocate memory
    ALLOCATE (NDMFS(5,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
    NDMFS=0
    ELEPAR(1,:)=0.0
    ELEPAR(2,:)=10000
!
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND1,IND2,RAD1,RAD2,KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NUMELE"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND1        
       NDMFS(3,N)=IND2     
!       ELEPAR(1,N)=RAD1
!       ELEPAR(2,N)=RAD2
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=IND1OLD
                 NDMFS(3,K)=IND2OLD
!                 ELEPAR(1,K)=R1OLD
!                 ELEPAR(2,K)=R2OLD
                 WRITE(IOUT,2250) K,NII,IND1OLD,IND2OLD,R1OLD,R2OLD               
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND1,IND2,RAD1,RAD2
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        IND1OLD=IND1
        IND2OLD=IND2
        R1OLD=RAD1
        R2OLD=RAD2
!
        READ(IIN,*) N,II,IND1,IND2,RAD1,RAD2,KN
    END DO
!---------------------------------------------------------------
!
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!
    IF(NumOFE.GT.0) THEN
       INTIND=.TRUE.
       ALLOCATE (NDINT(3,NumOFE),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
       CALL NewLine()
!---------------------------------------------------------------
       WRITE(IOUT,2000)
       NCOUNT=0
       READ(IIN,*) N, II, JJ, KK,KN
       DO WHILE(N.NE.0)
          IF(N.GT.NumOFE) THEN
             WRITE(*,*) "ERROR N .GT. NUMOFE"
             STOP
          END IF
          NCOUNT=NCOUNT+1
   
          NDINT(1,N)=II
          NDINT(2,N)=JJ          
          NDINT(3,N)=KK
!
          NOD1=NDMFS(1,II)
          NOD2=NDMFS(1,JJ)
          NOD3=NDMFS(1,KK)
!
          XY1(1)=X(NOD1)
          XY1(2)=Y(NOD1)
          XY2(1)=X(NOD2)
          XY2(2)=Y(NOD2)
          XY3(1)=X(NOD3)
          XY3(2)=Y(NOD3)
!
          XYC1=XY1-XY2
          XYC2=XY1-XY3
          XYC3=XY2-XY3
!
          L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
!
          L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
!
          L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
!
          R1=MAX(L1,L2)
          R2=MAX(L1,L3)
          R3=MAX(L2,L3)
          IF(ELEPAR(1,II).LT.R1) ELEPAR(1,II)=R1
          IF(ELEPAR(1,JJ).LT.R2) ELEPAR(1,JJ)=R2
          IF(ELEPAR(1,KK).LT.R3) ELEPAR(1,KK)=R3
          R1=MIN(L1,L2)
          R2=MIN(L1,L3)
          R3=MIN(L2,L3)
          IF(ELEPAR(2,II).GT.R1) ELEPAR(2,II)=R1
          IF(ELEPAR(2,JJ).GT.R2) ELEPAR(2,JJ)=R2
          IF(ELEPAR(2,KK).GT.R3) ELEPAR(2,KK)=R3
!
          IF (KNOLD.NE.0) THEN
             NUM=(N-NOLD)/KNOLD
             NUMN=NUM - 1
             IF (NUMN.GE.1) THEN
                 K=NOLD
                 NII=IOLD
                 NJJ=JOLD
                 NKK=KOLD
                 DO J=1,NUMN
                    K=K + KNOLD
   
                    NII=NII+KNOLD
                    NJJ=NJJ+KNOLD
                    NKK=NKK+KNOLD
   
                    NDINT(1,K)=NII
                    NDINT(2,K)=NJJ
                    NDINT(3,K)=NKK
!
                    NOD1=NDMFS(1,NII)
                    NOD2=NDMFS(1,NJJ)
                    NOD3=NDMFS(1,NKK)
          !
                    XY1(1)=X(NOD1)
                    XY1(2)=Y(NOD1)
                    XY2(1)=X(NOD2)
                    XY2(2)=Y(NOD2)
                    XY3(1)=X(NOD3)
                    XY3(2)=Y(NOD3)
          !
                    XYC1=XY1-XY2
                    XYC2=XY1-XY3
                    XYC3=XY2-XY3
          !
                    L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
          !
                    L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
          !
                    L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
          !
                    R1=MAX(L1,L2)
                    R2=MAX(L1,L3)
                    R3=MAX(L2,L3)
                    IF(ELEPAR(1,NII).LT.R1) ELEPAR(1,NII)=R1
                    IF(ELEPAR(1,NJJ).LT.R2) ELEPAR(1,NJJ)=R2
                    IF(ELEPAR(1,NKK).LT.R3) ELEPAR(1,NKK)=R3
!
                    R1=MIN(L1,L2)
                    R2=MIN(L1,L3)
                    R3=MIN(L2,L3)
                    IF(ELEPAR(2,NII).GT.R1) ELEPAR(2,NII)=R1
                    IF(ELEPAR(2,NJJ).GT.R2) ELEPAR(2,NJJ)=R2
                    IF(ELEPAR(2,NKK).GT.R3) ELEPAR(2,NKK)=R3
!
                    WRITE(IOUT,2050) K,NII,NJJ,NKK                  
                    CALL COLHT_OFE1_QUADS(NII,NJJ,NKK)
                    NCOUNT=NCOUNT+1
                 END DO
             END IF
           ENDIF
!
           WRITE(IOUT,2050) N,II,JJ,KK
           CALL COLHT_OFE1_QUADS(II,JJ,KK)
!     
           NOLD=N
           KNOLD=KN
           IOLD=II
           JOLD=JJ
           KOLD=KK
!
           READ(IIN,*) N,II,JJ,KK,KN
       END DO
!
       IF(NCOUNT.NE.NumOFE) THEN
          WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumOFE"
          STOP
       END IF
    ENDIF
!
    IF(NumCpl.GT.0) THEN
       INTIND=.TRUE.
       ALLOCATE (CPLINT(4,NumCPL),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
       CALL NewLine()
!---------------------------------------------------------------
       WRITE(IOUT,2010)
       NCOUNT=0
       READ(IIN,*) N, II, JJ, KK,LL, KN
       DO WHILE(N.NE.0)
          IF(N.GT.NumCPL) THEN
             WRITE(*,*) "ERROR N .GT. NumCPL"
             STOP
          END IF
          NCOUNT=NCOUNT+1
   
          CPLINT(1,N)=II
          CPLINT(2,N)=JJ          
          CPLINT(3,N)=KK
          CPLINT(4,N)=LL
!
          NOD1=NDMFS(1,II)
          NOD2=NDMFS(1,JJ)
          NOD3=NDMFS(1,KK)
!
          XY1(1)=X(NOD1)
          XY1(2)=Y(NOD1)
          XY2(1)=X(NOD2)
          XY2(2)=Y(NOD2)
          XY3(1)=X(NOD3)
          XY3(2)=Y(NOD3)
!
          XYC1=XY1-XY2
          XYC2=XY1-XY3
!
          L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
!
          L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
!
          R1=L3
          IF(ELEPAR(1,II).LT.R1) ELEPAR(1,II)=R1
          IF(ELEPAR(1,JJ).LT.R1) ELEPAR(1,JJ)=R1
          IF(ELEPAR(1,KK).LT.R1) ELEPAR(1,KK)=R1
          IF(ELEPAR(1,LL).LT.R1) ELEPAR(1,LL)=R1

          R1=L1
          IF(ELEPAR(2,II).GT.R1) ELEPAR(2,II)=R1
          IF(ELEPAR(2,JJ).GT.R1) ELEPAR(2,JJ)=R1
          IF(ELEPAR(2,KK).GT.R1) ELEPAR(2,KK)=R1
          IF(ELEPAR(2,LL).GT.R1) ELEPAR(2,LL)=R1
!
          IF (KNOLD.NE.0) THEN
             NUM=(N-NOLD)/KNOLD
             NUMN=NUM - 1
             IF (NUMN.GE.1) THEN
                 K=NOLD
                 NII=IOLD
                 NJJ=JOLD
                 NKK=KOLD
                 NLL=LOLD
                 DO J=1,NUMN
                    K=K + KNOLD
                    
                    NII=NII+KNOLD
                    NJJ=NJJ+KNOLD
                    NKK=NKK+KNOLD
                    NLL=NLL+KNOLD
  
                    CPLINT(1,K)=NII
                    CPLINT(2,K)=NJJ          
                    CPLINT(3,K)=NKK
                    CPLINT(4,K)=NLL

                    NOD1=NDMFS(1,NII)
                    NOD2=NDMFS(1,NJJ)
                    NOD3=NDMFS(1,NKK)

                    XY1(1)=X(NOD1)
                    XY1(2)=Y(NOD1)
                    XY2(1)=X(NOD2)
                    XY2(2)=Y(NOD2)
                    XY3(1)=X(NOD3)
                    XY3(2)=Y(NOD3)
          !
                    XYC1=XY1-XY2
                    XYC2=XY1-XY3
          !
                    L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
          !
                    L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
          !
                    R1=L3
                    IF(ELEPAR(1,NII).LT.R1) ELEPAR(1,NII)=R1
                    IF(ELEPAR(1,NJJ).LT.R1) ELEPAR(1,NJJ)=R1
                    IF(ELEPAR(1,NKK).LT.R1) ELEPAR(1,NKK)=R1
                    IF(ELEPAR(1,NLL).LT.R1) ELEPAR(1,NLL)=R1
          
                    R1=L1
                    IF(ELEPAR(2,NII).GT.R1) ELEPAR(2,NII)=R1
                    IF(ELEPAR(2,NJJ).GT.R1) ELEPAR(2,NJJ)=R1
                    IF(ELEPAR(2,NKK).GT.R1) ELEPAR(2,NKK)=R1
                    IF(ELEPAR(2,NLL).GT.R1) ELEPAR(2,NLL)=R1
!
                    WRITE(IOUT,2060) K,NII,NJJ,NKK,NLL                 
                    CALL COLHT_OFE1_COUPLE(NII,NJJ,NKK,NLL)
                    NCOUNT=NCOUNT+1
                 END DO
             END IF
           ENDIF
!
           WRITE(IOUT,2060) N,II,JJ,KK,LL
           CALL COLHT_OFE1_COUPLE(II,JJ,KK,LL)
!     
           NOLD=N
           KNOLD=KN
           IOLD=II
           JOLD=JJ
           KOLD=KK
           LOLD=LL
!
           READ(IIN,*) N,II,JJ,KK,LL,KN
       END DO
!
       IF(NCOUNT.NE.NumCPL) THEN
          WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumCPL"
          STOP
       ENDIF
    ENDIF
!
    IF(NumCpl2.GT.0) THEN
       INTIND=.TRUE.
       ALLOCATE (CPLINT(6,NumCpl2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!---------------------------------------------------------------
       CALL NewLine()
!---------------------------------------------------------------
       WRITE(IOUT,2020)
       NCOUNT=0
       READ(IIN,*) N, I1, I2, I3, I4, I5, I6, KN
       DO WHILE(N.NE.0)
          IF(N.GT.NumCpl2) THEN
             WRITE(*,*) "ERROR N .GT. NumCPL2"
             STOP
          END IF
          NCOUNT=NCOUNT+1
!   
          CPLINT(1,N)=I1
          CPLINT(2,N)=I2          
          CPLINT(3,N)=I3
          CPLINT(4,N)=I4
          CPLINT(5,N)=I5
          CPLINT(6,N)=I6
!
          NOD1=NDMFS(1,I1)
          NOD2=NDMFS(1,I2)
          NOD3=NDMFS(1,I3)        
!
          XY1(1)=X(NOD1)
          XY1(2)=Y(NOD1)
          XY2(1)=X(NOD2)
          XY2(2)=Y(NOD2)
          XY3(1)=X(NOD3)
          XY3(2)=Y(NOD3)
!
          XYC1=XY1-XY2
          XYC2=XY1-XY3
          XYC3=XY2-XY3
!
          L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
!
          L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
!
          L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
!
          R1=MAX(L1,L2)
          R2=MAX(L1,L3)
          R3=MAX(L2,L3)
          IF(ELEPAR(1,I1).LT.R1) ELEPAR(1,I1)=R1
          IF(ELEPAR(1,I2).LT.R2) ELEPAR(1,I2)=R2
          IF(ELEPAR(1,I3).LT.R3) ELEPAR(1,I3)=R3
!
          R1=MIN(L1,L2)
          R2=MIN(L1,L3)
          R3=MIN(L2,L3)
          IF(ELEPAR(2,I1).GT.R1) ELEPAR(2,I1)=R1
          IF(ELEPAR(2,I2).GT.R2) ELEPAR(2,I2)=R2
          IF(ELEPAR(2,I3).GT.R3) ELEPAR(2,I3)=R3
!
          ELEPAR(:,I4)=1.0D0
          ELEPAR(:,I5)=1.0D0
          ELEPAR(:,I6)=1.0D0
!
          IF (KNOLD.NE.0) THEN
             NUM=(N-NOLD)/KNOLD
             NUMN=NUM - 1
             IF (NUMN.GE.1) THEN
                 K=NOLD
                 NI1=I1OLD
                 NI2=I2OLD
                 NI3=I3OLD
                 NI4=I4OLD
                 NI5=I5OLD
                 NI6=I6OLD
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NI1=NI1+KNOLD
                    NI2=NI2+KNOLD
                    NI3=NI3+KNOLD
                    NI4=NI4+KNOLD
                    NI5=NI5+KNOLD
                    NI6=NI6+KNOLD
!   
                    CPLINT(1,K)=NI1
                    CPLINT(2,K)=NI2
                    CPLINT(3,K)=NI3
                    CPLINT(4,K)=NI4
                    CPLINT(5,K)=NI5
                    CPLINT(6,K)=NI6
!
                    NOD1=NDMFS(1,NI1)
                    NOD2=NDMFS(1,NI2)
                    NOD3=NDMFS(1,NI3)
          !
                    XY1(1)=X(NOD1)
                    XY1(2)=Y(NOD1)
                    XY2(1)=X(NOD2)
                    XY2(2)=Y(NOD2)
                    XY3(1)=X(NOD3)
                    XY3(2)=Y(NOD3)
          !
                    XYC1=XY1-XY2
                    XYC2=XY1-XY3
                    XYC3=XY2-XY3
          !
                    L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
          !
                    L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
          !
                    L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
          !
                    R1=MAX(L1,L2)
                    R2=MAX(L1,L3)
                    R3=MAX(L2,L3)
!
                    IF(ELEPAR(1,NI1).LT.R1) ELEPAR(1,NI1)=R1
                    IF(ELEPAR(1,NI2).LT.R2) ELEPAR(1,NI2)=R2
                    IF(ELEPAR(1,NI3).LT.R3) ELEPAR(1,NI3)=R3
!
                    R1=MIN(L1,L2)
                    R2=MIN(L1,L3)
                    R3=MIN(L2,L3)
!
                    IF(ELEPAR(2,NI1).GT.R1) ELEPAR(2,NI1)=R1
                    IF(ELEPAR(2,NI2).GT.R2) ELEPAR(2,NI2)=R2
                    IF(ELEPAR(2,NI3).GT.R3) ELEPAR(2,NI3)=R3
!
                    ELEPAR(:,NI4)=1.0D0
                    ELEPAR(:,NI5)=1.0D0
                    ELEPAR(:,NI6)=1.0D0
!
                    WRITE(IOUT,2070) K,NI1,NI2,NI3,NI4,NI5,NI6                 
                    CALL COLHT_OFE1_COUPLE2(NI1,NI2,NI3,NI4,NI5,NI6)
                    NCOUNT=NCOUNT+1
                 END DO
             END IF
           ENDIF
!
           WRITE(IOUT,2070) N,I1,I2,I3,I4,I5,I6
           CALL COLHT_OFE1_COUPLE2(I1,I2,I3,I4,I5,I6)
!     
           NOLD=N
           KNOLD=KN
           I1OLD=I1
           I2OLD=I2
           I3OLD=I3
           I4OLD=I4
           I5OLD=I5
           I6OLD=I6
!
           READ(IIN,*) N, I1, I2, I3, I4, I5, I6, KN
       END DO
!
       IF(NCOUNT.NE.NumCpl2) THEN
          WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumCpl2"
          STOP
       ENDIF
    ENDIF
!
    IF (.NOT.INTIND) THEN
       WRITE(*,*) "ERROR, OVERLAP REGIONS SHOULD BE GREATER THAN 0!"
       STOP
    ENDIF
!   
    ELEPAR(2,:)=ELEPAR(1,:)/ELEPAR(2,:)
    ELEPAR(1,:)=ELEPAR(2,:)**3*ELEPAR(1,:)
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE:QUADS            =',I5,/,&
              '               1: PLANE STRAIN          ',/,&
              '               2: PLANE STRESS          ',/,&       
              '    THE MATERIAL SET NUMBER            =',I5,/,&
              '    THE NUMBER OF OVERLAPPING ELEMENTS =',I5,/,&
              '    THE NUMBER OF OVERLAP REGIONS      =',I5,/,&
              '    THE NUMBER OF COUPLING REGIONS     =',I5,/,&
              '    THE NUMBER OF COUPLING REGIONS2     =',I5,/,&
              '    THE NUMBER OF INTEGRATION POINTS (X-DIRECTION) FOR THE OVERLAP REGIONS =',I5,/,&
              '    THE WEIGHT FUNCTION NUMBER         =',I5,/,&
              '    THE DEGREE PRECISION OF THE GAUSS INTEGRATION     =',I5)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND(2,N)     RADIUS(2,N) ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,2I5,2X,2ES13.3)

 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     ELE           ',/,&
              ' NUMBER-N      I        J       K             ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5)
 2010 FORMAT(/,&
              ' INTEGRATION   ELE      ELE      ELE     ELE     ',/,&
              ' NUMBER-N      I        J        K       L       ',/)
 2060 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5)
 2020 FORMAT(/,&
              ' INTEGRATION   ELE      ELE      ELE     ELE     ELE     ELE ',/,&
              ' NUMBER-N      1        2        3       4       5       6',/)      
 2070 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5)
!}}}
    END SUBROUTINE OFE1_QUADSD
!
    SUBROUTINE COLHT_OFE1_QUADS(II,JJ,KK)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, JJ, KK
!   IN-FUNCTION
    INTEGER :: NI,NJ,NK
    INTEGER :: IDD,LS   
    INTEGER :: I,J
    INTEGER :: KKK,ME
!       
    NI=NDMFS(1,II)
    NJ=NDMFS(1,JJ)
    NK=NDMFS(1,KK)
!
    LS=100000
!
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
    ENDDO
!------------------------------------------
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           DO I=1,NDOF(J,NI)
              KKK=IDD+I-1
              ME=KKK-LS
              IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
           ENDDO
       ENDIF
!
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NJ)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NK)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE1_QUADS
!
    SUBROUTINE COLHT_OFE1_COUPLE(II,JJ,KK,LL)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, JJ, KK,LL
!   IN-FUNCTION
    INTEGER :: NI,NJ,NK,NL
    INTEGER :: IDD,LS   
    INTEGER :: I,J
    INTEGER :: KKK,ME
!       
    NI=NDMFS(1,II)
    NJ=NDMFS(1,JJ)
    NK=NDMFS(1,KK)
    NL=NDMFS(1,LL)
!
    LS=100000
!
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NL)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
    ENDDO
!------------------------------------------
    DO J=1,2
       IDD=ID(J,NI)
       IF(IDD.NE.0) THEN
           DO I=1,NDOF(J,NI)
              KKK=IDD+I-1
              ME=KKK-LS
              IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
           ENDDO
       ENDIF
!
       IDD=ID(J,NJ)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NJ)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NK)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NK)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NL)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NL)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE1_COUPLE
!
    SUBROUTINE COLHT_OFE1_COUPLE2(I1,I2,I3,I4,I5,I6)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: I1,I2,I3,I4,I5,I6
!   IN-FUNCTION
    INTEGER :: NI1,NI2,NI3,NI4,NI5,NI6
    INTEGER :: NI,NJ,NK,NL
    INTEGER :: IDD,LS   
    INTEGER :: I,J
    INTEGER :: KKK,ME
!       
    NI1=NDMFS(1,I1)
    NI2=NDMFS(1,I2)
    NI3=NDMFS(1,I3)
    NI4=NDMFS(1,I4)
    NI5=NDMFS(1,I5)
    NI6=NDMFS(1,I6)
!
    LS=100000
!
    DO J=1,2
       IDD=ID(J,NI1)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NI2)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NI3)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NI4)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NI5)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
   !
       IDD=ID(J,NI6)
       IF(IDD.NE.0) THEN
           IF(IDD-LS.LT.0) LS=IDD
       ENDIF
    ENDDO
!------------------------------------------
    DO J=1,2
       IDD=ID(J,NI1)
       IF(IDD.NE.0) THEN
           DO I=1,NDOF(J,NI1)
              KKK=IDD+I-1
              ME=KKK-LS
              IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
           ENDDO
       ENDIF
!
       IDD=ID(J,NI2)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NI2)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NI3)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NI3)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NI4)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NI4)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NI5)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NI5)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
!
       IDD=ID(J,NI6)
       IF(IDD.NE.0) THEN
          DO I=1,NDOF(J,NI6)
             KKK=IDD+I-1
             ME=KKK-LS
             IF(ME.GT.MHT(KKK)) MHT(KKK) = ME
          ENDDO
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OFE1_COUPLE2
!
    SUBROUTINE COLHT_OVLP_QUADS(NEleAss,II,JJ,KK,NumI)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,II,JJ,KK
    INTEGER,INTENT(IN) :: NumI
!   IN-FUNCTION
    INTEGER :: NumLoopII,NumLoopJJ,NumLoopKK
    INTEGER :: LSMin
    INTEGER :: ELE
    INTEGER :: N,IDD,IDD2,ME
    INTEGER :: I,J,K
!
    LSMin=100000
!
    NumLoopII=NumLM(II)
    NumLoopJJ=NumLM(JJ)
    NumLoopKK=NumLM(KK)
!
    DO I=1,NumLoopII
       ELE=LM(I,II)
       N=NDMFS(1,ELE)
!
       DO J=1,2
          IDD=ID(J,N) 
          IF(IDD.NE.0) THEN
             IF(IDD-LSMin.LT.0)  LSMin=IDD
          ENDIF
       ENDDO
!
    ENDDO
!
    DO I=1,NumLoopJJ
       ELE=LM(I,JJ)
       N=NDMFS(1,ELE)
!
       DO J=1,2
          IDD=ID(J,N) 
          IF(IDD.NE.0) THEN
             IF(IDD-LSMin.LT.0)  LSMin=IDD
          ENDIF
       ENDDO
    ENDDO
!
    DO I=1,NumLoopKK
       ELE=LM(I,KK)
       N=NDMFS(1,ELE)
!
       DO J=1,2
          IDD=ID(J,N) 
          IF(IDD.NE.0) THEN
             IF(IDD-LSMin.LT.0)  LSMin=IDD
          ENDIF
       ENDDO
    ENDDO  
!
    DO I=1,NumLoopII
       ELE=LM(I,II)
       N=NDMFS(1,ELE)
       DO J=1,2
          IDD=ID(J,N)
          IF(IDD.NE.0) THEN
             DO K=1,NDOF(J,N)
                IDD2=IDD+K-1
                ME=IDD2-LSMin
                IF(ME.GT.MHT(IDD2)) MHT(IDD2) = ME
             ENDDO
          ENDIF
       ENDDO
    ENDDO 
!
    DO I=1,NumLoopJJ
       ELE=LM(I,JJ)
       N=NDMFS(1,ELE)
       DO J=1,2
          IDD=ID(J,N)
          IF(IDD.NE.0) THEN
             DO K=1,NDOF(J,N)
                IDD2=IDD+K-1
                ME=IDD2-LSMin
                IF(ME.GT.MHT(IDD2)) MHT(IDD2) = ME
             ENDDO
          ENDIF
       ENDDO
    ENDDO 
!
    DO I=1,NumLoopKK
       ELE=LM(I,KK)
       N=NDMFS(1,ELE)
       DO J=1,2
          IDD=ID(J,N)
          IF(IDD.NE.0) THEN
             DO K=1,NDOF(J,N)
                IDD2=IDD+K-1
                ME=IDD2-LSMin
                IF(ME.GT.MHT(IDD2)) MHT(IDD2) = ME
             ENDDO
          ENDIF
       ENDDO
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OVLP_QUADS
!
    SUBROUTINE OVLP_QUADS_ASSD(NumInt,NEleAss,NumEle)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM 
! *      NumInt   :  Number of Integration regions
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      NDINDTD :  
! *      LMD  :
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: NumInt, NEleAss, NumEle
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER :: ELE1,ELE2,ELE3,IND1,IND2,IND3
    INTEGER :: K,K1,K2,K3
    INTEGER :: IND
!
    DO I=1,NumEle
       NumLM(I)=NumLM(I)+1
       K=NumLM(I)
!
       IF(K.GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!",K,I
          STOP
       ENDIF
!
       LM(K,I)=I 
    ENDDO
!
    DO I=1,NumInt
       ELE1=NDINT(1,I)
       ELE2=NDINT(2,I)
       ELE3=NDINT(3,I)
       IND1=NDMFS(2,ELE1)
       IND2=NDMFS(2,ELE2)
       IND3=NDMFS(2,ELE3)
!
       IF (IND1.EQ.0) THEN      
           NumLM(ELE1)=NumLM(ELE1)+1
           K1=NumLM(ELE1)
           IF(K1.GT.NEleAss) THEN
              WRITE(*,*) "ERROR, NEleAss too small!",K1,ELE1
              STOP
           ENDIF
           LM(K1,ELE1)=ELE2
       ENDIF
!
       IF (IND2.EQ.0) THEN
           NumLM(ELE2)=NumLM(ELE2)+1
           K2=NumLM(ELE2)
           IF(K2.GT.NEleAss) THEN
              WRITE(*,*) "ERROR, NEleAss too small!",K2,ELE2
              STOP
           ENDIF
           LM(K2,ELE2)=ELE3
       ENDIF
!       
       IF (IND3.EQ.0) THEN
           NumLM(ELE3)=NumLM(ELE3)+1
           K3=NumLM(ELE3)
           IF(K3.GT.NEleAss) THEN
              WRITE(*,*) "ERROR, NEleAss too small!",K3,ELE3
              STOP
           ENDIF
           LM(K3,ELE3)=ELE1
       ENDIF
! 
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OVLP_QUADS_ASSD
!
    SUBROUTINE TRUSSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 1
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumEle
    INTEGER :: N,II,JJ,KN
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    REAL(8) :: AREA,AOLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(3) :: KWTRUSS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS(1) = 'MTYPE:'
    KWTRUSS(2) = 'NUMELE:'
    KWTRUSS(3) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
!   
    I = KeyWord(KWTRUSS,3)
    DO WHILE(I<=3)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS,3)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3)
!
    IF(NPAR(3).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
    NumEle=NPAR(3)
!
!-------- Allocate memory
!
        ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,6
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ, AREA, KN
        DO WHILE(N.NE.0)
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           XYZ(1,N)=X(II)
           XYZ(2,N)=Y(II)
           XYZ(3,N)=Z(II)
!
           XYZ(4,N)=X(JJ)
           XYZ(5,N)=Y(JJ)
           XYZ(6,N)=Z(JJ)
!          
           ELEPAR(1,N)=AREA
           DO I=1,3
             LM(I,N)=ID(I,II)
             LM(I+3,N)=ID(I,JJ)
           END DO
           CALL COLHT(6,LM(:,N))
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                  NJJ=JOLD
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+KNOLD
                    NJJ=NJJ+KNOLD
!
                    XYZ(1,K)=X(NII)
                    XYZ(2,K)=Y(NII)
                    XYZ(3,K)=Z(NII)
!
                    XYZ(4,K)=X(NJJ)
                    XYZ(5,K)=Y(NJJ)
                    XYZ(6,K)=Z(NJJ)
!
                    ELEPAR(1,K)=AOLD
!
                   DO I=1,3
                     LM(I,K)=ID(I,NII)
                     LM(I+3,K)=ID(I,NJJ)
                   END DO
                   CALL COLHT(6,LM(:,K))
                   WRITE(IOUT,2050) K,NII,NJJ,AOLD                   
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
            WRITE(IOUT,2050) N,II,JJ,AREA
!        
            NOLD=N
            KNOLD=KN
            IOLD=II
            JOLD=JJ
            AOLD=AREA
!
           READ(IIN,*) N,II,JJ,AREA,KN
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS     =',I5,/)
 2000 FORMAT(/,&
              ' ELEMENT     NODE     NODE     AREA   ',/,&
              ' NUMBER-N      I        J       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
!}}}
    END SUBROUTINE TRUSSD
!
    SUBROUTINE TRUSS_OVLPD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_OVLPD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 6
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMINT
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NEleAss
! *         NPAR(6) = NITG
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumInt,NumEle,NEleAss
    INTEGER :: N,M,II,JJ,KN
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    REAL(8) :: AREA,AOLD
    INTEGER :: IND,INDOLD
    REAL(8) :: RADIUS,ROLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(6) :: KWTRUSS_OVLP
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_OVLP(1) = 'MTYPE:'
    KWTRUSS_OVLP(2) = 'NUMINT:'
    kWTRUSS_OVLP(3) = 'NUMELE:'
    KWTRUSS_OVLP(4) = 'ELEASS:'
    KWTRUSS_OVLP(5) = 'NITG:'
    KWTRUSS_OVLP(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=0
    NPAR(6)=NITP_IN
!   
    I = KeyWord(KWTRUSS_OVLP,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_OVLP,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF INTEGRATION REGION .EQ.0"
        STOP
    END IF
    NumInt=NPAR(3)
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF ELEMENT NUMBER .EQ.0"
        STOP
    END IF
    NumEle=NPAR(4)
    NEleAss=NPAR(5)
!-------- Allocate memory
    ALLOCATE (NDINT(2,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(1,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (NDMFS(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
    ALLOCATE (WEIND(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
!---------------------------------------------------------------
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       INTPAR(1,I)=0.0D0
    END DO

    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
          WEIND(J,I)=0.0D0
       END DO
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!---------------------------------------------------------------
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND, RADIUS, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND         
       ELEPAR(1,N)=RADIUS
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
!
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=INDOLD
                 ELEPAR(1,K)=ROLD
                 WRITE(IOUT,2250) K,NII,INDOLD,ROLD                 
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND,RADIUS
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        INDOLD=IND
        ROLD=RADIUS
!
        READ(IIN,*) N,II,IND,RADIUS,KN
    END DO
!
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, AREA, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumInt) THEN
          WRITE(*,*) "ERROR N .GT. NumInt"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       INTPAR(1,N)=AREA
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD

                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 INTPAR(1,K)=AOLD
                 WRITE(IOUT,2050) K,NII,NJJ,AOLD                   
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2050) N,II,JJ,AREA
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        AOLD=AREA
!
        READ(IIN,*) N,II,JJ,AREA,KN
    END DO
!
    IF(NCOUNT.NE.NumInt) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumInt"
       STOP
    END IF
!--------------------------------------------------
    CALL OVLP_ASSD(2,NumInt,NEleAss,NumEle,NDINT(1,1),LM(1,1))
!   GENERATE WEIGHTS OF EACH NODE
    DO I=1,NumEle
       CALL WEIND_CAL(I,NEleAss,LM(1,I),WEIND(1,I))
    ENDDO        
!
!   GENERATE MHT
!
    DO I=1,NumInt
       II=NDINT(1,I)
       JJ=NDINT(2,I)
       CALL COLHT_OVLP(NEleAss,II,JJ)
    ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS     =',I5,/,&
              '    THE NUMBER OF TRUSS NODES        =',I5,/,&
              '    THE MAX ASSOCIATION NUMBER       =',I5,/,&
              '    THE GAUSS INTEGRATION NUMBER     =',I5)
 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     AREA   ',/,&
              ' NUMBER-N      I        J              ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND     RADIUS      ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
!}}}
    END SUBROUTINE TRUSS_OVLPD
!
    SUBROUTINE TRUSS_OVLPD2(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_OVLPD
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 6
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMINT
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NEleAss
! *         NPAR(6) = NITG
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumInt,NumEle,NEleAss
    INTEGER :: N,M,II,JJ,KN
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    REAL(8) :: AREA,AOLD
    INTEGER :: IND,INDOLD
    REAL(8) :: RADIUS,ROLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(6) :: KWTRUSS_OVLP
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_OVLP(1) = 'MTYPE:'
    KWTRUSS_OVLP(2) = 'NUMINT:'
    kWTRUSS_OVLP(3) = 'NUMELE:'
    KWTRUSS_OVLP(4) = 'ELEASS:'
    KWTRUSS_OVLP(5) = 'NITG:'
    KWTRUSS_OVLP(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=0
    NPAR(6)=NITP_IN
!   
    I = KeyWord(KWTRUSS_OVLP,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_OVLP,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF INTEGRATION REGION .EQ.0"
        STOP
    END IF
    NumInt=NPAR(3)
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF ELEMENT NUMBER .EQ.0"
        STOP
    END IF
    NumEle=NPAR(4)
    NEleAss=NPAR(5)
!-------- Allocate memory
    ALLOCATE (NDINT(2,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(1,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (NDMFS(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"         
!---------------------------------------------------------------
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       INTPAR(1,I)=0.0D0
    END DO

    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
       END DO
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!---------------------------------------------------------------
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND, RADIUS, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND         
       ELEPAR(1,N)=RADIUS
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
!
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=INDOLD
                 ELEPAR(1,K)=ROLD
                 WRITE(IOUT,2250) K,NII,INDOLD,ROLD                 
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND,RADIUS
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        INDOLD=IND
        ROLD=RADIUS
!
        READ(IIN,*) N,II,IND,RADIUS,KN
    END DO
!
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, AREA, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumInt) THEN
          WRITE(*,*) "ERROR N .GT. NumInt"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       INTPAR(1,N)=AREA
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD

                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 INTPAR(1,K)=AOLD
                 WRITE(IOUT,2050) K,NII,NJJ,AOLD                   
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2050) N,II,JJ,AREA
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        AOLD=AREA
!
        READ(IIN,*) N,II,JJ,AREA,KN
    END DO
!
    IF(NCOUNT.NE.NumInt) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumInt"
       STOP
    END IF
!--------------------------------------------------
    CALL OVLP_ASSD(2,NumInt,NEleAss,NumEle,NDINT(1,1),LM(1,1))
!   GENERATE MHT
!
    DO I=1,NumInt
       II=NDINT(1,I)
       JJ=NDINT(2,I)
       CALL COLHT_OVLP(NEleAss,II,JJ)
    ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS     =',I5,/,&
              '    THE NUMBER OF TRUSS NODES        =',I5,/,&
              '    THE MAX ASSOCIATION NUMBER       =',I5,/,&
              '    THE GAUSS INTEGRATION NUMBER     =',I5)
 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     AREA   ',/,&
              ' NUMBER-N      I        J              ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND     RADIUS      ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
!}}}
    END SUBROUTINE TRUSS_OVLPD2
!
    SUBROUTINE TRUSS_OVLPD3(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_OVLPD3
! *      READ AND STORE TRUSS DATA
! *      TRUSS
! *         NPAR(1) = 6
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMINT
! *         NPAR(4) = NUMELE
! *         NPAR(5) = NEleAss
! *         NPAR(6) = NITG
! *      INPUT:
! *         IIN, IOUT 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: IDD
    INTEGER :: NumInt,NumEle,NEleAss
    INTEGER :: N,M,II,JJ,KN,ME
    INTEGER :: KNOLD,NOLD,NUM,NUMN
    INTEGER :: IOLD,JOLD,NII,NJJ
    REAL(8) :: AREA,AOLD
    INTEGER :: IND,INDOLD
    REAL(8) :: RADIUS,ROLD
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(6) :: KWTRUSS_OVLP
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_OVLP(1) = 'MTYPE:'
    KWTRUSS_OVLP(2) = 'NUMINT:'
    kWTRUSS_OVLP(3) = 'NUMELE:'
    KWTRUSS_OVLP(4) = 'ELEASS:'
    KWTRUSS_OVLP(5) = 'NITG:'
    KWTRUSS_OVLP(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=0
    NPAR(5)=0
    NPAR(6)=NITP_IN
!   
    I = KeyWord(KWTRUSS_OVLP,6)
    DO WHILE(.TRUE.)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_OVLP,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6)
!
    IF(NPAR(3).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF INTEGRATION REGION .EQ.0"
        STOP
    END IF
    NumInt=NPAR(3)
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF ELEMENT NUMBER .EQ.0"
        STOP
    END IF
    NumEle=NPAR(4)
    NEleAss=NPAR(5)
!-------- Allocate memory
    ALLOCATE (NDINT(4,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(1,NumInt),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (NDMFS(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"         
    ALLOCATE (NumLM(NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!---------------------------------------------------------------
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       NDINT(3,I)=0
       NDINT(4,I)=0
!
       INTPAR(1,I)=0.0D0
    END DO

    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
       END DO
       NumLM(I)=0
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!---------------------------------------------------------------
    WRITE(IOUT,2200)
    NCOUNT=0
    READ(IIN,*) N, II, IND, RADIUS, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumEle) THEN
          WRITE(*,*) "ERROR N .GT. NumEle"
          STOP
       END IF
       NCOUNT=NCOUNT+1
!
       NDMFS(1,N)=II
       NDMFS(2,N)=IND         
       ELEPAR(1,N)=RADIUS
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              DO J=1,NUMN
                 K=K + KNOLD
!
                 NII=NII+KNOLD
!
                 NDMFS(1,K)=NII
                 NDMFS(2,K)=INDOLD
                 ELEPAR(1,K)=ROLD
                 WRITE(IOUT,2250) K,NII,INDOLD,ROLD                 
                 NCOUNT=NCOUNT+1
             END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2250) N,II,IND,RADIUS
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        INDOLD=IND
        ROLD=RADIUS
!
        READ(IIN,*) N,II,IND,RADIUS,KN
    END DO
!
    IF(NCOUNT.NE.NumEle) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle!"
       STOP
    ENDIF
!---------------------------------------------------------------
    CALL NewLine()
!---------------------------------------------------------------
    WRITE(IOUT,2000)
    NCOUNT=0
    READ(IIN,*) N, II, JJ, AREA, KN
    DO WHILE(N.NE.0)
       IF(N.GT.NumInt) THEN
          WRITE(*,*) "ERROR N .GT. NumInt"
          STOP
       END IF
       NCOUNT=NCOUNT+1

       NDINT(1,N)=II
       NDINT(2,N)=JJ          
       INTPAR(1,N)=AREA
!
       IF (KNOLD.NE.0) THEN
          NUM=(N-NOLD)/KNOLD
          NUMN=NUM - 1
          IF (NUMN.GE.1) THEN
              K=NOLD
              NII=IOLD
              NJJ=JOLD
              DO J=1,NUMN
                 K=K + KNOLD

                 NII=NII+KNOLD
                 NJJ=NJJ+KNOLD

                 NDINT(1,K)=NII
                 NDINT(2,K)=NJJ
                 INTPAR(1,K)=AOLD
                 WRITE(IOUT,2050) K,NII,NJJ,AOLD                   
                 NCOUNT=NCOUNT+1
              END DO
          END IF
        ENDIF
!
        WRITE(IOUT,2050) N,II,JJ,AREA
!     
        NOLD=N
        KNOLD=KN
        IOLD=II
        JOLD=JJ
        AOLD=AREA
!
        READ(IIN,*) N,II,JJ,AREA,KN
    END DO
!
    IF(NCOUNT.NE.NumInt) THEN
       WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumInt"
       STOP
    END IF
!--------------------------------------------------
    CALL OVLP_ASSD3(4,NumInt,NEleAss,NumEle,NDINT(1,1),LM(1,1))
!
!   GENERATE MHT
!
    DO I=1,NumInt
       II=NDINT(1,I)
       JJ=NDINT(2,I)
       CALL COLHT_OVLP3(NEleAss,II,JJ,I)
    ENDDO
!
    DO I=1,NumEle
       N=NDMFS(1,I)
       IDD=ID(1,N)
       IF (IDD.NE.0) THEN
           ME=MHT(IDD)
           DO J=2,NDOF(1,N)
              K=IDD+J-1
              MHT(K)=ME+J-1
           ENDDO
       ENDIF
    ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER          =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS     =',I5,/,&
              '    THE NUMBER OF TRUSS NODES        =',I5,/,&
              '    THE MAX ASSOCIATION NUMBER       =',I5,/,&
              '    THE GAUSS INTEGRATION NUMBER     =',I5)
 2000 FORMAT(/,&
              ' INTEGRATION   ELE      ELE     AREA   ',/,&
              ' NUMBER-N      I        J              ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
 2200 FORMAT(/,&
              ' ELEMENT     II    IND     RADIUS      ',/,&
              ' NUMBER-N                              ',/)
 2250 FORMAT(I5,6X,I5,4X,I5,2X,ES13.3)
!}}}
    END SUBROUTINE TRUSS_OVLPD3
!
    SUBROUTINE WEIND_CAL(II,NEleAss,LMM,WEINDD)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II, NEleAss
    INTEGER,DIMENSION(NEleAss) :: LMM
    REAL(8),DIMENSION(NEleAss) :: WEINDD
!   ASSISTANT VARIABLES
    INTEGER :: I
    REAL(8) :: WEIT_SUM
    INTEGER :: N,M,ELE
    REAL(8) :: RAD,S(3),SS
!
    WEIT_SUM=0.0D0
    N=NDMFS(1,II)
    RAD=ELEPAR(1,N)
    RAD=1.0/RAD
    DO I=1,NEleAss
       ELE=LMM(I)
       IF (ELE.NE.0) THEN
           IF(ELE.EQ.II)  WEINDD(I)=1.0D0
!          M=NDMFS(1,ELE)
!          SS=ABS(X(M)-X(N))*RAD
!          IF ((1.0-SS).GT.1E-12) THEN
!              S(1)=SS*SS
!              S(2)=S(1)*SS
!              S(3)=S(2)*SS
!              WEINDD(I)=1-6*S(1)+8*S(2)-3*S(3)
!          ENDIF
!          WEIT_SUM=WEIT_SUM+WEINDD(I)
       ENDIF
    ENDDO
!
!    DO I=1,NEleAss
!       WEINDD(I)=WEINDD(I)/WEIT_SUM
!    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE WEIND_CAL
!
    SUBROUTINE COLHT_OVLP(NEleAss,II,JJ)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,II,JJ
!   IN-FUNCTION
    INTEGER :: LS
    INTEGER :: ELE
    INTEGER :: N,IDD,ME
    INTEGER :: I,J,KK
!
    LS=100000
!
    DO I=1,NEleAss
       ELE=LM(I,II)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IDD=ID(1,N)
          
          IF(IDD.NE.0) THEN
             IF(IDD-LS.LT.0) LS=IDD
          ENDIF
       ENDIF
    ENDDO
!
    DO I=1,NEleAss
       ELE=LM(I,JJ)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IDD=ID(1,N)
          IF(IDD.NE.0) THEN
             IF(IDD-LS.LT.0) LS=IDD
          ENDIF
       ENDIF
    ENDDO
!
    DO I=1,NEleAss
       ELE=LM(I,II)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IDD=ID(1,N)
          IF(IDD.NE.0) THEN
             DO J=1,NDOF(1,N)
                KK=IDD+J-1
                ME=KK-LS
                IF(ME.GT.MHT(KK)) MHT(KK) = ME
             ENDDO
          ENDIF
       ENDIF
    ENDDO 
!
    DO I=1,NEleAss
       ELE=LM(I,JJ)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IDD=ID(1,N)
          IF(IDD.NE.0) THEN
             DO J=1,NDOF(1,N)
                KK=IDD+J-1
                ME=KK-LS
                IF(ME.GT.MHT(KK)) MHT(KK) = ME
             ENDDO
          ENDIF
       ENDIF
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OVLP
!
    SUBROUTINE COLHT_OVLP3(NEleAss,II,JJ,NumI)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,II,JJ
    INTEGER,INTENT(IN) :: NumI
!   IN-FUNCTION
    INTEGER :: NumLoopII,NumLoopJJ
    INTEGER :: LSMin,LSMax
    INTEGER :: ELE
    INTEGER :: N,IDD,IDD2,ME
    INTEGER :: I
!
    LSMin=100000
    LSMax=0
!
    NumLoopII=NumLM(II)
    NumLoopJJ=NumLM(JJ)
!
    DO I=1,NumLoopII
       ELE=LM(I,II)
       N=NDMFS(1,ELE)
!
       IDD=ID(1,N) 
       IDD2=IDD+NDOF(1,N)-1       
       IF(IDD.NE.0) THEN
          IF(IDD-LSMin.LT.0) LSMin=IDD
          IF(IDD2-LSMax.GT.0) LSMax=IDD2
       ENDIF
    ENDDO
!
    DO I=1,NumLoopJJ
       ELE=LM(I,JJ)
       N=NDMFS(1,ELE)
!
       IDD=ID(1,N)
       IDD2=IDD+NDOF(1,N)-1
       IF(IDD.NE.0) THEN
          IF(IDD-LSMin.LT.0) LSMin=IDD
          IF(IDD2-LSMax.GT.0) LSMax=IDD2
       ENDIF
    ENDDO
!
    NDINT(3,NumI)=LSMin
    NDINT(4,NumI)=LSMax    
!
    DO I=1,NumLoopII
       ELE=LM(I,II)
       N=NDMFS(1,ELE)
       IDD=ID(1,N)
       IF(IDD.NE.0) THEN
          ME=IDD-LSMin
          IF(ME.GT.MHT(IDD)) MHT(IDD) = ME
       ENDIF
    ENDDO 
!
    DO I=1,NumLoopJJ
       ELE=LM(I,JJ)
       N=NDMFS(1,ELE)
       IDD=ID(1,N)
       IF(IDD.NE.0) THEN
          ME=IDD-LSMin
          IF(ME.GT.MHT(IDD)) MHT(IDD) = ME
       ENDIF
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE COLHT_OVLP3
!
    SUBROUTINE OVLP_ASSD(NUM1,NumInt,NEleAss,NumEle,NDINTD,LMD)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM 
! *      NumInt   :  Number of Integration regions
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      NDINDTD :  
! *      LMD  :
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER,INTENT(IN) :: NUM1, NumInt, NEleAss, NumEle
    INTEGER,DIMENSION(NUM1,NumInt)   :: NDINTD
    INTEGER,DIMENSION(NEleAss,NumEle):: LMD
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER,DIMENSION(NumEle) :: KK
    INTEGER :: ELE1,ELE2
    INTEGER :: K,K1,K2
!
    DO I=1,NumEle
       KK(I)=1
    ENDDO
!
    DO I=1,NumInt
       ELE1=NDINTD(1,I)
       ELE2=NDINTD(2,I)
       K1=KK(ELE1)
       K2=KK(ELE2)
!
       LMD(K2,ELE2)=ELE1
       KK(ELE2)=KK(ELE2)+1
       LMD(K1,ELE1)=ELE2
       KK(ELE1)=KK(ELE1)+1
       IF(KK(ELE1).GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!"
          STOP
       ENDIF
       IF(KK(ELE2).GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!"
          STOP
       ENDIF
    ENDDO
!
    DO I=1,NumEle
       K=KK(I)
       LMD(K,I)=I 
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OVLP_ASSD
!
    SUBROUTINE OVLP_ASSD3(NUM1,NumInt,NEleAss,NumEle,NDINTD,LMD)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM 
! *      NumInt   :  Number of Integration regions
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      NDINDTD :  
! *      LMD  :
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: NUM1, NumInt, NEleAss, NumEle
    INTEGER, DIMENSION(NUM1,NumInt)   :: NDINTD
    INTEGER, DIMENSION(NEleAss,NumEle):: LMD
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER :: ELE1,ELE2
    INTEGER :: K,K1,K2
!
    DO I=1,NumInt
       ELE1=NDINTD(1,I)
       ELE2=NDINTD(2,I)
!
       NumLM(ELE1)=NumLM(ELE1)+1
       K1=NumLM(ELE1)
       NumLM(ELE2)=NumLM(ELE2)+1
       K2=NumLM(ELE2)
       IF(K1.GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!",K1
          STOP
       ENDIF
!
       IF(K2.GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!",K2
          STOP
       ENDIF
!
       LMD(K1,ELE1)=ELE2
       LMD(K2,ELE2)=ELE1
    ENDDO
!
    DO I=1,NumEle
       NumLM(I)=NumLM(I)+1
       K=NumLM(I)
!
       IF(K.GT.NEleAss) THEN
          WRITE(*,*) "ERROR, NEleAss too small!",K
          STOP
       ENDIF
!
       LMD(K,I)=I 
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OVLP_ASSD3
!
    SUBROUTINE QUADSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE QUADSD
! *      READ AND STORE TRUSS DATA
! *  QUADS
! *      NPAR(1) = 2
! *      NPAR(2) = ELEMENT TYPE
! *                0: AXISYMMETRIC
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      NPAR(3) = MTYPE
! *      NPAR(4) = NUMELE
! *      NPAR(5) = NINTG
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS
! *        NINTG    -  NUMBER OF NUMERICAL INTEGRATION POINTS
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle
    REAL(8) :: THIC
    INTEGER :: N,II,JJ,LL,MM
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(5) :: KWQUADS
!- - - - - - - - - - - - - - -  
!-     VARIABLE DEFINITION   -
!- - - - - - - - - - - - - - -
    KWQUADS(1) = 'ITYPE:'
    KWQUADS(2) = 'MTYPE:'
    KWQUADS(3) = 'NUMELE:'
    KWQUADS(4) = 'NINTG:'
    KWQUADS(5) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=0
    NPAR(3)=1
    NPAR(4)=0
    NPAR(5)=2
!  
    I = KeyWord(KWQUADS,5)
    DO WHILE(I<=5)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS,5)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
    NumEle=NPAR(4)
!
!-------- Allocate memory
!
        ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
        ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,8
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, JJ, LL, MM, THIC
        DO WHILE(N.NE.0)
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR N .GT. NumEle"
              STOP
           END IF
           NCOUNT=NCOUNT+1
!
           XYZ(1,N)=X(II)
           XYZ(2,N)=Y(II)
!
           XYZ(3,N)=X(JJ)
           XYZ(4,N)=Y(JJ)
!
           XYZ(5,N)=X(LL)
           XYZ(6,N)=Y(LL)
!
           XYZ(7,N)=X(MM)
           XYZ(8,N)=Y(MM)
!
           ELEPAR(1,N)=THIC
           DO I=1,2
             LM(I  ,N)=ID(I,II)
             LM(I+2,N)=ID(I,JJ)
             LM(I+4,N)=ID(I,LL)
             LM(I+6,N)=ID(I,MM)
           END DO
           CALL COLHT(8,LM(:,N))
           WRITE(IOUT,2050) N,II,JJ,LL,MM,THIC
!
           READ(IIN,*) N,II,JJ,LL,MM,THIC
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: QUADS          =',I5,/,&
              '               0: AXISYMMETRIC            ',/,&
              '               1: PLANE STRAIN            ',/,&
              '               2: PLANE STRESS            ',/,&                  
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF QUADS ELEMENTS      =',I5,/,&
              '    GAUSS NUMERICAL INTEGRATION ORDER =',I5,/)
 2000 FORMAT(/,&
              ' ELEMENT     NODE     NODE     NODE     NODE     THIC   ',/,&
              ' NUMBER-N      I        J        L        M       ',/)
 2050 FORMAT(I5,6X,I5,4X,I5,4X,I5,4X,I5,2X,ES10.3)
!}}}
    END SUBROUTINE QUADSD
!
    SUBROUTINE TRUSS_MFSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE TRUSS_MFSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS_MFS
! *         NPAR(1) = 3
! *         NPAR(2) = MTYPE
! *         NPAR(3) = NUMELE
! *         NPAR(4) = ELEASS
! *      INPUT:
! *        MTYPE    -  MATERIAL SET NUMBER
! *    	   NUMELE   -  THE NUMBER OF TRUSS ELEMENTS  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumEle
    INTEGER :: N,M,II,KN,IND
    INTEGER :: KNOLD,NOLD,NUM,NUMN,IND_OLD
    REAL(8) :: DIRN,DIRN_OLD
    INTEGER :: IOLD,NII
    REAL(8) :: RAD,ROLD
    INTEGER :: NCOUNT
    INTEGER :: NEleAss
    INTEGER :: INDSP,INDSP_OLD
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(6) :: KWTRUSS_MFS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWTRUSS_MFS(1) = 'MTYPE:'
    KWTRUSS_MFS(2) = 'NUMELE:'
    KWTRUSS_MFS(3) = 'ELEASS:'
    KWTRUSS_MFS(4) = 'NITP:'
    KWTRUSS_MFS(5) = 'NIDIV:'
    KWTRUSS_MFS(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=0
    NPAR(4)=NUMNP
!   
    I = KeyWord(KWTRUSS_MFS,6)
    DO WHILE(I<=6)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NITP=J
         CASE(5)
           J = GetInt()
           NIDIV=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN TRUSS_MFSD"
           STOP
       END SELECT
       I = KeyWord(KWTRUSS_MFS,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NITP,NIDIV
!
    IF(NPAR(3).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF TRUSS ELEMENT .EQ.0"
            STOP
    END IF
!
    NumEle  = NPAR(3)
    NEleAss = NPAR(4)
!
!-------- Allocate memory
!
        ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
        ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
        ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
! ----------------------------------------------------------------
        DO I=1,NumEle
           DO J=1,3
              NDMFS(J,I)=0
           ENDDO
           ELEPAR(1,I)=0.
           ELEPAR(2,I)=0.
!
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N, II, IND,RAD,DIRN,INDSP,KN
        DO WHILE(N.NE.0)
!           IF(II.GT.NumEle) THEN
!              WRITE(*,*) "ERROR N .GT. NumEle"
!              STOP
!           END IF
           NCOUNT=NCOUNT+1
!          
           ELEPAR(2,N)=RAD
           ELEPAR(1,N)=DIRN
           NDMFS(1,N)=II
           NDMFS(2,N)=IND
           NDMFS(3,N)=INDSP
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+KNOLD
                    NDMFS(1,K)=NII
                    NDMFS(2,K)=IND_OLD
                    NDMFS(3,K)=INDSP_OLD                 
!
                    ELEPAR(2,K)=ROLD
                    ELEPAR(1,K)=DIRN_OLD
!
                   WRITE(IOUT,2050) K,NII,IND_OLD,ROLD,DIRN_OLD,INDSP_OLD              
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
            ENDIF
            WRITE(IOUT,2050) N,II,IND,RAD,DIRN,INDSP
!        
            NOLD=N
            KNOLD=KN
            IOLD=II
            IND_OLD=IND
            ROLD=RAD
            DIRN_OLD=DIRN
            INDSP_OLD=INDSP
!
            READ(IIN,*) N, II,IND,RAD,DIRN,INDSP,KN
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT.NumEle"
              STOP
        END IF
        CALL MFS_TRUSS_ASSD(3,2,NumEle,NEleAss,NDMFS(1,1),ELEPAR(1,1),LM)
!       CONSTRUCT MHT
        DO N=1,NumEle
           DO I=1,NEleAss
             IF(LM(I,N).NE.0) THEN
              M=LM(I,N)
              CALL COLHT_MFS(NDMFS(1,N),NDMFS(1,M))
             ENDIF
           ENDDO
        ENDDO
!
    RETURN
!--------------------------------------------------
 1000 Format(/,&
              '    THE ELELMENT TYPE: TRUSS',/,&
              '    THE MATERIAL SET NUMBER                   =',I5,/,&
              '    THE NUMBER OF TRUSS ELEMENTS              =',I5,/,&
              '    THE MAXIMUM ELEMENT ASSOCIATION NUMBER    =',I5,/,&
              '    THE NUMBER GAUSS-LEGENDRE INTEGRATION POINTS     =',I5,/,&
              '    THE NUMBER OF INTEGRATION SECTIONS        =',I5)
 2000 FORMAT(/,&
              ' ELEMENT     NODE        IND       RAD        DIRN     INDSP',/,&
              ' NUMBER-N   NUMBER               ',/)
 2050 FORMAT(I5,6X,I5,6X,I5,2X,ES13.3,2X,ES10.3,2X,I5)
!}}}
    END SUBROUTINE TRUSS_MFSD
!
    SUBROUTINE COUPE1D(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE COUPE1D
! *      READ AND STORE TRUSS DATA
! *      TRUSS_MFS
! *         NPAR(1) = 4
! *      NPAR(2) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      NPAR(3) = MTYPE
! *      NPAR(4) = ELEASS
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle
    INTEGER :: NCOUNT 
    INTEGER :: N,I1,I2,I3,I4,IND1,IND2,IND3,IND4
    INTEGER :: N1,N2
    REAL(8) :: R1,R2,R3,R4
    INTEGER :: ALLocateStatus 
!
    CHARACTER(20),DIMENSION(6) :: KWQUADS_COUP
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_COUP(1) = 'ITYPE:'
    KWQUADS_COUP(2) = 'MTYPE:'
    KWQUADS_COUP(3) = 'NUMELE:'
    KWQUADS_COUP(4) = 'NITP_IN:'
    KWQUADS_COUP(5) = 'NIDIV_IN:'
    KWQUADS_COUP(6) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=1
    NPAR(4)=0
    NPAR(5)=NITP_IN
    NPAR(6)=NIDIV_IN
!
    I = KeyWord(KWQUADS_COUP,6)
    DO WHILE(I<=6)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. 3) THEN
               WRITE(*,*) "ERROR ITYPE .GT. 3"
               STOP
           ENDIF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN COUP1"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_COUP,6)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(5),NPAR(6)
!
    IF(NPAR(4).EQ.0) THEN
        WRITE(*,*) "ERROR, NUMBER OF COUP1 ELEMENT .EQ.0"
        STOP
    END IF
!
    NumEle  = NPAR(4)
!
!-------- Allocate memory
!
        ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
        ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
!
        DO I=1,NumEle
           DO J=1,8
              NDMFS(J,I)=0
           ENDDO
           DO J=1,4
              ELEPAR(J,I)=0.
           ENDDO
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N,I1,I2,I3,I4,IND1,IND2,IND3,IND4,R1,R2,R3,R4
        DO WHILE(N.NE.0)
           NCOUNT=NCOUNT+1
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR, COUPLING ELEMENT NUMBER .GT. NumEle"
              STOP
           END IF
!
           WRITE(IOUT,2050) N,I1,I2,I3,I4,IND1,IND2,IND3,IND4,R1,R2,R3,R4
!          
           ELEPAR(1,N)=R1
           ELEPAR(2,N)=R2
           ELEPAR(3,N)=R3
           ELEPAR(4,N)=R4
!
           NDMFS(1,N)=I1
           NDMFS(2,N)=I2
           NDMFS(3,N)=I3
           NDMFS(4,N)=I4
!
           NDMFS(5,N)=IND1
           NDMFS(6,N)=IND2
           NDMFS(7,N)=IND3
           NDMFS(8,N)=IND4
!
           DO I=1,4
              N1=NDMFS(I,N)
              DO J=I,4
                 N2=NDMFS(J,N)
                 CALL COLHT_MFS(N1,N2)
              ENDDO
           ENDDO
!
           READ(IIN,*) N,I1,I2,I3,I4,IND1,IND2,IND3,IND4,R1,R2,R3,R4
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
           WRITE(*,*) "ERROR, COUPLING ELEMENT NUMBER .LT. NumEle", NCOUNT
           STOP
        END IF
!
        RETURN
!
 1000 Format(/,&
              '    THE ELELMENT TYPE: QUADS          =',I5,/,&
              '               1: PLANE STRAIN            ',/,&
              '               2: PLANE STRESS            ',/,&                  
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF COUPLING ELEMENTS   =',I5,/,&
              '    THE NUMBER OF INTEGRATION POINTS  =',I5,/,&
              '    THE NUMBER OF INTEGRATION SECTION =',I5)
 2000 FORMAT(/,&
              ' ELEMENT       NODE(4,N)              IND(4,N)                              RAD(4,N) ',/,&
              ' NUMBER-N      NUMBER                 ',/)
 2050 FORMAT(I5,2X,4I5,2X,4I5,2X,4ES13.3)
!!}}}
    END SUBROUTINE COUPE1D
!
    SUBROUTINE QUADS_MFSD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE QUADS_MFSD
! *      READ AND STORE TRUSS DATA
! *      TRUSS_MFS
! *         NPAR(1) = 4
! *      NPAR(2) = ELEMENT TYPE
! *                1: PLANE STRAIN
! *                2: PLANE STRESS
! *      NPAR(3) = MTYPE
! *      NPAR(4) = NUMELE
! *      NPAR(5) = ELEASS
! *      INPUT:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    INTEGER :: NumEle,NEleAss
    INTEGER :: NCOUNT 
    INTEGER :: N,II,KN
    REAL(8) :: R_IN1,R_IN2,R_OUT1,R_OUT2
    INTEGER :: NOLD,IOLD,KNOLD
    INTEGER :: INDSP,INDSP_OLD
    INTEGER :: NUM,NUMN,K,NII
    INTEGER :: IND1,IND2,IND1O,IND2O
    REAL(8) :: R_IN_OLD1,R_IN_OLD2,R_OUT_OLD1,R_OUT_OLD2
    INTEGER :: ALLocateStatus 
    INTEGER :: M
    INTEGER :: DII
!
    CHARACTER(20),DIMENSION(7) :: KWQUADS_MFS
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWQUADS_MFS(1) = 'ITYPE:'
    KWQUADS_MFS(2) = 'MTYPE:'
    KWQUADS_MFS(3) = 'NUMELE:'
    KWQUADS_MFS(4) = 'ELEASS:'
    KWQUADS_MFS(5) = 'NITP_IN:'
    KWQUADS_MFS(6) = 'NIDIV_IN:'
    KWQUADS_MFS(7) = '\END'
!   SET DEFAULT VALUE
    NPAR(2)=1
    NPAR(3)=1
    NPAR(4)=1
    NPAR(5)=NUMNP
    NPAR(6)=NITP_IN
    NPAR(7)=NIDIV_IN
!
   I = KeyWord(KWQUADS_MFS,7)
    DO WHILE(I<=7)
       SELECT CASE(I)
         CASE(1)
           J = GetInt()
           IF (J .GT. 3) THEN
               WRITE(*,*) "ERROR ITYPE .GT. 3"
               STOP
           ENDIF
           NPAR(2)=J
         CASE(2)
           J = GetInt()
           IF (J .GT. NMAT) THEN
               WRITE(*,*) "ERROR MATERIAL SET NUMBER .GT. NMAT"
               STOP
           END IF
           NPAR(3)=J
         CASE(3)
           J = GetInt()
           NPAR(4)=J
         CASE(4)
           J = GetInt()
           NPAR(5)=J
         CASE(5)
           J = GetInt()
           NPAR(6)=J
         CASE(6)
           J = GetInt()
           NPAR(7)=J
         CASE(7)
           EXIT
         CASE DEFAULT
           WRITE(*,*) "ERROR, ERROR INPUT IN QUADS_MFSD"
           STOP
       END SELECT
       I = KeyWord(KWQUADS_MFS,7)
    END DO
!
    WRITE(IOUT,1000) NPAR(2),NPAR(3),NPAR(4),NPAR(6),NPAR(7),NPAR(5)
!
    IF(NPAR(4).EQ.0) THEN
            WRITE(*,*) "ERROR, NUMBER OF QUADS ELEMENT .EQ.0"
            STOP
    END IF
!
    NumEle  = NPAR(4)
    NEleAss = NPAR(5)
!
!-------- Allocate memory
!
        ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
        ALLOCATE (ELEPAR(6,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"        
        ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
        IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
        DO I=1,NumEle
           DO J=1,8
              NDMFS(J,I)=0
           ENDDO
!
           DO J=1,6
              ELEPAR(J,I)=0.
           ENDDO
!
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
! ----------------------------------------------------------------
        WRITE(IOUT,2000)
        NCOUNT=0
        READ(IIN,*) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP,KN
        DO WHILE(N.NE.0)
           NCOUNT=NCOUNT+1
           IF(N.GT.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .GT. NumEle"
              STOP
           END IF
!          
           ELEPAR(1,N)=R_IN1
           ELEPAR(2,N)=R_IN2
           ELEPAR(3,N)=R_OUT1
           ELEPAR(4,N)=R_OUT2
           NDMFS(1,N)=II
           NDMFS(2,N)=IND1
           NDMFS(3,N)=IND2
           NDMFS(6,N)=INDSP
!
           IF (KNOLD.NE.0) THEN
              NUM=(N-NOLD)/KNOLD
              NUMN=NUM - 1
              IF (NUMN.GE.1) THEN
                  K=NOLD
                  NII=IOLD
                  DII=(II-IOLD)/NUM
                 DO J=1,NUMN
                    K=K + KNOLD
!
                    NII=NII+DII
                    NDMFS(1,K)=NII
                    NDMFS(2,K)=IND1O
                    NDMFS(3,K)=IND2O 
                    NDMFS(6,K)=INDSP_OLD                   
!
                    ELEPAR(1,K)=R_IN_OLD1
                    ELEPAR(2,K)=R_IN_OLD2
                    ELEPAR(3,K)=R_OUT_OLD1
                    ELEPAR(4,K)=R_OUT_OLD2
!
                   WRITE(IOUT,2050) K,NII,IND1O,IND2O,R_IN_OLD1,R_IN_OLD2,R_OUT_OLD1,R_OUT_OLD2,INDSP_OLD             
                   NCOUNT=NCOUNT+1
                 END DO
              END IF
           ENDIF
           WRITE(IOUT,2050) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP
!        
           NOLD=N
           KNOLD=KN
           IOLD=II
           R_IN_OLD1=R_IN1
           R_IN_OLD2=R_IN2
           R_OUT_OLD1=R_OUT1
           R_OUT_OLD2=R_OUT2
           IND1O=IND1
           IND2O=IND2
           INDSP_OLD=INDSP
!
           READ(IIN,*) N,II,IND1,IND2,R_IN1,R_IN2,R_OUT1,R_OUT2,INDSP,KN
        END DO
!
        IF(NCOUNT.NE.NumEle) THEN
              WRITE(*,*) "ERROR, ELEMENT NUMBER .LT. NumEle", NCOUNT
              STOP
        END IF
!
        CALL MFS_2D_ASSD(8,6,NumEle,NEleAss,NDMFS(1,1),ELEPAR(1,1),LM(1,1))
!       CONSTRUCT MHT
        DO N=1,NumEle
           DO I=1,NEleAss
             M=LM(I,N)
             IF((M.NE.0).AND.(M.GE.N)) THEN
              CALL COLHT_MFS(NDMFS(1,N),NDMFS(1,M))
             ENDIF
           ENDDO
        ENDDO
!
        RETURN
!
 1000 Format(/,&
              '    THE ELELMENT TYPE: QUADS          =',I5,/,&
              '               1: PLANE STRAIN            ',/,&
              '               2: PLANE STRESS            ',/,&                  
              '    THE MATERIAL SET NUMBER           =',I5,/,&
              '    THE NUMBER OF QUADS ELEMENTS      =',I5,/,&
              '    THE NUMBER OF INTEGRATION POINTS  =',I5,/,&
              '    THE NUMBER OF INTEGRATION SECTION =',I5,/,&
              '    THE MAXIMUM ELEMENT ASSOCIATION NUMBER    =',I5)
 2000 FORMAT(/,&
              ' ELEMENT     NODE      IND(2,N)                            RAD(4,N)                       INDSP ',/,&
              ' NUMBER-N   NUMBER                 ',/)
 2050 FORMAT(I5,6X,I5,6X,2I5,2X,4ES13.3,2X,I5)
!}}}
    END SUBROUTINE QUADS_MFSD
!
    SUBROUTINE MFS_TRUSS_ASSD(NEle1,NEle2,NumEle,NEleAss,MFS_NODE,MFS_PAR,MFS_ASS)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM      
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      MFS_NODE :  
! *      MFS_PAR  :
! *      MFS_ASS  : 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER,INTENT(IN) :: NumEle, NEleAss, NEle1,NEle2
    INTEGER,DIMENSION(NEle1,NumEle)   :: MFS_NODE
    REAL(8),DIMENSION(NEle2,NumEle)   :: MFS_PAR
    INTEGER,DIMENSION(NEleAss,NumEle) :: MFS_ASS
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER,DIMENSION(NumEle) :: K
    INTEGER :: N,M
    INTEGER :: KK
    REAL(8) :: RAD1,RAD2
    REAL(8) :: DIST  
    LOGICAL :: EXPRN  
!
    KK=NEleAss+1
!
    DO I=1,NumEle
       K(I)=1
    ENDDO
    RAD1=0
    RAD2=0
    DO I=1,NumEle
       N=MFS_NODE(1,I)
       RAD1=MFS_PAR(NEle2,I) 
       DO J=I,NumEle              
          M=MFS_NODE(1,J)
          RAD2=MFS_PAR(NEle2,J)
          DIST=(RAD1+RAD2)
!
          DIST=ABS(X(N)-X(M))-DIST
!
          EXPRN=DIST>=(-1E-8) 
          IF(EXPRN) CYCLE 
!          
          MFS_ASS(K(I),I)=J
          K(I)=K(I)+1
          IF(I.NE.J) THEN
             MFS_ASS(K(J),J)=I
             K(J)=K(J)+1
          ENDIF
          IF(K(I).GT.KK) THEN
             WRITE(*,*) "ERROR, NEleAss too small!"
             STOP
          ENDIF
       ENDDO
    ENDDO
!
    DO I=1,NumEle
       IF(K(I).LE.2) THEN
          WRITE(*,*) "ERROR, ELEMENT DOES NOT OVERLAP!"
          STOP
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE MFS_TRUSS_ASSD
!
    SUBROUTINE MFS_2D_ASSD(NEle1,NEle2,NumEle,NEleAss,MFS_NODE,MFS_PAR,MFS_ASS)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE FOR CALCULATING NODE ASSOCIATION LM      
! *      NumEle   :  Number of elements
! *      EleAss   :  Maximum element association number
! *      MFS_NODE :  
! *      MFS_PAR  :
! *      MFS_ASS  : 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
!
    INTEGER,INTENT(IN) :: NumEle, NEleAss, NEle1,NEle2
    INTEGER,DIMENSION(NEle1,NumEle)   :: MFS_NODE
    REAL(8),DIMENSION(NEle2,NumEle)   :: MFS_PAR
    INTEGER,DIMENSION(NEleAss,NumEle) :: MFS_ASS
!
!   ASSISTANT VARIABLES
!
    INTEGER :: I,J
    INTEGER,DIMENSION(NumEle) :: K
    INTEGER :: N,M
    INTEGER :: KK
    REAL(8) :: RAD_OUT1,RAD_OUT2,RAD2_OUT1,RAD2_OUT2
    REAL(8) :: DIST1,DIST2
    LOGICAL :: EXPRN1,EXPRN2
!
    KK=NEleAss+1
!
    DO I=1,NumEle
       K(I)=1
    ENDDO
!
    DO I=1,NumEle
       RAD_OUT1=MFS_PAR(3,I)
       RAD_OUT2=MFS_PAR(4,I)
       N=MFS_NODE(1,I)
       DO J=I,NumEle              
          M=MFS_NODE(1,J)
!
          RAD2_OUT1=MFS_PAR(3,J)
          RAD2_OUT2=MFS_PAR(4,J)
!
          DIST1=RAD_OUT1+RAD2_OUT1
          DIST1=ABS(X(N)-X(M))-DIST1
          EXPRN1=DIST1.GT.(-1E-12)
          IF(EXPRN1) CYCLE
          DIST2=RAD_OUT2+RAD2_OUT2
          DIST2=ABS(Y(N)-Y(M))-DIST2
          EXPRN2=DIST2.GT.(-1E-12)
          IF(EXPRN2) CYCLE
!          
          MFS_ASS(K(I),I)=J
          K(I)=K(I)+1
!
          IF(I.NE.J) THEN
             MFS_ASS(K(J),J)=I
             K(J)=K(J)+1
          ENDIF
          IF(K(I).GT.KK) THEN
             WRITE(*,*) "ELleAss should .GE.",K(I)
             WRITE(*,*) "N AND M ARE:"
             WRITE(*,*) N,M
             WRITE(*,*) "X(N),Y(M),RADI1,RADI2=", X(N),Y(N),RAD_OUT1,RAD_OUT2
             WRITE(*,*) "X(N),Y(M),RADJ1,RADJ2=", X(M),Y(M),RAD2_OUT1,RAD2_OUT2
             write(*,*) "LM is:"
             WRITE(*,*) MFS_NODE(1,MFS_ASS(:,I))
             WRITE(*,*) "ERROR, NEleAss too small!"
             STOP
          ENDIF
       ENDDO
    ENDDO
!
    DO I=1,NumEle
       IF(K(I).LE.4) THEN
          WRITE(*,*) "ERROR, ELEMENT DOES NOT OVERLAP!"
          STOP
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE MFS_2D_ASSD 
!**************************************************
!
!   SUBROUTINES FOR SURFACE INPUT
!
!**************************************************
    SUBROUTINE SURFACE_IN(IIN,IOUT,ITEM2)
!{{{
!--------------------------------------------------
!   ALLOCATE SURFACE MEMORY AND SELECT SURFACE TYPE
!   NPAR()
!   POINT-POINT
!       NPAR(1) = 1
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,ITEM2
!
    INTEGER :: I,J,L,M
!
    INTEGER,PARAMETER :: KwSurf=1
    CHARACTER(20),DIMENSION(KwSurf) :: KwSurface
!
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KwSurface(1)='\Pt-Pt'
!
    WRITE(IOUT,2000)
!
    IF (NSPAIR<1E-7) THEN
            WRITE(*,*) " ERROR, NSPAIR EQUALS TO ZEROS"
            STOP
    ENDIF
!
    NCQ=0
    DO J=1,NSPAIR
       I = KeyWord(KwSurface,KwSurf)
       WRITE(IOUT,1000) KwSurface(I)
       SELECT CASE(I)
       CASE(1)
          NSFPAR=0
          NSFPAR(1)=1
          CALL PTPTD(IIN,IOUT)
          WRITE(ITEM2) (NSFPAR(L),L=1,7)
          WRITE(ITEM2) ((XYZC(L,M),L=1,3),M=1,NSFPAR(2)),((NSFC(L,M),L=1,1),M=1,NSFPAR(2)),((XYZT(L,M),L=1,3),M=1,NSFPAR(3)),&
          ((NSFT(L,M),L=1,1),M=1,NSFPAR(3))
          NCQ=NCQ+NSFPAR(2)
          ! Deallocate variables
          DEALLOCATE (NSFC) 
          DEALLOCATE (NSFT)
          DEALLOCATE (XYZC)
          DEALLOCATE (XYZT)      
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
    END DO
    RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A20)
 2000   FORMAT(//,' S U R F A C E   I N F O R M A T I O N',///)

!}}}
    END SUBROUTINE SURFACE_IN
!
    SUBROUTINE PTPTD(IIN,IOUT)
!{{{
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *      SUBROUTINE PTPTD
! *      READ AND STORE POINT-POINT SURFACE PAIR
! *      PTPTD
! *         NSFPAR(1) =   1
! *         NSFPAR(2) --- POINT NUMBER IN CONTRACTOR
! *         NSFPAR(3) --- POINT NUMBER IN TARGET 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT
!
    INTEGER :: I,J,K
!
    INTEGER :: N,NODE
    INTEGER :: NCOUNT
    INTEGER :: AllocateStatus
!
    CHARACTER(20),DIMENSION(3) :: KWPTPT
! - - - - - - - - - - - - - - -  
! -     VARIABLE DEFINITION   -
! - - - - - - - - - - - - - - -
    KWPTPT(1) = 'Contactor:'
    KWPTPT(2) = 'Target:'
    KWPTPT(3) = '\END'
!
    DO I=1,3
       J = KeyWord(KWPTPT,3)
       IF(J.LE.2) THEN
         K = GetInt()
         WRITE(IOUT,2000) KWPTPT(J),K
       ENDIF
         
       SELECT CASE(J)
       CASE(1)
          NSFPAR(2)=K
!
!-------- Allocate memory
!
          ALLOCATE (NSFC(1,NSFPAR(2)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
          ALLOCATE (XYZC(3,NSFPAR(2)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
          WRITE(IOUT,2050)
          N=0
          NODE=0
          NCOUNT=0
          READ(IIN,*) N,NODE
          DO WHILE(N.NE.0) 
             IF(NODE.GT.NUMNP) THEN
                    WRITE(*,*) "ERROR, NODE .GT. NUMNP IN SURFACE INPUT DATA"
                    STOP
             ENDIF
             NCOUNT=NCOUNT+1
             NSFC(1,N)=NODE
             XYZC(1,N)=X(NODE)
             XYZC(2,N)=Y(NODE)
             XYZC(3,N)=Z(NODE)
!
             WRITE(IOUT,3000) N,NODE
             READ(IIN,*) N,NODE

          END DO
          IF(NCOUNT.NE.NSFPAR(2)) THEN
                  WRITE(*,*) "ERROR, NCOUNT.NE.NSFPAR(2) IN SURFACE INPUT DATA"
                  STOP
          END IF
       CASE(2)
          NSFPAR(3)=K
!
!-------- Allocate memory
!
          ALLOCATE (NSFT(1,NSFPAR(3)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
          ALLOCATE (XYZT(3,NSFPAR(3)),STAT = ALLocateStatus)
          IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
          WRITE(IOUT,2050)
          N=0
          NODE=0
          NCOUNT=0
          READ(IIN,*) N,NODE
          DO WHILE(N.NE.0) 
             IF(NODE.GT.NUMNP) THEN
                    WRITE(*,*) "ERROR, NODE .GT. NUMNP IN SURFACE INPUT DATA"
                    STOP
             ENDIF
             NCOUNT=NCOUNT+1
             NSFT(1,N)=NODE
             XYZT(1,N)=X(NODE)
             XYZT(2,N)=Y(NODE)
             XYZT(3,N)=Z(NODE)
!
             WRITE(IOUT,3000) N,NODE
             READ(IIN,*) N,NODE
          END DO
          IF(NCOUNT.NE.NSFPAR(3)) THEN
                  WRITE(*,*) "ERROR, NCOUNT.NE.NSFPAR(3) IN SURFACE INPUT DATA"
                  STOP
          END IF
       CASE(3)
          EXIT
       CASE DEFAULT
          WRITE(*,*) "ERROR, INPUT ERROR IN ELEMENTD"
          STOP
       END SELECT 
!   
    END DO
    RETURN    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   FORMAT(/,A20)
 2000   FORMAT(/,A20,2X,I5)
 2050   FORMAT(/,&
              ' SURFACE     NODE ',/,&
              ' NUMBER-N      I  ',/)
 3000   FORMAT(I5,6X,I5)
!}}}
    END SUBROUTINE PTPTD
end module INPUTM
