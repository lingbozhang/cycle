! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Storage module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
module WAREHOUSE
    implicit none
    PUBLIC :: ALLOCATEM
    PUBLIC :: RELEASE
!{{{
!        CONSTANTS
         INTEGER :: NFILE
         CHARACTER(LEN=3) :: FILE_NAME
         REAL(8), PARAMETER :: Pi=3.1415926535898
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE GLOBALD
!--------------------------------------------------
         CHARACTER(256) :: HEAD='HEAD: TEST'
         INTEGER :: NUMNP   = 0
         INTEGER :: NLOAD   = 0
         INTEGER :: NBOUND  = 0
         INTEGER :: NMAT    = 0  
         INTEGER :: NBGROUP = 0       
         INTEGER :: NEGROUP = 0
         INTEGER :: NSPAIR  = 0
         INTEGER :: DYNAMIC = 0
         INTEGER :: CONTACT = 0
         INTEGER :: NPOR_MAX= 0
         INTEGER :: PENALTY = 0
! NUMBER OF INTEGRATION POINTS                            
         INTEGER :: NITP_BC=40
         INTEGER :: NIDIV_BC=1
         INTEGER :: NITP_IN=60
         INTEGER :: NIDIV_IN=4
!
         INTEGER :: NEQ     = 0
         INTEGER :: NCQ     = 0
!-------------------------------
!    PARAMETERS
!-------------------------------
!    FOR PENALTY METHOD
!-------------------------------
         REAL(8),PARAMETER :: BETA=1E6
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE NODED (FEM)
!--------------------------------------------------
         INTEGER,DIMENSION(:,:),  ALLOCATABLE :: ID
         REAL(8),DIMENSION(:),    ALLOCATABLE :: X
         REAL(8),DIMENSION(:),    ALLOCATABLE :: Y
         REAL(8),DIMENSION(:),    ALLOCATABLE :: Z
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE MATERIALD 
!--------------------------------------------------
!   CONSTITUTIVE MATRIX 
         REAL(8),DIMENSION(:),  ALLOCATABLE :: KCS
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: K3CS
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE NODED (MFS)
!        THE DEGREE OF FREEDOM
!--------------------------------------------------       
         INTEGER,DIMENSION(:,:),    ALLOCATABLE :: NDOF       
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE LOADD
!--------------------------------------------------
         REAL(8),DIMENSION(:), ALLOCATABLE :: R
         INTEGER,DIMENSION(:), ALLOCATABLE :: NRF
!
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: RNAT
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: RESS
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE BOUNDD
!--------------------------------------------------
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE INITLD
!--------------------------------------------------
         INTEGER :: IND_INITIAL=0
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: UI
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: VI
         REAL(8),DIMENSION(:),   ALLOCATABLE :: U
         REAL(8),DIMENSION(:),   ALLOCATABLE :: V
         REAL(8),DIMENSION(:),   ALLOCATABLE :: ACC
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE MATERIALD
!--------------------------------------------------
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: MATP
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE BOUNDARY
!--------------------------------------------------
         INTEGER,DIMENSION(:,:), ALLOCATABLE :: BCED
         INTEGER,DIMENSION(:,:), ALLOCATABLE :: BCIND
         REAL(8),DIMENSION(:,:), ALLOCATABLE :: BCPAR        
!**************************************************
!   WAREHOUSE FOR SUBROUTINE ELEMENT_IN
!**************************************************
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
         INTEGER,DIMENSION(12) :: NPAR
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE ELEMENTS (FEM)
!--------------------------------------------------
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: XYZ
         INTEGER,DIMENSION(:,:),   ALLOCATABLE :: LM
         INTEGER,DIMENSION(:),     ALLOCATABLE :: NumLM
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: ELEPAR     
!-------------------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE ELEMENTS (MFS)
!        STORING NODE NUMBER OF THE METHODS OF FINITE SPHERES
!        MFS_ITP :  INTEGRATION POINT
!        NITP    :  NUMBER OF INTEGRATION POINTS
!-------------------------------------------------------------
         INTEGER  :: NITP=20
         INTEGER  :: NIDIV=2
         INTEGER,DIMENSION(:,:),   ALLOCATABLE :: NDMFS
! MID-POINT INTEGRATION METHODS FOR THE BOUNDARY
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: XABSC_WEIG_BC
! Gauss-Legendre METHOD FOR THE INNER SPHERE
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: XABSC_WEIG_IN
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: ACOOR_IN
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: STIFF_IJ
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: OFE1_QUADS_S
!--------------------------------------------------
!   WAREHOUSE FOR NEW OVERLAPPING ELEMENTS
!--------------------------------------------------
         INTEGER,DIMENSION(:,:),   ALLOCATABLE :: NDINT
         INTEGER,DIMENSION(:,:),   ALLOCATABLE :: CPLINT
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: INTPAR
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: WEIND
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: DHI_OVLP_TRUSS
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: HI_OVLP_TRUSS
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: H_OVLP_TRUSS
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: DH_OVLP_QUADS
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: H_OVLP_QUADS
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: Rmatrix_all
         REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: Qmatrix_all
         REAL(8),DIMENSION(:,:),   ALLOCATABLE :: Wmatrix_all
!        BC INDICATOR
         LOGICAL :: OFE1BCONOFF=.FALSE.
         LOGICAL :: QUADSBCONOFF=.FALSE.
         LOGICAL :: FEMCBCONOFF=.FALSE.
!**************************************************
!   WAREHOUSE FOR SUBROUTINE ELEMENT_IN
!**************************************************
!
!**************************************************
!   WAREHOUSE FOR SUBROUTINE SURFACE_IN
!**************************************************
!   POINT
!       NSFPAR(1) = 1
!       NSFPAR(2) = NUMSF
!--------------------------------------------------
        INTEGER,DIMENSION(7) :: NSFPAR
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE SURFACES
!   NSFC:     NODE NUMBER OF THE CONTACTOR
!   XYZC:     XYZ COORDINATES OF THE CONTACTORS
!   NSFT:     NODE NUMBER OF THE TARGET
!   XYZT:     XYZ COORDINATES OF THE TARGETS
!--------------------------------------------------
        INTEGER :: IND_SURFACE=0
        INTEGER,DIMENSION(:,:),   ALLOCATABLE :: NSFC
        REAL(8),DIMENSION(:,:),   ALLOCATABLE :: XYZC
        INTEGER,DIMENSION(:,:),   ALLOCATABLE :: NSFT
        REAL(8),DIMENSION(:,:),   ALLOCATABLE :: XYZT      
!**************************************************
!   WAREHOUSE FOR SUBROUTINE SURFACE_IN
!**************************************************
!
!**************************************************
!   WAREHOUSE FOR SUBROUTINE GENERATE -- CONTACT
!**************************************************
         INTEGER :: NWK_CONT=0
!**************************************************
!   WAREHOUSE FOR SUBROUTINE SOLVER
!**************************************************
!
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE SOLVER_COLSOL
!--------------------------------------------------
        INTEGER,DIMENSION(:), ALLOCATABLE :: MHT
        INTEGER,DIMENSION(:), ALLOCATABLE :: MHT_ND
!   MAXIMUM POSSIBLE HALF-BANDWIDTH IN CONTACT MATRIX
        INTEGER                           :: MHTC        
        INTEGER*8,DIMENSION(:), ALLOCATABLE :: MAXA_COLSOL
        INTEGER*8,DIMENSION(:), ALLOCATABLE :: MAXA_ND
        INTEGER*8 :: NWK_ND=0
        INTEGER*8 :: NWK_COLSOL=0
        REAL(8),DIMENSION(:), ALLOCATABLE :: K_COLSOL
        LOGICAL,DIMENSION(:), ALLOCATABLE :: BAND_FLAG
        REAL(8),DIMENSION(:), ALLOCATABLE :: KMATRIX
        REAL(8),DIMENSION(:), ALLOCATABLE :: M_COLSOL
        REAL(8),DIMENSION(:), ALLOCATABLE :: C_COLSOL
        INTEGER                           :: HALF_BAND
!
        REAL(8),DIMENSION(:), ALLOCATABLE :: DISP        
!--------------------------------------------------
!
!**************************************************
!   WAREHOUSE FOR SUBROUTINE ODE_SOLVER
!**************************************************
!       TPARA(1) -- TSTART
!       TPARA(2) -- TEND
!       TPARA(3) -- DTIME
!       INDO_ODE: INDICATOR FOR ODE SOLVER
!                 1  --  NEWMARK
!                 2  --  BATHE
!--------------------------------------------------
        INTEGER :: IND_ODE=0
        REAL(8),DIMENSION(3) :: TPARA=0
        REAL(8) :: TIME=0
        REAL(8) :: DTIME=0
        REAL(8),DIMENSION(:), ALLOCATABLE :: KODE
        REAL(8),DIMENSION(:), ALLOCATABLE :: RODE
        REAL(8),DIMENSION(:), ALLOCATABLE :: UODE
        REAL(8),DIMENSION(:), ALLOCATABLE :: VODE
        REAL(8),DIMENSION(:), ALLOCATABLE :: AODE
!
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE OUTPUT
!--------------------------------------------------
        REAL(8) :: ESTRAIN      ! STRAIN ENERGY
        REAL(8),DIMENSION(:), ALLOCATABLE :: HESTRAIN
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE UPLOAD
!--------------------------------------------------
!        REAL(8) :: K=-10000    ! TANGENT VALUE FOR LINEAR LOAD FUNCTION
!}}}
        CONTAINS
!
      SUBROUTINE ALLOCATEM()
!{{{
!--------------------------------------------------
!       ALLOCATE MEMORY
!--------------------------------------------------
        IMPLICIT NONE
        INTEGER :: I,NN
        INTEGER :: ALLocateStatus
!**************************************************
!   WAREHOUSE FOR SUBROUTINE GENERATE
!**************************************************
! Build MAXA, Kmatrix, Mmatrix, CMATRIX
!{{{
!--- ALLOCATE VARIABLES
!    ALLOCATE (MAXA(NEQ+NCQ+1),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!!--------------------------------------------------
!    DO I=1,NEQ+NCQ+1
!       MAXA(I)=0
!    END DO
!!--------------------------------------------------
!!   SET DEFAULT VALUE
!    MAXA(1)=1
!    MAXA(2)=2
!    IF ((NEQ+NCQ).GT.1) THEN
!       DO I=2,NEQ+NCQ
!          MAXA(I+1)=MAXA(I)+I
!       END DO
!    END IF
!!--- ALLOCATE VARIABLES
!    NWK=MAXA(NEQ+1)-MAXA(1)
!    NWK_CONT=MAXA(NEQ+NCQ+1)-MAXA(1)
!!
!    ALLOCATE (KMATRIX(NWK_CONT),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!    ALLOCATE (MMATRIX(NWK_CONT),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!    ALLOCATE (CMATRIX(NWK_CONT),STAT = ALLocateStatus)
!    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!!
!    DO I=1,NWK_CONT
!       KMatrix(I)=0
!       MMatrix(I)=0
!       CMatrix(I)=0
!    END DO    
!}}}
!**************************************************
!   WAREHOUSE FOR SUBROUTINE ODE_SOLVER
!**************************************************
!    IF(DYNAMIC.EQ.1) THEN
!      ALLOCATE (KODE(NWK_CONT),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!      ALLOCATE (RODE(NEQ+NCQ),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"         
!      ALLOCATE (UODE(NEQ+NCQ),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!      ALLOCATE (VODE(NEQ+NCQ),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!      ALLOCATE (AODE(NEQ+NCQ),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!!
!    DO I=1,NWK_CONT
!       KODE(I)=0
!    END DO
!    DO I=1,NEQ+NCQ
!       RODE(I) = 0
!       UODE(I) = 0 
!       VODE(I) = 0
!       AODE(I) = 0
!    END DO
!   ENDIF
!}}}
      END SUBROUTINE ALLOCATEM
!
      SUBROUTINE RELEASE()
!{{{
!--------------------------------------------------
!       RELEASE MEMORY
!--------------------------------------------------
        IMPLICIT NONE
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE NODED
!--------------------------------------------------
         DEALLOCATE (ID)
         DEALLOCATE (X)
         DEALLOCATE (Y)
         DEALLOCATE (Z)
         DEALLOCATE (NDOF)
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE LOADD
!--------------------------------------------------
         DEALLOCATE (R)
         DEALLOCATE (NRF)
         DEALLOCATE (RNAT)
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE BOUNDD
!--------------------------------------------------
         DEALLOCATE (RESS)
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE BOUNDARY
!--------------------------------------------------

!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE INITLD
!--------------------------------------------------
       IF(IND_INITIAL.EQ.1) THEN
           DEALLOCATE (UI)
           DEALLOCATE (VI)
           DEALLOCATE (U)
           DEALLOCATE (V)
           DEALLOCATE (ACC)
       END IF
!--------------------------------------------------
!   DEALLOCATE WAREHOUSE FOR SUBROUTINE SURFACE
!--------------------------------------------------
       IF(IND_SURFACE.EQ.1) THEN
          DEALLOCATE (NSFC)
          DEALLOCATE (XYZC)
          DEALLOCATE (NSFT)
          DEALLOCATE (XYZT)
       END IF
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE MATERIALD
!--------------------------------------------------
         DEALLOCATE (MATP) 
!**************************************************
!   WAREHOUSE FOR SUBROUTINE GENERATE
!**************************************************
!         DEALLOCATE(MAXA)
!         DEALLOCATE(KMATRIX)
!         DEALLOCATE(MMATRIX)
!         DEALLOCATE(CMATRIX)
!--------------------------------------------------
!   WAREHOUSE FOR SUBROUTINE SOLVER_COLSOL
!--------------------------------------------------
          DEALLOCATE(MHT)
          IF (ALLOCATED(MHT_ND)) DEALLOCATE(MHT_ND)
          IF (ALLOCATED(DISP)) DEALLOCATE(DISP)
          IF (ALLOCATED(MAXA_COLSOL)) DEALLOCATE(MAXA_COLSOL)
          IF (ALLOCATED(MAXA_ND)) DEALLOCATE(MAXA_ND)
          IF (ALLOCATED(K_COLSOL)) DEALLOCATE(K_COLSOL)
          IF (ALLOCATED(BAND_FLAG)) DEALLOCATE(BAND_FLAG)
          IF (ALLOCATED(KMATRIX)) DEALLOCATE(KMATRIX)
          IF (ALLOCATED(HESTRAIN)) DEALLOCATE(HESTRAIN)
!         DEALLOCATE(M_COLSOL)
!         DEALLOCATE(C_COLSOL)
!*************************************************
!   WAREHOUSE FOR SUBROUTINE ODE_SOLVER
!**************************************************
       IF(DYNAMIC.EQ.1) THEN
         DEALLOCATE(KODE)
         DEALLOCATE(RODE)
         DEALLOCATE(UODE)
         DEALLOCATE(VODE)
         DEALLOCATE(AODE)
       END IF
!
         RETURN         
!}}}
      END SUBROUTINE RELEASE
end module WAREHOUSE
