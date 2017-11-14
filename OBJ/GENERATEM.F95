! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Generation module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

module GENERATEM
    USE WAREHOUSE
    USE SOLVERM
    USE GAUSSM
    USE OPERATIONM
    USE SORTM
    USE PKEAST,ONLY : KEAST_ORDER_NUM,KEAST_RULE
    implicit none
    PUBLIC  GENERATE
!   GENERATE LOAD VECTOR(MFS)
    PRIVATE  MFS_LINE_G
    PRIVATE  MFS_LINE_STIFF_WEAK
    PRIVATE  MFS_LINE_STIFF_PENALTY
    PRIVATE  MFS_LINE_ADDBAN
!   GENERATE STIFF MATRIX(FEM)
    PRIVATE TRUSS_G
    PRIVATE QUADS_G
    PRIVATE QUADS_LINE_G
!   GENERATE STIFF MATRIX(COUP)
    PRIVATE COUP1_G
    PRIVATE COUP1_STIFF
    PRIVATE COUP1_GAUSSG
    PRIVATE COUP1_ADDBAN
!   GENERATE STIFF MATRIX(MFS)
    PRIVATE TRUSS_MFS_G
    PRIVATE TRUSS_MFS_STIFF
    PRIVATE TRUSS_MFS_BODY
!
    PRIVATE QUADS_MFS_G
    PRIVATE QUADS_MFS_STIFF_V2
    PRIVATE QUADS_MFS_GAUSSG_V2
    PRIVATE QUADS_MFS_ADDBAN_V2
!   GENERATE STIFF MATRIX (NOVLP)
    PRIVATE TRUSS_OVLP_G
    PRIVATE TRUSS_OVLP2_G
    PRIVATE TRUSS_OVLP_STIFF
        PRIVATE TRUSS_OVLP_BODY
        PRIVATE TRUSS_OVLP_STIFF_IJ
        PRIVATE TRUSS_OVLP_ADDBAN
        PRIVATE TRUSS_OVLP_DHI
        PRIVATE TRUSS_OVLP2_DHI
        PRIVATE WEI_CAL
        PRIVATE DHI_TRUSS
        PRIVATE DHI_TRUSS2
    PRIVATE TRUSS_OVLP3_G
        PRIVATE LeastSQ
        PRIVATE TRUSS_OVLP_STIFF3
          PRIVATE H_TRUSS
    PRIVATE QUADS_OVLP_G
        PRIVATE OVLP_LINE_G
        PRIVATE QUADS_OVLP_STIFF
        PRIVATE DH_QUADS
        PUBLIC  WEI_CAL_QUADS
        PRIVATE H_QUADS
        PRIVATE WEI_CAL_QUADS2
        PRIVATE FBODY_OVLP
!   GENERATE STIFF MATRIX (OFE)
    PRIVATE OFE_TRUSS_G
    PRIVATE OFE_TRUSS_STIFF
    PRIVATE OFE_TRUSS_ADDBAN
!   GENERATE STIFF MATRIX (OFE1)
    PRIVATE OFE1_TRUSS_G
    PRIVATE OFE1_TRUSS_STIFF
!   GENERATE STIFF MATRIX (OFE1_QUADS)
    PRIVATE OFE1_QUADS_G
    PRIVATE OFE1_QUADS_STIFF
    PRIVATE OFE1_COUPLE_STIFF
    PRIVATE OFE1_CPL2_STIFF
!
    PRIVATE OFE1_LINE_G
    PRIVATE OFE1B_LINE_G
    PRIVATE OFE1_QUADS_ADDBAN
    PRIVATE OFE1_COUPLE_ADDBAN
    PRIVATE OFE1_CPL2_ADDBAN
!   GENERATE STIFF MATRIX (FEMC_QUADS)
    PRIVATE FEMC_QUADS_G
    PRIVATE FEMC_QUADS_STIFF
    PRIVATE FEMC_LINE_G
!   GENERATE CONSTITUTIVE MATRIX
    PRIVATE ADDBAN
    PRIVATE ADDBAN_ZLB
    PRIVATE OUTPUTKM
!   HEX8
    PRIVATE HEX8_G
    PRIVATE HEX8D
!   OFE_3D
    PRIVATE OFE_TETRA_G
    PRIVATE OFE_SURFACE
    PRIVATE OFE_TETRA_STIFF
    PRIVATE OFE_TETRA_ADDBAN
    PRIVATE OLR_3D_ADDBAN
!
    CONTAINS
!
    SUBROUTINE GENERATE(ADK,ADM,ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE 
!     SELECT GENERATOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM,ADELE,ADBC   
!
    INTEGER :: I,J,BcType,EleType
!   ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
!-------------------------------------------------------------
    REWIND(ADELE)
    REWIND(ADBC)
!
    DO I=1,NEGROUP
       READ(ADELE) (NPAR(J),J=1,12)
       EleType=NPAR(1)
       SELECT CASE(EleType)
       CASE(1)
            CALL TRUSS_G(ADELE)
       CASE(2)
            WRITE(*,*) "Generating the stiff matrix and load vector/QUADS_G" 
            CALL QUADS_G(ADELE,ADBC)    
       CASE(3)
            CALL TRUSS_MFS_G(ADELE)
       CASE(4)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... QUAD_MFS" 
            CALL QUADS_MFS_G(ADELE,ADBC)
       CASE(5)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... COUP1_G" 
            CALL COUP1_G(ADELE)
       CASE(6)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... Truss_OVLP"
            CALL TRUSS_OVLP_G(ADELE)
       CASE(7)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... Truss_OVLP2"
            CALL TRUSS_OVLP2_G(ADELE)
       CASE(8)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... Truss_OVLP3"
            CALL TRUSS_OVLP3_G(ADELE)
       CASE(9)
           WRITE(*,*) "Generating the stiff matrix and load vector ......... QUADS_OVLP"
            CALL QUADS_OVLP_G(ADELE,ADBC)
       CASE(10)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... OFE_TRUSS"
            CALL OFE_TRUSS_G(ADELE)
       CASE(11)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... OFE1_TRUSS"
            CALL OFE1_TRUSS_G(ADELE)
       CASE(12)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... OFE1_QUADS"
            CALL OFE1_QUADS_G(ADELE,ADBC)
       CASE(13)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... FEMC_QUADS"
            CALL FEMC_QUADS_G(ADELE,ADBC) 
       CASE(14)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... OFE_3D"
            CALL OFE_TETRA_G(ADELE,ADBC) 
       CASE(16)
            WRITE(*,*) "Generating the stiff matrix and load vector ......... HEX8"
            CALL HEX8_G(ADELE,ADBC) 
       CASE DEFAULT
            WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE ELEMENT PROCEDURE"
            STOP
       END SELECT       
    END DO
!    CALL OUTPUTKM(ADK,ADM)
!    stop
    RETURN
!}}}
    END SUBROUTINE GENERATE
!
    SUBROUTINE HEX8_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   VARIABLES FOR K & M
    INTEGER :: N
    REAL(8) :: S(24,24)
    INTEGER :: NINTP
    INTEGER :: NN
    REAL(8),DIMENSION(:,:,:),   ALLOCATABLE :: NI
    REAL(8),DIMENSION(:),ALLOCATABLE :: WEI
!   ASSISTANT VARIABLE FOR K*U
    INTEGER :: I,J,K
    INTEGER :: Mtype,NumEle,NINTG  
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: ALLocateStatus
!   READ 
    Mtype  = NPAR(2)
    NumEle = NPAR(3)
    NINTG  = NPAR(4)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (XYZ(24,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(24,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!--------------------------------------------------
    XYZ=0
    LM=0
    NINTP=NINTG*NINTG*NINTG
!--------------------------------------------------
    READ(ADELE) ((XYZ(I,J),I=1,24),J=1,NumEle),((LM(I,J),I=1,24),J=1,NumEle)
!   GENERATE
    ALLOCATE (XABSC_WEIG_BC(NINTG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"     
    CALL GAULEG(NINTG,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   CALCULATE Pre-set Value
    ALLOCATE (NI(NINTP,3,8),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY LI ***"
    ALLOCATE (WEI(NINTP),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY LI ***" 
!
    NN=0
    DO I=1,NINTG
       DO J=1,NINTG
          DO K=1,NINTG
             NN=NN+1
             WEI=XABSC_WEIG_BC(I,2)*XABSC_WEIG_BC(J,2)*XABSC_WEIG_BC(K,2)
             NI(NN,1,1)=-(1.0D0-XABSC_WEIG_BC(J,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,2,1)=-(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,3,1)=-(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(J,1))
!
             NI(NN,1,2)=(1.0D0-XABSC_WEIG_BC(J,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,2,2)=-(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,3,2)=-(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(J,1))
!
             NI(NN,1,3)=(1.0D0+XABSC_WEIG_BC(J,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,2,3)=(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,3,3)=-(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(J,1))
!
             NI(NN,1,4)=-(1.0D0+XABSC_WEIG_BC(J,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,2,4)=(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(K,1))
             NI(NN,3,4)=-(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(J,1))
!
             NI(NN,1,5)=-(1.0D0-XABSC_WEIG_BC(J,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,2,5)=-(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,3,5)=(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(J,1))
!
             NI(NN,1,6)=(1.0D0-XABSC_WEIG_BC(J,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,2,6)=-(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,3,6)=(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0-XABSC_WEIG_BC(J,1))
!
             NI(NN,1,7)=(1.0D0+XABSC_WEIG_BC(J,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,2,7)=(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,3,7)=(1.0D0+XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(J,1))
!
             NI(NN,1,8)=-(1.0D0+XABSC_WEIG_BC(J,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,2,8)=(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(K,1))
             NI(NN,3,8)=(1.0D0-XABSC_WEIG_BC(I,1))*(1.0D0+XABSC_WEIG_BC(J,1))
          ENDDO
       ENDDO
    ENDDO
    NI=NI/8.0D0
!   CONSTITUTIVE MATRIX
    CALL G_KCS(3,0,E,PR)
!   STIFFNESS MATRIX
    DO N=1,NumEle
       S=0.0D0
       CALL HEX8D(N,NINTP,S,NI,WEI)     
!
       CALL ADDBAN_ZLB(24,S,LM(:,N),1)  
    END DO
!
    DEALLOCATE(XYZ)
    DEALLOCATE(LM)
    DEALLOCATE(NI)
    DEALLOCATE(WEI)
    DEALLOCATE(XABSC_WEIG_BC)
    DEALLOCATE(K3CS)
    RETURN
!}}}
    END SUBROUTINE HEX8_G
!
    SUBROUTINE HEX8D(N,NINTP,S,NI,WEI)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N,NINTP
    REAL(8) :: S(24,24),NI(NINTP,3,8),WEI(NINTP)
! IN-FUNCTION VARIABLES
    REAL(8),DIMENSION(3,3) :: D,DINV
    REAL(8),DIMENSION(6,24) :: B_M
    REAL(8) :: DD
    INTEGER :: I,J,N1,N2,N3
!
    B_M=0.0D0
    DO I=1,NINTP
! Calculate Jacobian matrix
       D(1,1)=SUM(XYZ(1:8,N)*NI(I,1,:))
       D(1,2)=SUM(XYZ(9:16,N)*NI(I,1,:))
       D(1,3)=SUM(XYZ(17:24,N)*NI(I,1,:))
!
       D(2,1)=SUM(XYZ(1:8,N)*NI(I,2,:))
       D(2,2)=SUM(XYZ(9:16,N)*NI(I,2,:))
       D(2,3)=SUM(XYZ(17:24,N)*NI(I,2,:))
!
       D(3,1)=SUM(XYZ(1:8,N)*NI(I,3,:))
       D(3,2)=SUM(XYZ(9:16,N)*NI(I,3,:))
       D(3,3)=SUM(XYZ(17:24,N)*NI(I,3,:))
!
       DINV(1,1)=D(2,2)*D(3,3)-D(2,3)*D(3,2)
       DINV(2,2)=D(3,3)*D(1,1)-D(3,1)*D(1,3)
       DINV(3,3)=D(1,1)*D(2,2)-D(1,2)*D(2,1)
       DINV(1,2)=D(3,2)*D(1,3)-D(1,2)*D(3,3)
       DINV(2,3)=D(1,3)*D(2,1)-D(2,3)*D(1,1)
       DINV(3,1)=D(2,1)*D(3,2)-D(3,1)*D(2,2)
       DINV(2,1)=D(2,3)*D(3,1)-D(2,1)*D(3,3)
       DINV(3,2)=D(3,1)*D(1,2)-D(3,2)*D(1,1)
       DINV(1,3)=D(1,2)*D(2,3)-D(1,3)*D(2,2)
       DD=D(1,1)*DINV(1,1)+D(1,2)*DINV(2,1)+D(1,3)*DINV(3,1)
!            
       DINV=DINV/DD
!      Build B_M matrix
       DO J=1,8
          N1=3*(J-1)+1
          B_M(1,N1)=SUM(DINV(1,:)*NI(I,:,J))
          N2=3*(J-1)+2
          B_M(2,N2)=SUM(DINV(2,:)*NI(I,:,J))
          N3=3*J
          B_M(3,N3)=SUM(DINV(3,:)*NI(I,:,J))
!
          B_M(4,N1)=SUM(DINV(2,:)*NI(I,:,J))
          B_M(4,N2)=SUM(DINV(1,:)*NI(I,:,J))
!
          B_M(5,N2)=SUM(DINV(3,:)*NI(I,:,J))
          B_M(5,N3)=SUM(DINV(2,:)*NI(I,:,J))
!
          B_M(6,N1)=SUM(DINV(3,:)*NI(I,:,J))
          B_M(6,N3)=SUM(DINV(1,:)*NI(I,:,J))
       ENDDO
       S=S+MATMUL(TRANSPOSE(B_M),MATMUL(K3CS,B_M))*DD*WEI(I)       
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE HEX8D
!
    SUBROUTINE OFE_TETRA_G(ADELE,ADBC)
!{{{
!----------------------------------------------------------------
!     PURPOSE:
!           CALCULATE THE HALF-BANDWIDTH WITH 
!           ZERO ENTRIES ELIMINATED
!           UPDATE HALF_BAND
!     OUTPUT:
!           HALF_BAND  : EFFECTIVE HALF BANDWIDTH LOGICAL MATRIX
!-----------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!IN FUNCTION VARIABLES
    INTEGER :: I
    INTEGER :: MTYPE,NumELE,NumOLR,NGDOF
    INTEGER :: NumBc
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: N,L,M
    INTEGER :: RULE
    INTEGER,DIMENSION(8):: ELE_NM,ND_NM
    INTEGER :: KK,HT
    INTEGER :: ALLocateStatus
!
    MTYPE =NPAR(2)
    NumEle=NPAR(3)
    NumOLR=NPAR(4)
    NGDOF =NPAR(5)
!  SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)       
!   CALCULATE THE NUMBER OF INTEGRATION POINTS
    SELECT CASE(NGDOF)
    CASE(1)
        RULE=1
        NITP_BC=1
    CASE(2)
        RULE=2
        NITP_BC=3
    CASE(3)
        RULE=4
        NITP_BC=6
    CASE(4)
        RULE=6
        NITP_BC=7
    CASE(5)
        RULE=7
        NITP_BC=9
    CASE(6)
        RULE=8
        NITP_BC=12
    CASE(7)
        RULE=9
        NITP_BC=13
    CASE(8)
        RULE=10
        NITP_BC=16
    CASE DEFAULT
        WRITE(*,*) "ERROR,GENERATIONM/OFE_TETRA_G"
        STOP
    END SELECT
!
!-------- Allocate memory
!
    ALLOCATE (NDMFS(7,NumELE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(5,NumELE),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    READ(ADELE) ((NDMFS(L,M),L=1,7),M=1,NPAR(3)),((ELEPAR(L,M),L=1,5),M=1,NPAR(3))
    ALLOCATE (NDINT(4,NumOLR),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    READ(ADELE) ((NDINT(L,M),L=1,4),M=1,NPAR(4))
!
!   CALCULATE BOUNDARY
    READ(ADBC) (NPAR(I),I=1,7)
    NumBc=NPAR(2)
    ALLOCATE (BCED(3,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(6,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(3,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"      
    READ(ADBC) ((BCED(L,M),L=1,3),M=1,NumBc),((BCIND(L,M),L=1,6),M=1,NumBc),((BCPAR(L,M),L=1,3),M=1,NumBc)
!
    CALL OFE_SURFACE(NGDOF,NumBc)
!
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR)
!------- Allocate memory
!
    CALL KEAST_ORDER_NUM(RULE,NITP_IN)
!   CALL MATERIAL MODULE
    CALL G_KCS(3,0,E,PR)
!    DO I=1,NUMELE
!       WRITE(*,*) NDMFS(:,I)
!       WRITE(*,*) ELEPAR(:,I)
!    ENDDO
!
!   GENERATING STIFF MATRIX
    CALL OFE_TETRA_STIFF(NUMOLR,RULE)
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(NDINT)
    DEALLOCATE(K3CS)
    RETURN
!}}}
    END SUBROUTINE OFE_TETRA_G
!{{{
    SUBROUTINE OFE_SURFACE(NGDOF,NumBc)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NGDOF,NumBc
! IN-FUNCTION VARIABLES
!------------------------------------------
    REAL(8),DIMENSION(NITP_BC)   :: LILJ,LJLK
    REAL(8),DIMENSION(NITP_BC,3) :: AREA_COOR
    REAL(8),DIMENSION(NITP_BC)   :: GAUSS_WEI
!------------------------------------------
    INTEGER              :: I,J,K,L
    INTEGER              :: Ele1,Ele2,Ele3,N1,N2,N3
    REAL(8)              :: RADI,RADJ,RADK
    REAL(8),DIMENSION(3) :: XY1,XY2,XY3,XYC1,XYC2,XYC3,XYCM1,XYCM2,XYCM3
    REAL(8)              :: L1,L2,L3,L1_HF,L2_HF,L3_HF,LM1,LM2,LM3,AREA2
    REAL(8),DIMENSION(3) :: S1,S2,S3,S4,WI,WJ,WK
    REAL(8),DIMENSION(2) :: PHI1,PHI2,PHI3,PHI4,PHI5,PHI6,PHI7,PHI8,PHI9,CONS1,CONS2,CONS3,CONS4,CONS5,CONS6
    REAL(8),DIMENSION(2) :: C1,C2,C3,C4,C5,C6
    REAL(8),DIMENSION(NITP_BC) :: XTMP,YTMP,ZTMP,X_BAR,Y_BAR,Z_BAR
    REAL(8),DIMENSION(NITP_BC) :: HIW,HJW,HKW
    REAL(8),DIMENSION(NITP_BC) :: HIWC2,HIWC1,HJWC2,HJWC1
!
    INTEGER :: INDBCI,INDBCJ,INDBCK
    INTEGER :: NDOFI,NDOFJ,NDOFK
    INTEGER :: NCOUNT
    REAL(8) :: CONST1
    REAL(8),DIMENSION(NITP_BC)  :: CONSTX,CONSTY,CONSTZ
!
    REAL(8),DIMENSION(3,3) :: XYZ
    REAL(8),DIMENSION(3)   :: AXYZ,BXYZ
    LOGICAL :: LG_SURFX,LG_SURFY,LG_SURFZ
    REAL(8) :: NX,NY,NZ,NORM
    REAL(8),DIMENSION(NITP_BC,3) :: LOAD
    INTEGER :: NUM1,NUM2,JJ2,JJ1,J1
    REAL(8),DIMENSION(10,3) :: SI,SJ,SK
!
    CALL GAUTRIANGLE_OFE1(NGDOF,NITP_BC,AREA_COOR,GAUSS_WEI)
!   CALCULATE Pre-Set Value
    LILJ=AREA_COOR(:,1)*AREA_COOR(:,2)
    LJLK=AREA_COOR(:,2)*AREA_COOR(:,3)
!
    DO I=1,NumBc
       Ele1=BCED(1,I)
       Ele2=BCED(2,I)
       Ele3=BCED(3,I)
!
       LG_SURFX=BCIND(1,I).EQ.1
       LG_SURFY=BCIND(2,I).EQ.1
       LG_SURFZ=BCIND(3,I).EQ.1
       IF ((.NOT. LG_SURFX) .AND. (.NOT. LG_SURFY) .AND. (.NOT. LG_SURFZ)) CYCLE
!
       IF (LG_SURFX) THEN
          NDMFS(5,Ele1)=1
          NDMFS(5,Ele2)=1
          NDMFS(5,Ele3)=1
       ENDIF
!
       IF (LG_SURFY) THEN
          NDMFS(6,Ele1)=1
          NDMFS(6,Ele2)=1
          NDMFS(6,Ele3)=1
       ENDIF
!
       IF (LG_SURFZ) THEN
          NDMFS(7,Ele1)=1
          NDMFS(7,Ele2)=1
          NDMFS(7,Ele3)=1
       ENDIF
!
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
!      N1
       XYZ(1,1)=X(N1)
       XYZ(2,1)=Y(N1)
       XYZ(3,1)=Z(N1)
!      N2
       XYZ(1,2)=X(N2)
       XYZ(2,2)=Y(N2)
       XYZ(3,2)=Z(N2)
!      N3
       XYZ(1,3)=X(N3)
       XYZ(2,3)=Y(N3)
       XYZ(3,3)=Z(N3)
!      
       AXYZ=XYZ(:,2)-XYZ(:,1)
       BXYZ=XYZ(:,3)-XYZ(:,1)
!
       NORM=0.0D0
       NX=AXYZ(2)*BXYZ(3)-AXYZ(3)*BXYZ(2)
       IF (ABS(NX)<1E-12) THEN
           NX=0
       ELSE
           NORM=NX*NX
       ENDIF
       NY=AXYZ(3)*BXYZ(1)-AXYZ(1)*BXYZ(3)
       IF (ABS(NY)<1E-12) THEN
           NY=0
       ELSE
           NORM=NORM+NY*NY
       ENDIF
       NZ=AXYZ(1)*BXYZ(2)-AXYZ(2)*BXYZ(1)
       IF (ABS(NZ)<1E-12) THEN
           NZ=0
       ELSE
           NORM=NORM+NZ*NZ
       ENDIF
!
       NORM=SQRT(NORM)
       IF (ABS(NX)<1E-12) THEN
          NX=0
       ELSE
          NX=NX/NORM
       ENDIF
       IF (ABS(NY)<1E-12) THEN
          NY=0
       ELSE
          NY=NY/NORM
       ENDIF
       IF (ABS(NZ)<1E-12) THEN
          NZ=0
       ELSE
          NZ=NZ/NORM
       ENDIF
!
       ELEPAR(3,Ele1)=NX
       ELEPAR(4,Ele1)=NY
       ELEPAR(5,Ele1)=NZ
!      
       ELEPAR(3,Ele2)=NX
       ELEPAR(4,Ele2)=NY
       ELEPAR(5,Ele2)=NZ
!      
       ELEPAR(3,Ele3)=NX
       ELEPAR(4,Ele3)=NY
       ELEPAR(5,Ele3)=NZ
    ENDDO
!
    DO I=1,NumBc     
       IF ((BCIND(1,I).EQ.1).AND.(BCIND(2,I).EQ.1).AND.(BCIND(3,I).EQ.1)) CYCLE
!
       Ele1=BCED(1,I)
       Ele2=BCED(2,I)
       Ele3=BCED(3,I)
! 
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
!
       NDOFI=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFJ=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFK=NDOF(1,N3)+NDMFS(2,Ele3)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
       RADK=ELEPAR(1,Ele3)
!
       INDBCI=NDMFS(5,Ele1)
       INDBCJ=NDMFS(5,Ele2)
       INDBCK=NDMFS(5,Ele3)
!
       XY1(1)=X(N1)
       XY1(2)=Y(N1)
       XY1(3)=Z(N1)
       XY2(1)=X(N2)
       XY2(2)=Y(N2)
       XY2(3)=Z(N2)
       XY3(1)=X(N3)
       XY3(2)=Y(N3)
       XY3(3)=Z(N3)
! -----------------------------------
       XYC1=XY1-XY2
       XYC2=XY1-XY3
       XYC3=XY2-XY3
       XYCM1=XYC1+XYC2
       XYCM2=XYC3-XYC1
       XYCM3=XYC3+XYC2
! -----------------------------------
       L1=SQRT(DOT_PRODUCT(XYC1,XYC1))
       L1_HF=L1/2.0D0
       L3=SQRT(DOT_PRODUCT(XYC2,XYC2))
       L3_HF=L3/2.0D0
       L2=SQRT(DOT_PRODUCT(XYC3,XYC3))
       L2_HF=L2/2.0D0
       LM1=SQRT(DOT_PRODUCT(XYCM1,XYCM1))/2.0D0
       LM2=SQRT(DOT_PRODUCT(XYCM2,XYCM2))/2.0D0
       LM3=SQRT(DOT_PRODUCT(XYCM3,XYCM3))/2.0D0
! -----------------------------------   42 FLOP
       AXYZ=XY2-XY1
       BXYZ=XY3-XY1
!
       AREA2=(AXYZ(2)*BXYZ(3)-AXYZ(3)*BXYZ(2))*(AXYZ(2)*BXYZ(3)-AXYZ(3)*BXYZ(2))
       AREA2=AREA2+(AXYZ(3)*BXYZ(1)-AXYZ(1)*BXYZ(3))*(AXYZ(3)*BXYZ(1)-AXYZ(1)*BXYZ(3))
       AREA2=AREA2+(AXYZ(1)*BXYZ(2)-AXYZ(2)*BXYZ(1))*(AXYZ(1)*BXYZ(2)-AXYZ(2)*BXYZ(1))
       AREA2=SQRT(AREA2)
! -----------------------------------   WI
       S1=(/L1_HF,LM1,L3_HF/)/RADI 
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WI=1.0D0-S2+S3-S4
! -----------------------------------   WJ
       S1=(/L2_HF,LM2,L1_HF/)/RADJ
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WJ=1.0D0-S2+S3-S4 
! -----------------------------------   WK
       S1=(/L3_HF,LM3,L2_HF/)/RADK
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WK=1.0D0-S2+S3-S4 
!
       IF(WI(2).LT.0.0D0) WI(2)=0.0D0
       IF(WJ(2).LT.0.0D0) WJ(2)=0.0D0
       IF(WK(2).LT.0.0D0) WK(2)=0.0D0
! -----------------------------------   PHI_I
       PHI1(1)=0.0D0
       PHI2(1)=0.0D0
       PHI3(1)=WI(1)/(WJ(3)+WI(1))
       PHI4(1)=WI(3)/(WK(1)+WI(3))
       PHI5(1)=WI(2)/(WI(2)+WJ(1)+WK(3))+WI(3)/(WI(3)+WJ(2)+WK(1))+WI(1)/(WI(1)+WJ(3)+WK(2))
       PHI6(1)=1.0D0
       PHI7(1)=1.0D0
! -----------------------------------   PHI_J
       PHI1(2)=0.0D0
       PHI2(2)=0.0D0
       PHI3(2)=WJ(1)/(WK(3)+WJ(1))
       PHI4(2)=WJ(3)/(WI(1)+WJ(3))
       PHI5(2)=WJ(2)/(WJ(2)+WK(1)+WI(3))+WJ(3)/(WJ(3)+WK(2)+WI(1))+WJ(1)/(WJ(1)+WK(3)+WI(2))
       PHI6(2)=1.0D0
       PHI7(2)=1.0D0
! -----------------------------------   197+46=243 FLOP
! -----------------------------------   C1 C2 C3 C4 C5 C6
       CONS1=-1.0D0
       CONS2=2.0D0*PHI4
       CONS3=2.0D0*(PHI3-PHI5)
       CONS4=-1.0D0+4.0D0*PHI4
       CONS5=CONS2
       CONS6=1.0D0+CONS2
!
       C1=2.0D0*(CONS3+CONS5)
       C2=2.0D0*(CONS2+CONS3)
       C3=4.0D0*(-CONS2+PHI5)
       C4=CONS4
       C5=2.0D0-4.0D0*PHI4
       C6=0.0D0 
! -----------------------------------   243+48 = 291 FLOP
!------------------------------------------  XTMP,YTMP 
       XTMP=AREA_COOR(:,1)*XY1(1)+AREA_COOR(:,2)*XY2(1)+AREA_COOR(:,3)*XY3(1)
       YTMP=AREA_COOR(:,1)*XY1(2)+AREA_COOR(:,2)*XY2(2)+AREA_COOR(:,3)*XY3(2)
       ZTMP=AREA_COOR(:,1)*XY1(3)+AREA_COOR(:,2)*XY2(3)+AREA_COOR(:,3)*XY3(3)
!------------------------------------------ HIW,HIWC1,HIWC2, HJW, HJWC1,HJWC2
       HIWC1=AREA_COOR(:,2)*C1(1)+AREA_COOR(:,1)*C2(1)+C3(1)
       HIWC2=C4(1)+AREA_COOR(:,1)*C5(1)
       HIW  =LILJ*HIWC1+AREA_COOR(:,1)*HIWC2 
!
       HJWC1=AREA_COOR(:,3)*C1(2)+AREA_COOR(:,2)*C2(2)+C3(2)
       HJWC2=C4(2)+AREA_COOR(:,2)*C5(2)
       HJW  =LJLK*HJWC1+AREA_COOR(:,2)*HJWC2
!
       HKW  =1.0D0-HIW-HJW
!
!       HIW=AREA_COOR(:,1)
!       HJW=AREA_COOR(:,2)
!       HKW=AREA_COOR(:,3)
!
       LOAD(:,1)=BCPAR(1,I)*GAUSS_WEI
       LOAD(:,2)=BCPAR(2,I)*GAUSS_WEI
       LOAD(:,3)=BCPAR(3,I)*GAUSS_WEI
!
       CONST1=AREA2/2.0D0
       SI=0.0D0
       SJ=0.0D0
       SK=0.0D0
! **************************** 
!{
       IF(INDBCI.EQ.1) THEN
          NX=ELEPAR(3,Ele1)
          NY=ELEPAR(4,Ele1)
          NZ=ELEPAR(5,Ele1)
!    
          CONSTX=(XTMP-X(N1))/RADI
          CONSTY=(YTMP-Y(N1))/RADI
          CONSTZ=(ZTMP-Z(N1))/RADI
!  
          X_BAR=NX*CONSTX+NY*CONSTY+NZ*CONSTZ
          Y_BAR=-NY*CONSTX+NX*CONSTY
          Z_BAR=-NX*NZ*CONSTX-NY*NZ*CONSTY+(NX*NX+NY*NY)*CONSTZ
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1      
              SI(1,3)=SUM(HIW*LOAD(:,3))*CONST1  
!
              SI(2,1)=SUM(HIW*X_BAR*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*X_BAR*LOAD(:,2))*CONST1      
              SI(2,3)=SUM(HIW*X_BAR*LOAD(:,3))*CONST1  
!
              SI(3,1)=SUM(HIW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SI(3,3)=SUM(HIW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SI(4,1)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SI(4,3)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SI(5,1)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SI(5,2)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SI(5,3)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,3))*CONST1         
       ELSE
          X_BAR=(XTMP-X(N1))/RADI
          Y_BAR=(YTMP-Y(N1))/RADI
          Z_BAR=(ZTMP-Z(N1))/RADI
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1      
              SI(1,3)=SUM(HIW*LOAD(:,3))*CONST1  
!
              SI(2,1)=SUM(HIW*X_BAR*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*X_BAR*LOAD(:,2))*CONST1      
              SI(2,3)=SUM(HIW*X_BAR*LOAD(:,3))*CONST1  
!
              SI(3,1)=SUM(HIW*Y_BAR*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*Y_BAR*LOAD(:,2))*CONST1      
              SI(3,3)=SUM(HIW*Y_BAR*LOAD(:,3))*CONST1  
!
              SI(4,1)=SUM(HIW*Z_BAR*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*Z_BAR*LOAD(:,2))*CONST1      
              SI(4,3)=SUM(HIW*Z_BAR*LOAD(:,3))*CONST1              
!
              SI(5,1)=SUM(HIW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SI(5,2)=SUM(HIW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SI(5,3)=SUM(HIW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SI(6,1)=SUM(HIW*Y_BAR*Y_BAR*LOAD(:,1))*CONST1
              SI(6,2)=SUM(HIW*Y_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SI(6,3)=SUM(HIW*Y_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SI(7,1)=SUM(HIW*Z_BAR*Z_BAR*LOAD(:,1))*CONST1
              SI(7,2)=SUM(HIW*Z_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SI(7,3)=SUM(HIW*Z_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SI(8,1)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SI(8,2)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SI(8,3)=SUM(HIW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SI(9,1)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SI(9,2)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SI(9,3)=SUM(HIW*X_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SI(10,1)=SUM(HIW*Y_BAR*Z_BAR*LOAD(:,1))*CONST1
              SI(10,2)=SUM(HIW*Y_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SI(10,3)=SUM(HIW*Y_BAR*Z_BAR*LOAD(:,3))*CONST1  
       ENDIF
!}}}
!{{{
       IF(INDBCJ.EQ.1) THEN
          NX=ELEPAR(3,Ele2)
          NY=ELEPAR(4,Ele2)
          NZ=ELEPAR(5,Ele2)
!    
          CONSTX=(XTMP-X(N2))/RADJ
          CONSTY=(YTMP-Y(N2))/RADJ
          CONSTZ=(ZTMP-Z(N2))/RADJ
!  
          X_BAR=NX*CONSTX+NY*CONSTY+NZ*CONSTZ
          Y_BAR=-NY*CONSTX+NX*CONSTY
          Z_BAR=-NX*NZ*CONSTX-NY*NZ*CONSTY+(NX*NX+NY*NY)*CONSTZ
!
              SJ(1,1)=SUM(HJW*LOAD(:,1))*CONST1
              SJ(1,2)=SUM(HJW*LOAD(:,2))*CONST1      
              SJ(1,3)=SUM(HJW*LOAD(:,3))*CONST1  
!
              SJ(2,1)=SUM(HJW*X_BAR*LOAD(:,1))*CONST1
              SJ(2,2)=SUM(HJW*X_BAR*LOAD(:,2))*CONST1      
              SJ(2,3)=SUM(HJW*X_BAR*LOAD(:,3))*CONST1  
!
              SJ(3,1)=SUM(HJW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SJ(3,2)=SUM(HJW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SJ(3,3)=SUM(HJW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SJ(4,1)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SJ(4,2)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SJ(4,3)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SJ(5,1)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SJ(5,2)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SJ(5,3)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,3))*CONST1  
      
       ELSE
          X_BAR=(XTMP-X(N2))/RADJ
          Y_BAR=(YTMP-Y(N2))/RADJ
          Z_BAR=(ZTMP-Z(N2))/RADJ
              SJ(1,1)=SUM(HJW*LOAD(:,1))*CONST1
              SJ(1,2)=SUM(HJW*LOAD(:,2))*CONST1      
              SJ(1,3)=SUM(HJW*LOAD(:,3))*CONST1  
!
              SJ(2,1)=SUM(HJW*X_BAR*LOAD(:,1))*CONST1
              SJ(2,2)=SUM(HJW*X_BAR*LOAD(:,2))*CONST1      
              SJ(2,3)=SUM(HJW*X_BAR*LOAD(:,3))*CONST1  
!
              SJ(3,1)=SUM(HJW*Y_BAR*LOAD(:,1))*CONST1
              SJ(3,2)=SUM(HJW*Y_BAR*LOAD(:,2))*CONST1      
              SJ(3,3)=SUM(HJW*Y_BAR*LOAD(:,3))*CONST1  
!
              SJ(4,1)=SUM(HJW*Z_BAR*LOAD(:,1))*CONST1
              SJ(4,2)=SUM(HJW*Z_BAR*LOAD(:,2))*CONST1      
              SJ(4,3)=SUM(HJW*Z_BAR*LOAD(:,3))*CONST1              
!
              SJ(5,1)=SUM(HJW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SJ(5,2)=SUM(HJW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SJ(5,3)=SUM(HJW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SJ(6,1)=SUM(HJW*Y_BAR*Y_BAR*LOAD(:,1))*CONST1
              SJ(6,2)=SUM(HJW*Y_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SJ(6,3)=SUM(HJW*Y_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SJ(7,1)=SUM(HJW*Z_BAR*Z_BAR*LOAD(:,1))*CONST1
              SJ(7,2)=SUM(HJW*Z_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SJ(7,3)=SUM(HJW*Z_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SJ(8,1)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SJ(8,2)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SJ(8,3)=SUM(HJW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SJ(9,1)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SJ(9,2)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SJ(9,3)=SUM(HJW*X_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SJ(10,1)=SUM(HJW*Y_BAR*Z_BAR*LOAD(:,1))*CONST1
              SJ(10,2)=SUM(HJW*Y_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SJ(10,3)=SUM(HJW*Y_BAR*Z_BAR*LOAD(:,3))*CONST1  

       ENDIF
!}}}
! **************************** dHKmdx,dHKmdy
!{{{
       IF(INDBCK.EQ.1) THEN
          NX=ELEPAR(3,Ele3)
          NY=ELEPAR(4,Ele3)
          NZ=ELEPAR(5,Ele3)
!    
          CONSTX=(XTMP-X(N3))/RADK
          CONSTY=(YTMP-Y(N3))/RADK
          CONSTZ=(ZTMP-Z(N3))/RADK
!  
          X_BAR=NX*CONSTX+NY*CONSTY+NZ*CONSTZ
          Y_BAR=-NY*CONSTX+NX*CONSTY
          Z_BAR=-NX*NZ*CONSTX-NY*NZ*CONSTY+(NX*NX+NY*NY)*CONSTZ
!
              SK(1,1)=SUM(HKW*LOAD(:,1))*CONST1
              SK(1,2)=SUM(HKW*LOAD(:,2))*CONST1      
              SK(1,3)=SUM(HKW*LOAD(:,3))*CONST1  
!
              SK(2,1)=SUM(HKW*X_BAR*LOAD(:,1))*CONST1
              SK(2,2)=SUM(HKW*X_BAR*LOAD(:,2))*CONST1      
              SK(2,3)=SUM(HKW*X_BAR*LOAD(:,3))*CONST1  
!
              SK(3,1)=SUM(HKW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SK(3,2)=SUM(HKW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SK(3,3)=SUM(HKW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SK(4,1)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SK(4,2)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SK(4,3)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SK(5,1)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SK(5,2)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SK(5,3)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,3))*CONST1  
      
       ELSE
          X_BAR=(XTMP-X(N3))/RADK
          Y_BAR=(YTMP-Y(N3))/RADK
          Z_BAR=(ZTMP-Z(N3))/RADK
              SK(1,1)=SUM(HKW*LOAD(:,1))*CONST1
              SK(1,2)=SUM(HKW*LOAD(:,2))*CONST1      
              SK(1,3)=SUM(HKW*LOAD(:,3))*CONST1  
!
              SK(2,1)=SUM(HKW*X_BAR*LOAD(:,1))*CONST1
              SK(2,2)=SUM(HKW*X_BAR*LOAD(:,2))*CONST1      
              SK(2,3)=SUM(HKW*X_BAR*LOAD(:,3))*CONST1  
!
              SK(3,1)=SUM(HKW*Y_BAR*LOAD(:,1))*CONST1
              SK(3,2)=SUM(HKW*Y_BAR*LOAD(:,2))*CONST1      
              SK(3,3)=SUM(HKW*Y_BAR*LOAD(:,3))*CONST1  
!
              SK(4,1)=SUM(HKW*Z_BAR*LOAD(:,1))*CONST1
              SK(4,2)=SUM(HKW*Z_BAR*LOAD(:,2))*CONST1      
              SK(4,3)=SUM(HKW*Z_BAR*LOAD(:,3))*CONST1              
!
              SK(5,1)=SUM(HKW*X_BAR*X_BAR*LOAD(:,1))*CONST1
              SK(5,2)=SUM(HKW*X_BAR*X_BAR*LOAD(:,2))*CONST1      
              SK(5,3)=SUM(HKW*X_BAR*X_BAR*LOAD(:,3))*CONST1  
!
              SK(6,1)=SUM(HKW*Y_BAR*Y_BAR*LOAD(:,1))*CONST1
              SK(6,2)=SUM(HKW*Y_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SK(6,3)=SUM(HKW*Y_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SK(7,1)=SUM(HKW*Z_BAR*Z_BAR*LOAD(:,1))*CONST1
              SK(7,2)=SUM(HKW*Z_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SK(7,3)=SUM(HKW*Z_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SK(8,1)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,1))*CONST1
              SK(8,2)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,2))*CONST1      
              SK(8,3)=SUM(HKW*X_BAR*Y_BAR*LOAD(:,3))*CONST1  
!
              SK(9,1)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,1))*CONST1
              SK(9,2)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SK(9,3)=SUM(HKW*X_BAR*Z_BAR*LOAD(:,3))*CONST1  
!
              SK(10,1)=SUM(HKW*Y_BAR*Z_BAR*LOAD(:,1))*CONST1
              SK(10,2)=SUM(HKW*Y_BAR*Z_BAR*LOAD(:,2))*CONST1      
              SK(10,3)=SUM(HKW*Y_BAR*Z_BAR*LOAD(:,3))*CONST1  
       ENDIF
!}}}
       DO J=1,3
          IF(BCIND(J,I).EQ.1) CYCLE
!
          J1=1+J
!
          NUM1=ID(J,N1)
          NUM2=NUM1+NDOF(J,N1)-1
!
          JJ1=1+NDMFS(J1,Ele1)
          JJ2=NDOF(J,N1)+NDMFS(J1,Ele1)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SI(JJ1:JJ2,J)
!--------------------------------------------------
          NUM1=ID(J,N2)
          NUM2=NUM1+NDOF(J,N2)-1
!
          JJ1=1+NDMFS(J1,Ele2)
          JJ2=NDOF(J,N2)+NDMFS(J1,Ele2)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SJ(JJ1:JJ2,J)
!--------------------------------------------------
          NUM1=ID(J,N3)
          NUM2=NUM1+NDOF(J,N3)-1
!
          JJ1=1+NDMFS(J1,Ele3)
          JJ2=NDOF(J,N3)+NDMFS(J1,Ele3)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SK(JJ1:JJ2,J)

       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OFE_SURFACE
!
    SUBROUTINE OFE_TETRA_STIFF(NUMOLR,RULE)
!{{{
!---------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NUMOLR,RULE
!---------------------------------------------
    REAL(8),DIMENSION(NITP_IN,4) :: FI,DFIX,DFIY,DFIZ
    REAL(8) :: DETJ
    REAL(8),DIMENSION(3,NITP_IN) :: AREA_COORD
    REAL(8),DIMENSION(NITP_IN)   :: WEI
    REAL(8),DIMENSION(NITP_IN,4) :: HI,HI_ASS
    REAL(8),DIMENSION(NITP_IN,10):: HI2
!
    INTEGER :: I,J,K,L
    INTEGER :: ELE1,ELE2,ELE3,ELE4,N1,N2,N3,N4
    REAL(8),DIMENSION(4) :: RAD
    REAL(8),DIMENSION(3,10) :: XYZ
!
    INTEGER,DIMENSION(5) :: NDOFI
    INTEGER,DIMENSION(4) :: INDBC,ELEI,NOD
    REAL(8),DIMENSION(NITP_IN,10) :: POLY,DXPOLY,DYPOLY,DZPOLY
    REAL(8),DIMENSION(NITP_IN) :: XTMP,YTMP,ZTMP,X_BAR,Y_BAR,Z_BAR,X_CONST,Y_CONST,Z_CONST
    REAL(8) :: NX,NY,NZ
!
    INTEGER :: IS1,IS2,IS3,IE1,IE2,IE3
    REAL(8),DIMENSION(NITP_IN,12*NPOR_MAX,6) :: Bmat
    REAL(8),DIMENSION(12*NPOR_MAX,6) :: Bmat2
    REAL(8),DIMENSION(12*NPOR_MAX,12*NPOR_MAX) :: S
    INTEGER,DIMENSION(12*NPOR_MAX) :: LM
    INTEGER :: NN,IDD,KK,L0
!-----------------------------------------------
    CALL KEAST_RULE(RULE,NITP_IN,AREA_COORD,WEI)
!   HI
    HI(:,1)=1.0D0-AREA_COORD(1,:)-AREA_COORD(2,:)-AREA_COORD(3,:)
    HI(:,2)=AREA_COORD(1,:)
    HI(:,3)=AREA_COORD(2,:)
    HI(:,4)=AREA_COORD(3,:)
!   HI_ASS
    HI_ASS=4.0D0*HI-1.0D0  
!   HI2
    HI2=HI*(2.0D0*HI-1.0D0)
!
    HI2(:,5)=4.0D0*HI(:,1)*HI(:,2)
    HI2(:,6)=4.0D0*HI(:,2)*HI(:,3)
    HI2(:,7)=4.0D0*HI(:,3)*HI(:,1)
    HI2(:,8)=4.0D0*HI(:,1)*HI(:,4)
    HI2(:,9)=4.0D0*HI(:,2)*HI(:,4)
    HI2(:,10)=4.0D0*HI(:,3)*HI(:,4)
!------------------------------------------   
    DO I=1,NUMOLR
!   ELE,N,RAD
       ELE1=NDINT(1,I)
       ELE2=NDINT(2,I)
       ELE3=NDINT(3,I)
       ELE4=NDINT(4,I)
!
       ELEI(1)=ELE1
       ELEI(2)=ELE2
       ELEI(3)=ELE3
       ELEI(4)=ELE4
!
       N1=NDMFS(1,ELE1)
       N2=NDMFS(1,ELE2)
       N3=NDMFS(1,ELE3)
       N4=NDMFS(1,ELE4)
!
       NOD(1)=N1
       NOD(2)=N2
       NOD(3)=N3
       NOD(4)=N4
!
       NDOFI(1)=0
       NDOFI(2)=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFI(3)=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFI(4)=NDOF(1,N3)+NDMFS(2,Ele3)
       NDOFI(5)=NDOF(1,N4)+NDMFS(2,Ele4)
!
       INDBC(1)=NDMFS(5,Ele1)
       INDBC(2)=NDMFS(5,Ele2)
       INDBC(3)=NDMFS(5,Ele3)
       INDBC(4)=NDMFS(5,Ele4)
!      N1
       XYZ(1,1)=X(N1)
       XYZ(2,1)=Y(N1)
       XYZ(3,1)=Z(N1)
!      N2
       XYZ(1,2)=X(N2)
       XYZ(2,2)=Y(N2)
       XYZ(3,2)=Z(N2)
!      N3
       XYZ(1,3)=X(N3)
       XYZ(2,3)=Y(N3)
       XYZ(3,3)=Z(N3)
!      N4
       XYZ(1,4)=X(N4)
       XYZ(2,4)=Y(N4)
       XYZ(3,4)=Z(N4)
!      
       XYZ(:,5)=(XYZ(:,1)+XYZ(:,2))/2.0D0
       XYZ(:,6)=(XYZ(:,2)+XYZ(:,3))/2.0D0
       XYZ(:,7)=(XYZ(:,1)+XYZ(:,3))/2.0D0
       XYZ(:,8)=(XYZ(:,1)+XYZ(:,4))/2.0D0
       XYZ(:,9)=(XYZ(:,2)+XYZ(:,4))/2.0D0
       XYZ(:,10)=(XYZ(:,3)+XYZ(:,4))/2.0D0
!    
       RAD(1)=ELEPAR(1,ELE1)
       RAD(2)=ELEPAR(1,ELE2)
       RAD(3)=ELEPAR(1,ELE3)
       RAD(4)=ELEPAR(1,ELE4)
!
       CALL Get_ShapeFun_Derivative(XYZ,RAD,HI,HI_ASS,HI2,FI,DFIX,DFIY,DFIZ,DETJ)
!
       XTMP=HI(:,1)*XYZ(1,1)+HI(:,2)*XYZ(1,2)+HI(:,3)*XYZ(1,3)+HI(:,4)*XYZ(1,4)
       YTMP=HI(:,1)*XYZ(2,1)+HI(:,2)*XYZ(2,2)+HI(:,3)*XYZ(2,3)+HI(:,4)*XYZ(2,4)
       ZTMP=HI(:,1)*XYZ(3,1)+HI(:,2)*XYZ(3,2)+HI(:,3)*XYZ(3,3)+HI(:,4)*XYZ(3,4)
!
       IS1=1
!
       S=0.0D0
       Bmat=0.0D0
       Bmat2=0.0D0
       DO J=1,4
          IF (INDBC(J)==1) THEN
              X_CONST=(XTMP-XYZ(1,J))/RAD(J)
              Y_CONST=(YTMP-XYZ(2,J))/RAD(J)
              Z_CONST=(ZTMP-XYZ(3,J))/RAD(J)
!
              NX=ELEPAR(3,ELEI(J))
              NY=ELEPAR(4,ELEI(J))
              NZ=ELEPAR(5,ELEI(J))
!
              X_BAR=NX*X_CONST+NY*Y_CONST+NZ*Z_CONST
              Y_BAR=-NY*X_CONST+NX*Y_CONST
              Z_BAR=-NX*NZ*X_CONST-NY*NZ*Y_CONST+(NX*NX+NY*NY)*Z_CONST
!
              POLY(:,1) =1.0D0
              POLY(:,2) =X_BAR
              POLY(:,3) =X_BAR*X_BAR
              POLY(:,4) =X_BAR*Y_BAR
              POLY(:,5) =X_BAR*Z_BAR
!
              DXPOLY(:,1)=0.0D0
              DXPOLY(:,2)=NX/RAD(J)
              DXPOLY(:,3)=2.0D0*X_BAR*DXPOLY(:,2)
              DXPOLY(:,4)=DXPOLY(:,2)*Y_BAR-X_BAR*NY/RAD(J)
              DXPOLY(:,5)=DXPOLY(:,2)*Z_BAR-X_BAR*NX*NZ/RAD(J)
              DXPOLY=DXPOLY*DETJ
!
              DYPOLY(:,1)=0.0D0
              DYPOLY(:,2)=NY/RAD(J)
              DYPOLY(:,3)=2.0D0*X_BAR*DYPOLY(:,2)
              DYPOLY(:,4)=DYPOLY(:,2)*Y_BAR+X_BAR*NX/RAD(J)
              DYPOLY(:,5)=DYPOLY(:,2)*Z_BAR-X_BAR*NY*NZ/RAD(J)
              DYPOLY=DYPOLY*DETJ
!
              DZPOLY(:,1)=0.0D0
              DZPOLY(:,2)=NZ/RAD(J)
              DZPOLY(:,3)=2.0D0*X_BAR*DZPOLY(:,2)
              DZPOLY(:,4)=DZPOLY(:,2)*Y_BAR
              DZPOLY(:,5)=DZPOLY(:,2)*Z_BAR+X_BAR*(NX*NX+NY*NY)/RAD(J)
              DZPOLY=DZPOLY*DETJ
          ELSE
              X_BAR=(XTMP-XYZ(1,J))/RAD(J)
              Y_BAR=(YTMP-XYZ(2,J))/RAD(J)
              Z_BAR=(ZTMP-XYZ(3,J))/RAD(J)
!             
              POLY(:,1) =1.0D0
              POLY(:,2) =X_BAR
              POLY(:,3) =Y_BAR
              POLY(:,4) =Z_BAR
              POLY(:,5) =X_BAR*X_BAR
              POLY(:,6) =Y_BAR*Y_BAR
              POLY(:,7) =Z_BAR*Z_BAR
              POLY(:,8) =X_BAR*Y_BAR
              POLY(:,9) =X_BAR*Z_BAR
              POLY(:,10)=Y_BAR*Z_BAR          
!             
              DXPOLY(:,1) =0.0D0
              DXPOLY(:,2) =1.0D0/RAD(J)
              DXPOLY(:,3) =0.0D0
              DXPOLY(:,4) =0.0D0
              DXPOLY(:,5) =2.0D0*X_BAR*DXPOLY(:,2)
              DXPOLY(:,6) =0.0D0
              DXPOLY(:,7) =0.0D0
              DXPOLY(:,8) =Y_BAR*DXPOLY(:,2)
              DXPOLY(:,9) =Z_BAR*DXPOLY(:,2)
              DXPOLY(:,10)=0.0D0
              DXPOLY=DXPOLY*DETJ
!             
              DYPOLY(:,1) =0.0D0
              DYPOLY(:,2) =0.0D0
              DYPOLY(:,3) =1.0D0/RAD(J)
              DYPOLY(:,4) =0.0D0
              DYPOLY(:,5) =0.0D0
              DYPOLY(:,6) =2.0D0*Y_BAR*DYPOLY(:,3)
              DYPOLY(:,7) =0.0D0
              DYPOLY(:,8) =X_BAR*DYPOLY(:,3)
              DYPOLY(:,9) =0.0D0
              DYPOLY(:,10)=Z_BAR*DYPOLY(:,3)
              DYPOLY=DYPOLY*DETJ
!             
              DZPOLY(:,1) =0.0D0
              DZPOLY(:,2) =0.0D0
              DZPOLY(:,3) =0.0D0
              DZPOLY(:,4) =1.0D0/RAD(J)
              DZPOLY(:,5) =0.0D0
              DZPOLY(:,6) =0.0D0
              DZPOLY(:,7) =2.0D0*Z_BAR*DZPOLY(:,4)
              DZPOLY(:,8) =0.0D0
              DZPOLY(:,9) =X_BAR*DZPOLY(:,4)
              DZPOLY(:,10)=Y_BAR*DZPOLY(:,4)
              DZPOLY=DZPOLY*DETJ
          ENDIF
          IS1=IS1+3*NDOFI(J)
          IE1=IS1+NDOFI(J+1)-1
          IS2=IS1+NDOFI(J+1)
          IE2=IS2+NDOFI(J+1)-1
          IS3=IS2+NDOFI(J+1)
          IE3=IS3+NDOFI(J+1)-1

          DO K=1,NDOFI(J+1)
             Bmat(:,IS1+K-1,1)=FI(:,J)*DXPOLY(:,K)+DFIX(:,J)*POLY(:,K)
             Bmat(:,IS2+K-1,2)=FI(:,J)*DYPOLY(:,K)+DFIY(:,J)*POLY(:,K)
             Bmat(:,IS3+K-1,3)=FI(:,J)*DZPOLY(:,K)+DFIZ(:,J)*POLY(:,K)
!
             Bmat(:,IS1+K-1,4)=FI(:,J)*DYPOLY(:,K)+DFIY(:,J)*POLY(:,K)
             Bmat(:,IS2+K-1,4)=FI(:,J)*DXPOLY(:,K)+DFIX(:,J)*POLY(:,K)       
!
             Bmat(:,IS2+K-1,5)=FI(:,J)*DZPOLY(:,K)+DFIZ(:,J)*POLY(:,K)
             Bmat(:,IS3+K-1,5)=FI(:,J)*DYPOLY(:,K)+DFIY(:,J)*POLY(:,K)       
!
             Bmat(:,IS1+K-1,6)=FI(:,J)*DZPOLY(:,K)+DFIZ(:,J)*POLY(:,K)
             Bmat(:,IS3+K-1,6)=FI(:,J)*DXPOLY(:,K)+DFIX(:,J)*POLY(:,K)
          ENDDO
       ENDDO 
!
       DO J=1,NITP_IN
          Bmat2=Bmat(J,1:12*NPOR_MAX,1:6)
          S=S+MATMUL(Bmat2,MATMUL(K3CS,TRANSPOSE(Bmat2)))*WEI(J) 
       ENDDO
       S=S/DETJ/6.0D0
!------------------------------
!      ASSEMBLE PROCESS
!------------------------------
!      CONSTRUCT LM
!------------------------------
!{{{
       LM=0
       IS1=1
       DO J=1,4
          IS1=IS1+3*NDOFI(J)
          NN=NOD(J) 
          DO K=1,3
             IS2=IS1+(K-1)*NDOFI(J+1)
             IDD=ID(K,NN)
             IF (IDD.NE.0) THEN
                L0=NDMFS(1+K,ELEI(J))+1
                DO L=L0,NDOFI(J+1)
                   IS3=IS2+L-1
                   KK=IDD+L-L0
                   LM(IS3)=KK
                  ! write(*,*) IS3,KK
                ENDDO
             ENDIF
          ENDDO
       ENDDO
!}}}
       CALL OFE_TETRA_ADDBAN(12*NPOR_MAX,LM,S)
    ENDDO
!
    RETURN 
!}}}
    END SUBROUTINE OFE_TETRA_STIFF
!
    SUBROUTINE Get_ShapeFun_Derivative(XYZ,RAD,HI,HI_ASS,HI2,FI,DFIX,DFIY,DFIZ,DETJ)
!{{{
!------------------------------------------
       IMPLICIT NONE
       REAL(8),DIMENSION(3,10) :: XYZ
       REAL(8),DIMENSION(4)   :: RAD
       REAL(8),DIMENSION(NITP_IN,4)  :: HI,HI_ASS
       REAL(8),DIMENSION(NITP_IN,10) :: HI2
       REAL(8),DIMENSION(NITP_IN,4)  :: FI,DFIX,DFIY,DFIZ
       REAL(8) :: DETJ
! IN-FUNCTION VARIABLES
       INTEGER :: I,J
       REAL(8),DIMENSION(10,4):: WI
       REAL(8),DIMENSION(10,4,4) :: PHI
       REAL(8),DIMENSION(NITP_IN,4,4) :: HI2_PHI,DXHI2_PHI,DYHI2_PHI,DZHI2_PHI
       REAL(8),DIMENSION(NITP_IN,4)  :: DHIX,DHIY,DHIZ
       REAL(8),DIMENSION(NITP_IN,10) :: DHIX2,DHIY2,DHIZ2
       REAL(8) :: J11,J12,J13,J14,J21,J22,J23,J24,J31,J32,J33,J34
       REAL(8),DIMENSION(4,3) :: JINV
       REAL(8),DIMENSION(3)   :: S
       REAL(8) :: S1,S2,S3,S4
! ------------------------------------------
!      WI
!------------------------------------------
       DO I=1,4
          DO J=1,10
             S=XYZ(:,I)-XYZ(:,J)
             S2=DOT_PRODUCT(S,S)
             S1=SQRT(S2)/RAD(I)
             S2=S1*S1
             S3=S1*S2
             S4=S2*S2
!
             WI(J,I)=1.0D0-6.0D0*S2+8.0D0*S3-3.0D0*S4
          ENDDO 
       ENDDO
!------------------------------------------
       PHI=0.0D0
!------ PH1 --------------- 
       PHI(1,1,1) = 1.0D0
       PHI(5,1,1) = WI(5,1)/(WI(5,1)+WI(5,2))
       PHI(6,1,1) = WI(6,1)/(WI(6,1)+WI(6,2)+WI(6,3))
       PHI(7,1,1) = WI(7,1)/(WI(7,1)+WI(7,3))
       PHI(8,1,1) = WI(8,1)/(WI(8,1)+WI(8,4))
       PHI(9,1,1) = WI(9,1)/(WI(9,1)+WI(9,2)+WI(9,4))
       PHI(10,1,1)= WI(10,1)/(WI(10,1)+WI(10,3)+WI(10,4))
!
       PHI(1,2,1) = 1.0D0
       PHI(5,2,1) = WI(5,1)/(WI(5,1)+WI(5,2))
       PHI(7,2,1) = WI(7,1)/(WI(7,1)+WI(7,2)+WI(7,3))
       PHI(8,2,1) = WI(8,1)/(WI(8,1)+WI(8,2)+WI(8,4))
!
       PHI(1,3,1) = 1.0D0
       PHI(5,3,1) = WI(5,1)/(WI(5,1)+WI(5,2)+WI(5,3))
       PHI(7,3,1) = WI(7,1)/(WI(7,1)+WI(7,3))
       PHI(8,3,1) = WI(8,1)/(WI(8,1)+WI(8,3)+WI(8,4))
!
       PHI(1,4,1) = 1.0D0
       PHI(5,4,1) = WI(5,1)/(WI(5,1)+WI(5,2)+WI(5,4))
       PHI(7,4,1) = WI(7,1)/(WI(7,1)+WI(7,3)+WI(7,4))
       PHI(8,4,1) = WI(8,1)/(WI(8,1)+WI(8,4))
!
!       PHI(:,:,1)=0.0D0
!       PHI(:,1,1)=1.0D0
!------ PH2 --------------- 
       PHI(2,2,2) = 1.0D0
       PHI(5,2,2) = WI(5,2)/(WI(5,2)+WI(5,1))
       PHI(6,2,2) = WI(6,2)/(WI(6,2)+WI(6,3))
       PHI(7,2,2) = WI(7,2)/(WI(7,2)+WI(7,3)+WI(7,1))
       PHI(8,2,2) = WI(8,2)/(WI(8,2)+WI(8,4)+WI(8,1))
       PHI(9,2,2) = WI(9,2)/(WI(9,2)+WI(9,4))
       PHI(10,2,2)= WI(10,2)/(WI(10,2)+WI(10,3)+WI(10,4))
!
       PHI(2,1,2) = 1.0D0
       PHI(5,1,2) = WI(5,2)/(WI(5,2)+WI(5,1))
       PHI(6,1,2) = WI(6,2)/(WI(6,2)+WI(6,1)+WI(6,3))
       PHI(9,1,2) = WI(9,2)/(WI(9,2)+WI(9,1)+WI(9,4))
!
       PHI(2,3,2) = 1.0D0
       PHI(5,3,2) = WI(5,2)/(WI(5,2)+WI(5,3)+WI(5,1))
       PHI(6,3,2) = WI(6,2)/(WI(6,2)+WI(6,3))
       PHI(9,3,2) = WI(9,2)/(WI(9,2)+WI(9,3)+WI(9,4))
!
       PHI(2,4,2) = 1.0D0
       PHI(5,4,2) = WI(5,2)/(WI(5,2)+WI(5,4)+WI(5,1))
       PHI(6,4,2) = WI(6,2)/(WI(6,2)+WI(6,4)+WI(6,3))
       PHI(9,4,2) = WI(9,2)/(WI(9,2)+WI(9,4))
!
!       PHI(:,:,2)=0.0D0
!       PHI(:,2,2)=1.0D0
!------ PH3 --------------- 
       PHI(3,3,3) = 1.0D0
       PHI(5,3,3) = WI(5,3)/(WI(5,3)+WI(5,1)+WI(5,2))
       PHI(6,3,3) = WI(6,3)/(WI(6,3)+WI(6,2))
       PHI(7,3,3) = WI(7,3)/(WI(7,3)+WI(7,1))
       PHI(8,3,3) = WI(8,3)/(WI(8,3)+WI(8,4)+WI(8,1))
       PHI(9,3,3) = WI(9,3)/(WI(9,3)+WI(9,4)+WI(9,2))
       PHI(10,3,3)= WI(10,3)/(WI(10,3)+WI(10,4))
!
       PHI(3,1,3) = 1.0D0
       PHI(6,1,3) = WI(6,3)/(WI(6,3)+WI(6,1)+WI(6,2))
       PHI(7,1,3) = WI(7,3)/(WI(7,3)+WI(7,1))
       PHI(10,1,3)= WI(10,3)/(WI(10,3)+WI(10,1)+WI(10,4))
!
       PHI(3,2,3) = 1.0D0
       PHI(6,2,3) = WI(6,3)/(WI(6,3)+WI(6,2))
       PHI(7,2,3) = WI(7,3)/(WI(7,3)+WI(7,2)+WI(7,1))
       PHI(10,2,3)= WI(10,3)/(WI(10,3)+WI(10,2)+WI(10,4))
!
       PHI(3,4,3) = 1.0D0
       PHI(6,4,3) = WI(6,3)/(WI(6,3)+WI(6,4)+WI(6,2))
       PHI(7,4,3) = WI(7,3)/(WI(7,3)+WI(7,4)+WI(7,1))
       PHI(10,4,3)= WI(10,3)/(WI(10,3)+WI(10,4))
!!
!       PHI(:,:,3)=0.0D0
!       PHI(:,3,3)=1.0D0
!------ PH4 --------------- 
!       PHI(:,:,4)=0.0D0
!       PHI(:,4,4)=1.0D0
       PHI(:,1,4)=1.0D0-PHI(:,1,1)-PHI(:,1,2)-PHI(:,1,3)
       PHI(:,2,4)=1.0D0-PHI(:,2,1)-PHI(:,2,2)-PHI(:,2,3)
       PHI(:,3,4)=1.0D0-PHI(:,3,1)-PHI(:,3,2)-PHI(:,3,3)
       PHI(:,4,4)=1.0D0-PHI(:,4,1)-PHI(:,4,2)-PHI(:,4,3)
!---------------------------------------------------------
!
!------------------------------------------
       J11=XYZ(1,1)
       J12=XYZ(1,2)
       J13=XYZ(1,3)
       J14=XYZ(1,4)
       J21=XYZ(2,1)
       J22=XYZ(2,2)
       J23=XYZ(2,3)
       J24=XYZ(2,4)
       J31=XYZ(3,1)
       J32=XYZ(3,2)
       J33=XYZ(3,3)
       J34=XYZ(3,4)
!
       DETJ=(J13*J22-J12*J23+J14*J23-J14*J22+J12*J24-J13*J24)*J31
       DETJ=DETJ-(J13*J21-J11*J23+J14*J23-J14*J21+J11*J24-J13*J24)*J32
       DETJ=DETJ+(J12*J21-J11*J22+J14*J22-J14*J21+J11*J24-J12*J24)*J33
       DETJ=DETJ-(J12*J21-J11*J22+J13*J22-J13*J21+J11*J23-J12*J23)*J34
! 
       JINV(1,1)=J22*(J34-J33)-J23*(J34-J32)+J24*(J33-J32)
       JINV(1,2)=-J12*(J34-J33)+J13*(J34-J32)-J14*(J33-J32)
       JINV(1,3)=J12*(J24-J23)-J13*(J24-J22)+J14*(J23-J22)
       JINV(2,1)=-J21*(J34-J33)+J23*(J34-J31)-J24*(J33-J31)
       JINV(2,2)=J11*(J34-J33)-J13*(J34-J31)+J14*(J33-J31)
       JINV(2,3)=-J11*(J24-J23)+J13*(J24-J21)-J14*(J23-J21)
       JINV(3,1)=J21*(J34-J32)-J22*(J34-J31)+J24*(J32-J31)
       JINV(3,2)=-J11*(J34-J32)+J12*(J34-J31)-J14*(J32-J31)
       JINV(3,3)=J11*(J24-J22)-J12*(J24-J21)+J14*(J22-J21)
       JINV(4,1)=-J21*(J33-J32)+J22*(J33-J31)-J23*(J32-J31)
       JINV(4,2)=J11*(J33-J32)-J12*(J33-J31)+J13*(J32-J31)
       JINV(4,3)=-J11*(J23-J22)+J12*(J23-J21)-J13*(J22-J21)
!      CALCULATE DHIX,DHIY
       DHIX(:,1)=JINV(1,1)
       DHIX(:,2)=JINV(2,1)
       DHIX(:,3)=JINV(3,1)
       DHIX(:,4)=JINV(4,1)
!
       DHIY(:,1)=JINV(1,2)
       DHIY(:,2)=JINV(2,2)
       DHIY(:,3)=JINV(3,2)
       DHIY(:,4)=JINV(4,2)
!
       DHIZ(:,1)=JINV(1,3)
       DHIZ(:,2)=JINV(2,3)
       DHIZ(:,3)=JINV(3,3)
       DHIZ(:,4)=JINV(4,3)
!      CALCULATE DHIX2,DHIY2 (node 1-4)
       DHIX2(:,1)=HI_ASS(:,1)*JINV(1,1)
       DHIX2(:,2)=HI_ASS(:,2)*JINV(2,1)
       DHIX2(:,3)=HI_ASS(:,3)*JINV(3,1)
       DHIX2(:,4)=HI_ASS(:,4)*JINV(4,1)
!
       DHIY2(:,1)=HI_ASS(:,1)*JINV(1,2)
       DHIY2(:,2)=HI_ASS(:,2)*JINV(2,2)
       DHIY2(:,3)=HI_ASS(:,3)*JINV(3,2)
       DHIY2(:,4)=HI_ASS(:,4)*JINV(4,2)
!
       DHIZ2(:,1)=HI_ASS(:,1)*JINV(1,3)
       DHIZ2(:,2)=HI_ASS(:,2)*JINV(2,3)
       DHIZ2(:,3)=HI_ASS(:,3)*JINV(3,3)
       DHIZ2(:,4)=HI_ASS(:,4)*JINV(4,3)
!      CALCULATE DHIX2,DHIY2 (node 5-10)
       DHIX2(:,5) =4.0D0*(JINV(1,1)*HI(:,2)+JINV(2,1)*HI(:,1))
       DHIX2(:,6) =4.0D0*(JINV(2,1)*HI(:,3)+JINV(3,1)*HI(:,2))
       DHIX2(:,7) =4.0D0*(JINV(3,1)*HI(:,1)+JINV(1,1)*HI(:,3))
       DHIX2(:,8) =4.0D0*(JINV(1,1)*HI(:,4)+JINV(4,1)*HI(:,1))
       DHIX2(:,9) =4.0D0*(JINV(2,1)*HI(:,4)+JINV(4,1)*HI(:,2))
       DHIX2(:,10)=4.0D0*(JINV(3,1)*HI(:,4)+JINV(4,1)*HI(:,3))
!
       DHIY2(:,5) =4.0D0*(JINV(1,2)*HI(:,2)+JINV(2,2)*HI(:,1))
       DHIY2(:,6) =4.0D0*(JINV(2,2)*HI(:,3)+JINV(3,2)*HI(:,2))
       DHIY2(:,7) =4.0D0*(JINV(3,2)*HI(:,1)+JINV(1,2)*HI(:,3))
       DHIY2(:,8) =4.0D0*(JINV(1,2)*HI(:,4)+JINV(4,2)*HI(:,1))
       DHIY2(:,9) =4.0D0*(JINV(2,2)*HI(:,4)+JINV(4,2)*HI(:,2))
       DHIY2(:,10)=4.0D0*(JINV(3,2)*HI(:,4)+JINV(4,2)*HI(:,3))
!
       DHIZ2(:,5) =4.0D0*(JINV(1,3)*HI(:,2)+JINV(2,3)*HI(:,1))
       DHIZ2(:,6) =4.0D0*(JINV(2,3)*HI(:,3)+JINV(3,3)*HI(:,2))
       DHIZ2(:,7) =4.0D0*(JINV(3,3)*HI(:,1)+JINV(1,3)*HI(:,3))
       DHIZ2(:,8) =4.0D0*(JINV(1,3)*HI(:,4)+JINV(4,3)*HI(:,1))
       DHIZ2(:,9) =4.0D0*(JINV(2,3)*HI(:,4)+JINV(4,3)*HI(:,2))
       DHIZ2(:,10)=4.0D0*(JINV(3,3)*HI(:,4)+JINV(4,3)*HI(:,3))
!  
       DO I=1,4
          HI2_PHI(:,:,I)  =MATMUL(HI2(:,:),PHI(:,:,I))
          DXHI2_PHI(:,:,I)=MATMUL(DHIX2(:,:),PHI(:,:,I))
          DYHI2_PHI(:,:,I)=MATMUL(DHIY2(:,:),PHI(:,:,I)) 
          DZHI2_PHI(:,:,I)=MATMUL(DHIZ2(:,:),PHI(:,:,I)) 
       ENDDO      
!
       DO I=1,4
          FI(:,I)  =SUM(HI(:,:)*HI2_PHI(:,:,I),DIM=2)
          DFIX(:,I)=SUM(DHIX(:,:)*HI2_PHI(:,:,I),DIM=2)+SUM(HI(:,:)*DXHI2_PHI(:,:,I),DIM=2)
          DFIY(:,I)=SUM(DHIY(:,:)*HI2_PHI(:,:,I),DIM=2)+SUM(HI(:,:)*DYHI2_PHI(:,:,I),DIM=2)
          DFIZ(:,I)=SUM(DHIZ(:,:)*HI2_PHI(:,:,I),DIM=2)+SUM(HI(:,:)*DZHI2_PHI(:,:,I),DIM=2)
       ENDDO
!
       RETURN 
!}}}
    END SUBROUTINE Get_ShapeFun_Derivative
!
    SUBROUTINE OFE_TETRA_ADDBAN(ND,LMM,S)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble upper triangular element stiffness into             .
! .   compacted global stiffness                                      .
! .                                                                   .
! .         A = GLOBAL STIFFNESS                                      .
! .         S = ELEMENT STIFFNESS                                     .
! .         ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS              .
! .                                                                   .
! .                   S(1)        S(2)        S(3)        . . .       .
! .         S   =                 S(ND+1)     S(ND+2)     . . .       .
! .                                           S(2*ND)     . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! .                   A(1)        A(3)        A(6)        . . .       .
! .         A   =                 A(2)        A(5)        . . .       .
! .                                           A(4)        . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: ND
   INTEGER,DIMENSION(ND),INTENT(IN) :: LMM
   REAL(8),DIMENSION(ND,ND),INTENT(IN) :: S
!
   INTEGER :: I,J,L,M,K
   INTEGER :: II,JJ
!
   DO I=1,ND
      L=LMM(I)
      IF(L.GT.0) THEN
        DO J=1,ND
           M=LMM(J)
           IF(M.GT.0) THEN
              K=M-L
              IF(K.GE.0) THEN
                K=MAXA_COLSOL(M)+K
                K_COLSOL(K)=K_COLSOL(K)+S(I,J)
                KMATRIX(K)=KMATRIX(K)+S(I,J)
              END IF
           ENDIF
        ENDDO
      ENDIF
   ENDDO

!
   RETURN
!}}}
    END SUBROUTINE OFE_TETRA_ADDBAN
!
    SUBROUTINE OLR_3D_ADDBAN(ND_NM,NN)
!{{{
    IMPLICIT NONE
    INTEGER :: NN
    INTEGER,DIMENSION(NN),INTENT(INOUT)  :: ND_NM
!   IN-FUNCTION VARIABLES
    INTEGER :: I,I1,I2,I3,J,J1,J2,J3
    INTEGER :: K_TRI,KK
    INTEGER :: NDI,NDJ
    
!   SORT THE ARRAY 
    CALL SORT(ND_NM,NN)
!   LOOP COLUMN
    DO I=1,NN
       NDI=ND_NM(I)
       DO I1=1,3
          DO I2=1,NDOF(I1,NDI)
             I3=ID(I1,NDI)+I2-1
             K_TRI=MAXA_COLSOL(I3)
!   LOOP ROW 
             DO J=1,I-1
                NDJ=ND_NM(J)
                DO J1=1,3
                   DO J2=1,NDOF(J1,NDJ)
                      J3=ID(J1,NDJ)+J2-1
!                      write(*,*) NDI,NDJ,I3,J3,"***"
                      KK=K_TRI+I3-J3
                      K_COLSOL(KK)=1
                   ENDDO
                ENDDO
             ENDDO
!
                NDJ=ND_NM(I)
                DO J1=1,I1-1
                   DO J2=1,NDOF(J1,NDJ)
                      J3=ID(J1,NDJ)+J2-1
!                      write(*,*) NDI,NDJ,I3,J3,"***"
                      KK=K_TRI+I3-J3
                      K_COLSOL(KK)=1
                   ENDDO
                ENDDO
                J1=I1
                   DO J2=1,I2
                      J3=ID(J1,NDJ)+J2-1
!                      write(*,*) NDI,NDJ,I3,J3,"***"
                      KK=K_TRI+I3-J3
                      K_COLSOL(KK)=1
                   ENDDO
          ENDDO
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE  OLR_3D_ADDBAN
!
    SUBROUTINE TRUSS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
    INTEGER :: NumEle,Mtype   
!   VARIABLES FOR K & M
    INTEGER :: N,L,KL,K
    REAL(8) :: XL2,XL,D(3),XX,YY
    REAL(8) :: ST(6)
    REAL(8) :: S(21)
    REAL(8) :: SMASS(21)
!   ASSISTANT VARIABLE FOR K*U
    REAL(8) :: FORCE
!   ASSISTANT VARIABLES
    INTEGER :: I,J
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   READ 
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,6
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
!--------------------------------------------------
    READ(ADELE) ((XYZ(I,J),I=1,6),J=1,NumEle),((LM(I,J),I=1,6),J=1,NumEle),(ELEPAR(1,I),I=1,NumEle)
!   GENERATE
!   STIFFNESS MATRIX
    DO N=1,NumEle
       AREA=ELEPAR(1,N)
       XL2=0.
       DO L=1,3
          D(L)=XYZ(L,N) - XYZ(L+3,N)
          XL2=XL2 + D(L)*D(L)
       END DO
       XL=SQRT(XL2)
       XX=E*AREA*XL
       DO L=1,3
          ST(L)=D(L)/XL2
          ST(L+3)=-ST(L)
       END DO
!
       KL=0
       DO L=1,6
          YY=ST(L)*XX
          DO K=L,6
             KL=KL + 1
             S(KL)=ST(K)*YY                
          END DO
       END DO  
!  
       CALL ADDBAN(S(1),21,LM(1,N),6,1)
!   MASS MATRIX: CONSISTENT MASS 
       DO L=1,21
          SMASS(L)=0
       END DO
       SMASS(1)  = AREA*DEN*XL/3
       SMASS(4)  = AREA*DEN*XL/6
       SMASS(7)  = AREA*DEN*XL/3
       SMASS(10) = AREA*DEN*XL/6
       SMASS(12) = AREA*DEN*XL/3
       SMASS(15) = AREA*DEN*XL/6
       SMASS(16) = AREA*DEN*XL/3
       SMASS(19) = AREA*DEN*XL/3
       SMASS(21) = AREA*DEN*XL/3
       CALL ADDBAN(SMASS(1),21,LM(1,N),6,2)
!   CALCULATE K*U
       IF(DYNAMIC.EQ.1) THEN
         FORCE=0        
         DO L=1,3
            I=LM(L,N)
            IF ((I.GT.0).AND.(U(I).NE.0)) FORCE = FORCE + ST(L)*U(I)*E*AREA
            J=LM(L+3,N)
            IF ((J.GT.0).AND.(U(J).NE.0)) FORCE = FORCE + ST(L+3)*U(J)*E*AREA
         END DO 
         IF(FORCE.NE.0) THEN     
           DO L=1,3
              I=LM(L,N)
              IF(I.GT.0) RODE(I)=RODE(I)+FORCE*D(L)/XL
              J=LM(L+3,N)
              IF(J.GT.0) RODE(J)=RODE(J)-FORCE*D(L)/XL
           END DO
         END IF
       END IF
    END DO
!
    DEALLOCATE(XYZ)
    DEALLOCATE(LM)
    DEALLOCATE(ELEPAR)
    RETURN
!}}}
    END SUBROUTINE TRUSS_G
!
    SUBROUTINE QUADS_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
    INTEGER :: ITYPE,Mtype,NumEle,NINTG  
!   VARIABLES FOR K & M
    INTEGER :: N
    REAL(8) :: S(8,8)=0
    REAL(8) :: SMASS(8,8)=0
!   ASSISTANT VARIABLE FOR K*U
    REAL(8) :: FORCE
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,KK
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: ALLocateStatus
!
    REAL(8), DIMENSION(2,4) :: XX
!   READ 
    ITYPE  = NPAR(2)
    Mtype  = NPAR(3)
    NumEle = NPAR(4)
    NINTG  = NPAR(5)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!**************************
!   CONSTRUCT LOAD VECTOR 
!**************************！    
    IF (QUADSBCONOFF) THEN
       NITP_BC=2
       ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"     
!
       CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!
!   Initialization LOAD
!
       READ(ADBC) (NPAR(J),J=1,7)
       IF(NPAR(1).NE.3) THEN
          WRITE(*,*) "ERROR, QUADS_LINE_G"
          STOP
       ENDIF
!
       CALL QUADS_LINE_G(ADBC,THIC)
! -----------------------------------
       DEALLOCATE(XABSC_WEIG_BC)  
    ENDIF
!   ALLOCATE VARIABLES
    ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"  
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,8
              XYZ(J,I)=0
              LM(J,I)=0
           END DO
           ELEPAR(1,I)=0.
        END DO
!--------------------------------------------------
    READ(ADELE) ((XYZ(I,J),I=1,8),J=1,NumEle),((LM(I,J),I=1,8),J=1,NumEle),(ELEPAR(1,I),I=1,NumEle)
!   GENERATE
!   STIFFNESS MATRIX
    DO N=1,NumEle
!      THIC=ELEPAR(1,N)
!      write(*,*) THIC
       DO I=1,8
          DO J=1,8
             S(I,J)=0.
             SMASS(I,J)=0.
          ENDDO
       ENDDO
       DO I=1,4
          XX(1,I)=XYZ(2*I-1,N)
          XX(2,I)=XYZ(2*I,N)
       ENDDO
       CALL QUADS(NumEle,ITYPE,NINTG,THIC,E,PR,XX,S,SMASS,DEN)     
       CALL ADDBAN_ZLB(8,S,LM(:,N),1)  
!   MASS MATRIX: LUMPED MASS  
!       DO I=1,8
!          DO J=1,8
!             IF (J.NE.I) THEN
!               SMASS(I,I)=SMASS(I,I)+SMASS(I,J)
!               SMASS(I,J)=0
!             ENDIF
!          ENDDO
!       ENDDO
!       CALL ADDBAN_ZLB(8,SMASS,LM(:,N),2)
!   CALCULATE K*U
!       IF(DYNAMIC.EQ.1) THEN
!          DO I=1,8
!             FORCE=0
!             KK=LM(I,N)
!             IF(KK.GT.0) THEN
!               DO J=1,8
!                 K=LM(J,N)
!                 IF(K.GT.0) THEN
!                    FORCE=FORCE+S(I,J)*U(K)
!                 ENDIF
!               ENDDO
!               RODE(KK)=RODE(KK)+FORCE
!             ENDIF
!          ENDDO 
!       END IF
    END DO
!
    DEALLOCATE(XYZ)
    DEALLOCATE(LM)
    DEALLOCATE(ELEPAR)
    RETURN
!}}}
    END SUBROUTINE QUADS_G
!
    SUBROUTINE QUADS_LINE_G(ADBC,THIC)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC
    REAL(8),INTENT(IN) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    INTEGER :: N1,N2
    REAL(8) :: LENGTH,NX,NY,NXX,NYY
    REAL(8) :: XM,XL,YM,YL
    REAL(8),DIMENSION(2)         :: FLOAD
    INTEGER,DIMENSION(2)         :: FUNC
    REAL(8),DIMENSION(NITP_BC)   :: HI,XTMP,YTMP
    REAL(8),DIMENSION(2)         :: LOADD
    REAL(8),DIMENSION(NITP_BC,2) :: LOAD
    REAL(8),DIMENSION(1,2)       :: SI,SJ
    INTEGER :: NDOFF
    REAL(8) :: CONST1
    INTEGER :: NUM1,NUM2,JJ2
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCED ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCIND ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCPAR ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
!   CALCULATE LOAD VECTOR
!
    DO I=1,NumBc
       N1=BCED(1,I)
       N2=BCED(2,I)
!
       XM=(X(N1)+X(N2))/2.0D0
       XL=(X(N2)-X(N1))/2.0D0
       YM=(Y(N1)+Y(N2))/2.0D0
       YL=(Y(N2)-Y(N1))/2.0D0
!
       LENGTH=XL*XL+YL*YL
       LENGTH=2.0D0*SQRT(LENGTH)
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1) =BCIND(3,I)
       FUNC(2) =BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!----- XTMP YTMP -----------------------
       XTMP=XM+XL*XABSC_WEIG_BC(:,1)
       YTMP=YM+YL*XABSC_WEIG_BC(:,1)
!----- HIW  HJW ------------------------
       DO J=1,NITP_BC
          HI(J)=(-XABSC_WEIG_BC(J,1)+1.0D0)/2.0D0
!-----   LOAD  ---------------------------
          LOADD= F_LOAD(FLOAD,FUNC,NX,NY,XTMP(J),YTMP(J),N1,N2)
          LOAD(J,1)=LOADD(1)*XABSC_WEIG_BC(J,2)
          LOAD(J,2)=LOADD(2)*XABSC_WEIG_BC(J,2)
       ENDDO
!
       CONST1=LENGTH*THIC/2.0D0
! ------ PmI --------------
       SI(1,1)=SUM(HI*LOAD(:,1))*CONST1
       SI(1,2)=SUM(HI*LOAD(:,2))*CONST1
! ------- PmJ --------------------
       SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
       SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!   
       DO J=1,2
          IF(BCIND(J,I).EQ.1) CYCLE
!
          NUM1=ID(J,N1)
          R(NUM1)=R(NUM1)+SI(1,J)
!--------------------------------------------------
          NUM1=ID(J,N2)
          R(NUM1)=R(NUM1)+SJ(1,J)
       ENDDO
!
    ENDDO      
!------------------------------------------
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE QUADS_LINE_G  
!
    SUBROUTINE QUADS(NEL,ITYPE,NINT,THIC,YM,PR,XX,S,SMASS,DEN)              
!{{{                                                                      
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO CALCULATE ISOPARAMETRIC QUADRILATERAL ELEMENT STIFFNESS . 
! .        MATRIX FOR AXISYMMETRIC, PLANE STRESS, AND PLANE STRAIN    . 
! .        CONDITIONS                                                 . 
! .                                                                   . 
! .  - - INPUT VARIABLES - -                                          . 
! .        NEL       = NUMBER OF ELEMENT                              . 
! .        ITYPE     = ELEMENT TYPE                                   . 
! .                        EQ.0 = AXISYMMETRIC                        . 
! .                        EQ.1 = PLANE STRAIN                        . 
! .                        EQ.2 = PLANE STRESS                        . 
! .        NINT      = GAUSS NUMERICAL INTEGRATION ORDER              . 
! .        THIC      = THICKNESS OF ELEMENT                           . 
! .        YM        = YOUNG'S MODULUS                                . 
! .        PR        = POISSON'S RATIO                                . 
! .        XX(2,4)   = ELEMENT NODE COORDINATES                       . 
! .        S(8,8)    = STORAGE FOR STIFFNESS MATRIX                   . 
! .                                                                   . 
! .  - - OUTPUT - -                                                   . 
! .        S(8,8)    = CALCULATED STIFFNESS MATRIX                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT NONE                               
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . 
! .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . 
! .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . 
! .   SINGLE PRECISION ARITHMETIC.                                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      INTEGER :: NEL,ITYPE,NINT
      REAL(8) :: DEN
      REAL(8) :: YM,PR,XBAR,THIC
      REAL(8) :: F,G,H
      REAL(8) :: A,RI,SI,DET,WT,STIFF
      INTEGER :: I,J,K,IST,LX,LY,L,M
      REAL(8) :: D(4,4),B(4,8),HH(8,8),XX(2,4),S(8,8),SMASS(8,8),XG(4,4),WGT(4,4),DB(4)     
!                                                                       
!     MATRIX XG STORES GAUSS - LEGENDRE SAMPLING POINTS                 
!                                                                       
      DATA XG/   0.D0,   0.D0,   0.D0,   0.D0,   -.5773502691896D0,&    
      .5773502691896D0,   0.D0,   0.D0,   -.7745966692415D0,   0.D0,&   
      .7745966692415D0,   0.D0,   -.8611363115941D0,&                   
      -.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /       
!                                                                       
!     MATRIX WGT STORES GAUSS - LEGENDRE WEIGHTING FACTORS              
!                                                                       
      DATA WGT /  2.D0,   0.D0,   0.D0,   0.D0,   1.D0,   1.D0,&        
      0.D0,   0.D0,   .5555555555556D0,   .8888888888889D0,&            
      .5555555555556D0,   0.D0,   .3478548451375D0,   .6521451548625D0,&
      .6521451548625D0,   .3478548451375D0 /                            
!                                                                       
!     O B T A I N  S T R E S S - S T R A I N  L A W                     
!                                                                       
      F=YM/(1.+PR)                                                      
      G=F*PR/(1.-2.*PR)                                                 
      H=F + G                                                           
!                                                                       
!     PLANE STRAIN ANALYSIS                                             
!                                                                       
      D(1,1)=H                                                          
      D(1,2)=G                                                          
      D(1,3)=0.                                                         
      D(2,1)=G                                                          
      D(2,2)=H                                                          
      D(2,3)=0.                                                         
      D(3,1)=0.                                                         
      D(3,2)=0.                                                         
      D(3,3)=F/2.                                                       
      IF (ITYPE.EQ.1) THEN                                              
      THIC=1.                                                           
      GO TO 20                                                          
      ENDIF                                                             
!                                                                       
!     AXISYMMETRIC ANALYSIS                                             
!                                                                       
      D(1,4)=G                                                          
      D(2,4)=G                                                          
      D(3,4)=0.                                                         
      D(4,1)=G                                                          
      D(4,2)=G                                                          
      D(4,3)=0.                                                         
      D(4,4)=H                                                          
      IF (ITYPE.EQ.0) GO TO 20                                          
!                                                                       
!     FOR PLANE STRESS ANALYSIS CONDENSE STRESS-STRAIN MATRIX           
!                                                                       
      DO 10 I=1,3                                                       
      A=D(I,4)/D(4,4)                                                   
      DO 10 J=I,3                                                       
      D(I,J)=D(I,J) - D(4,J)*A                                          
   10 D(J,I)=D(I,J)                                                     
!                                                                       
!     C A L C U L A T E  E L E M E N T  S T I F F N E S S               
!                                                                       
   20 DO 30 I=1,8                                                       
      DO 30 J=1,8                                                       
      S(I,J)=0.
   30 SMASS(I,J)=0. 
!   
      IST=3                                                             
      IF (ITYPE.EQ.0) IST=4                                             
      DO 80 LX=1,NINT                                                   
      RI=XG(LX,NINT)                                                    
      DO 80 LY=1,NINT                                                   
      SI=XG(LY,NINT)                                                    
!                                                                       
!     EVALUATE DERIVATIVE OPERATOR B AND THE JACOBIAN DETERMINANT DET   
!                                                                       
      CALL STDM (XX,B,HH,DET,RI,SI,XBAR,NEL,ITYPE)                         
!                                                                       
!     ADD CONTRIBUTION TO ELEMENT STIFFNESS                             
!                                                                       
      IF (ITYPE.GT.0) XBAR=THIC                                         
      WT=WGT(LX,NINT)*WGT(LY,NINT)*XBAR*DET                             
      DO 70 J=1,8                                                       
      DO 40 K=1,IST                                                     
      DB(K)=0.0                                                         
      DO 40 L=1,IST                                                     
   40 DB(K)=DB(K) + D(K,L)*B(L,J)                                       
      DO 60 I=J,8                                                       
      STIFF=0.0                                                         
      DO 50 L=1,IST                                                     
   50 STIFF=STIFF + B(L,I)*DB(L)                                        
   60 S(I,J)=S(I,J) + STIFF*WT                                          
   70 CONTINUE
       DO L=1,8
         DO M=1,8
            SMASS(L,M)=SMASS(L,M)+DEN*HH(L,M)*WT
         ENDDO
       ENDDO
   80 CONTINUE                                                          
!
      DO J=1,8                                                       
         DO I=J,8                                                       
           S(J,I)=S(I,J)                                                     
         ENDDO
      ENDDO
!   LUMPED MASS MATRIX
      RETURN                                                            
!                                                                       
      END                                                               
      SUBROUTINE STDM (XX,B,HH,DET,R,S,XBAR,NEL,ITYPE)                     
!                                                                       
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .     TO EVALUATE THE STRAIN-DISPLACEMENT TRANSFORMATION MATRIX B   . 
! .     AT POINT (R,S) FOR A QUADRILATERAL ELEMENT                    . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT NONE
      INTEGER :: NEL,ITYPE
      INTEGER :: I,J,K2,K
      REAL(8) :: R,S,RP,SP,RM,SM
      REAL(8) :: DET,XBAR,DUM     
      REAL(8) :: XX(2,4),B(4,8),HH(8,8),H(4),P(2,4),XJ(2,2),XJI(2,2)             
!                                                                       
      RP = 1.0 + R                                                      
      SP = 1.0 + S                                                      
      RM = 1.0 - R                                                      
      SM = 1.0 - S                                                      
!                                                                       
!     INTERPOLATION FUNCTIONS                                           
!                                                                       
      H(1) = 0.25* RP* SP                                               
      H(2) = 0.25* RM* SP                                               
      H(3) = 0.25* RM* SM                                               
      H(4) = 0.25* RP* SM                                               
!
      HH=0.
      DO I=1,8
         If(mod(I,2) == 0) then
           J=I/2
           HH(I,2)=H(J)*H(1)
           HH(I,4)=H(J)*H(2)
           HH(I,6)=H(J)*H(3)
           HH(I,8)=H(J)*H(4)
         ELSE
           J=(I+1)/2
           HH(I,1)=H(J)*H(1)
           HH(I,3)=H(J)*H(2)
           HH(I,5)=H(J)*H(3)
           HH(I,7)=H(J)*H(4)
         ENDIF
      ENDDO
! 
!     NATURAL COORDINATE DERIVATIVES OF THE INTERPOLATION FUNCTIONS     
!                                                                       
!        1. WITH RESPECT TO R                                           
!                                                                       
      P(1,1) = 0.25* SP                                                 
      P(1,2) = - P(1,1)                                                 
      P(1,3) = - 0.25* SM                                               
      P(1,4) = - P(1,3)                                                 
!                                                                       
!        2. WITH RESPECT TO S                                           
!                                                                       
      P(2,1) = 0.25* RP                                                 
      P(2,2) = 0.25* RM                                                 
      P(2,3) = - P(2,2)                                                 
      P(2,4) = - P(2,1)                                                 
!                                                                       
!     EVALUATE THE JACOBIAN MATRIX AT POINT (R,S)                       
!                                                                       
   10 DO 30 I=1,2                                                       
      DO 30 J=1,2                                                       
      DUM = 0.0                                                         
      DO 20 K=1,4                                                       
   20 DUM=DUM + P(I,K)*XX(J,K)                                          
   30 XJ(I,J)=DUM                                                       
!                                                                       
!     COMPUTE THE DETERMINANT OF THE JACOBIAN MATRIX AT POINT (R,S)     
!                                                                       
      DET = XJ(1,1)* XJ(2,2) - XJ(2,1)* XJ(1,2)                         
      IF (DET.GT.0.00000001) GO TO 40                                   
      WRITE (*,2000) NEL                                                
      GO TO 800                                                         
!                                                                       
!     COMPUTE INVERSE OF THE JACOBIAN MATRIX                            
!                                                                       
   40 DUM=1./DET                                                        
      XJI(1,1) = XJ(2,2)* DUM                                           
      XJI(1,2) =-XJ(1,2)* DUM                                           
      XJI(2,1) =-XJ(2,1)* DUM                                           
      XJI(2,2) = XJ(1,1)* DUM                                           
!                                                                       
!     EVALUATE GLOBAL DERIVATIVE OPERATOR B                             
!                                                                       
      K2=0                                                              
      DO 60 K=1,4                                                       
      K2=K2 + 2                                                         
      B(1,K2-1) = 0.                                                    
      B(1,K2  ) = 0.                                                    
      B(2,K2-1) = 0.                                                    
      B(2,K2  ) = 0.                                                    
      DO 50 I=1,2                                                       
      B(1,K2-1) = B(1,K2-1) + XJI(1,I) * P(I,K)                         
   50 B(2,K2  ) = B(2,K2  ) + XJI(2,I) * P(I,K)                         
      B(3,K2  ) = B(1,K2-1)                                             
   60 B(3,K2-1) = B(2,K2  )                                             
!                                                                       
!     IN CASE OF PLANE STRAIN OR PLANE STRESS ANALYSIS DO NOT INCLUDE   
!     THE NORMAL STRAIN COMPONENT                                       
!                                                                       
      IF (ITYPE.GT.0) GO TO 900                                         
!                                                                       
!     COMPUTE THE RADIUS AT POINT (R,S)                                 
!                                                                       
      XBAR=0.0                                                          
      DO 70 K=1,4                                                       
   70 XBAR=XBAR + H(K)*XX(1,K)                                          
!                                                                       
!     EVALUATE THE HOOP STRAIN-DISPLACEMENT RELATION                    
!                                                                       
      IF (XBAR.GT.0.00000001) GO TO 90                                  
!                                                                       
!     FOR THE CASE OF ZERO RADIUS EQUATE RADIAL TO HOOP STRAIN          
!                                                                       
      DO 80 K=1,8                                                       
   80 B(4,K)=B(1,K)                                                     
      GO TO 900                                                         
!                                                                       
!     NON-ZERO RADIUS                                                   
!                                                                       
   90 DUM=1./XBAR                                                       
      K2=0                                                              
      DO 100 K=1,4                                                      
      K2=K2 + 2                                                         
      B(4,K2  ) = 0.                                                    
  100 B(4,K2-1) = H(K)*DUM                                              
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 2000 FORMAT (//,' *** ERROR *** ',&                                     
         ' ZERO OR NEGATIVE JACOBIAN DETERMINANT FOR ELEMENT (',I8,')')
!}}}                                                                     
    END 
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
!  
!     SUBROUTINES FOR TRUSS_MFS
!
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!{{{
!
    SUBROUTINE TRUSS_MFS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE      
!   ASSISTANT VARIABLES
    INTEGER :: N,M
    INTEGER :: I,J,K,L,II,JJ
    INTEGER :: NumEle,Mtype,NEleAss
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   ASSISTANT VARIABLES
    INTEGER :: NUM,NUM1
    INTEGER :: INDSP
    INTEGER :: KI
    REAL(8) :: Radius
    REAL(8) :: S        
    REAL(8) :: HI(NPOR_MAX),DHI(NPOR_MAX)
    REAL(8) :: XX,YY
    INTEGER :: NPOR
!   READ 
    NEleAss= NPAR(4)
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
    AREA = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
    ALLOCATE (XABSC_WEIG_IN(NITP,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
        DO I=1,NumEle
           DO J=1,3
              NDMFS(J,I)=0
           ENDDO
           ELEPAR(1,I)=0.
           ELEPAR(2,I)=0.
           DO J=1,NEleAss
              LM(J,I)=0
           ENDDO
        END DO
!--------------------------------------------------
    READ(ADELE) ((NDMFS(I,J),I=1,3),J=1,NumEle),((ELEPAR(J,I),J=1,2),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE
    CALL GAULEG(NITP,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    ALLOCATE (STIFF_IJ(1,NPOR_MAX*NPOR_MAX),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
!   STIFFNESS MATRIX
    DO I=1,NumEle
!
!   GENERATE LOAD VECTOR
!
!   Neumann boundary
!
       INDSP=NDMFS(3,I)
       IF(NDMFS(2,I).EQ.2) THEN
          N=NDMFS(1,I)
          KI=ID(1,N)
          XX=X(N)
          HI=TRUSS_MFS_HIHJ(NEleAss,I,NDOF(1,N),XX,INDSP)
          DO JJ=1,NDOF(1,N)
             NUM=KI+JJ-1
             NPOR=JJ
             R(NUM)=RNAT(1,N)*HI(NPOR)*AREA
          ENDDO
       ENDIF
!       
!   Dirichlet boundary
!
!       IF(NDMFS(2,I).EQ.1) THEN
!          N=NDMFS(1,I)
!          KI=ID(1,N)
!          IF(KI.NE.0) THEN
!             XX=X(N)
!             DHI=TRUSS_MFS_HIHJ_V2(NEleAss,I,NDOF(1,N)+1,XX,INDSP)
!             DO JJ=1,NDOF(1,N)
!                NPOR=JJ+1
!                NUM=KI+JJ-1
!                R(NUM) = RESS(1,N)*DHI(NPOR)*AREA*E
!             ENDDO
!          ENDIF
!       ENDIF
!
!   BODY FORCE
!
!       CALL TRUSS_MFS_BODY(I,NEleAss,AREA,INDSP)           
!
!   CALCULATE STIFFNESS MATRIX
!
       DO J=1,NEleAss
          K=LM(J,I)
          IF (K.GE.I) THEN        
             CALL TRUSS_MFS_STIFF(I,K,E,DEN,AREA,NEleAss)
          ENDIF
       ENDDO 
!
    ENDDO              
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(STIFF_IJ)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_G
!
    SUBROUTINE TRUSS_MFS_STIFF(II,JJ,E,DEN,AREA,NEleAss)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE STIFFNESS MATRIX OF TRUSS_MFS                         .
! .          A(NWK)       :  GLOBAL STIFFNESS MATRIX                  .
! .          AMAXA(NEQ+1) :  ADDRESS MATRIX                           .
! .          N,M          :  NODE NUMBER                              .
! .          AREA1,RAD1,AREA2,RAD2                                    .
! .          E            :  YOUNG'S MODULUS                          .
! .          DEN          :  DENSITY                                  .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ
    REAL(8)            :: E,DEN,AREA
    INTEGER            :: NEleAss
!   ASSISTANT VARIABLES
    INTEGER            :: N,M
    INTEGER            :: KI,KJ,KK2
!
    INTEGER            :: I,J,K,K1,L,LL
    INTEGER            :: INDK,INDL
    REAL(8)            :: XM,XL,XL0,XTMP,WEIGHT
    INTEGER            :: IND1,IND2
    LOGICAL            :: IND1_X,IND2_X,EXPN1,EXPN2
    INTEGER            :: NPOR1,NPOR2
    REAL(8)            :: HIHJ(NPOR_MAX*NPOR_MAX)
    REAL(8)            :: HI(NPOR_MAX),HJ(NPOR_MAX),DHI(NPOR_MAX),DHJ(NPOR_MAX)
    REAL(8)            :: S2,CONST2,CONST1,CONST3,DIRNI,DIRNJ
    REAL(8)            :: STIFF(NPOR_MAX*NPOR_MAX)
    INTEGER            :: NUM
    INTEGER            :: INDSPI,INDSPJ
    REAL(8)            :: AREA2
!
    N       =  NDMFS(1,II)
    M       =  NDMFS(1,JJ)
    INDSPI  =  NDMFS(3,II)
    INDSPJ  =  NDMFS(3,JJ)
!
    IND1 = NDMFS(2,II)
    IND2 = NDMFS(2,JJ)
    IND1_X = (IND1.EQ.1)
    IND2_X = (IND2.EQ.1)
!
    EXPN1=(IND1.NE.0).AND.(IND2.NE.0)
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    L=NPOR1*NPOR2 
    DO I=1,L
      HIHJ(I)=0.0D0
    ENDDO 
!
    XL0 = ELEPAR(2,II)
    XL = XL0/NIDIV
    DIRNI=ELEPAR(1,II)
    CONST2=1.0/ELEPAR(2,JJ)
    DO K1=1,NIDIV
       XM = X(N)-XL0+XL*(2*K1-1)
       DO I=1,NITP
          XTMP=XM+XL*XABSC_WEIG_IN(I,1)
          AREA2=1
          S2=ABS(XTMP-X(M))*CONST2
          IF(S2.GT.1) CYCLE
          WEIGHT=XABSC_WEIG_IN(I,2)
          IF(EXPN1) THEN
             CONST3=(XTMP-X(N))*DIRNI
             EXPN2=CONST3.LT.0.0
             IF(EXPN2) CYCLE
   !
             DHI=TRUSS_MFS_HIHJ_V2(NEleAss,II,NPOR1,XTMP,INDSPI)
             DHJ=TRUSS_MFS_HIHJ_V2(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
   !
             DO J=1,NPOR1
               DO K=1,NPOR2
                  LL=(J-1)*NPOR2+K
                  HIHJ(LL)=HIHJ(LL)+DHI(J)*DHJ(K)*WEIGHT*AREA2
               ENDDO
             ENDDO             
          ELSE
   !
             DHI=TRUSS_MFS_HIHJ_V2(NEleAss,II,NPOR1,XTMP,INDSPI)
             DHJ=TRUSS_MFS_HIHJ_V2(NEleAss,JJ,NPOR2,XTMP,INDSPJ)
   !
             DO J=1,NPOR1
               DO K=1,NPOR2
                  LL=(J-1)*NPOR2+K
                  HIHJ(LL)=HIHJ(LL)+DHI(J)*DHJ(K)*WEIGHT*AREA2
               ENDDO
             ENDDO    
          ENDIF
       ENDDO
    ENDDO
!
    CONST1=E*AREA*XL
    DO I=1,L
       STIFF_IJ(1,I)=HIHJ(I)*CONST1
    ENDDO                     
!
!   UPDATE LOAD VECTORS
!
    IF(IND1_X) THEN
        IF(ID(1,M).NE.0) THEN
           DO I=1,NDOF(1,M)
              IF(IND2_X) THEN
                 K=I+1
              ELSE
                 K=I
              ENDIF
              LL=K
              NUM=ID(1,M)+I-1
              R(NUM)=R(NUM)-STIFF_IJ(1,LL)*RESS(1,N)
           ENDDO           
        ENDIF
    ENDIF
!
    IF(IND2_X.AND.(N.NE.M)) THEN
       DO I=1,NDOF(1,N)
          LL=(I-1)*NPOR2+1
          NUM=ID(1,N)+I-1
          R(NUM)=R(NUM)-STIFF_IJ(1,LL)*RESS(1,M)
       ENDDO
    ENDIF 
!
!   ASSEMBLY
!
    K=1
    INDK=1
    IF(IND1_X) THEN
       K=2
       INDK=2
    ENDIF
    DO WHILE(K.LE.NPOR1)
       L=1
       INDL=1
       IF(IND2_X) THEN
          L=2
          INDL=2
       ENDIF
       DO WHILE(L.LE.NPOR2)
          KI=ID(1,N)+K-INDK
          KJ=ID(1,M)+L-INDL
          KK2=(K-1)*NPOR2+L
          IF((N.NE.M).OR.(KJ.GE.KI)) THEN
             IF(KJ.LE.KI) THEN
                LL=MAXA_COLSOL(KI)+ABS(KJ-KI)
             ELSE
                LL=MAXA_COLSOL(KJ)+ABS(KJ-KI)
             ENDIF
             K_COLSOL(LL)=K_COLSOL(LL)+STIFF_IJ(1,KK2)
             KMATRIX(LL)=KMATRIX(LL)+STIFF_IJ(1,KK2)
          ENDIF
         L=L+1
       ENDDO
       K=K+1
    ENDDO          
! 
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_STIFF
!
    SUBROUTINE TRUSS_MFS_BODY(II,NEleAss,AREA,INDSP)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE  BODY FORCE CONTRIBUTION TO LOAD VECTOR               .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER :: II
    INTEGER :: NEleAss
    REAL(8) :: AREA
    INTEGER :: INDSP
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: N,NUM,NPOR,NPORMAX
    REAL(8) :: HI(NPOR_MAX),STIFF(NPOR_MAX)
    REAL(8) :: XM,XL,XL0,XTMP,WEIGHT
    INTEGER :: IND1
    INTEGER :: INC
    LOGICAL :: EXPN2,EXPN3
    REAL(8) :: CONST1,CONST2,DIRN
!
       N   = NDMFS(1,II)
       IND1= NDMFS(2,II)
       IF(IND1.EQ.1) THEN
          NPORMAX=NDOF(1,N)+1
          INC=1
       ELSE
          NPORMAX=NDOF(1,N)
          INC=0
       ENDIF
! 
       DO I=1,NPOR_MAX
          STIFF(I)=0.0D0
       ENDDO
!
       NUM=ID(1,N)
       IF(NUM.NE.0) THEN
         XL0= ELEPAR(2,II)
         XL = XL0/NIDIV
         DIRN=ELEPAR(1,II)
         DO L=1,NIDIV
            XM = X(N)-XL0+XL*(2*L-1)
            DO J=1,NITP
               XTMP = XM+XL*XABSC_WEIG_IN(J,1)
               WEIGHT = XABSC_WEIG_IN(J,2)
  !
               EXPN2=(IND1.NE.0)
               IF(EXPN2) THEN
                  CONST2=(XTMP-X(N))*DIRN
                  EXPN3=CONST2.LT.0.0
                  IF(EXPN3) CYCLE
                  HI=TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XTMP,INDSP)
                  DO K=1,NDOF(1,N)
                     NPOR=K+INC
                     STIFF(K)=STIFF(K)+HI(NPOR)*WEIGHT
                  ENDDO
               ELSE
                  HI=TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XTMP,INDSP)
                  DO K=1,NDOF(1,N)
                     NPOR=K+INC
                     STIFF(K)=STIFF(K)+HI(NPOR)*WEIGHT
                  ENDDO
               ENDIF
  !               
            ENDDO
         ENDDO
       ENDIF
!
       CONST1=AREA*XL
       DO I=1,NDOF(1,N)
          NUM   = ID(1,N)+I-1
          R(NUM)= R(NUM) + STIFF(I)*CONST1
       ENDDO                                 
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_MFS_BODY
!
!}}}
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
    SUBROUTINE QUADS_MFS_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      MASS      MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: NEleAss,NumEle,Mtype,ITYPE
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: BcType
!   ASSISTANT VARIABLES
    INTEGER :: KI
    INTEGER :: ALLocateStatus
    REAL(8) :: TT1,TT2
!   READ 
    NEleAss= NPAR(5)
    NumEle = NPAR(4)
    Mtype  = NPAR(3)
    ITYPE  = NPAR(2)
    NITP_IN= NPAR(6)
    NIDIV_IN=NPAR(7)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(6,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
    ALLOCATE (XABSC_WEIG_IN(NITP_IN,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
    DO I=1,NumEle
       DO J=1,6
          NDMFS(J,I)=0
          ELEPAR(J,I)=0.0
       ENDDO
       NDMFS(7,I)=0
       NDMFS(8,I)=0
    ENDDO
!
    READ(ADELE) ((NDMFS(I,J),I=1,8),J=1,NumEle),((ELEPAR(J,I),J=1,6),I=1,NumEle),((LM(J,I),J=1,NEleAss),I=1,NumEle)
!   GENERATE CONSTITUTIVE MATRIX
    CALL G_KCS(2,ITYPE,E,PR)
!   GENERATE LOAD VECTOR  
    CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   CALL CPU_TIME(TT1)
    CALL GAULEG(NITP_IN,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2)) 
!   CALL CPU_TIME(TT2)
!   WRITE(*,*) TT2-TT1
!
    READ(ADBC) (NPAR(J),J=1,7)
!
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,4),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    BcType=NPAR(1)
    SELECT CASE(BcType)
    CASE(1)
        CALL MFS_LINE_G(ADBC,NEleAss,THIC)
    CASE DEFAULT
        WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
        STOP
    END SELECT
!   GENERATE STIFFNESS MATRIX    
!    DO I=1,NumEle
!       WRITE(*,*) NDMFS(7,I),NDMFS(8,I)
!       WRITE(*,*) "******"
!       WRITE(*,*) ELEPAR(5,I),ELEPAR(6,I)
!    ENDDO
!
    DO I=1,NumEle
!   CALCULATE STIFFNESS MATRIX
       DO J=1,NEleAss
          K=LM(J,I)
          IF (K.GE.I) THEN
!   CALCULATE ENTRIES OF STIFFNESS MATRIX  
             CALL QUADS_MFS_STIFF_V2(I,K,NEleAss,THIC)
          ENDIF
       ENDDO              
    ENDDO
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM)
    DEALLOCATE(XABSC_WEIG_BC)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(STIFF_IJ)
    DEALLOCATE(KCS)
!    
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_G   
!
    SUBROUTINE COUP1_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: N
    INTEGER :: II,JJ
    INTEGER :: NumEle,Mtype,ITYPE
    REAL(8) :: E,PR,DEN,THIC
!   ASSISTANT VARIABLES
!    INTEGER :: KI
    INTEGER :: ALLocateStatus
!   READ 
    NumEle = NPAR(4)
    Mtype  = NPAR(3)
    ITYPE  = NPAR(2)
    NITP_IN= NPAR(5)
    NIDIV_IN=NPAR(6)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   ALLOCATE VARIABLES
    ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    ALLOCATE (XABSC_WEIG_IN(NITP_IN,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
    DO I=1,NumEle
       DO J=1,8
          NDMFS(J,I)=0
       ENDDO
       DO J=1,4
          ELEPAR(J,I)=0.0
       ENDDO
    ENDDO
!
    READ(ADELE) ((NDMFS(I,J),I=1,8),J=1,NumEle),((ELEPAR(J,I),J=1,4),I=1,NumEle)
!   GENERATE CONSTITUTIVE MATRIX
    CALL G_KCS(2,ITYPE,E,PR)
!   GENERATE LOAD VECTOR  
    CALL GAULEG(NITP_IN,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2)) 
!
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,4),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   GENERATE STIFFNESS MATRIX   
    DO N=1,NumEle
       DO J=1,4
          DO K=J,4
             CALL COUP1_STIFF(J,K,N,THIC)
          ENDDO
       ENDDO   
    END DO
!
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(STIFF_IJ)
    DEALLOCATE(KCS)
!    
    RETURN
!}}}
    END SUBROUTINE COUP1_G 
!
    SUBROUTINE TRUSS_OVLP_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumInt,NumEle,NEleAss,NITG
    REAL(8) :: E
    INTEGER :: I,J,K,L
    INTEGER :: N,ELE,IND,IDD
    INTEGER :: ALLocateStatus
!
    Mtype  = NPAR(2)
    NumInt = NPAR(3)
    NumEle = NPAR(4)
    NEleAss= NPAR(5)
    NITG   = NPAR(6) 
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
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
!
    ALLOCATE (XABSC_WEIG_IN(NITG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       INTPAR(1,I)=0.0D0
    END DO
!
    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
          WEIND(J,I)=0.0D0
       END DO
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!   READ ELEMENT DATA
    READ(ADELE) ((NDMFS(I,J),I=1,2),J=1,NPAR(4)),(ELEPAR(1,I),I=1,NPAR(4)),((LM(I,J),I=1,NPAR(5)),J=1,NPAR(4)),&
    ((WEIND(I,J),I=1,NPAR(5)),J=1,NPAR(4)),((NDINT(I,J),I=1,2),J=1,NPAR(3)),(INTPAR(1,I),I=1,NPAR(3))    
!
    CALL GAULEG(NITG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    ALLOCATE (DHI_OVLP_TRUSS(NITG,NPOR_MAX,2*NEleAss),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (HI_OVLP_TRUSS(NITG,NPOR_MAX,2*NEleAss),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,1),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!   CALCULATE THE LOAD VECTOR
    DO I=1,NumEle
       IND=NDMFS(2,I)
       IF(IND.EQ.2) THEN
          N=NDMFS(1,I)
          IDD=ID(1,N)
          R(IDD)=1.0D0  
       ENDIF
    ENDDO
!
!   GENERATE STIFFNESS MATRIX   
!
    DO N=1,NumInt
       DO I=1,2*NEleAss
          DO J=1,NPOR_MAX
             DO K=1,NITG
                DHI_OVLP_TRUSS(K,J,I)=0.0D0
                HI_OVLP_TRUSS(K,J,I)=0.0D0
             ENDDO
          ENDDO
       ENDDO
       CALL TRUSS_OVLP_STIFF(NITG,NEleAss,N,E)
    END DO
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM) 
!
    DEALLOCATE(WEIND)
    DEALLOCATE(DHI_OVLP_TRUSS)
    DEALLOCATE(HI_OVLP_TRUSS)
    DEALLOCATE(STIFF_IJ)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_G
!
    SUBROUTINE OFE_TRUSS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumOlp,NITG,NWEI,NDIV
    REAL(8) :: E,AREA
    INTEGER :: I,J
    INTEGER :: NI,NJ,INDI,INDJ,IDD      
    REAL(8) :: FLOAD
    INTEGER :: ALLocateStatus
!
    Mtype  = NPAR(2)
    NumOlp = NPAR(3)
    NITG   = NPAR(4)
    NWEI   = NPAR(5)
    NDIV   = NPAR(6)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    AREA = MATP(4,Mtype)
!-------- Allocate memory
    ALLOCATE (NDINT(4,NumOlp),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(2,NumOlp),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
    ALLOCATE (XABSC_WEIG_IN(NITG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   INITIALIZATION
!   READ ELEMENT DATA
    READ(ADELE) ((NDINT(I,J),I=1,4),J=1,NumOlp),((INTPAR(I,J),I=1,2),J=1,NumOlp)
!
    CALL GAULEG(NITG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!   CALCULATE THE LOAD VECTOR
!    FLOAD=1000D0
!    DO I=1,NumOlp
!       NI  = NDINT(1,I)
!       NJ  = NDINT(2,I)
!       INDI= NDINT(3,I)
!       INDJ= NDINT(4,I)
!       IF(INDI.EQ.2) THEN
!          IDD=ID(1,NI)
!          R(IDD)=FLOAD
!       ENDIF
!       IF(INDJ.EQ.2) THEN
!          IDD=ID(1,NJ)
!          R(IDD)=FLOAD
!       ENDIF
!    ENDDO
!
!   GENERATE STIFFNESS MATRIX   
!
    CALL OFE_TRUSS_STIFF(NumOlp,NITG,NDIV,NWEI,E,AREA)
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(XABSC_WEIG_IN)
!
    RETURN
!}}}
    END SUBROUTINE OFE_TRUSS_G
!
    SUBROUTINE OFE1_TRUSS_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumOlp,NITG,NWEI,NDIV
    REAL(8) :: E,AREA
    INTEGER :: I,J
    INTEGER :: NI,NJ,INDI,INDJ,IDD      
    REAL(8) :: FLOAD
    INTEGER :: ALLocateStatus
!
    Mtype  = NPAR(2)
    NumOlp = NPAR(3)
    NITG   = NPAR(4)
    NWEI   = NPAR(5)
    NDIV   = NPAR(6)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    AREA = MATP(4,Mtype)
!-------- Allocate memory
    ALLOCATE (NDINT(4,NumOlp),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (INTPAR(2,NumOlp),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
    ALLOCATE (XABSC_WEIG_IN(NITG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   INITIALIZATION
!   READ ELEMENT DATA
    READ(ADELE) ((NDINT(I,J),I=1,4),J=1,NumOlp),((INTPAR(I,J),I=1,2),J=1,NumOlp)
!
    CALL GAULEG(NITG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!   CALCULATE THE LOAD VECTOR
!    FLOAD=1000D0
!    DO I=1,NumOlp
!       NI  = NDINT(1,I)
!       NJ  = NDINT(2,I)
!       INDI= NDINT(3,I)
!       INDJ= NDINT(4,I)
!       IF(INDI.EQ.2) THEN
!          IDD=ID(1,NI)
!          R(IDD)=FLOAD
!       ENDIF
!       IF(INDJ.EQ.2) THEN
!          IDD=ID(1,NJ)
!          R(IDD)=FLOAD
!       ENDIF
!    ENDDO
!
!   GENERATE STIFFNESS MATRIX   
!
    CALL OFE1_TRUSS_STIFF(NumOlp,NITG,NDIV,NWEI,E,AREA)
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(XABSC_WEIG_IN)
!
    RETURN
!}}}
    END SUBROUTINE OFE1_TRUSS_G
!
    SUBROUTINE OFE_TRUSS_STIFF(NumOlp,NITG,NDIV,NWEI,E,AREA) 
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NumOlp,NITG,NWEI,NDIV
    REAL(8),INTENT(IN) :: E
    REAL(8),INTENT(INOUT):: AREA
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J,JJ,K,KK
    INTEGER :: NI,NJ
    INTEGER :: INDI,INDJ
    REAL(8) :: XI,XJ,XM1,XM2,RADI,RADJ,XTMP,WEIGHT,XM
    REAL(8) :: LENGTH,XL
    REAL(8) :: DHDX(NPOR_MAX,2),STIFF(2*NPOR_MAX,2*NPOR_MAX)
    REAL(8) :: H(2),DH(2),WEI(2),DWEI(2),PHI(2),DPHI(2),WEISUM,DWEISUM
    REAL(8) :: P(NPOR_MAX,2),DP(NPOR_MAX,2)
    REAL(8) :: S,SS(3)
!
    DO I=1,NumOlp
       NI=NDINT(1,I)
       NJ=NDINT(2,I)
       INDI=NDINT(3,I)
       INDJ=NDINT(4,I)
       XI=X(NI)
       XJ=X(NJ) 
       RADI=INTPAR(1,I)
       RADJ=INTPAR(2,I)
       LENGTH=ABS(XI-XJ)/2.0D0
!
       XM=(XI+XJ)/2.0D0
!
       IF(RADI.LT.(2*LENGTH-1e-12)) THEN
          WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! --- RADIUS IS TOO SMALL!"
          STOP
       ENDIF
!
       IF(RADJ.LT.(2*LENGTH-1e-12)) THEN
          WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! --- RADIUS IS TOO SMALL!"
          STOP
       ENDIF
!
       STIFF=0.0D0
!
!      SECTION 1
!
!{{{
       DO JJ=1,NDIV
          XL = LENGTH/NDIV/2.0D0
          XM1= XI+(2*JJ-1)*XL
          DO J=1,NITG
             XTMP=XM1+XL*XABSC_WEIG_IN(J,1)
             WEIGHT=XABSC_WEIG_IN(J,2)
             AREA=1.0D0+XTMP
   !
             H(1)=1.0D0-(XTMP-XI)/LENGTH
             DH(1)=-1.0D0/LENGTH
             H(2)=1.0D0-(XM-XTMP)/LENGTH
             DH(2)=-DH(1)
             SELECT CASE(NWEI)
             CASE(1)
   !         CALCULATE WEI_I
                 S=ABS(XTMP-XI)/RADI
                 SS(1)=S*S
                 SS(2)=SS(1)*S
                 SS(3)=SS(2)*S
                 WEI(1)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
                 DWEI(1)=(-12.0D0*S+24.0D0*SS(1)-12.0D0*SS(2))/RADI
   !         CALCULATE WEI_J
                 S=ABS(XTMP-XJ)/RADJ
                 SS(1)=S*S
                 SS(2)=SS(1)*S
                 SS(3)=SS(2)*S
                 WEI(2)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
                 DWEI(2)=(12.0D0*S-24.0D0*SS(1)+12.0D0*SS(2))/RADJ
             CASE(2)
   !         CALCULATE WEI_I
                 S=XTMP-XI
                 SS(1)=(RADI-S)/RADI/S
                 WEI(1)=SS(1)*SS(1)
                 DWEI(1)=2*SS(1)*(-1.0D0/RADI/S-SS(1)/S)
   !         CALCULATE WEI_J
                 S=XJ-XTMP
                 SS(1)=(RADJ-S)/RADJ/S
                 WEI(2)=SS(1)*SS(1)
                 DWEI(2)=2*SS(1)*(1.0D0/RADJ/S+SS(1)/S)
             CASE(3)
   !         CALCULATE WEI_I
                 S=XTMP-XI
                 WEI(1)=1.0D0/S
                 DWEI(1)=-1.0D0/S/S
   !         CALCULATE WEI_J
                 S=XJ-XTMP
                 WEI(2)=1.0D0/S
                 DWEI(2)=1.0D0/S/S
             CASE DEFAULT
                 WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
                 STOP
             END SELECT
   !         CALCULATE PHI
             WEISUM=WEI(1)+WEI(2)
             DWEISUM=DWEI(1)+DWEI(2)
             PHI(1)=WEI(1)/WEISUM
             PHI(2)=WEI(2)/WEISUM
             DPHI(1)=DWEI(1)/WEISUM-PHI(1)/WEISUM*DWEISUM
             DPHI(2)=DWEI(2)/WEISUM-PHI(2)/WEISUM*DWEISUM
             SELECT CASE(NPOR_MAX)
             CASE(1)
                 DHDX(1,1)=(DH(1)+DH(2)*PHI(1)+H(2)*DPHI(1))
                 DHDX(1,2)=(DH(2)*PHI(2)+H(2)*DPHI(2))
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(1,2)=STIFF(1,2)+DHDX(1,1)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,2)=STIFF(2,2)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
             CASE(2)
   !         CALCULATE P
                 P(1,1)=(XTMP-XI)/RADI
                 P(1,2)=(XTMP-XJ)/RADJ
   !          
                 DP(1,1)=1.0D0/RADI
                 DP(1,2)=1.0D0/RADJ
   !
                 DHDX(1,1)=(DH(1)+DH(2)*PHI(1)+H(2)*DPHI(1))
                 DHDX(1,2)=(DH(2)*PHI(2)+H(2)*DPHI(2))
   !
                 DHDX(2,1)=DHDX(1,1)*P(1,1)+(H(1)+H(2)*PHI(1))*DP(1,1)
                 DHDX(2,2)=DHDX(1,2)*P(1,2)+(H(2)*PHI(2))*DP(1,2)              
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(3,3)=STIFF(3,3)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,3)=STIFF(2,3)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,3)=STIFF(1,3)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(4,4)=STIFF(4,4)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(3,4)=STIFF(3,4)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,4)=STIFF(2,4)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,4)=STIFF(1,4)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT
             CASE(3)
   !         CALCULATE P
                 P(1,1)=(XTMP-XI)/RADI
                 P(2,1)=P(1,1)*P(1,1)
                 P(1,2)=(XTMP-XJ)/RADJ
                 P(2,2)=P(1,2)*P(1,2)
   !          
                 DP(1,1)=1.0D0/RADI
                 DP(2,1)=2.0D0*P(1,1)/RADI
                 DP(1,2)=1.0D0/RADJ
                 DP(2,2)=2.0D0*P(1,2)/RADJ
   !
                 DHDX(1,1)=(DH(1)+DH(2)*PHI(1)+H(2)*DPHI(1))
                 DHDX(1,2)=(DH(2)*PHI(2)+H(2)*DPHI(2))
   !
                 DHDX(2,1)=DHDX(1,1)*P(1,1)+(H(1)+H(2)*PHI(1))*DP(1,1)
                 DHDX(2,2)=DHDX(1,2)*P(1,2)+(H(2)*PHI(2))*DP(1,2)
   !
                 DHDX(3,1)=DHDX(1,1)*P(2,1)+(H(1)+H(2)*PHI(1))*DP(2,1)
                 DHDX(3,2)=DHDX(1,2)*P(2,2)+(H(2)*PHI(2))*DP(2,2)
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(3,3)=STIFF(3,3)+DHDX(3,1)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,3)=STIFF(2,3)+DHDX(3,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,3)=STIFF(1,3)+DHDX(3,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(4,4)=STIFF(4,4)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,4)=STIFF(3,4)+DHDX(1,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,4)=STIFF(2,4)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,4)=STIFF(1,4)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(5,5)=STIFF(5,5)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(4,5)=STIFF(4,5)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,5)=STIFF(3,5)+DHDX(2,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,5)=STIFF(2,5)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,5)=STIFF(1,5)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT            
   !
                 STIFF(6,6)=STIFF(6,6)+DHDX(3,2)*E*DHDX(3,2)*AREA*WEIGHT
                 STIFF(5,6)=STIFF(5,6)+DHDX(3,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(4,6)=STIFF(4,6)+DHDX(3,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,6)=STIFF(3,6)+DHDX(3,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,6)=STIFF(2,6)+DHDX(3,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,6)=STIFF(1,6)+DHDX(3,2)*E*DHDX(1,1)*AREA*WEIGHT
             CASE DEFAULT
                 WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
                 STOP
             END SELECT     
          ENDDO
       ENDDO
!}}}
!
!      SECTION 2
!
!!{{{
       DO JJ=1,NDIV
          XL = LENGTH/NDIV/2.0D0
          XM2= XM+(2*JJ-1)*XL
          DO J=1,NITG
             XTMP=XM2+XL*XABSC_WEIG_IN(J,1)
             WEIGHT=XABSC_WEIG_IN(J,2)
             AREA=1.0D0+XTMP
   !
             H(1)=1.0D0-(XTMP-XM)/LENGTH
             DH(1)=-1.0D0/LENGTH
             H(2)=1.0D0-(XJ-XTMP)/LENGTH
             DH(2)=-DH(1)
   !         CALCULATE WEI_I
             SELECT CASE(NWEI)
             CASE(1)
                 S=ABS(XTMP-XI)/RADI
                 SS(1)=S*S
                 SS(2)=SS(1)*S
                 SS(3)=SS(2)*S
                 WEI(1)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
                 DWEI(1)=(-12.0D0*S+24.0D0*SS(1)-12.0D0*SS(2))/RADI
   !         CALCULATE WEI_J
                 S=ABS(XTMP-XJ)/RADJ
                 SS(1)=S*S
                 SS(2)=SS(1)*S
                 SS(3)=SS(2)*S
                 WEI(2)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
                 DWEI(2)=(12.0D0*S-24.0D0*SS(1)+12.0D0*SS(2))/RADJ
             CASE(2)
   !         CALCULATE WEI_I
                 S=XTMP-XI
                 SS(1)=(RADI-S)/RADI/S
                 WEI(1)=SS(1)*SS(1)
                 DWEI(1)=2*SS(1)*(-1.0D0/RADI/S-SS(1)/S)
   !         CALCULATE WEI_J
                 S=XJ-XTMP
                 SS(1)=(RADJ-S)/RADJ/S
                 WEI(2)=SS(1)*SS(1)
                 DWEI(2)=2*SS(1)*(1.0D0/RADJ/S+SS(1)/S)
             CASE(3)
   !         CALCULATE WEI_I
                 S=XTMP-XI
                 WEI(1)=1.0D0/S
                 DWEI(1)=-1.0D0/S/S
   !         CALCULATE WEI_J
                 S=XJ-XTMP
                 WEI(2)=1.0D0/S
                 DWEI(2)=1.0D0/S/S
             CASE DEFAULT
                 WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
                 STOP
             END SELECT
   !         CALCULATE PHI
             WEISUM=WEI(1)+WEI(2)
             DWEISUM=DWEI(1)+DWEI(2)
             PHI(1)=WEI(1)/WEISUM
             PHI(2)=WEI(2)/WEISUM
             DPHI(1)=DWEI(1)/WEISUM-PHI(1)/WEISUM*DWEISUM
             DPHI(2)=DWEI(2)/WEISUM-PHI(2)/WEISUM*DWEISUM
             SELECT CASE(NPOR_MAX)
             CASE(1)
                 DHDX(1,1)=(DH(1)*PHI(1)+H(1)*DPHI(1))
                 DHDX(1,2)=(DH(2)+DH(1)*PHI(2)+H(1)*DPHI(2))
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(1,2)=STIFF(1,2)+DHDX(1,1)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,2)=STIFF(2,2)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
             CASE(2)
   !         CALCULATE P
                 P(1,1)=(XTMP-XI)/RADI
                 P(1,2)=(XTMP-XJ)/RADJ
   !          
                 DP(1,1)=1.0D0/RADI
                 DP(1,2)=1.0D0/RADJ
   !
                 DHDX(1,1)=(DH(1)*PHI(1)+H(1)*DPHI(1))
                 DHDX(1,2)=(DH(2)+DH(1)*PHI(2)+H(1)*DPHI(2))
   !
                 DHDX(2,1)=DHDX(1,1)*P(1,1)+(H(1)*PHI(1))*DP(1,1)
                 DHDX(2,2)=DHDX(1,2)*P(1,2)+(H(1)*PHI(2)+H(2))*DP(1,2)             
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(3,3)=STIFF(3,3)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,3)=STIFF(2,3)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,3)=STIFF(1,3)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(4,4)=STIFF(4,4)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(3,4)=STIFF(3,4)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(2,4)=STIFF(2,4)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,4)=STIFF(1,4)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT
             CASE(3)
   !         CALCULATE P
                 P(1,1)=(XTMP-XI)/RADI
                 P(2,1)=P(1,1)*P(1,1)
                 P(1,2)=(XTMP-XJ)/RADJ
                 P(2,2)=P(1,2)*P(1,2)
   !          
                 DP(1,1)=1.0D0/RADI
                 DP(2,1)=2.0D0*P(1,1)/RADI
                 DP(1,2)=1.0D0/RADJ
                 DP(2,2)=2.0D0*P(1,2)/RADJ
   !
                 DHDX(1,1)=(DH(1)*PHI(1)+H(1)*DPHI(1))
                 DHDX(1,2)=(DH(2)+DH(1)*PHI(2)+H(1)*DPHI(2))
   !
                 DHDX(2,1)=DHDX(1,1)*P(1,1)+(H(1)*PHI(1))*DP(1,1)
                 DHDX(2,2)=DHDX(1,2)*P(1,2)+(H(1)*PHI(2)+H(2))*DP(1,2) 
   !
                 DHDX(3,1)=DHDX(1,1)*P(2,1)+(H(1)*PHI(1))*DP(2,1)
                 DHDX(3,2)=DHDX(1,2)*P(2,2)+(H(1)*PHI(2)+H(2))*DP(2,2)
   !
                 STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(3,3)=STIFF(3,3)+DHDX(3,1)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,3)=STIFF(2,3)+DHDX(3,1)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,3)=STIFF(1,3)+DHDX(3,1)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(4,4)=STIFF(4,4)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,4)=STIFF(3,4)+DHDX(1,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,4)=STIFF(2,4)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,4)=STIFF(1,4)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
   !
                 STIFF(5,5)=STIFF(5,5)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(4,5)=STIFF(4,5)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,5)=STIFF(3,5)+DHDX(2,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,5)=STIFF(2,5)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,5)=STIFF(1,5)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT            
   !
                 STIFF(6,6)=STIFF(6,6)+DHDX(3,2)*E*DHDX(3,2)*AREA*WEIGHT
                 STIFF(5,6)=STIFF(5,6)+DHDX(3,2)*E*DHDX(2,2)*AREA*WEIGHT
                 STIFF(4,6)=STIFF(4,6)+DHDX(3,2)*E*DHDX(1,2)*AREA*WEIGHT
                 STIFF(3,6)=STIFF(3,6)+DHDX(3,2)*E*DHDX(3,1)*AREA*WEIGHT
                 STIFF(2,6)=STIFF(2,6)+DHDX(3,2)*E*DHDX(2,1)*AREA*WEIGHT
                 STIFF(1,6)=STIFF(1,6)+DHDX(3,2)*E*DHDX(1,1)*AREA*WEIGHT
             CASE DEFAULT
                 WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
                 STOP
             END SELECT     
          ENDDO
       ENDDO
!!}}}
!
       CALL OFE_TRUSS_ADDBAN(I,NI,NJ,INDI,INDJ,STIFF,XL)
!
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OFE_TRUSS_STIFF
!
    SUBROUTINE OFE1_TRUSS_STIFF(NumOlp,NITG,NDIV,NWEI,E,AREA) 
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NumOlp,NITG,NWEI,NDIV
    REAL(8),INTENT(IN) :: E
    REAL(8),INTENT(INOUT):: AREA
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J,JJ,K,KK
    INTEGER :: NI,NJ
    INTEGER :: INDI,INDJ
    REAL(8) :: XI,XJ,RADI,RADJ,XTMP,WEIGHT,XM
    REAL(8) :: LENGTH,XL
    REAL(8) :: DHDX(NPOR_MAX,2),STIFF(2*NPOR_MAX,2*NPOR_MAX)
    REAL(8) :: H(2),PHI(2),DPHI(2)
    REAL(8) :: P(NPOR_MAX,2),DP(NPOR_MAX,2)
    REAL(8) :: S,SS(3)
!   SHAPE FUNCTIONS
    REAL(8) :: NSP_2(3),DNSP_2(3),NSP_1(2),DNSP_1(2)
    REAL(8) :: X_XI(3),XI_XJ(3)
!   WEIGHT FUNCTIONS
    REAL(8) :: WEI(3,2),WEISUM(3)

!
    DO I=1,NumOlp
       NI=NDINT(1,I)
       NJ=NDINT(2,I)
       INDI=NDINT(3,I)
       INDJ=NDINT(4,I)
       XI=X(NI)
       XJ=X(NJ) 
       RADI=INTPAR(1,I)
       RADJ=INTPAR(2,I)
       LENGTH=ABS(XI-XJ)
       XL=LENGTH/2.0D0
!
       XM=(XI+XJ)/2.0D0
!
       IF(RADI.LT.(LENGTH-1e-12)) THEN
          WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! --- RADIUS IS TOO SMALL!"
          STOP
       ENDIF
!
       IF(RADJ.LT.(LENGTH-1e-12)) THEN
          WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! --- RADIUS IS TOO SMALL!"
          STOP
       ENDIF
!
       STIFF=0.0D0
!
       XI_XJ(1)=XI-XM
       XI_XJ(2)=XI-XJ
       XI_XJ(3)=XM-XJ
!      CALCULATE WEIGHT FUNCTION
       SELECT CASE(NWEI)
       CASE(1)
!      CALCULATE WEI_I
           WEI(1,1)=1.0D0
           S=ABS(XI_XJ(1))/RADI
           SS(1)=S*S
           SS(2)=SS(1)*S
           SS(3)=SS(2)*S
           WEI(2,1)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
!
           S=ABS(XI_XJ(2))/RADI
           SS(1)=S*S
           SS(2)=SS(1)*S
           SS(3)=SS(2)*S
           WEI(3,1)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
!      CALCULATE WEI_J
           S=ABS(XI_XJ(2))/RADJ
           SS(1)=S*S
           SS(2)=SS(1)*S
           SS(3)=SS(2)*S
           WEI(1,2)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
!
           S=ABS(XI_XJ(3))/RADJ
           SS(1)=S*S
           SS(2)=SS(1)*S
           SS(3)=SS(2)*S
           WEI(2,2)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
!
           WEI(3,2)=1.0D0
       CASE DEFAULT
           WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
           STOP
       END SELECT
!
       WEISUM(1)=WEI(1,1)+WEI(1,2)
       WEISUM(2)=WEI(2,1)+WEI(2,2)
       WEISUM(3)=WEI(3,1)+WEI(3,2)
!
       DO J=1,3
          WEI(J,1)=WEI(J,1)/WEISUM(J)
          WEI(J,2)=WEI(J,2)/WEISUM(J)
       ENDDO
!
       DNSP_1(1)=1.0D0/XI_XJ(2)
       DNSP_1(2)=-DNSP_1(1)
       DO J=1,NITG
          XTMP=XM+XL*XABSC_WEIG_IN(J,1)
          WEIGHT=XABSC_WEIG_IN(J,2)
          AREA=1.0D0+XTMP
!         CALCULATE SHAPE FUNCTION
          X_XI(1)=XTMP-XI
          X_XI(2)=XTMP-XM
          X_XI(3)=XTMP-XJ
!
          NSP_2(1)=X_XI(2)*X_XI(3)/XI_XJ(1)/XI_XJ(2)
          NSP_2(2)=-X_XI(1)*X_XI(3)/XI_XJ(1)/XI_XJ(3)
          NSP_2(3)=X_XI(1)*X_XI(2)/XI_XJ(2)/XI_XJ(3)
!
          DNSP_2(1)=(2.0D0*XTMP-(XM+XJ))/XI_XJ(1)/XI_XJ(2)
          DNSP_2(2)=-(2.0D0*XTMP-(XI+XJ))/XI_XJ(1)/XI_XJ(3)
          DNSP_2(3)=(2.0D0*XTMP-(XI+XM))/XI_XJ(2)/XI_XJ(3)
!
          NSP_1(1)=X_XI(3)/XI_XJ(2)
          NSP_1(2)=-X_XI(1)/XI_XJ(2)
!
          SELECT CASE(NPOR_MAX)
          CASE(1)
              DHDX(1,1)=DNSP_1(1)*(NSP_2(1)+NSP_2(2)*WEI(2,1)+NSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(1)*(DNSP_2(1)+DNSP_2(2)*WEI(2,1)+DNSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+DNSP_1(2)*(NSP_2(1)*WEI(1,1)+NSP_2(2)*WEI(2,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(2)*(DNSP_2(1)*WEI(1,1)+DNSP_2(2)*WEI(2,1))
!
              DHDX(1,2)=DNSP_1(1)*(NSP_2(2)*WEI(2,2)+NSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+NSP_1(1)*(DNSP_2(2)*WEI(2,2)+DNSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+DNSP_1(2)*(NSP_2(1)*WEI(1,2)+NSP_2(2)*WEI(2,2)+NSP_2(3))
              DHDX(1,2)=DHDX(1,2)+NSP_1(2)*(DNSP_2(1)*WEI(1,2)+DNSP_2(2)*WEI(2,2)+DNSP_2(3))
!
              STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(1,2)=STIFF(1,2)+DHDX(1,1)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(2,2)=STIFF(2,2)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
          CASE(2)
!         CALCULATE P
              P(1,1)=(XTMP-XI)/RADI
              P(1,2)=(XTMP-XJ)/RADJ
!          
              DP(1,1)=1.0D0/RADI
              DP(1,2)=1.0D0/RADJ
!
              H(1)=NSP_1(1)*(NSP_2(1)+NSP_2(2)*WEI(2,1)+NSP_2(3)*WEI(3,1))
              H(1)=H(1)+NSP_1(2)*(NSP_2(1)*WEI(1,1)+NSP_2(2)*WEI(2,1))
!
              H(2)=NSP_1(1)*(NSP_2(2)*WEI(2,2)+NSP_2(3)*WEI(3,2))
              H(2)=H(2)+NSP_1(2)*(NSP_2(1)*WEI(1,2)+NSP_2(2)*WEI(2,2)+NSP_2(3))
!
              DHDX(1,1)=DNSP_1(1)*(NSP_2(1)+NSP_2(2)*WEI(2,1)+NSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(1)*(DNSP_2(1)+DNSP_2(2)*WEI(2,1)+DNSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+DNSP_1(2)*(NSP_2(1)*WEI(1,1)+NSP_2(2)*WEI(2,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(2)*(DNSP_2(1)*WEI(1,1)+DNSP_2(2)*WEI(2,1))
!
              DHDX(1,2)=DNSP_1(1)*(NSP_2(2)*WEI(2,2)+NSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+NSP_1(1)*(DNSP_2(2)*WEI(2,2)+DNSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+DNSP_1(2)*(NSP_2(1)*WEI(1,2)+NSP_2(2)*WEI(2,2)+NSP_2(3))
              DHDX(1,2)=DHDX(1,2)+NSP_1(2)*(DNSP_2(1)*WEI(1,2)+DNSP_2(2)*WEI(2,2)+DNSP_2(3))
!
              DHDX(2,1)=DHDX(1,1)*P(1,1)+H(1)*DP(1,1)
              DHDX(2,2)=DHDX(1,2)*P(1,2)+H(2)*DP(1,2)              
!
              STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(3,3)=STIFF(3,3)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(2,3)=STIFF(2,3)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,3)=STIFF(1,3)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(4,4)=STIFF(4,4)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
              STIFF(3,4)=STIFF(3,4)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(2,4)=STIFF(2,4)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,4)=STIFF(1,4)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT
          CASE(3)
!         CALCULATE P
              P(1,1)=(XTMP-XI)/RADI
              P(2,1)=P(1,1)*P(1,1)
              P(1,2)=(XTMP-XJ)/RADJ
              P(2,2)=P(1,2)*P(1,2)
!          
              DP(1,1)=1.0D0/RADI
              DP(2,1)=2.0D0*P(1,1)/RADI
              DP(1,2)=1.0D0/RADJ
              DP(2,2)=2.0D0*P(1,2)/RADJ
!
              H(1)=NSP_1(1)*(NSP_2(1)+NSP_2(2)*WEI(2,1)+NSP_2(3)*WEI(3,1))
              H(1)=H(1)+NSP_1(2)*(NSP_2(1)*WEI(1,1)+NSP_2(2)*WEI(2,1))
!
              H(2)=NSP_1(1)*(NSP_2(2)*WEI(2,2)+NSP_2(3)*WEI(3,2))
              H(2)=H(2)+NSP_1(2)*(NSP_2(1)*WEI(1,2)+NSP_2(2)*WEI(2,2)+NSP_2(3))
!
              DHDX(1,1)=DNSP_1(1)*(NSP_2(1)+NSP_2(2)*WEI(2,1)+NSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(1)*(DNSP_2(1)+DNSP_2(2)*WEI(2,1)+DNSP_2(3)*WEI(3,1))
              DHDX(1,1)=DHDX(1,1)+DNSP_1(2)*(NSP_2(1)*WEI(1,1)+NSP_2(2)*WEI(2,1))
              DHDX(1,1)=DHDX(1,1)+NSP_1(2)*(DNSP_2(1)*WEI(1,1)+DNSP_2(2)*WEI(2,1))
!
              DHDX(1,2)=DNSP_1(1)*(NSP_2(2)*WEI(2,2)+NSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+NSP_1(1)*(DNSP_2(2)*WEI(2,2)+DNSP_2(3)*WEI(3,2))
              DHDX(1,2)=DHDX(1,2)+DNSP_1(2)*(NSP_2(1)*WEI(1,2)+NSP_2(2)*WEI(2,2)+NSP_2(3))
              DHDX(1,2)=DHDX(1,2)+NSP_1(2)*(DNSP_2(1)*WEI(1,2)+DNSP_2(2)*WEI(2,2)+DNSP_2(3))
!
              DHDX(2,1)=DHDX(1,1)*P(1,1)+H(1)*DP(1,1)
              DHDX(2,2)=DHDX(1,2)*P(1,2)+H(2)*DP(1,2)
!
              DHDX(3,1)=DHDX(1,1)*P(2,1)+H(1)*DP(2,1)
              DHDX(3,2)=DHDX(1,2)*P(2,2)+H(2)*DP(2,2)
!
              STIFF(1,1)=STIFF(1,1)+DHDX(1,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(2,2)=STIFF(2,2)+DHDX(2,1)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,2)=STIFF(1,2)+DHDX(2,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(3,3)=STIFF(3,3)+DHDX(3,1)*E*DHDX(3,1)*AREA*WEIGHT
              STIFF(2,3)=STIFF(2,3)+DHDX(3,1)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,3)=STIFF(1,3)+DHDX(3,1)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(4,4)=STIFF(4,4)+DHDX(1,2)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(3,4)=STIFF(3,4)+DHDX(1,2)*E*DHDX(3,1)*AREA*WEIGHT
              STIFF(2,4)=STIFF(2,4)+DHDX(1,2)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,4)=STIFF(1,4)+DHDX(1,2)*E*DHDX(1,1)*AREA*WEIGHT
!
              STIFF(5,5)=STIFF(5,5)+DHDX(2,2)*E*DHDX(2,2)*AREA*WEIGHT
              STIFF(4,5)=STIFF(4,5)+DHDX(2,2)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(3,5)=STIFF(3,5)+DHDX(2,2)*E*DHDX(3,1)*AREA*WEIGHT
              STIFF(2,5)=STIFF(2,5)+DHDX(2,2)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,5)=STIFF(1,5)+DHDX(2,2)*E*DHDX(1,1)*AREA*WEIGHT            
!
              STIFF(6,6)=STIFF(6,6)+DHDX(3,2)*E*DHDX(3,2)*AREA*WEIGHT
              STIFF(5,6)=STIFF(5,6)+DHDX(3,2)*E*DHDX(2,2)*AREA*WEIGHT
              STIFF(4,6)=STIFF(4,6)+DHDX(3,2)*E*DHDX(1,2)*AREA*WEIGHT
              STIFF(3,6)=STIFF(3,6)+DHDX(3,2)*E*DHDX(3,1)*AREA*WEIGHT
              STIFF(2,6)=STIFF(2,6)+DHDX(3,2)*E*DHDX(2,1)*AREA*WEIGHT
              STIFF(1,6)=STIFF(1,6)+DHDX(3,2)*E*DHDX(1,1)*AREA*WEIGHT
          CASE DEFAULT
              WRITE(*,*) "ERROR, /GENERATEM/OEL_TRUSS_STIFF! "
              STOP
          END SELECT     
       ENDDO
       CALL OFE_TRUSS_ADDBAN(I,NI,NJ,INDI,INDJ,STIFF,XL)        
!
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OFE1_TRUSS_STIFF
!
    SUBROUTINE OFE_TRUSS_ADDBAN(INTI,NI,NJ,INDI,INDJ,STIFF,XL)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: INTI,NI,NJ,INDI,INDJ
    REAL(8),INTENT(IN) :: STIFF(2*NPOR_MAX,2*NPOR_MAX),XL
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J,I1,J1
    INTEGER :: NNI,NNJ
    INTEGER :: NUMIJ(2)
    INTEGER :: NDOFIJ(2),SFI,SFJ,KI,KJ
    INTEGER :: LL1
!
    NUMIJ=0
    IF(INDI.EQ.1) NUMIJ(1)=1
    IF(INDJ.EQ.1) NUMIJ(2)=1
!
!
    NDOFIJ(1)=NDOF(1,NI)
    NDOFIJ(2)=NDOF(1,NJ)
!
    DO I=1,2
       NNI=NDINT(I,INTI)
       DO I1=1,NDOFIJ(I)
          KI=ID(1,NNI)+I1-1
          SFI=(I-1)*NPOR_MAX+NUMIJ(I)+I1
          DO J=I,2
             NNJ=NDINT(J,INTI)
             DO J1=1,NDOFIJ(J)
                KJ=ID(1,NNJ)+J1-1
                SFJ=(J-1)*NPOR_MAX+NUMIJ(J)+J1
!
                IF(SFI.GT.SFJ) CYCLE
!
                IF(KJ.LE.KI) THEN
                   LL1=MAXA_COLSOL(KI)+ABS(KJ-KI)
                ELSE
                   LL1=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                ENDIF
!
                K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(SFI,SFJ)*XL
                KMATRIX(LL1) =KMATRIX(LL1) +STIFF(SFI,SFJ)*XL
             ENDDO                 
          ENDDO
       ENDDO
    ENDDO
!
    RETURN
!
!}}}
    END SUBROUTINE OFE_TRUSS_ADDBAN
!
    SUBROUTINE TRUSS_OVLP2_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumInt,NumEle,NEleAss,NITG
    REAL(8) :: E
    INTEGER :: I,J,K,L
    INTEGER :: N,ELE,IND,IDD
    INTEGER :: ALLocateStatus
!
    Mtype  = NPAR(2)
    NumInt = NPAR(3)
    NumEle = NPAR(4)
    NEleAss= NPAR(5)
    NITG   = NPAR(6) 
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
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
!
    ALLOCATE (XABSC_WEIG_IN(NITG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       INTPAR(1,I)=0.0D0
    END DO
!
    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
       END DO
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!   READ ELEMENT DATA
    READ(ADELE) ((NDMFS(I,J),I=1,2),J=1,NPAR(4)),(ELEPAR(1,I),I=1,NPAR(4)),((LM(I,J),I=1,NPAR(5)),J=1,NPAR(4)),&
    ((NDINT(I,J),I=1,2),J=1,NPAR(3)),(INTPAR(1,I),I=1,NPAR(3))    
!
    CALL GAULEG(NITG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    ALLOCATE (DHI_OVLP_TRUSS(NITG,NPOR_MAX,2*NEleAss),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (HI_OVLP_TRUSS(NITG,NPOR_MAX,2*NEleAss),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,1),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!   CALCULATE THE LOAD VECTOR
    DO I=1,NumEle
       IND=NDMFS(2,I)
       IF(IND.EQ.2) THEN
          N=NDMFS(1,I)
          IDD=ID(1,N)
          R(IDD)=1000  
       ENDIF
    ENDDO
!
!   GENERATE STIFFNESS MATRIX   
!
    DO N=1,NumInt
       DO I=1,2*NEleAss
          DO J=1,NPOR_MAX
             DO K=1,NITG
                DHI_OVLP_TRUSS(K,J,I)=0.0D0
                HI_OVLP_TRUSS(K,J,I)=0.0D0
             ENDDO
          ENDDO
       ENDDO
       CALL TRUSS_OVLP_STIFF(NITG,NEleAss,N,E)
    END DO
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM) 
!
    DEALLOCATE(DHI_OVLP_TRUSS)
    DEALLOCATE(HI_OVLP_TRUSS)
    DEALLOCATE(STIFF_IJ)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP2_G
!
    SUBROUTINE TRUSS_OVLP3_G(ADELE)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumInt,NumEle,NEleAss,NITG
    REAL(8) :: E
    INTEGER :: I,J,K,L,M
    INTEGER :: N,ELE,IND,IDD
    INTEGER :: IDDMin,IDDMax
    INTEGER :: ALLocateStatus
!
    Mtype  = NPAR(2)
    NumInt = NPAR(3)
    NumEle = NPAR(4)
    NEleAss= NPAR(5)
    NITG   = NPAR(6) 
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
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
!
    ALLOCATE (XABSC_WEIG_IN(NITG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!   
    DO I=1,NumInt
       NDINT(1,I)=0
       NDINT(2,I)=0
       NDINT(3,I)=0
       NDINT(4,I)=0
       INTPAR(1,I)=0.0D0
    END DO
!
    DO I=1,NumEle
       DO J=1,NEleAss
          LM(J,I)=0
       END DO
       NumLM(I)=0
       ELEPAR(1,I)=0.0D0
       NDMFS(1,I)=0
       NDMFS(2,I)=0
    END DO
!   READ ELEMENT DATA
    READ(ADELE) ((NDMFS(I,J),I=1,2),J=1,NPAR(4)),(ELEPAR(1,I),I=1,NPAR(4)),((LM(I,J),I=1,NPAR(5)),J=1,NPAR(4)),&
    ((NDINT(I,J),I=1,4),J=1,NPAR(3)),(INTPAR(1,I),I=1,NPAR(3)),(NumLM(I),I=1,NPAR(4))   
!
    CALL GAULEG(NITG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    ALLOCATE (Rmatrix_all(3,3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (Qmatrix_all(3,NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (Wmatrix_all(NEleAss,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
!   CALCULATE Rmatrix, Qmatrix, Wmatrix
!
    DO I=1,NumEle
       IND=NDMFS(2,I)
       M=NumLM(I)
       CALL LeastSQ(IND,3,M,I,Rmatrix_all(1,1,I),Qmatrix_all(1,1,I),Wmatrix_all(1,I)) 
    ENDDO
!    DO I=1,NumEle
!        DO J=1,3
!            WRITE(*,*) Rmatrix_all(J,:,I)
!        ENDDO
!        WRITE(*,*) " QMATRIX"
!        DO J=1,3
!            WRITE(*,*) Qmatrix_all(J,:,I)
!        ENDDO
!        WRITE(*,*) " WMATRIX"
!        WRITE(*,*) Wmatrix_all(:,I)
!    ENDDO
!   CALCULATE THE LOAD VECTOR
    DO I=1,NumEle
       IND=NDMFS(2,I)
       IF(IND.EQ.2) THEN
          N=NDMFS(1,I)
          IDD=ID(1,N)
          R(IDD)=1  
       ENDIF
    ENDDO
!
!   GENERATE STIFFNESS MATRIX   
!
    ALLOCATE (H_OVLP_TRUSS(2,NITG,NEQ),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
    DO N=1,NumInt
       IDDMin=NDINT(3,N)
       IDDMax=NDINT(4,N)
       DO I=IDDMin,IDDMax
          DO J=1,NITG
             H_OVLP_TRUSS(1,J,I)=0.0D0
             H_OVLP_TRUSS(2,J,I)=0.0D0
          ENDDO
       ENDDO
       CALL TRUSS_OVLP_STIFF3(NITG,NEleAss,N,E)
    END DO
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM) 
    DEALLOCATE(NumLM)
!
    DEALLOCATE(H_OVLP_TRUSS)
    DEALLOCATE(Rmatrix_all)
    DEALLOCATE(Qmatrix_all)
    DEALLOCATE(Wmatrix_all)
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP3_G
!
    SUBROUTINE QUADS_OVLP_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumInt,NumEle,NEleAss,NGDOF,ITYPE
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: BcType
!
    INTEGER :: I,J,K,L,M
    INTEGER :: N,ELE,IND,IDD
    INTEGER :: EleMin,EleMax
    INTEGER :: ALLocateStatus
!   READ
    Mtype  = NPAR(2)
    NumInt = NPAR(3)
    NumEle = NPAR(4)
    NEleAss= NPAR(5)
    NGDOF  = NPAR(6) 
    ITYPE  = NPAR(7)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
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
!   READ ELEMENT DATA
    READ(ADELE) ((NDMFS(I,J),I=1,6),J=1,NPAR(4)),((ELEPAR(I,J),I=1,4),J=1,NPAR(4)),((LM(I,J),I=1,NPAR(5)),J=1,NPAR(4)),&
    ((NDINT(I,J),I=1,3),J=1,NPAR(3)),(INTPAR(1,I),I=1,NPAR(3)),(NumLM(I),I=1,NPAR(4))
!   CALCULATE GAUSS POINT
    SELECT CASE(NGDOF)
    CASE(2)
       NITP_IN=3
       NITP_BC=2
    CASE(3)
       NITP_IN=6
       NITP_BC=2
    CASE(4)
       NITP_IN=7
       NITP_BC=3
    CASE(5)
       NITP_IN=9
       NITP_BC=3
    CASE(6)
       NITP_IN=12
       NITP_BC=4
    CASE(7)
       NITP_IN=13
       NITP_BC=4
    CASE(8)
       NITP_IN=16
       NITP_BC=4
    CASE(9)
       NITP_IN=19
       NITP_BC=4
    CASE(10)
       NITP_IN=25
       NITP_BC=4
    CASE DEFAULT
       WRITE(*,*) "ERROR, /GENERATEM/QUADS_OVLP_G!"
       STOP
    END SELECT 
!
    ALLOCATE (XABSC_WEIG_IN(NITP_IN,3),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (ACOOR_IN(NITP_IN,3),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
    ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   GENERATE CONSTITUTIVE MATRIX
    CALL G_KCS(2,ITYPE,E,PR)
!   CALCULATE THE LOAD VECTOR
    ALLOCATE (H_OVLP_QUADS(NPOR_MAX,NITP_BC,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!   INITIALIZATION
    H_OVLP_QUADS=0.0D0
!
    READ(ADBC) (NPAR(J),J=1,7)
    BcType=NPAR(1)
    SELECT CASE(BcType)
    CASE(1)
        CALL OVLP_LINE_G(ADBC,NEleAss,THIC)
    CASE DEFAULT
        WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
        STOP
    END SELECT
    DEALLOCATE(H_OVLP_QUADS)
!
!   GENERATE STIFFNESS MATRIX   
!
    ALLOCATE (DH_OVLP_QUADS(2*NPOR_MAX,NITP_IN,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (H_OVLP_QUADS(NPOR_MAX,NITP_IN,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!   INITIALIZATION
    DH_OVLP_QUADS=0.0D0
    H_OVLP_QUADS=0.0D0
! 
    ALLOCATE (STIFF_IJ(NPOR_MAX*NPOR_MAX,4),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!
    DO N=1,NumInt
       CALL QUADS_OVLP_STIFF(NGDOF,NEleAss,NumEle,N,THIC)
    END DO
!
    DEALLOCATE(NDINT)
    DEALLOCATE(INTPAR)
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(LM) 
    DEALLOCATE(NumLM)
! **************************
    DEALLOCATE(XABSC_WEIG_BC)
    DEALLOCATE(XABSC_WEIG_IN)
    DEALLOCATE(ACOOR_IN)
    DEALLOCATE(KCS)
! **************************
    DEALLOCATE(DH_OVLP_QUADS)
    DEALLOCATE(H_OVLP_QUADS)
    DEALLOCATE(STIFF_IJ)
!
    RETURN
!}}}
    END SUBROUTINE QUADS_OVLP_G
!
    SUBROUTINE OFE1_QUADS_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumOFE,NumEle,NGDOF,ITYPE,NWEI,NumCpl,NINTG,NumCpl2
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: I,J
!
    REAL(8),DIMENSION(:,:),   ALLOCATABLE :: LI
    REAL(8) :: L1,L2
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    Mtype  = NPAR(2)
    NumOFE = NPAR(3)
    NumEle = NPAR(4)
    NGDOF  = NPAR(5) 
    ITYPE  = NPAR(6)
    NWEI   = NPAR(7)
    NumCpl = NPAR(8)
    NINTG  = NPAR(9)
    NumCPL2= NPAR(10)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   CALCULATE GAUSS POINT
    SELECT CASE(NGDOF)
    CASE(1)
       NITP_IN=1
       NITP_BC=1
    CASE(2)
       NITP_IN=3
       NITP_BC=2
    CASE(3)
       NITP_IN=6
       NITP_BC=2
    CASE(4)
       NITP_IN=7
       NITP_BC=3
    CASE(5)
       NITP_IN=9
       NITP_BC=3
    CASE(6)
       NITP_IN=12
       NITP_BC=4
    CASE(7)
       NITP_IN=13
       NITP_BC=4
    CASE(8)
       NITP_IN=16
       NITP_BC=5
    CASE(9)
       NITP_IN=19
       NITP_BC=5
    CASE(10)
       NITP_IN=25
       NITP_BC=6
    CASE DEFAULT
       WRITE(*,*) "ERROR, /GENERATEM/QUADS_OVLP_G!"
       STOP
    END SELECT 
!-------- Allocate memory -------
    ALLOCATE (NDMFS(5,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDMFS***"
    ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ELEPAR***"
    READ(ADELE) ((NDMFS(I,J),I=1,5),J=1,NumEle),((ELEPAR(I,J),I=1,4),J=1,NumEle)
!**************************
!   CONSTRUCT LOAD VECTOR 
!**************************！    
    IF(OFE1BCONOFF) THEN
       ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"     
!
       CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   CALCULATE Pre-set Value
       ALLOCATE (LI(6,NITP_BC),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY LI ***" 
       LI(5,:)=(1.0D0-XABSC_WEIG_BC(:,1))/2.0D0
       LI(6,:)=1.0D0-LI(5,:)
!  ------ LI -------------------------
       LI(1,:)=2.0D0*LI(5,:)*LI(5,:)*(LI(5,:)-0.5D0)
       LI(2,:)=4.0D0*LI(5,:)*LI(6,:)
       LI(3,:)=2.0D0*LI(5,:)*LI(6,:)*(LI(6,:)-0.5D0)
       LI(4,:)=2.0D0*LI(5,:)*LI(6,:)*(LI(5,:)-0.5D0)
!
!   Initialization LOAD
!
       READ(ADBC) (NPAR(J),J=1,7)
       IF(NPAR(1).NE.2) THEN
          WRITE(*,*) "ERROR, OFE1_LINE_G"
          STOP
       ENDIF
!
       CALL OFE1_LINE_G(ADBC,THIC,LI)
! -----------------------------------
!
       DEALLOCATE(XABSC_WEIG_BC)  
       DEALLOCATE(LI)
!
    ENDIF
! -----------------------------------
!   GENERATE CONSTITUTIVE MATRIX 
    CALL G_KCS(2,ITYPE,E,PR)
!   READ ELEMENT DATA
    IF(NumOFE.GT.0) THEN    
       ALLOCATE (NDINT(3,NumOFE),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
       READ(ADELE) ((NDINT(I,J),I=1,3),J=1,NumOFE)
!*******************************
!   CONSTRUCT STIFFNESS MATRIX
!*******************************
       CALL OFE1_QUADS_STIFF(NUMOFE,NGDOF,THIC)
!
       DEALLOCATE(NDINT)
    ENDIF
!
    IF(NumCPL.GT.0) THEN
       ALLOCATE (CPLINT(4,NumCPL),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
       READ(ADELE) ((CPLINT(I,J),I=1,4),J=1,NumCPL)
!   Calculate contributions from coupling regions       
       CALL OFE1_COUPLE_STIFF(NUMCPL,NINTG,THIC)
!**************************
       DEALLOCATE(CPLINT)
    ENDIF   
!
    IF(NumCpl2.GT.0) THEN
       ALLOCATE (CPLINT(6,NumCpl2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
       READ(ADELE) ((CPLINT(I,J),I=1,6),J=1,NumCpl2)
!   Calculate contributions from coupling regions       
       CALL OFE1_CPL2_STIFF(NUMCPL2,NGDOF,THIC)
!**************************
       DEALLOCATE(CPLINT)
    ENDIF 
! **************************
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(KCS)
!
    RETURN
!}}}
    END SUBROUTINE OFE1_QUADS_G
!
    SUBROUTINE FEMC_QUADS_G(ADELE,ADBC)
!{{{
!--------------------------------------------------
!   GENERATE AND STORE
!      STIFFNESS MATRIX
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADELE,ADBC
!   ASSISTANT VARIABLES
    INTEGER :: Mtype,NumOFE,NumEle,NGDOF,ITYPE
    REAL(8) :: E,PR,DEN,THIC
    INTEGER :: I,J
    REAL(8),DIMENSION(:,:),   ALLOCATABLE :: LI
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    Mtype  = NPAR(2)
    NumOFE = NPAR(3)
    NumEle = NPAR(4)
    NGDOF  = NPAR(5) 
    ITYPE  = NPAR(6)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    PR   = MATP(2,Mtype)
    DEN  = MATP(3,Mtype)
    THIC = MATP(4,Mtype)
!   CALCULATE GAUSS POINT
    SELECT CASE(NGDOF)
    CASE(2)
       NITP_IN=3
       NITP_BC=2
    CASE(3)
       NITP_IN=6
       NITP_BC=2
    CASE(4)
       NITP_IN=7
       NITP_BC=3
    CASE(5)
       NITP_IN=9
       NITP_BC=3
    CASE(6)
       NITP_IN=12
       NITP_BC=4
    CASE(7)
       NITP_IN=13
       NITP_BC=4
    CASE(8)
       NITP_IN=16
       NITP_BC=5
    CASE(9)
       NITP_IN=19
       NITP_BC=5
    CASE(10)
       NITP_IN=25
       NITP_BC=6
    CASE DEFAULT
       WRITE(*,*) "ERROR, /GENERATEM/QUADS_OVLP_G!"
       STOP
    END SELECT 
!-------- Allocate memory -------
    ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDMFS***"
    ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ELEPAR***"
    READ(ADELE) ((NDMFS(I,J),I=1,3),J=1,NumEle),(ELEPAR(1,J),J=1,NumEle)
!**************************
!   CONSTRUCT LOAD VECTOR 
!**************************！    
    IF(FEMCBCONOFF) THEN
       ALLOCATE (XABSC_WEIG_BC(NITP_BC,2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"     
       CALL GAULEG(NITP_BC,XABSC_WEIG_BC(1,1),XABSC_WEIG_BC(1,2))
!   CALCULATE Pre-set Value
       ALLOCATE (LI(NITP_BC,2),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY LI ***" 
       LI(:,1)=(1.0D0-XABSC_WEIG_BC(:,1))/2.0D0
       LI(:,2)=1.0D0-LI(:,1)
!
!   Initialization LOAD
!
       READ(ADBC) (NPAR(J),J=1,7)
       IF(NPAR(1).NE.4) THEN
          WRITE(*,*) "ERROR, FEMC_LINE_G"
          STOP
       ENDIF
! 
       CALL FEMC_LINE_G(ADBC,THIC,LI(:,1),LI(:,2))
!
       DEALLOCATE(XABSC_WEIG_BC)  
       DEALLOCATE(LI)
    ENDIF
! -----------------------------------
!   GENERATE CONSTITUTIVE MATRIX 
    CALL G_KCS(2,ITYPE,E,PR)
!   READ ELEMENT DATA
    IF(NumOFE.GT.0) THEN    
       ALLOCATE (NDINT(3,NumOFE),STAT = ALLocateStatus)
       IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
       READ(ADELE) ((NDINT(I,J),I=1,3),J=1,NumOFE)
!*******************************
!   CONSTRUCT STIFFNESS MATRIX
!*******************************
       CALL FEMC_QUADS_STIFF(NUMOFE,NGDOF,THIC)
       DEALLOCATE(NDINT)
    ENDIF
! **************************
    DEALLOCATE(NDMFS)
    DEALLOCATE(ELEPAR)
    DEALLOCATE(KCS)
!
    RETURN
!}}}
    END SUBROUTINE FEMC_QUADS_G
!
    SUBROUTINE FEMC_QUADS_STIFF(NUMOFE,NGDOF,THIC)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NUMOFE,NGDOF
    REAL(8),INTENT(IN) :: THIC
! IN-FUNCTION VARIABLES
!------------------------------------------
    REAL(8),DIMENSION(NITP_IN,3) :: AREA_COOR
    REAL(8),DIMENSION(NITP_IN)   :: GAUSS_WEI
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,6) :: STIFF
!------------------------------------------
    INTEGER              :: I,J,K,L
    INTEGER              :: Ele1,Ele2,Ele3,N1,N2,N3
    REAL(8)              :: RADI,RADJ,RADK
    REAL(8),DIMENSION(2) :: XY1,XY2,XY3,XYC1,XYC2,XYC3
    REAL(8),DIMENSION(NITP_IN) :: XTMP,YTMP,X_BAR,Y_BAR
!
    REAL(8),DIMENSION(NITP_IN) :: HIW,dHIWdx,dHIWdy,HIWC1,HIWC2,HIWC3
    REAL(8),DIMENSION(NITP_IN) :: HJW,dHJWdx,dHJWdy,HJWC1,HJWC2,HJWC3
    REAL(8),DIMENSION(NITP_IN) :: HKW,dHKWdx,dHKWdy
!
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHImdx,dHImdy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHJmdx,dHJmdy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHKmdx,dHKmdy
!
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHImdxW,dHImdyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHJmdxW,dHJmdyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHKmdxW,dHKmdyW
!
    REAL(8) :: dLIdx,dLIdy,dLJdx,dLJdy,dLKdx,dLKdy
!
    INTEGER :: NDOFI,NDOFJ,NDOFK
    INTEGER :: NCOUNT
    REAL(8) :: CONST
    REAL(8) :: AREA2
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK
!
    REAL(8) :: DHIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: KCS1,KCS2,KCS3
!
    CALL GAUTRIANGLE_OFE1(NGDOF,NITP_IN,AREA_COOR,GAUSS_WEI)
!
    DO I=1,NumOFE       
       Ele1=NDINT(1,I)
       Ele2=NDINT(2,I)
       Ele3=NDINT(3,I)
! 
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
!
       NDOFI=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFJ=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFK=NDOF(1,N3)+NDMFS(2,Ele3)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
       RADK=ELEPAR(1,Ele3)
!
       XY1(1)=X(N1)
       XY1(2)=Y(N1)
       XY2(1)=X(N2)
       XY2(2)=Y(N2)
       XY3(1)=X(N3)
       XY3(2)=Y(N3)
! -----------------------------------
       XYC1=XY1-XY2
       XYC2=XY1-XY3
       XYC3=XY2-XY3
! -----------------------------------   
       XTMP=AREA_COOR(:,1)*XY1(1)+AREA_COOR(:,2)*XY2(1)+AREA_COOR(:,3)*XY3(1)
       YTMP=AREA_COOR(:,1)*XY1(2)+AREA_COOR(:,2)*XY2(2)+AREA_COOR(:,3)*XY3(2)
! -----------------------------------
       AREA2=X(N1)*XYC3(2)-X(N2)*XYC2(2)+X(N3)*XYC1(2)
!
       dLIdx=XYC3(2)/AREA2 
       dLIdy=-XYC3(1)/AREA2
       dLJdx=-XYC2(2)/AREA2
       dLJdy=XYC2(1)/AREA2
       dLKdx=XYC1(2)/AREA2
       dLKdy=-XYC1(1)/AREA2
!
       HIW(:)=AREA_COOR(:,1)
       HJW(:)=AREA_COOR(:,2)
       HKW(:)=AREA_COOR(:,3)
!
       dHIWdx(:)=dLIdx
       dHIWdy(:)=dLIdy
       dHJWdx(:)=dLJdx
       dHJWdy(:)=dLJdy
       dHKWdx(:)=dLKdx
       dHKWdy(:)=dLKdy
!----------------------------
       dHImdx=0.0D0
       dHImdy=0.0D0
       dHJmdx=0.0D0
       dHJmdy=0.0D0
       dHKmdx=0.0D0
       dHKmdy=0.0D0
! **************************** dHImdx,dHImdy
!
       X_BAR=(XTMP-X(N1))/RADI
       Y_BAR=(YTMP-Y(N1))/RADI
       SELECT CASE(NDOFI)
       CASE(1)
          dHImdx(:,1)=dHIWdx
          dHImdy(:,1)=dHIWdy       
!!
!          dHImdx(:,1)=dHIWdx*AREA_COOR(:,1)+HIW*dLIdx
!          dHImdy(:,1)=dHIWdy*AREA_COOR(:,1)+HIW*dLIdy
!!
!          dHJmdx(:,1)=dHIWdx*AREA_COOR(:,2)+HIW*dLJdx
!          dHJmdy(:,1)=dHIWdy*AREA_COOR(:,2)+HIW*dLJdy 
!!
!          dHKmdx(:,1)=dHIWdx*AREA_COOR(:,3)+HIW*dLKdx
!          dHKmdy(:,1)=dHIWdy*AREA_COOR(:,3)+HIW*dLKdy 
       CASE(3)
          HIWC1=HIW/RADI
!
          dHImdx(:,1)=dHIWdx
          dHImdy(:,1)=dHIWdy
!
          dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
          dHImdy(:,2)=dHIWdy*X_BAR
!
          dHImdx(:,3)=dHIWdx*Y_BAR
          dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1           
       CASE(4)
          HIWC1=HIW/RADI
!
          dHImdx(:,1)=dHIWdx
          dHImdy(:,1)=dHIWdy
!
          dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
          dHImdy(:,2)=dHIWdy*X_BAR
!
          dHImdx(:,3)=dHIWdx*Y_BAR
          dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1
!
          dHImdx(:,4)=dHImdx(:,2)*Y_BAR
          dHImdy(:,4)=dHImdy(:,3)*X_BAR          
       CASE(6)
          HIWC1=HIW/RADI
!
          dHImdx(:,1)=dHIWdx
          dHImdy(:,1)=dHIWdy
!
          dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
          dHImdy(:,2)=dHIWdy*X_BAR
!
          dHImdx(:,3)=dHIWdx*Y_BAR
          dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1
!
          dHImdx(:,4)=dHImdx(:,2)*Y_BAR
          dHImdy(:,4)=dHImdy(:,3)*X_BAR
!
          dHImdx(:,5)=(dHImdx(:,2)+HIWC1)*X_BAR
          dHImdy(:,5)=dHImdy(:,2)*X_BAR
!
          dHImdx(:,6)=dHImdx(:,3)*Y_BAR
          dHImdy(:,6)=(dHImdy(:,3)+HIWC1)*Y_BAR              
       CASE DEFAULT
          WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
          STOP
       END SELECT 
! **************************** dHJmdx,dHJmdy
       X_BAR=(XTMP-X(N2))/RADJ
       Y_BAR=(YTMP-Y(N2))/RADJ
       SELECT CASE(NDOFJ)
       CASE(1)
          dHJmdx(:,1)=dHJWdx
          dHJmdy(:,1)=dHJWdy
!
!          dHImdx(:,1)=dHImdx(:,1)+dHJWdx*AREA_COOR(:,1)+HJW*dLIdx
!          dHImdy(:,1)=dHImdy(:,1)+dHJWdy*AREA_COOR(:,1)+HJW*dLIdy
!!
!          dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx*AREA_COOR(:,2)+HJW*dLJdx
!          dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy*AREA_COOR(:,2)+HJW*dLJdy 
!!
!          dHKmdx(:,1)=dHKmdx(:,1)+dHJWdx*AREA_COOR(:,3)+HJW*dLKdx
!          dHKmdy(:,1)=dHKmdy(:,1)+dHJWdy*AREA_COOR(:,3)+HJW*dLKdy 
       CASE(3)
          HJWC1=HJW/RADJ

          dHJmdx(:,1)=dHJWdx
          dHJmdy(:,1)=dHJWdy

          dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
          dHJmdy(:,2)=dHJWdy*X_BAR

          dHJmdx(:,3)=dHJWdx*Y_BAR
          dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1           
       CASE(4)
          HJWC1=HJW/RADJ

          dHJmdx(:,1)=dHJWdx
          dHJmdy(:,1)=dHJWdy

          dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
          dHJmdy(:,2)=dHJWdy*X_BAR

          dHJmdx(:,3)=dHJWdx*Y_BAR
          dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1 

          dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
          dHJmdy(:,4)=dHJmdy(:,3)*X_BAR          
       CASE(6)
          HJWC1=HJW/RADJ

          dHJmdx(:,1)=dHJWdx
          dHJmdy(:,1)=dHJWdy

          dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
          dHJmdy(:,2)=dHJWdy*X_BAR

          dHJmdx(:,3)=dHJWdx*Y_BAR
          dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1

          dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
          dHJmdy(:,4)=dHJmdy(:,3)*X_BAR

          dHJmdx(:,5)=(dHJmdx(:,2)+HJWC1)*X_BAR
          dHJmdy(:,5)=dHJmdy(:,2)*X_BAR

          dHJmdx(:,6)=dHJmdx(:,3)*Y_BAR
          dHJmdy(:,6)=(dHJmdy(:,3)+HJWC1)*Y_BAR 
       CASE DEFAULT
          WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
          STOP
       END SELECT 
! **************************** dHKmdx,dHKmdy
       X_BAR=(XTMP-X(N3))/RADK
       Y_BAR=(YTMP-Y(N3))/RADK
       SELECT CASE(NDOFK)
       CASE(1)
          dHKmdx(:,1)=dHKWdx
          dHKmdy(:,1)=dHKWdy
!
!          dHImdx(:,1)=dHImdx(:,1)+dHKWdx*AREA_COOR(:,1)+HKW*dLIdx
!          dHImdy(:,1)=dHImdy(:,1)+dHKWdy*AREA_COOR(:,1)+HKW*dLIdy
!!
!          dHJmdx(:,1)=dHJmdx(:,1)+dHKWdx*AREA_COOR(:,2)+HKW*dLJdx
!          dHJmdy(:,1)=dHJmdy(:,1)+dHKWdy*AREA_COOR(:,2)+HKW*dLJdy 
!!
!          dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx*AREA_COOR(:,3)+HKW*dLKdx
!          dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy*AREA_COOR(:,3)+HKW*dLKdy 
       CASE(3)
          HJWC1=HKW/RADK

          dHKmdx(:,1)=dHKWdx
          dHKmdy(:,1)=dHKWdy

          dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
          dHKmdy(:,2)=dHKWdy*X_BAR

          dHKmdx(:,3)=dHKWdx*Y_BAR
          dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1           
       CASE(4)
          HJWC1=HKW/RADK

          dHKmdx(:,1)=dHKWdx
          dHKmdy(:,1)=dHKWdy

          dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
          dHKmdy(:,2)=dHKWdy*X_BAR

          dHKmdx(:,3)=dHKWdx*Y_BAR
          dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1  

          dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
          dHKmdy(:,4)=dHKmdy(:,3)*X_BAR          
       CASE(6)
          HJWC1=HKW/RADK

          dHKmdx(:,1)=dHKWdx
          dHKmdy(:,1)=dHKWdy

          dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
          dHKmdy(:,2)=dHKWdy*X_BAR

          dHKmdx(:,3)=dHKWdx*Y_BAR
          dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1

          dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
          dHKmdy(:,4)=dHKmdy(:,3)*X_BAR

          dHKmdx(:,5)=(dHKmdx(:,2)+HJWC1)*X_BAR
          dHKmdy(:,5)=dHKmdy(:,2)*X_BAR

          dHKmdx(:,6)=dHKmdx(:,3)*Y_BAR
          dHKmdy(:,6)=(dHKmdy(:,3)+HJWC1)*Y_BAR 
       CASE DEFAULT
          WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
          STOP
       END SELECT 
!------------------------------------------ 297+NITP_IN*82+NITP_IN*39=297+NITP_IN*121 FLOP
       CONST=AREA2*THIC/2.0D0
       KCS1=KCS(1)*CONST
       KCS2=KCS(2)*CONST
       KCS3=KCS(3)*CONST
!------------------------------------------ 297+NITP_IN*121+5=302+NITP_IN*121 FLOP
!------------------------------------------ STIFFII,STIFFJJ,STIFFKK
!      STIFFII
       NCOUNT=0
       DO J=1,NDOFI
          dHImdxW(:,J)=dHImdx(:,J)*GAUSS_WEI(:)
          dHImdyW(:,J)=dHImdy(:,J)*GAUSS_WEI(:)
          ADI(J)=NCOUNT
          DO K=J,NDOFI
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHImdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHImdy(:,K))  
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,1)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,1)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,1)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,1)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))           
!      STIFFJJ
       NCOUNT=0
       DO J=1,NDOFJ
          dHJmdxW(:,J)=dHJmdx(:,J)*GAUSS_WEI(:)
          dHJmdyW(:,J)=dHJmdy(:,J)*GAUSS_WEI(:)
          ADJ(J)=NCOUNT
          DO K=J,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,2)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,2)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,2)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,2)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!
!      STIFFKK
       NCOUNT=0
       DO J=1,NDOFK
          dHKmdxW(:,J)=dHKmdx(:,J)*GAUSS_WEI(:)
          dHKmdyW(:,J)=dHKmdy(:,J)*GAUSS_WEI(:)
          ADK(J)=NCOUNT
          DO K=J,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHKmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHKmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHKmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHKmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,3)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,3)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,3)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,3)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IJ
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,4)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,4)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,4)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,4)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IK
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,5)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,5)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,5)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,5)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF JK
       NCOUNT=0
       DO J=1,NDOFJ
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,6)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,6)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,6)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,6)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))    
!------------------------------------------ 302+NITP_IN*121+NITP_IN*648-78=224+NITP*769 FLOP
       CALL OFE1_QUADS_ADDBAN(Ele1,Ele2,Ele3,N1,N2,N3,NDOFI,NDOFJ,NDOFK,ADI,ADJ,ADK,STIFF)
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE FEMC_QUADS_STIFF
!
    SUBROUTINE FEMC_LINE_G(ADBC,THIC,HIW,HJW)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC
    REAL(8),INTENT(IN) :: THIC
    REAL(8),INTENT(IN) :: HIW(NITP_BC),HJW(NITP_BC)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    INTEGER :: Ele1,Ele2,N1,N2
    REAL(8) :: LENGTH_HF,LENGTH,NX,NY
    REAL(8) :: XM,XL,YM,YL
    REAL(8) :: RADI,RADJ
    REAL(8),DIMENSION(2) :: FLOAD
    INTEGER,DIMENSION(2) :: FUNC
    REAL(8),DIMENSION(NITP_BC) :: XTMP,YTMP,X_BAR,Y_BAR
    REAL(8),DIMENSION(2) :: LOADD
    REAL(8),DIMENSION(NITP_BC,2) :: LOAD
    REAL(8),DIMENSION(NPOR_MAX,2):: SI,SJ
    REAL(8),DIMENSION(NITP_BC,NPOR_MAX) :: PmI,PmJ
    INTEGER :: NDOFF
    REAL(8) :: CONST1
    INTEGER :: NUM1,NUM2,JJ1,JJ2,J1
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCED ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCIND ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCPAR ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
!   CALCULATE LOAD VECTOR
!
    DO I=1,NumBc
       ELe1=BCED(1,I)
       ELe2=BCED(2,I)
!
       N1=NDMFS(1,ELe1)
       N2=NDMFS(1,ELe2)
!
       XM=(X(N1)+X(N2))/2.0D0
       XL=(X(N2)-X(N1))/2.0D0
       YM=(Y(N1)+Y(N2))/2.0D0
       YL=(Y(N2)-Y(N1))/2.0D0
!
       LENGTH=XL*XL+YL*YL
       LENGTH_HF=SQRT(LENGTH)
       LENGTH=LENGTH_HF*2.0D0
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
!----- XTMP YTMP -----------------------
       XTMP=XM+XL*XABSC_WEIG_BC(:,1)
       YTMP=YM+YL*XABSC_WEIG_BC(:,1)
!----- HIW  HJW ------------------------
       DO J=1,NITP_BC
!-----   LOAD  ---------------------------
          LOADD= F_LOAD(FLOAD,FUNC,NX,NY,XTMP(J),YTMP(J),N1,N2)
          LOAD(J,1)=LOADD(1)*XABSC_WEIG_BC(J,2)
          LOAD(J,2)=LOADD(2)*XABSC_WEIG_BC(J,2)
       ENDDO
!
       CONST1=LENGTH*THIC/2.0D0
! ------ PmI --------------
       NDOFF=NDOF(1,N1)+NDMFS(2,Ele1)
!
       PmI=1.0D0
       X_BAR=(XTMP-X(N1))/RADI
       Y_BAR=(YTMP-Y(N1))/RADI
       SELECT CASE(NDOFF)
       CASE(1)
           SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
           SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
       CASE(2)
           PmI(:,2)=X_BAR

           SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
           SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1

           SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
           SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
       CASE(3)
           PmI(:,2)=X_BAR
           PmI(:,3)=Y_BAR

           SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
           SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1

           SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
           SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1

           SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
           SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1              
       CASE(4)
           PmI(:,2)=X_BAR
           PmI(:,3)=Y_BAR
           PmI(:,4)=PmI(:,2)*Y_BAR

           SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
           SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1

           SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
           SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1

           SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
           SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1

           SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
           SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1
       CASE(6)
           PmI(:,2)=X_BAR
           PmI(:,3)=Y_BAR
           PmI(:,4)=PmI(:,2)*Y_BAR
           PmI(:,5)=PmI(:,2)*X_BAR
           PmI(:,6)=PmI(:,3)*Y_BAR

           SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
           SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1

           SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
           SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1

           SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
           SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1

           SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
           SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1

           SI(5,1)=SUM(HIW*PmI(:,5)*LOAD(:,1))*CONST1
           SI(5,2)=SUM(HIW*PmI(:,5)*LOAD(:,2))*CONST1

           SI(6,1)=SUM(HIW*PmI(:,6)*LOAD(:,1))*CONST1
           SI(6,2)=SUM(HIW*PmI(:,6)*LOAD(:,2))*CONST1
       CASE DEFAULT
          WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
          STOP
       END SELECT
! ------- PmJ --------------------
       NDOFF=NDOF(1,N2)+NDMFS(2,Ele2)
!
       PmJ=1.0D0
       X_BAR=(XTMP-X(N2))/RADJ
       Y_BAR=(YTMP-Y(N2))/RADJ
       SELECT CASE(NDOFF)
       CASE(1)
           SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
           SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
       CASE(2)
           PmJ(:,2)=X_BAR

           SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
           SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)

           SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
           SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
       CASE(3)
           PmJ(:,2)=X_BAR
           PmJ(:,3)=Y_BAR

           SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
           SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)

           SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
           SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1

           SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
           SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
       CASE(4)
           PmJ(:,2)=X_BAR
           PmJ(:,3)=Y_BAR
           PmJ(:,4)=PmJ(:,2)*Y_BAR

           SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
           SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)

           SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
           SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1

           SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
           SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1

           SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
           SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
       CASE(6)
           PmJ(:,2)=X_BAR
           PmJ(:,3)=Y_BAR
           PmJ(:,4)=PmJ(:,2)*Y_BAR
           PmJ(:,5)=PmJ(:,2)*X_BAR
           PmJ(:,6)=PmJ(:,3)*Y_BAR

           SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
           SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)

           SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
           SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1

           SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
           SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1

           SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
           SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1

           SJ(5,1)=SUM(HJW*LOAD(:,1)*PmJ(:,5))*CONST1
           SJ(5,2)=SUM(HJW*LOAD(:,2)*PmJ(:,5))*CONST1

           SJ(6,1)=SUM(HJW*LOAD(:,1)*PmJ(:,6))*CONST1
           SJ(6,2)=SUM(HJW*LOAD(:,2)*PmJ(:,6))*CONST1
       CASE DEFAULT
          WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
          STOP
       END SELECT
!   
       DO J=1,2
          J1=1+J
!
          NUM1=ID(J,N1)
          NUM2=NUM1+NDOF(J,N1)-1
!
          JJ1=1+NDMFS(J1,Ele1)
          JJ2=NDOF(J,N1)+NDMFS(J1,Ele1)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SI(JJ1:JJ2,J)
!--------------------------------------------------
          NUM1=ID(J,N2)
          NUM2=NUM1+NDOF(J,N2)-1
!
          JJ1=1+NDMFS(J1,Ele2)
          JJ2=NDOF(J,N2)+NDMFS(J1,Ele2)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SJ(JJ1:JJ2,J)
       ENDDO
!
    ENDDO      
!------------------------------------------
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE FEMC_LINE_G   
!
    SUBROUTINE OFE1_COUPLE_STIFF(NUMCPL,NINTG,THIC)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NUMCPL,NINTG
    REAL(8),INTENT(IN) :: THIC
!--------------------------------------!
!   IN-FUNCTION VARIABLES
!------------------------------------------
    INTEGER :: ALLocateStatus
    INTEGER :: I,J,K,L,NCOUNT
    INTEGER :: NS,NE
    REAL(8),DIMENSION(NINTG*NINTG) :: HI,HJ,HK,HL,HC1
    REAL(8),DIMENSION(NINTG*NINTG) :: dHIdx,dHIdy,dHJdx,dHJdy,dHKdx,dHKdy,dHLdx,dHLdy
    REAL(8),DIMENSION(NINTG*NINTG) :: XTMP,YTMP,X_BAR,Y_BAR
    REAL(8),DIMENSION(NINTG*NINTG) :: XTMPC,YTMPC
    REAL(8),DIMENSION(NINTG*NINTG) :: dHIdxC,dHIdyC,dHJdxC,dHJdyC,dHKdxC,dHKdyC
    REAL(8),DIMENSION(NINTG*NINTG) :: GAUSS_WEI
    REAL(8),DIMENSION(NINTG*NINTG,NPOR_MAX) :: dHImdx,dHImdy,dHJmdx,dHJmdy,dHKmdx,dHKmdy,dHLmdx,dHLmdy
!
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4) :: DHIHJ
    REAL(8) ::  CONST,KCS1,KCS2,KCS3
!
    INTEGER :: Ele1,Ele2,Ele3,Ele4
    INTEGER :: N1,N2,N3,N4
    INTEGER :: NDOFI,NDOFJ,NDOFK,NDOFL
    REAL(8) :: RADI,RADJ,RADK,RADL
    REAL(8) :: LENGTH,LENGTH2,XM,YM
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,10) :: STIFF
!
    REAL(8),DIMENSION(NINTG*NINTG,NPOR_MAX) :: dHImdxW,dHImdyW
    REAL(8),DIMENSION(NINTG*NINTG,NPOR_MAX) :: dHJmdxW,dHJmdyW
    REAL(8),DIMENSION(NINTG*NINTG,NPOR_MAX) :: dHKmdxW,dHKmdyW
    REAL(8),DIMENSION(NINTG*NINTG,NPOR_MAX) :: dHLmdxW,dHLmdyW
!
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK,ADL
!
    ALLOCATE (XABSC_WEIG_IN(NINTG,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"   
!    
    CALL GAULEG(NINTG,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    DO I=1,NINTG
       NS=(I-1)*NINTG+1
       NE=I*NINTG
!
       GAUSS_WEI(NS:NE)=XABSC_WEIG_IN(:,2)*XABSC_WEIG_IN(I,2)
!
       HI(NS:NE)=(-XABSC_WEIG_IN(:,1)+1.0D0)*(-XABSC_WEIG_IN(I,1)+1.0D0)/4.0D0
       HJ(NS:NE)=(XABSC_WEIG_IN(:,1)+1.0D0)*(-XABSC_WEIG_IN(I,1)+1.0D0)/4.0D0
       HK(NS:NE)=(XABSC_WEIG_IN(:,1)+1.0D0)*(XABSC_WEIG_IN(I,1)+1.0D0)/4.0D0
! ------------------------------------------------
       dHIdxC(NS:NE)=(XABSC_WEIG_IN(I,1)-1.0D0)/2.0D0
       dHIdyC(NS:NE)=(XABSC_WEIG_IN(:,1)-1.0D0)/2.0D0
! ------------------------------------------------
       dHJdxC(NS:NE)=(-XABSC_WEIG_IN(I,1)+1.0D0)/2.0D0
       dHJdyC(NS:NE)=-(XABSC_WEIG_IN(:,1)+1.0D0)/2.0D0
! ------------------------------------------------
       dHKdxC(NS:NE)=(XABSC_WEIG_IN(I,1)+1.0D0)/2.0D0
       dHKdyC(NS:NE)=(XABSC_WEIG_IN(:,1)+1.0D0)/2.0D0
!------------------------------------------
       XTMPC(NS:NE)=XABSC_WEIG_IN(:,1)
       YTMPC(NS:NE)=XABSC_WEIG_IN(I,1)
    ENDDO
    HL=1.0D0-HI-HJ-HK
!
    DO I=1,NUMCPL
       Ele1=CPLINT(1,I)
       Ele2=CPLINT(2,I)
       Ele3=CPLINT(3,I)
       Ele4=CPLINT(4,I)
! 
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
       N4=NDMFS(1,Ele4)
!
       XM=(X(N1)+X(N2)+X(N3)+X(N4))/4.0D0
       YM=(Y(N1)+Y(N2)+Y(N3)+Y(N4))/4.0D0
!
       NDOFI=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFJ=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFK=NDOF(1,N3)+NDMFS(2,Ele3)
       NDOFL=NDOF(1,N4)+NDMFS(2,Ele4)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
       RADK=ELEPAR(1,Ele3)
       RADL=ELEPAR(1,Ele4)
!
       LENGTH2=SQRT((X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2)))
       LENGTH=LENGTH2/2.0D0
! ------------------------------------------------
       dHIdx=dHIdxC/LENGTH2
       dHIdy=dHIdyC/LENGTH2
! ------------------------------------------------
       dHJdx=dHJdxC/LENGTH2
       dHJdy=dHJdyC/LENGTH2
! ------------------------------------------------
       dHKdx=dHKdxC/LENGTH2
       dHKdy=dHKdyC/LENGTH2
! ------------------------------------------------
       dHLdx=-dHIdx-dHJdx-dHKdx
       dHLdy=-dHIdy-dHJdy-dHKdy
! ------------------------------------------------
       XTMP=XM+XTMPC*LENGTH
       YTMP=YM+YTMPC*LENGTH
!      Initialization
       dHJmdx=0.0D0
       dHJmdy=0.0D0
!
       dHKmdx=0.0D0
       dHKmdy=0.0D0
!
       dHLmdx=0.0D0
       dHLmdy=0.0D0     
!------------------------------------------
!      dHImdx,dHImdy
!------------------------------------------
      X_BAR=(XTMP-X(N1))/RADI
      Y_BAR=(YTMP-Y(N1))/RADI
      SELECT CASE(NDOFI)
      CASE(1)
         dHImdx(:,1)=2.0D0*dHIdx*HI
         dHImdy(:,1)=2.0D0*dHIdy*HI
!
         dHJmdx(:,1)=dHIdx*HJ+HI*dHJdx
         dHJmdy(:,1)=dHIdy*HJ+HI*dHJdy
!
         dHKmdx(:,1)=dHIdx*HK+HI*dHKdx
         dHKmdy(:,1)=dHIdy*HK+HI*dHKdy
!
         dHLmdx(:,1)=dHIdx*HL+HI*dHLdx
         dHLmdy(:,1)=dHIdy*HL+HI*dHLdy
      CASE(3)
         HC1=HI/RADI
!
         dHImdx(:,1)=dHIdx
         dHImdy(:,1)=dHIdy
!
         dHImdx(:,2)=dHIdx*X_BAR+HC1
         dHImdy(:,2)=dHIdy*X_BAR
!
         dHImdx(:,3)=dHIdx*Y_BAR
         dHImdy(:,3)=dHIdy*Y_BAR+HC1           
      CASE(4)
         HC1=HI/RADI
!
         dHImdx(:,1)=dHIdx
         dHImdy(:,1)=dHIdy
!
         dHImdx(:,2)=dHIdx*X_BAR+HC1
         dHImdy(:,2)=dHIdy*X_BAR
!
         dHImdx(:,3)=dHIdx*Y_BAR
         dHImdy(:,3)=dHIdy*Y_BAR+HC1
!
         dHImdx(:,4)=dHImdx(:,2)*Y_BAR
         dHImdy(:,4)=dHImdy(:,3)*X_BAR          
      CASE(6)
         HC1=HI/RADI
!
         dHImdx(:,1)=dHIdx
         dHImdy(:,1)=dHIdy
!
         dHImdx(:,2)=dHIdx*X_BAR+HC1
         dHImdy(:,2)=dHIdy*X_BAR
!
         dHImdx(:,3)=dHIdx*Y_BAR
         dHImdy(:,3)=dHIdy*Y_BAR+HC1          
!
         dHImdx(:,4)=dHImdx(:,2)*Y_BAR
         dHImdy(:,4)=dHImdy(:,3)*X_BAR
!
         dHImdx(:,5)=(dHImdx(:,2)+HC1)*X_BAR
         dHImdy(:,5)=dHImdy(:,2)*X_BAR
!
         dHImdx(:,6)=dHImdx(:,3)*Y_BAR
         dHImdy(:,6)=(dHImdy(:,3)+HC1)*Y_BAR              
      CASE DEFAULT
         WRITE(*,*) "ERROR, GENERATEM/OFE1_COUPLE_STIFF!"
         STOP
      END SELECT 
!---------------------------------------
!      dHJmdx,dHJmdy
!---------------------------------------
      X_BAR=(XTMP-X(N2))/RADJ
      Y_BAR=(YTMP-Y(N2))/RADJ
      SELECT CASE(NDOFJ)
      CASE(1)
         dHImdx(:,1)=dHImdx(:,1)+dHJdx*HI+HJ*dHIdx
         dHImdy(:,1)=dHImdy(:,1)+dHJdy*HI+HJ*dHIdy
!
         dHJmdx(:,1)=dHJmdx(:,1)+2.0D0*dHJdx*HJ
         dHJmdy(:,1)=dHJmdy(:,1)+2.0D0*dHJdy*HJ
!
         dHKmdx(:,1)=dHKmdx(:,1)+dHJdx*HK+HJ*dHKdx
         dHKmdy(:,1)=dHKmdy(:,1)+dHJdy*HK+HJ*dHKdy
!
         dHLmdx(:,1)=dHLmdx(:,1)+dHJdx*HL+HJ*dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+dHJdy*HL+HJ*dHLdy
      CASE(3)
         HC1=HJ/RADJ
!
         dHJmdx(:,1)=dHJmdx(:,1)+dHJdx
         dHJmdy(:,1)=dHJmdy(:,1)+dHJdy
!
         dHJmdx(:,2)=dHJdx*X_BAR+HC1
         dHJmdy(:,2)=dHJdy*X_BAR
!
         dHJmdx(:,3)=dHJdx*Y_BAR
         dHJmdy(:,3)=dHJdy*Y_BAR+HC1           
      CASE(4)
         HC1=HJ/RADJ
!
         dHJmdx(:,1)=dHJmdx(:,1)+dHJdx
         dHJmdy(:,1)=dHJmdy(:,1)+dHJdy
!
         dHJmdx(:,2)=dHJdx*X_BAR+HC1
         dHJmdy(:,2)=dHJdy*X_BAR
!
         dHJmdx(:,3)=dHJdx*Y_BAR
         dHJmdy(:,3)=dHJdy*Y_BAR+HC1 
!
         dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
         dHJmdy(:,4)=dHJmdy(:,3)*X_BAR          
      CASE(6)
         HC1=HJ/RADJ
!
         dHJmdx(:,1)=dHJmdx(:,1)+dHJdx
         dHJmdy(:,1)=dHJmdy(:,1)+dHJdy
!
         dHJmdx(:,2)=dHJdx*X_BAR+HC1
         dHJmdy(:,2)=dHJdy*X_BAR
!
         dHJmdx(:,3)=dHJdx*Y_BAR
         dHJmdy(:,3)=dHJdy*Y_BAR+HC1     
!
         dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
         dHJmdy(:,4)=dHJmdy(:,3)*X_BAR
!
         dHJmdx(:,5)=(dHJmdx(:,2)+HC1)*X_BAR
         dHJmdy(:,5)=dHJmdy(:,2)*X_BAR
!
         dHJmdx(:,6)=dHJmdx(:,3)*Y_BAR
         dHJmdy(:,6)=(dHJmdy(:,3)+HC1)*Y_BAR 
      CASE DEFAULT
         WRITE(*,*) "ERROR, GENERATEM/OFE1_COUPLE_STIFF!"
         STOP
      END SELECT 
!---------------------------------------
!      dHKmdx,dHKmdy
!---------------------------------------
      X_BAR=(XTMP-X(N3))/RADK
      Y_BAR=(YTMP-Y(N3))/RADK
      SELECT CASE(NDOFK)
      CASE(1)
         dHImdx(:,1)=dHImdx(:,1)+dHKdx*HI+HK*dHIdx
         dHImdy(:,1)=dHImdy(:,1)+dHKdy*HI+HK*dHIdy
!
         dHJmdx(:,1)=dHJmdx(:,1)+dHKdx*HJ+HK*dHJdx
         dHJmdy(:,1)=dHJmdy(:,1)+dHKdy*HJ+HK*dHJdy
!
         dHKmdx(:,1)=dHKmdx(:,1)+2.0D0*dHKdx*HK
         dHKmdy(:,1)=dHKmdy(:,1)+2.0D0*dHKdy*HK
!
         dHLmdx(:,1)=dHLmdx(:,1)+dHKdx*HL+HK*dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+dHKdy*HL+HK*dHLdy
      CASE(3)
         HC1=HK/RADK
!
         dHKmdx(:,1)=dHKmdx(:,1)+dHKdx
         dHKmdy(:,1)=dHKmdy(:,1)+dHKdy
!
         dHKmdx(:,2)=dHKdx*X_BAR+HC1
         dHKmdy(:,2)=dHKdy*X_BAR
!
         dHKmdx(:,3)=dHKdx*Y_BAR
         dHKmdy(:,3)=dHKdy*Y_BAR+HC1           
      CASE(4)
         HC1=HK/RADK
!
         dHKmdx(:,1)=dHKmdx(:,1)+dHKdx
         dHKmdy(:,1)=dHKmdy(:,1)+dHKdy
!
         dHKmdx(:,2)=dHKdx*X_BAR+HC1
         dHKmdy(:,2)=dHKdy*X_BAR
!
         dHKmdx(:,3)=dHKdx*Y_BAR
         dHKmdy(:,3)=dHKdy*Y_BAR+HC1 
!
         dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
         dHKmdy(:,4)=dHKmdy(:,3)*X_BAR          
      CASE(6)
         HC1=HK/RADK
!
         dHKmdx(:,1)=dHKmdx(:,1)+dHKdx
         dHKmdy(:,1)=dHKmdy(:,1)+dHKdy
!
         dHKmdx(:,2)=dHKdx*X_BAR+HC1
         dHKmdy(:,2)=dHKdy*X_BAR
!
         dHKmdx(:,3)=dHKdx*Y_BAR
         dHKmdy(:,3)=dHKdy*Y_BAR+HC1     
!
         dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
         dHKmdy(:,4)=dHKmdy(:,3)*X_BAR
!
         dHKmdx(:,5)=(dHKmdx(:,2)+HC1)*X_BAR
         dHKmdy(:,5)=dHKmdy(:,2)*X_BAR
!
         dHKmdx(:,6)=dHKmdx(:,3)*Y_BAR
         dHKmdy(:,6)=(dHKmdy(:,3)+HC1)*Y_BAR 
      CASE DEFAULT
         WRITE(*,*) "ERROR, GENERATEM/OFE1_COUPLE_STIFF!"
         STOP
      END SELECT 
!---------------------------------------
!      dHLmdx,dHLmdy
!---------------------------------------
      X_BAR=(XTMP-X(N4))/RADL
      Y_BAR=(YTMP-Y(N4))/RADL
      SELECT CASE(NDOFL)
      CASE(1)
         dHImdx(:,1)=dHImdx(:,1)+dHLdx*HI+HL*dHIdx
         dHImdy(:,1)=dHImdy(:,1)+dHLdy*HI+HL*dHIdy
!
         dHJmdx(:,1)=dHJmdx(:,1)+dHLdx*HJ+HL*dHJdx
         dHJmdy(:,1)=dHJmdy(:,1)+dHLdy*HJ+HL*dHJdy
!
         dHKmdx(:,1)=dHKmdx(:,1)+dHLdx*HK+HL*dHKdx
         dHKmdy(:,1)=dHKmdy(:,1)+dHLdy*HK+HL*dHKdy
!
         dHLmdx(:,1)=dHLmdx(:,1)+2.0D0*HL*dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+2.0D0*HL*dHLdy
      CASE(3)
         HC1=HL/RADL
!
         dHLmdx(:,1)=dHLmdx(:,1)+dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+dHLdy
!
         dHLmdx(:,2)=dHLdx*X_BAR+HC1
         dHLmdy(:,2)=dHLdy*X_BAR
!
         dHLmdx(:,3)=dHLdx*Y_BAR
         dHLmdy(:,3)=dHLdy*Y_BAR+HC1           
      CASE(4)
         HC1=HL/RADL
!
         dHLmdx(:,1)=dHLmdx(:,1)+dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+dHLdy
!
         dHLmdx(:,2)=dHLdx*X_BAR+HC1
         dHLmdy(:,2)=dHLdy*X_BAR
!
         dHLmdx(:,3)=dHLdx*Y_BAR
         dHLmdy(:,3)=dHLdy*Y_BAR+HC1 
!
         dHLmdx(:,4)=dHLmdx(:,2)*Y_BAR
         dHLmdy(:,4)=dHLmdy(:,3)*X_BAR          
      CASE(6)
         HC1=HL/RADL
!
         dHLmdx(:,1)=dHLmdx(:,1)+dHLdx
         dHLmdy(:,1)=dHLmdy(:,1)+dHLdy
!
         dHLmdx(:,2)=dHLdx*X_BAR+HC1
         dHLmdy(:,2)=dHLdy*X_BAR
!
         dHLmdx(:,3)=dHLdx*Y_BAR
         dHLmdy(:,3)=dHLdy*Y_BAR+HC1     
!
         dHLmdx(:,4)=dHLmdx(:,2)*Y_BAR
         dHLmdy(:,4)=dHLmdy(:,3)*X_BAR
!
         dHLmdx(:,5)=(dHLmdx(:,2)+HC1)*X_BAR
         dHLmdy(:,5)=dHLmdy(:,2)*X_BAR
!
         dHLmdx(:,6)=dHLmdx(:,3)*Y_BAR
         dHLmdy(:,6)=(dHLmdy(:,3)+HC1)*Y_BAR 
      CASE DEFAULT
         WRITE(*,*) "ERROR, GENERATEM/OFE1_COUPLE_STIFF!"
         STOP
      END SELECT 
      
!------------------------------------------
       CONST=LENGTH*LENGTH*THIC
       KCS1=KCS(1)*CONST
       KCS2=KCS(2)*CONST
       KCS3=KCS(3)*CONST
!      STIFFII
       NCOUNT=0
       DO J=1,NDOFI
          dHImdxW(:,J)=dHImdx(:,J)*GAUSS_WEI(:)
          dHImdyW(:,J)=dHImdy(:,J)*GAUSS_WEI(:)
          ADI(J)=NCOUNT
          DO K=J,NDOFI
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHImdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHImdy(:,K))  
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,1)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,1)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,1)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,1)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFFJJ
       NCOUNT=0
       DO J=1,NDOFJ
          dHJmdxW(:,J)=dHJmdx(:,J)*GAUSS_WEI(:)
          dHJmdyW(:,J)=dHJmdy(:,J)*GAUSS_WEI(:)
          ADJ(J)=NCOUNT
          DO K=J,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,2)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,2)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,2)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,2)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!
!      STIFFKK
!
       NCOUNT=0
       DO J=1,NDOFK
          dHKmdxW(:,J)=dHKmdx(:,J)*GAUSS_WEI(:)
          dHKmdyW(:,J)=dHKmdy(:,J)*GAUSS_WEI(:)
          ADK(J)=NCOUNT
          DO K=J,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHKmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHKmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHKmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHKmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,3)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,3)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,3)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,3)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!
!      STIFFLL
       NCOUNT=0
       DO J=1,NDOFL
          dHLmdxW(:,J)=dHLmdx(:,J)*GAUSS_WEI(:)
          dHLmdyW(:,J)=dHLmdy(:,J)*GAUSS_WEI(:)
          ADL(J)=NCOUNT
          DO K=J,NDOFL
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHLmdxW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHLmdxW(:,J)*dHLmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHLmdyW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHLmdyW(:,J)*dHLmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,4)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,4)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,4)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,4)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IJ
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,5)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,5)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,5)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,5)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IK
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,6)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,6)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,6)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,6)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IL
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFL
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHLmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHLmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,7)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,7)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,7)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,7)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF JK
       NCOUNT=0
       DO J=1,NDOFJ
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,8)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,8)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,8)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,8)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF JL
       NCOUNT=0
       DO J=1,NDOFJ
          DO K=1,NDOFL
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHLmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHLmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,9)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,9)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,9)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,9)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF KL
       NCOUNT=0
       DO J=1,NDOFK
          DO K=1,NDOFL
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHKmdxW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHKmdxW(:,J)*dHLmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHKmdyW(:,J)*dHLmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHKmdyW(:,J)*dHLmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,10)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,10)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,10)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,10)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!
       CALL OFE1_COUPLE_ADDBAN(Ele1,Ele2,Ele3,Ele4,N1,N2,N3,N4,NDOFI,NDOFJ,NDOFK,NDOFL,ADI,ADJ,ADK,ADL,STIFF)
    ENDDO
!
    DEALLOCATE(XABSC_WEIG_IN)  
!
    RETURN
!}}}
    END SUBROUTINE OFE1_COUPLE_STIFF 
!
    SUBROUTINE OFE1_CPL2_STIFF(NUMCPL2,NGDOF,THIC)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NUMCPL2,NGDOF
    REAL(8),INTENT(IN) :: THIC
!--------------------------------------!
!   IN-FUNCTION VARIABLES
!------------------------------------------
    INTEGER :: I,J,K
    INTEGER :: Ele1,Ele2,Ele3,Ele4,Ele5,Ele6
    INTEGER :: N1,N2,N3,NDOFI,NDOFJ,NDOFK
    REAL(8) :: RADI,RADJ,RADK
    REAL(8),DIMENSION(2) :: XY1,XY2,XY3,XYC1,XYC2,XYC3,XYCM1,XYCM2,XYCM3
    REAL(8) :: AREA2
    REAL(8) :: L1,L2,L3,L1_HF,L2_HF,L3_HF,LM1,LM2,LM3
    REAL(8),DIMENSION(3) :: S1,S2,S3,S4,WI,WJ,WK
    REAL(8),DIMENSION(2) :: PHI1,PHI2,PHI3,PHI4,PHI5,PHI6,PHI7,PHI8,PHI9,CONS1,CONS2,CONS3,CONS4,CONS5,CONS6
    REAL(8),DIMENSION(2) :: C1,C2,C3,C4,C5,C6
    REAL(8),DIMENSION(NITP_IN) :: LILJ,LJLK
    REAL(8),DIMENSION(NITP_IN) :: dLILJdx,dLILJdy,dLJLKdx,dLJLKdy
    REAL(8),DIMENSION(NITP_IN) :: HIW,dHIWdx,dHIWdy,HIWC1,HIWC2,HIWC3
    REAL(8),DIMENSION(NITP_IN) :: HJW,dHJWdx,dHJWdy,HJWC1,HJWC2,HJWC3
    REAL(8),DIMENSION(NITP_IN) :: HKW,dHKWdx,dHKWdy
    REAL(8) :: CONST
!
    REAL(8),DIMENSION(NITP_IN,3) :: AREA_COOR
    REAL(8),DIMENSION(NITP_IN)   :: GAUSS_WEI
    REAL(8),DIMENSION(NITP_IN)   :: NSP1,NSP2,NSP3,NSP4,NSP5,NSP6
    REAL(8),DIMENSION(NITP_IN)   :: DXNSP1,DXNSP2,DXNSP3,DXNSP4,DXNSP5,DXNSP6
    REAL(8),DIMENSION(NITP_IN)   :: DYNSP1,DYNSP2,DYNSP3,DYNSP4,DYNSP5,DYNSP6
    REAL(8)                      :: dLIdx,dLIdy,dLJdx,dLJdy,dLKdx,dLKdy
    REAL(8),DIMENSION(NITP_IN)   :: XTMP,YTMP,X_BAR,Y_BAR
!
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI1dx,dHI1dy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI2dx,dHI2dy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI3dx,dHI3dy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI1dxW,dHI1dyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI2dxW,dHI2dyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHI3dxW,dHI3dyW
!
    REAL(8),DIMENSION(NITP_IN) :: dHI4dx,dHI4dy
    REAL(8),DIMENSION(NITP_IN) :: dHI5dx,dHI5dy
    REAL(8),DIMENSION(NITP_IN) :: dHI6dx,dHI6dy
!
    REAL(8) :: DHIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: KCS1,KCS2,KCS3
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,21) :: STIFF
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK
    INTEGER :: NCOUNT,L
!
    CALL GAUTRIANGLE_OFE1(NGDOF,NITP_IN,AREA_COOR,GAUSS_WEI)
!   CALCULATE Pre-Set Value
    LILJ=AREA_COOR(:,1)*AREA_COOR(:,2)
    LJLK=AREA_COOR(:,2)*AREA_COOR(:,3)
!
    NSP1=2.0D0*AREA_COOR(:,1)*(AREA_COOR(:,1)-0.50D0)
    NSP2=2.0D0*AREA_COOR(:,2)*(AREA_COOR(:,2)-0.50D0)
    NSP3=2.0D0*AREA_COOR(:,3)*(AREA_COOR(:,3)-0.50D0)
    NSP4=4.0D0*AREA_COOR(:,1)*AREA_COOR(:,2)
    NSP5=4.0D0*AREA_COOR(:,3)*AREA_COOR(:,2)
    NSP6=4.0D0*AREA_COOR(:,3)*AREA_COOR(:,1)
!
    DO I=1,NUMCPL2
       Ele1=CPLINT(1,I)
       Ele2=CPLINT(2,I)
       Ele3=CPLINT(3,I)
       Ele4=CPLINT(4,I)
       Ele5=CPLINT(5,I)
       Ele6=CPLINT(6,I)
!
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
!
       NDOFI=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFJ=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFK=NDOF(1,N3)+NDMFS(2,Ele3)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
       RADK=ELEPAR(1,Ele3)
!
       XY1(1)=X(N1)
       XY1(2)=Y(N1)
       XY2(1)=X(N2)
       XY2(2)=Y(N2)
       XY3(1)=X(N3)
       XY3(2)=Y(N3)
! -----------------------------------
       XYC1=XY1-XY2
       XYC2=XY1-XY3
       XYC3=XY2-XY3
       XYCM1=XYC1+XYC2
       XYCM2=XYC3-XYC1
       XYCM3=XYC3+XYC2
! -----------------------------------
       L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
       L1_HF=L1/2.0D0
       L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
       L3_HF=L3/2.0D0
       L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
       L2_HF=L2/2.0D0
       LM1=SQRT(XYCM1(1)*XYCM1(1)+XYCM1(2)*XYCM1(2))/2.0D0
       LM2=SQRT(XYCM2(1)*XYCM2(1)+XYCM2(2)*XYCM2(2))/2.0D0
       LM3=SQRT(XYCM3(1)*XYCM3(1)+XYCM3(2)*XYCM3(2))/2.0D0
!
       AREA2=X(N1)*XYC3(2)-X(N2)*XYC2(2)+X(N3)*XYC1(2)
! -----------------------------------   WI
       S1=(/L1_HF,LM1,L3_HF/)/RADI 
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WI=1.0D0-S2+S3-S4
! -----------------------------------   WJ
       S1=(/L2_HF,LM2,L1_HF/)/RADJ
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WJ=1.0D0-S2+S3-S4 
! -----------------------------------   WK
       S1=(/L3_HF,LM3,L2_HF/)/RADK
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WK=1.0D0-S2+S3-S4
! -----------------------------------   PHI_I
       PHI1(1)=0.0D0
       PHI2(1)=0.0D0
       PHI3(1)=WI(1)/(WJ(3)+WI(1))
       PHI4(1)=WI(3)/(WK(1)+WI(3))
       PHI5(1)=WI(2)/(WI(2)+WJ(1)+WK(3))+WI(3)/(WI(3)+WJ(2)+WK(1))+WI(1)/(WI(1)+WJ(3)+WK(2))
       PHI6(1)=1.0D0
       PHI7(1)=1.0D0
! -----------------------------------   PHI_J
       PHI1(2)=0.0D0
       PHI2(2)=0.0D0
       PHI3(2)=WJ(1)/(WK(3)+WJ(1))
       PHI4(2)=WJ(3)/(WI(1)+WJ(3))
       PHI5(2)=WJ(2)/(WJ(2)+WK(1)+WI(3))+WJ(3)/(WJ(3)+WK(2)+WI(1))+WJ(1)/(WJ(1)+WK(3)+WI(2))
       PHI6(2)=1.0D0
       PHI7(2)=1.0D0
! -----------------------------------   197+46=243 FLOP
! -----------------------------------   C1 C2 C3 C4 C5 C6
       CONS1=-1.0D0
       CONS2=2.0D0*PHI4
       CONS3=2.0D0*(PHI3-PHI5)
       CONS4=-1.0D0+4.0D0*PHI4
       CONS5=CONS2
       CONS6=1.0D0+CONS2
!
       C1=2.0D0*(CONS3+CONS5)
       C2=2.0D0*(CONS2+CONS3)
       C3=4.0D0*(-CONS2+PHI5)
       C4=CONS4
       C5=2.0D0-4.0D0*PHI4
       C6=0.0D0 
!------------------------------------------  XTMP,YTMP 
       XTMP=AREA_COOR(:,1)*XY1(1)+AREA_COOR(:,2)*XY2(1)+AREA_COOR(:,3)*XY3(1)
       YTMP=AREA_COOR(:,1)*XY1(2)+AREA_COOR(:,2)*XY2(2)+AREA_COOR(:,3)*XY3(2) 
!------------------------------------------ HIW,HIWC1,HIWC2, HJW, HJWC1,HJWC2
       HIWC1=AREA_COOR(:,2)*C1(1)+AREA_COOR(:,1)*C2(1)+C3(1)
       HIWC2=C4(1)+AREA_COOR(:,1)*C5(1)
       HIW  =LILJ*HIWC1+AREA_COOR(:,1)*HIWC2 
!
       HJWC1=AREA_COOR(:,3)*C1(2)+AREA_COOR(:,2)*C2(2)+C3(2)
       HJWC2=C4(2)+AREA_COOR(:,2)*C5(2)
       HJW  =LJLK*HJWC1+AREA_COOR(:,2)*HJWC2
!
       HKW  =1.0D0-HIW-HJW
!-------------------------------------------
       dLIdx=XYC3(2)/AREA2 
       dLIdy=-XYC3(1)/AREA2
       dLJdx=-XYC2(2)/AREA2
       dLJdy=XYC2(1)/AREA2
       dLKdx=XYC1(2)/AREA2
       dLKdy=-XYC1(1)/AREA2
!------------------------------------------ dHIWdx,dHIWdy,dHJWdx,dHJWdy,dLILJdx,dLILJdy,dLJLKdx,dLJLKdy
       dLILJdx=dLIdx*AREA_COOR(:,2)+dLJdx*AREA_COOR(:,1)
       dLILJdy=dLIdy*AREA_COOR(:,2)+dLJdy*AREA_COOR(:,1)
!
       dLJLKdx=dLJdx*AREA_COOR(:,3)+dLKdx*AREA_COOR(:,2)
       dLJLKdy=dLJdy*AREA_COOR(:,3)+dLKdy*AREA_COOR(:,2)
!
       HIWC3=AREA_COOR(:,1)*C5(1)
       HJWC3=AREA_COOR(:,2)*C5(2)
!
       dHIWdx=dLILJdx*HIWC1+dLIdx*HIWC2+LILJ*(dLJdx*C1(1)+dLIdx*C2(1))+dLIdx*HIWC3
       dHIWdy=dLILJdy*HIWC1+dLIdy*HIWC2+LILJ*(dLJdy*C1(1)+dLIdy*C2(1))+dLIdy*HIWC3        
!
       dHJWdx=dLJLKdx*HJWC1+dLJdx*HJWC2+LJLK*(dLKdx*C1(2)+dLJdx*C2(2))+dLJdx*HJWC3
       dHJWdy=dLJLKdy*HJWC1+dLJdy*HJWC2+LJLK*(dLKdy*C1(2)+dLJdy*C2(2))+dLJdy*HJWC3
!
       dHKWdx=-dHIWdx-dHJWdx
       dHKWdy=-dHIWdy-dHJWdy          
!
       DXNSP1=(4.0D0*AREA_COOR(:,1)-1.0D0)*dLIdx
       DYNSP1=(4.0D0*AREA_COOR(:,1)-1.0D0)*dLIdy
       DXNSP2=(4.0D0*AREA_COOR(:,2)-1.0D0)*dLJdx
       DYNSP2=(4.0D0*AREA_COOR(:,2)-1.0D0)*dLJdy
       DXNSP3=(4.0D0*AREA_COOR(:,3)-1.0D0)*dLKdx
       DYNSP3=(4.0D0*AREA_COOR(:,3)-1.0D0)*dLKdy
       DXNSP4=4.0D0*(AREA_COOR(:,1)*dLJdx+dLIdx*AREA_COOR(:,2))
       DYNSP4=4.0D0*(AREA_COOR(:,1)*dLJdy+dLIdy*AREA_COOR(:,2))
       DXNSP5=4.0D0*(AREA_COOR(:,3)*dLJdx+dLKdx*AREA_COOR(:,2))
       DYNSP5=4.0D0*(AREA_COOR(:,3)*dLJdy+dLKdy*AREA_COOR(:,2))
       DXNSP6=4.0D0*(AREA_COOR(:,3)*dLIdx+dLKdx*AREA_COOR(:,1))
       DYNSP6=4.0D0*(AREA_COOR(:,3)*dLIdy+dLKdy*AREA_COOR(:,1))
!-----------------------------
!       HIW=AREA_COOR(:,1)
!       HJW=AREA_COOR(:,2)
!       HKW=AREA_COOR(:,3)
!       dHIWdx=dLIdx 
!       dHIWdy=dLIdy
!       dHJWdx=dLJdx
!       dHJWdy=dLJdy
!       dHKWdx=dLKdx
!       dHKWdy=dLKdy
! **************************** INITIALIZATION
       dHI1dx=0.0D0
       dHI1dy=0.0D0
       dHI2dx=0.0D0
       dHI2dy=0.0D0
       dHI3dx=0.0D0
       dHI3dy=0.0D0
       dHI4dx=0.0D0
       dHI4dy=0.0D0
       dHI5dx=0.0D0
       dHI5dy=0.0D0
       dHI6dx=0.0D0
       dHI6dy=0.0D0       
! **************************** dHI1mdx,dHI1mdy
       IF(NDOFI.EQ.1) THEN
          dHI1dx(:,1)=dHIWdx*NSP1+HIW*DXNSP1
          dHI1dy(:,1)=dHIWdy*NSP1+HIW*DYNSP1
!
          dHI2dx(:,1)=dHIWdx*NSP2+HIW*DXNSP2
          dHI2dy(:,1)=dHIWdy*NSP2+HIW*DYNSP2
!
          dHI3dx(:,1)=dHIWdx*NSP3+HIW*DXNSP3
          dHI3dy(:,1)=dHIWdy*NSP3+HIW*DYNSP3
!
          dHI4dx(:)=dHIWdx*NSP4+HIW*DXNSP4
          dHI4dy(:)=dHIWdy*NSP4+HIW*DYNSP4
!
          dHI5dx(:)=dHIWdx*NSP5+HIW*DXNSP5
          dHI5dy(:)=dHIWdy*NSP5+HIW*DYNSP5
!
          dHI6dx(:)=dHIWdx*NSP6+HIW*DXNSP6
          dHI6dy(:)=dHIWdy*NSP6+HIW*DYNSP6          
       ELSE
          X_BAR=(XTMP-X(N1))/RADI
          Y_BAR=(YTMP-Y(N1))/RADI
!
          SELECT CASE(NDOFI)
          CASE(3)
              dHI1dx(:,1)=dHIWdx
              dHI1dy(:,1)=dHIWdy
!
              dHI1dx(:,2)=dHIWdx*X_BAR+HIW/RADI
              dHI1dy(:,2)=dHIWdy*X_BAR
!
              dHI1dx(:,3)=dHIWdx*Y_BAR
              dHI1dy(:,3)=dHIWdy*Y_BAR+HIW/RADI              
          CASE(4)
              dHI1dx(:,1)=dHIWdx
              dHI1dy(:,1)=dHIWdy
!
              dHI1dx(:,2)=dHIWdx*X_BAR+HIW/RADI
              dHI1dy(:,2)=dHIWdy*X_BAR
!
              dHI1dx(:,3)=dHIWdx*Y_BAR
              dHI1dy(:,3)=dHIWdy*Y_BAR+HIW/RADI  
!
              dHI1dx(:,4)=dHI1dx(:,2)*Y_BAR
              dHI1dy(:,4)=dHI1dy(:,3)*X_BAR
          CASE(6)
              dHI1dx(:,1)=dHIWdx
              dHI1dy(:,1)=dHIWdy
!
              dHI1dx(:,2)=dHIWdx*X_BAR+HIW/RADI
              dHI1dy(:,2)=dHIWdy*X_BAR
!
              dHI1dx(:,3)=dHIWdx*Y_BAR
              dHI1dy(:,3)=dHIWdy*Y_BAR+HIW/RADI  
!
              dHI1dx(:,4)=dHI1dx(:,2)*Y_BAR
              dHI1dy(:,4)=dHI1dy(:,3)*X_BAR
!    
              dHI1dx(:,5)=(dHI1dx(:,2)+HIW/RADI)*X_BAR
              dHI1dy(:,5)=dHI1dy(:,2)*X_BAR
!
              dHI1dx(:,6)=dHI1dx(:,3)*Y_BAR
              dHI1dy(:,6)=(dHI1dy(:,3)+HIW/RADI)*Y_BAR   
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_CPL2_STIFF!"
             STOP
          END SELECT 
       ENDIF
! **************************** dHI2mdx,dHI2mdy
       IF(NDOFJ.EQ.1) THEN
          dHI1dx(:,1)=dHI1dx(:,1)+dHJWdx*NSP1+HJW*DXNSP1
          dHI1dy(:,1)=dHI1dy(:,1)+dHJWdy*NSP1+HJW*DYNSP1
!
          dHI2dx(:,1)=dHI2dx(:,1)+dHJWdx*NSP2+HJW*DXNSP2
          dHI2dy(:,1)=dHI2dy(:,1)+dHJWdy*NSP2+HJW*DYNSP2
!
          dHI3dx(:,1)=dHI3dx(:,1)+dHJWdx*NSP3+HJW*DXNSP3
          dHI3dy(:,1)=dHI3dy(:,1)+dHJWdy*NSP3+HJW*DYNSP3
!
          dHI4dx(:)=dHI4dx(:)+dHJWdx*NSP4+HJW*DXNSP4
          dHI4dy(:)=dHI4dy(:)+dHJWdy*NSP4+HJW*DYNSP4
!
          dHI5dx(:)=dHI5dx(:)+dHJWdx*NSP5+HJW*DXNSP5
          dHI5dy(:)=dHI5dy(:)+dHJWdy*NSP5+HJW*DYNSP5
!
          dHI6dx(:)=dHI6dx(:)+dHJWdx*NSP6+HJW*DXNSP6
          dHI6dy(:)=dHI6dy(:)+dHJWdy*NSP6+HJW*DYNSP6          
       ELSE
          X_BAR=(XTMP-X(N2))/RADJ
          Y_BAR=(YTMP-Y(N2))/RADJ
          SELECT CASE(NDOFJ)
          CASE(3)
              dHI2dx(:,1)=dHI2dx(:,1)+dHJWdx
              dHI2dy(:,1)=dHI2dy(:,1)+dHJWdy
!
              dHI2dx(:,2)=dHJWdx*X_BAR+HJW/RADJ
              dHI2dy(:,2)=dHJWdy*X_BAR
!
              dHI2dx(:,3)=dHJWdx*Y_BAR
              dHI2dy(:,3)=dHJWdy*Y_BAR+HJW/RADJ              
          CASE(4)
              dHI2dx(:,1)=dHI2dx(:,1)+dHJWdx
              dHI2dy(:,1)=dHI2dy(:,1)+dHJWdy
!
              dHI2dx(:,2)=dHJWdx*X_BAR+HJW/RADJ
              dHI2dy(:,2)=dHJWdy*X_BAR
!
              dHI2dx(:,3)=dHJWdx*Y_BAR
              dHI2dy(:,3)=dHJWdy*Y_BAR+HJW/RADJ  
!
              dHI2dx(:,4)=dHI2dx(:,2)*Y_BAR
              dHI2dy(:,4)=dHI2dy(:,3)*X_BAR
          CASE(6)
              dHI2dx(:,1)=dHI2dx(:,1)+dHJWdx
              dHI2dy(:,1)=dHI2dy(:,1)+dHJWdy
!
              dHI2dx(:,2)=dHJWdx*X_BAR+HJW/RADJ
              dHI2dy(:,2)=dHJWdy*X_BAR
!
              dHI2dx(:,3)=dHJWdx*Y_BAR
              dHI2dy(:,3)=dHJWdy*Y_BAR+HJW/RADJ    
!
              dHI2dx(:,4)=dHI2dx(:,2)*Y_BAR
              dHI2dy(:,4)=dHI2dy(:,3)*X_BAR
!    
              dHI2dx(:,5)=(dHI2dx(:,2)+HJW/RADJ)*X_BAR
              dHI2dy(:,5)=dHI2dy(:,2)*X_BAR
!     
              dHI2dx(:,6)=dHI2dx(:,3)*Y_BAR
              dHI2dy(:,6)=(dHI2dy(:,3)+HJW/RADJ)*Y_BAR   
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_CPL2_STIFF!"
             STOP
          END SELECT 
       ENDIF
! **************************** dHI3mdx,dHI3mdy
       IF(NDOFK.EQ.1) THEN
          dHI1dx(:,1)=dHI1dx(:,1)+dHKWdx*NSP1+HKW*DXNSP1
          dHI1dy(:,1)=dHI1dy(:,1)+dHKWdy*NSP1+HKW*DYNSP1
!
          dHI2dx(:,1)=dHI2dx(:,1)+dHKWdx*NSP2+HKW*DXNSP2
          dHI2dy(:,1)=dHI2dy(:,1)+dHKWdy*NSP2+HKW*DYNSP2
!
          dHI3dx(:,1)=dHI3dx(:,1)+dHKWdx*NSP3+HKW*DXNSP3
          dHI3dy(:,1)=dHI3dy(:,1)+dHKWdy*NSP3+HKW*DYNSP3
!
          dHI4dx(:)=dHI4dx(:)+dHKWdx*NSP4+HKW*DXNSP4
          dHI4dy(:)=dHI4dy(:)+dHKWdy*NSP4+HKW*DYNSP4
!
          dHI5dx(:)=dHI5dx(:)+dHKWdx*NSP5+HKW*DXNSP5
          dHI5dy(:)=dHI5dy(:)+dHKWdy*NSP5+HKW*DYNSP5
!
          dHI6dx(:)=dHI6dx(:)+dHKWdx*NSP6+HKW*DXNSP6
          dHI6dy(:)=dHI6dy(:)+dHKWdy*NSP6+HKW*DYNSP6          
       ELSE
          X_BAR=(XTMP-X(N3))/RADK
          Y_BAR=(YTMP-Y(N3))/RADK
          SELECT CASE(NDOFK)
          CASE(3)
              dHI3dx(:,1)=dHI3dx(:,1)+dHKWdx
              dHI3dy(:,1)=dHI3dy(:,1)+dHKWdy
!
              dHI3dx(:,2)=dHKWdx*X_BAR+HKW/RADK
              dHI3dy(:,2)=dHKWdy*X_BAR
!
              dHI3dx(:,3)=dHKWdx*Y_BAR
              dHI3dy(:,3)=dHKWdy*Y_BAR+HKW/RADK            
          CASE(4)
              dHI3dx(:,1)=dHI3dx(:,1)+dHKWdx
              dHI3dy(:,1)=dHI3dy(:,1)+dHKWdy
!
              dHI3dx(:,2)=dHKWdx*X_BAR+HKW/RADK
              dHI3dy(:,2)=dHKWdy*X_BAR
!
              dHI3dx(:,3)=dHKWdx*Y_BAR
              dHI3dy(:,3)=dHKWdy*Y_BAR+HKW/RADK  
!
              dHI3dx(:,4)=dHI3dx(:,2)*Y_BAR
              dHI3dy(:,4)=dHI3dy(:,3)*X_BAR
          CASE(6)
              dHI3dx(:,1)=dHI3dx(:,1)+dHKWdx
              dHI3dy(:,1)=dHI3dy(:,1)+dHKWdy
!
              dHI3dx(:,2)=dHKWdx*X_BAR+HKW/RADK
              dHI3dy(:,2)=dHKWdy*X_BAR
!
              dHI3dx(:,3)=dHKWdx*Y_BAR
              dHI3dy(:,3)=dHKWdy*Y_BAR+HKW/RADK 
!
              dHI3dx(:,4)=dHI3dx(:,2)*Y_BAR
              dHI3dy(:,4)=dHI3dy(:,3)*X_BAR
!    
              dHI3dx(:,5)=(dHI3dx(:,2)+HKW/RADK)*X_BAR
              dHI3dy(:,5)=dHI3dy(:,2)*X_BAR
!     
              dHI3dx(:,6)=dHI3dx(:,3)*Y_BAR
              dHI3dy(:,6)=(dHI3dy(:,3)+HKW/RADK)*Y_BAR   
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_CPL2_STIFF!"
             STOP
          END SELECT 
       ENDIF
!------------------------------------------
       CONST=AREA2*THIC/2.0D0
       KCS1=KCS(1)*CONST
       KCS2=KCS(2)*CONST
       KCS3=KCS(3)*CONST
!------------------------------------------ STIFFII,STIFFJJ,STIFFKK
!      STIFF11
       NCOUNT=0
       DO J=1,NDOFI
          dHI1dxW(:,J)=dHI1dx(:,J)*GAUSS_WEI(:)
          dHI1dyW(:,J)=dHI1dy(:,J)*GAUSS_WEI(:)
          ADI(J)=NCOUNT
          DO K=J,NDOFI
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI1dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI1dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI1dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI1dy(:,K))  
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,1)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,1)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,1)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,1)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))           
!      STIFF22
       NCOUNT=0
       DO J=1,NDOFJ
          dHI2dxW(:,J)=dHI2dx(:,J)*GAUSS_WEI(:)
          dHI2dyW(:,J)=dHI2dy(:,J)*GAUSS_WEI(:)
          ADJ(J)=NCOUNT
          DO K=J,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI2dxW(:,J)*dHI2dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI2dxW(:,J)*dHI2dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI2dyW(:,J)*dHI2dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI2dyW(:,J)*dHI2dy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,2)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,2)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,2)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,2)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF33
       NCOUNT=0
       DO J=1,NDOFK
          dHI3dxW(:,J)=dHI3dx(:,J)*GAUSS_WEI(:)
          dHI3dyW(:,J)=dHI3dy(:,J)*GAUSS_WEI(:)
          ADK(J)=NCOUNT
          DO K=J,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI3dxW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI3dxW(:,J)*dHI3dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI3dyW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI3dyW(:,J)*dHI3dy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,3)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,3)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,3)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,3)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF44
       DHIHJ(1,1)=SUM(dHI4dx(:)*dHI4dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI4dx(:)*dHI4dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI4dy(:)*dHI4dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI4dy(:)*dHI4dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,4)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,4)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,4)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,4)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1)) 
!      STIFF55
       DHIHJ(1,1)=SUM(dHI5dx(:)*dHI5dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI5dx(:)*dHI5dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI5dy(:)*dHI5dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI5dy(:)*dHI5dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,5)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,5)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,5)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,5)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1)) 
!      STIFF66
       DHIHJ(1,1)=SUM(dHI6dx(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI6dx(:)*dHI6dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI6dy(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI6dy(:)*dHI6dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,6)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,6)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,6)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,6)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1)) 
!      STIFF 12
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI2dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI2dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI2dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI2dy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,7)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,7)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,7)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,7)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 13
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI3dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI3dy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,8)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,8)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,8)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,8)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 14
       NCOUNT=0
       DO J=1,NDOFI
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI4dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI4dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,9)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,9)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,9)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,9)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 15
       NCOUNT=0
       DO J=1,NDOFI
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI5dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI5dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,10)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,10)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,10)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,10)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 16
       NCOUNT=0
       DO J=1,NDOFI
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI1dxW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI1dxW(:,J)*dHI6dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI1dyW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI1dyW(:,J)*dHI6dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,11)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,11)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,11)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,11)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 23
       NCOUNT=0
       DO J=1,NDOFJ
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHI2dxW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHI2dxW(:,J)*dHI3dy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHI2dyW(:,J)*dHI3dx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHI2dyW(:,J)*dHI3dy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,12)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,12)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,12)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,12)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 24
       NCOUNT=0
       DO J=1,NDOFJ
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI2dxW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI2dxW(:,J)*dHI4dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI2dyW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI2dyW(:,J)*dHI4dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,13)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,13)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,13)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,13)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 25
       NCOUNT=0
       DO J=1,NDOFJ
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI2dxW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI2dxW(:,J)*dHI5dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI2dyW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI2dyW(:,J)*dHI5dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,14)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,14)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,14)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,14)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 26
       NCOUNT=0
       DO J=1,NDOFJ
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI2dxW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI2dxW(:,J)*dHI6dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI2dyW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI2dyW(:,J)*dHI6dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,15)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,15)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,15)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,15)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 34
       NCOUNT=0
       DO J=1,NDOFK
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI3dxW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI3dxW(:,J)*dHI4dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI3dyW(:,J)*dHI4dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI3dyW(:,J)*dHI4dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,16)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,16)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,16)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,16)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF 35
       NCOUNT=0
       DO J=1,NDOFK
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI3dxW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI3dxW(:,J)*dHI5dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI3dyW(:,J)*dHI5dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI3dyW(:,J)*dHI5dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,17)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,17)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,17)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,17)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))
!      STIFF 36
       NCOUNT=0
       DO J=1,NDOFK
          NCOUNT=NCOUNT+1
          DHIHJ(NCOUNT,1)=SUM(dHI3dxW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,2)=SUM(dHI3dxW(:,J)*dHI6dy(:))
          DHIHJ(NCOUNT,3)=SUM(dHI3dyW(:,J)*dHI6dx(:))
          DHIHJ(NCOUNT,4)=SUM(dHI3dyW(:,J)*dHI6dy(:))           
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,18)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,18)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,18)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,18)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))
!      STIFF 45
       DHIHJ(1,1)=SUM(dHI4dx(:)*dHI5dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI4dx(:)*dHI5dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI4dy(:)*dHI5dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI4dy(:)*dHI5dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,19)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,19)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,19)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,19)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1))
!      STIFF 46
       DHIHJ(1,1)=SUM(dHI4dx(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI4dx(:)*dHI6dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI4dy(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI4dy(:)*dHI6dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,20)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,20)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,20)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,20)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1))
!      STIFF 56
       DHIHJ(1,1)=SUM(dHI5dx(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,2)=SUM(dHI5dx(:)*dHI6dy(:)*GAUSS_WEI(:))
       DHIHJ(1,3)=SUM(dHI5dy(:)*dHI6dx(:)*GAUSS_WEI(:))
       DHIHJ(1,4)=SUM(dHI5dy(:)*dHI6dy(:)*GAUSS_WEI(:))           
!
       STIFF(1,1,21)=(KCS1*DHIHJ(1,1)+KCS3*DHIHJ(1,4))
       STIFF(1,2,21)=(KCS2*DHIHJ(1,2)+KCS3*DHIHJ(1,3))
       STIFF(1,3,21)=(KCS2*DHIHJ(1,3)+KCS3*DHIHJ(1,2))
       STIFF(1,4,21)=(KCS1*DHIHJ(1,4)+KCS3*DHIHJ(1,1))
!------------------------------------------ 
       CALL OFE1_CPL2_ADDBAN(Ele1,Ele2,Ele3,Ele4,Ele5,Ele6,N1,N2,N3,NDOFI,NDOFJ,NDOFK,ADI,ADJ,ADK,STIFF)
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OFE1_CPL2_STIFF 
!
    SUBROUTINE OFE1_CPL2_ADDBAN(Ele1,Ele2,Ele3,Ele4,Ele5,Ele6,N1,N2,N3,NDOFI,NDOFJ,NDOFK,ADI,ADJ,ADK,STIFF)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: Ele1,Ele2,Ele3,Ele4,Ele5,Ele6
    INTEGER,INTENT(IN)          :: N1,N2,N3,NDOFI,NDOFJ,NDOFK
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,21) :: STIFF
!   IN-FUNCTION VARIABLES
    INTEGER :: N4,N5,N6
    INTEGER,DIMENSION(2) :: NDOFIXY,NDOFJXY,NDOFKXY
    INTEGER,DIMENSION(2) :: IDI,IDJ,IDK
    INTEGER :: I,J,M,N
    INTEGER*8 :: LL1
    INTEGER :: KI,KJ,KK2,KK1
    LOGICAL :: LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10,LG11,LG12,LG13,LG14,LG15
!
    N4=NDMFS(1,Ele4)
    N5=NDMFS(1,Ele5)
    N6=NDMFS(1,Ele6)
!
    NDOFIXY(1)=1+NDMFS(2,Ele1)
    NDOFIXY(2)=1+NDMFS(3,Ele1)
    IDI(1)=ID(1,N1)
    IDI(2)=ID(2,N1)
! 
    NDOFJXY(1)=1+NDMFS(2,Ele2)
    NDOFJXY(2)=1+NDMFS(3,Ele2)
    IDJ(1)=ID(1,N2)
    IDJ(2)=ID(2,N2)
!
    NDOFKXY(1)=1+NDMFS(2,Ele3)
    NDOFKXY(2)=1+NDMFS(3,Ele3)
    IDK(1)=ID(1,N3)
    IDK(2)=ID(2,N3)
!*******************************************************************
    LG1=N1.LT.N2
    LG2=N1.LT.N3
    LG3=N2.LT.N3
    LG4=N1.LT.N4
    LG5=N1.LT.N5
    LG6=N1.LT.N6
    LG7=N2.LT.N4
    LG8=N2.LT.N5
    LG9=N2.LT.N6
    LG10=N3.LT.N4
    LG11=N3.LT.N5
    LG12=N3.LT.N6
    LG13=N4.LT.N5
    LG14=N4.LT.N6
    LG15=N5.LT.N6    
!*******************************************************************
! --------------------------------------------   STIFF11
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(1)+N-NDOFIXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(2),NDOFI
       KJ=IDI(2)+M-NDOFIXY(2)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       DO N=NDOFIXY(2),NDOFI
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADI(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,1)
          ELSE
             KK2=ADI(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,1)
          ENDIF
       ENDDO
    ENDDO
! --------------------------------------------   STIFF22
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(1)+N-NDOFJXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(2),NDOFJ
       KJ=IDJ(2)+M-NDOFJXY(2)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       DO N=NDOFJXY(2),NDOFJ
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADJ(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,2)
          ELSE
             KK2=ADJ(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,2)
          ENDIF
       ENDDO
    ENDDO
! --------------------------------------------   STIFF33
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(1)+N-NDOFKXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(2),NDOFK
       KJ=IDK(2)+M-NDOFKXY(2)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       DO N=NDOFKXY(2),NDOFK
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADK(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,3)
          ELSE
             KK2=ADK(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,3)
          ENDIF
       ENDDO
    ENDDO
! --------------------------------------------   STIFF44
    IF (ID(1,N4).NE.0) THEN
       KI=ID(1,N4)
       LL1=MAXA_COLSOL(KI)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,1,4)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,1,4)
    ENDIF
!   
    IF (ID(2,N4).NE.0) THEN
       KJ=ID(2,N4)
       LL1=MAXA_COLSOL(KJ)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,4,4)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,4,4)
    ENDIF
!
    IF ((ID(1,N4).NE.0).AND.(ID(2,N4).NE.0)) THEN
       LL1=MAXA_COLSOL(KJ)+KJ-KI
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,3,4)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,3,4)
    ENDIF 
! --------------------------------------------   STIFF55
    IF (ID(1,N5).NE.0) THEN
       KI=ID(1,N5)
       LL1=MAXA_COLSOL(KI)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,1,5)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,1,5)
    ENDIF
!   
    IF (ID(2,N5).NE.0) THEN
       KJ=ID(2,N5)
       LL1=MAXA_COLSOL(KJ)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,4,5)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,4,5)
    ENDIF
!
    IF ((ID(1,N5).NE.0).AND.(ID(2,N5).NE.0)) THEN
       LL1=MAXA_COLSOL(KJ)+KJ-KI
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,3,5)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,3,5)
    ENDIF 
! --------------------------------------------   STIFF66
    IF (ID(1,N6).NE.0) THEN
       KI=ID(1,N6)
       LL1=MAXA_COLSOL(KI)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,1,6)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,1,6)
    ENDIF
!   
    IF (ID(2,N6).NE.0) THEN
       KJ=ID(2,N6)
       LL1=MAXA_COLSOL(KJ)
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,4,6)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,4,6)
    ENDIF
!
    IF ((ID(1,N6).NE.0).AND.(ID(2,N6).NE.0)) THEN
       LL1=MAXA_COLSOL(KJ)+KJ-KI
!
       K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,3,6)
       KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,3,6)
    ENDIF 
! --------------------------------------------   STIFF12 STIFF13 STIFF23
    KK1=0
    DO I=1,2
       DO J=1,2
          KK1=KK1+1
!   STIFF12
          IF(LG1) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,7)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,7)                
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,7)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,7)                
                ENDDO
             ENDDO 
          ENDIF
!   STIFF13
          IF(LG2) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,8) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,8)               
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,8) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,8)               
                ENDDO
             ENDDO
          ENDIF
!   STIFF14
          IF(LG4) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,9)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,9)                
                ENDIF
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,9)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,9)                
                ENDIF
             ENDDO 
          ENDIF
!   STIFF15
          IF(LG5) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,10)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,10)                
                ENDIF
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,10)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,10)                
                ENDIF
             ENDDO 
          ENDIF
!   STIFF16
          IF(LG6) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,11)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,11)                
                ENDIF
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,11)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,11)                
                ENDIF
             ENDDO 
          ENDIF
!   STIFF23
          IF(LG3) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,12) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,12)               
                ENDDO
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,12) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,12)               
                ENDDO
             ENDDO
          ENDIF
!   STIFF24
          IF(LG7) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,13) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,13)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,13) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,13)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF25
          IF(LG8) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,14) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,14)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,14) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,14)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF26
          IF(LG9) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,15) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,15)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,15) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,15)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF34
          IF(LG10) THEN
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,16) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,16)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N4).NE.0) THEN
                   KI=ID(J,N4)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,16) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,16)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF35
          IF(LG11) THEN
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,17) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,17)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,17) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,17)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF36
          IF(LG12) THEN
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,18) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,18)               
                ENDIF
             ENDDO
          ELSE
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(M,KK1,18) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(M,KK1,18)               
                ENDIF
             ENDDO
          ENDIF
!   STIFF45
          IF(LG13) THEN
             IF(ID(I,N4).NE.0) THEN
                KJ=ID(I,N4)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,19) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,19)               
                ENDIF
             ENDIF
          ELSE
             IF(ID(I,N4).NE.0) THEN
                KJ=ID(I,N4)
                IF(ID(J,N5).NE.0) THEN
                   KI=ID(J,N5)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,19) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,19)               
                ENDIF
             ENDIF
          ENDIF
!   STIFF46
          IF(LG14) THEN
             IF(ID(I,N4).NE.0) THEN
                KJ=ID(I,N4)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,20) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,20)               
                ENDIF
             ENDIF
          ELSE
             IF(ID(I,N4).NE.0) THEN
                KJ=ID(I,N4)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,20) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,20)               
                ENDIF
             ENDIF
          ENDIF
!   STIFF56
          IF(LG15) THEN
             IF(ID(I,N5).NE.0) THEN
                KJ=ID(I,N5)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,21) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,21)               
                ENDIF
             ENDIF
          ELSE
             IF(ID(I,N5).NE.0) THEN
                KJ=ID(I,N5)
                IF(ID(J,N6).NE.0) THEN
                   KI=ID(J,N6)
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(1,KK1,21) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(1,KK1,21)               
                ENDIF
             ENDIF
          ENDIF
!
       ENDDO
    ENDDO
!*******************************************************************
!
    RETURN
!}}}
    END SUBROUTINE  OFE1_CPL2_ADDBAN
!
    SUBROUTINE OFE1_COUPLE_ADDBAN(Ele1,Ele2,Ele3,Ele4,N1,N2,N3,N4,NDOFI,NDOFJ,NDOFK,NDOFL,ADI,ADJ,ADK,ADL,STIFF)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: Ele1,Ele2,Ele3,Ele4
    INTEGER,INTENT(IN) :: N1,N2,N3,N4
    INTEGER,INTENT(IN) :: NDOFI,NDOFJ,NDOFK,NDOFL
!
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK,ADL
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,10) :: STIFF
!   IN-FUNCTION VARIABLES
    INTEGER,DIMENSION(2) :: NDOFIXY,NDOFJXY,NDOFKXY,NDOFLXY
    INTEGER,DIMENSION(2) :: IDI,IDJ,IDK,IDL
    INTEGER :: I,J,M,N
    INTEGER*8 :: LL1
    INTEGER :: KI,KJ,KK2,KK1
    LOGICAL :: LOG1,LOG2,LOG3,LOG4,LOG5,LOG6
!
    NDOFIXY(1)=1+NDMFS(2,Ele1)
    NDOFIXY(2)=1+NDMFS(3,Ele1)
    IDI(1)=ID(1,N1)
    IDI(2)=ID(2,N1)
! 
    NDOFJXY(1)=1+NDMFS(2,Ele2)
    NDOFJXY(2)=1+NDMFS(3,Ele2)
    IDJ(1)=ID(1,N2)
    IDJ(2)=ID(2,N2)
!
    NDOFKXY(1)=1+NDMFS(2,Ele3)
    NDOFKXY(2)=1+NDMFS(3,Ele3)
    IDK(1)=ID(1,N3)
    IDK(2)=ID(2,N3)
!
    NDOFLXY(1)=1+NDMFS(2,Ele4)
    NDOFLXY(2)=1+NDMFS(3,Ele4)
    IDL(1)=ID(1,N4)
    IDL(2)=ID(2,N4)
!*******************************************************************
    LOG1=N1.LT.N2
    LOG2=N1.LT.N3
    LOG3=N1.LT.N4
    LOG4=N2.LT.N3
    LOG5=N2.LT.N4
    LOG6=N3.LT.N4
!*******************************************************************
!   STIFFII
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(1)+N-NDOFIXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(2),NDOFI
       KJ=IDI(2)+M-NDOFIXY(2)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       DO N=NDOFIXY(2),NDOFI
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADI(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,1)
          ELSE
             KK2=ADI(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,1)
          ENDIF
       ENDDO
    ENDDO
!   STIFFJJ
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(1)+N-NDOFJXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(2),NDOFJ
       KJ=IDJ(2)+M-NDOFJXY(2)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       DO N=NDOFJXY(2),NDOFJ
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADJ(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,2)
          ELSE
             KK2=ADJ(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,2)
          ENDIF
       ENDDO
    ENDDO
!   STIFFKK
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(1)+N-NDOFKXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(2),NDOFK
       KJ=IDK(2)+M-NDOFKXY(2)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       DO N=NDOFKXY(2),NDOFK
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADK(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,3)
          ELSE
             KK2=ADK(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,3)
          ENDIF
       ENDDO
    ENDDO
!
!   STIFFLL
    DO M=NDOFLXY(1),NDOFL
       KJ=IDL(1)+M-NDOFLXY(1)
       KK2=ADL(M)
       DO N=M,NDOFL
          KK2=KK2+1
          KI=IDL(1)+N-NDOFLXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,4)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,4)
       ENDDO
    ENDDO
!
    DO M=NDOFLXY(2),NDOFL
       KJ=IDL(2)+M-NDOFLXY(2)
       KK2=ADL(M)
       DO N=M,NDOFL
          KK2=KK2+1
          KI=IDL(2)+N-NDOFLXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,4)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,4)
       ENDDO
    ENDDO
!
    DO M=NDOFLXY(1),NDOFL
       KJ=IDL(1)+M-NDOFLXY(1)
       DO N=NDOFLXY(2),NDOFL
          KI=IDL(2)+N-NDOFLXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADL(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,4)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,4)
          ELSE
             KK2=ADL(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,4)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,4)
          ENDIF
       ENDDO
    ENDDO
!*******************************************************************
    KK1=0
    DO I=1,2
       DO J=1,2
          KK1=KK1+1
!   STIFFIJ
          IF(LOG1) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,5)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,5)                
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,5)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,5)                
                ENDDO
             ENDDO 
          ENDIF
!   STIFFIK
          IF(LOG2) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,6) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,6)               
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,6) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,6)               
                ENDDO
             ENDDO
          ENDIF
!   STIFFIL
          IF(LOG3) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,7) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,7)               
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,7) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,7)               
                ENDDO
             ENDDO
          ENDIF
!   STIFFJK
          IF(LOG4) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,8) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,8)               
                ENDDO
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,8) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,8)               
                ENDDO
             ENDDO
          ENDIF
!   STIFFJL
          IF(LOG5) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,9) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,9)               
                ENDDO
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,9) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,9)               
                ENDDO
             ENDDO
          ENDIF
!   STIFFKL
          IF(LOG6) THEN
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,10) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,10)               
                ENDDO
             ENDDO
          ELSE
             DO M=NDOFKXY(I),NDOFK
                KJ=IDK(I)+M-NDOFKXY(I)
                DO N=NDOFLXY(J),NDOFL
                   KI=IDL(J)+N-NDOFLXY(J)
                   KK2=(M-1)*NDOFL+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,10) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,10)               
                ENDDO
             ENDDO
          ENDIF
!           
       ENDDO
     ENDDO
!*******************************************************************
!
    RETURN
!}}}
    END SUBROUTINE OFE1_COUPLE_ADDBAN
!
    SUBROUTINE OFE1_QUADS_STIFF(NUMOFE,NGDOF,THIC)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NUMOFE,NGDOF
    REAL(8),INTENT(IN) :: THIC
! IN-FUNCTION VARIABLES
!------------------------------------------
    REAL(8),DIMENSION(NITP_IN)   :: LILJ,LJLK
    REAL(8),DIMENSION(NITP_IN,3) :: AREA_COOR
    REAL(8),DIMENSION(NITP_IN)   :: GAUSS_WEI
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,6) :: STIFF
!------------------------------------------
    INTEGER              :: I,J,K,L
    INTEGER              :: Ele1,Ele2,Ele3,N1,N2,N3
    REAL(8)              :: RADI,RADJ,RADK
    REAL(8),DIMENSION(2) :: XY1,XY2,XY3,XYC1,XYC2,XYC3,XYCM1,XYCM2,XYCM3
    REAL(8)              :: L1,L2,L3,L1_HF,L2_HF,L3_HF,LM1,LM2,LM3,AREA2
    REAL(8),DIMENSION(3) :: S1,S2,S3,S4,WI,WJ,WK
    REAL(8),DIMENSION(2) :: PHI1,PHI2,PHI3,PHI4,PHI5,PHI6,PHI7,PHI8,PHI9,CONS1,CONS2,CONS3,CONS4,CONS5,CONS6
    REAL(8),DIMENSION(2) :: C1,C2,C3,C4,C5,C6
    REAL(8),DIMENSION(NITP_IN) :: XTMP,YTMP,X_BAR,Y_BAR
!
    REAL(8),DIMENSION(NITP_IN) :: HIW,dHIWdx,dHIWdy,HIWC1,HIWC2,HIWC3
    REAL(8),DIMENSION(NITP_IN) :: HJW,dHJWdx,dHJWdy,HJWC1,HJWC2,HJWC3
    REAL(8),DIMENSION(NITP_IN) :: HKW,dHKWdx,dHKWdy
!
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHImdx,dHImdy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHJmdx,dHJmdy
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHKmdx,dHKmdy
!
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHImdxW,dHImdyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHJmdxW,dHJmdyW
    REAL(8),DIMENSION(NITP_IN,NPOR_MAX) :: dHKmdxW,dHKmdyW
!
    REAL(8),DIMENSION(NITP_IN) :: dLILJdx,dLILJdy,dLJLKdx,dLJLKdy
    REAL(8)                    :: dLIdx,dLIdy,dLJdx,dLJdy,dLKdx,dLKdy
    INTEGER :: INDBCI,INDBCJ,INDBCK
    REAL(8) :: NXX,NYY
    INTEGER :: NDOFI,NDOFJ,NDOFK
    INTEGER :: NCOUNT
    REAL(8) :: CONST
    REAL(8),DIMENSION(NITP_IN)  :: CONSTX,CONSTY
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK
!
    REAL(8) :: DHIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: KCS1,KCS2,KCS3
!
    CALL GAUTRIANGLE_OFE1(NGDOF,NITP_IN,AREA_COOR,GAUSS_WEI)
!   CALCULATE Pre-Set Value
    LILJ=AREA_COOR(:,1)*AREA_COOR(:,2)
    LJLK=AREA_COOR(:,2)*AREA_COOR(:,3)
!
    DO I=1,NumOFE       
       Ele1=NDINT(1,I)
       Ele2=NDINT(2,I)
       Ele3=NDINT(3,I)
! 
       N1=NDMFS(1,Ele1)
       N2=NDMFS(1,Ele2)
       N3=NDMFS(1,Ele3)
!
       NDOFI=NDOF(1,N1)+NDMFS(2,Ele1)
       NDOFJ=NDOF(1,N2)+NDMFS(2,Ele2)
       NDOFK=NDOF(1,N3)+NDMFS(2,Ele3)
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
       RADK=ELEPAR(1,Ele3)
!
       INDBCI=NDMFS(4,Ele1)
       INDBCJ=NDMFS(4,Ele2)
       INDBCK=NDMFS(4,Ele3)
!
       XY1(1)=X(N1)
       XY1(2)=Y(N1)
       XY2(1)=X(N2)
       XY2(2)=Y(N2)
       XY3(1)=X(N3)
       XY3(2)=Y(N3)
! -----------------------------------
       XYC1=XY1-XY2
       XYC2=XY1-XY3
       XYC3=XY2-XY3
       XYCM1=XYC1+XYC2
       XYCM2=XYC3-XYC1
       XYCM3=XYC3+XYC2
! -----------------------------------
       L1=SQRT(XYC1(1)*XYC1(1)+XYC1(2)*XYC1(2))
       L1_HF=L1/2.0D0
       L3=SQRT(XYC2(1)*XYC2(1)+XYC2(2)*XYC2(2))
       L3_HF=L3/2.0D0
       L2=SQRT(XYC3(1)*XYC3(1)+XYC3(2)*XYC3(2))
       L2_HF=L2/2.0D0
       LM1=SQRT(XYCM1(1)*XYCM1(1)+XYCM1(2)*XYCM1(2))/2.0D0
       LM2=SQRT(XYCM2(1)*XYCM2(1)+XYCM2(2)*XYCM2(2))/2.0D0
       LM3=SQRT(XYCM3(1)*XYCM3(1)+XYCM3(2)*XYCM3(2))/2.0D0
! -----------------------------------   42 FLOP
       AREA2=X(N1)*XYC3(2)-X(N2)*XYC2(2)+X(N3)*XYC1(2)
! -----------------------------------   WI
!       S1=(/L1,L3,L1_HF,LM1,L3_HF/)/RADI 
!       S2=S1*S1
!       S3=S2*S1
!       S4=S3*S1
!!
!       S2=S2*6.0D0
!       S3=S3*8.0D0
!       S4=S4*3.0D0
!!
!       WI=1.0D0-S2+S3-S4
!! -----------------------------------   WJ
!       S1=(/L2,L1,L2_HF,LM2,L1_HF/)/RADJ
!       S2=S1*S1
!       S3=S2*S1
!       S4=S3*S1
!!
!       S2=S2*6.0D0
!       S3=S3*8.0D0
!       S4=S4*3.0D0
!!
!       WJ=1.0D0-S2+S3-S4 
!! -----------------------------------   WK
!       S1=(/L3,L2,L3_HF,LM3,L2_HF/)/RADK
!       S2=S1*S1
!       S3=S2*S1
!       S4=S3*S1
!!
!       S2=S2*6.0D0
!       S3=S3*8.0D0
!       S4=S4*3.0D0
!!
!       WK=1.0D0-S2+S3-S4 
!
       S1=(/L1_HF,LM1,L3_HF/)/RADI 
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WI=1.0D0-S2+S3-S4
! -----------------------------------   WJ
       S1=(/L2_HF,LM2,L1_HF/)/RADJ
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WJ=1.0D0-S2+S3-S4 
! -----------------------------------   WK
       S1=(/L3_HF,LM3,L2_HF/)/RADK
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WK=1.0D0-S2+S3-S4 
!
       IF(WI(2).LT.0.0D0) WI(2)=0.0D0
       IF(WJ(2).LT.0.0D0) WJ(2)=0.0D0
       IF(WK(2).LT.0.0D0) WK(2)=0.0D0
! -----------------------------------   42+155=197 FLOP
! -----------------------------------   PHI_I
! ----------------------------------    PHI_J
!      IF ((INDBCI.NE.0).OR.(INDBCJ.NE.0).OR.(INDBCK.NE.0)) THEN
!          SELECT CASE(INDBCI)
!          CASE(0)
!                PHI1(1)=WI(1)/(WI(1)+1.0D0)
!                PHI2(1)=WI(2)/(WI(2)+1.0D0)
!                PHI3(1)=WI(3)/(WJ(5)+WI(3))
!                PHI4(1)=WI(5)/(WK(3)+WI(5))
!                PHI5(1)=WI(4)/(WI(4)+WJ(3)+WK(5))
!                PHI6(1)=0.0D0
!                PHI7(1)=0.0D0
!                PHI8(1)=0.0D0
!                PHI9(1)=0.0D0
!    !
!                PHI1(2)=0.0D0
!                PHI2(2)=0.0D0
!                PHI3(2)=0.0D0
!                PHI4(2)=0.0D0
!                PHI5(2)=WJ(3)/(WJ(3)+WK(5)+WI(4))
!                PHI6(2)=0.0D0
!                PHI7(2)=1.0D0/(WI(1)+1.0D0)            
!                PHI8(2)=0.0D0
!                PHI9(2)=1-PHI3(1)
!          CASE DEFAULT
!                PHI1(1)=1.0D0
!                PHI2(1)=1.0D0
!                PHI3(1)=1.0D0
!                PHI4(1)=1.0D0
!                PHI5(1)=1.0D0
!                PHI6(1)=0.0D0
!                PHI7(1)=0.0D0
!                PHI8(1)=0.0D0
!                PHI9(1)=0.0D0
!    !
!                PHI1(2)=0.0D0
!                PHI2(2)=0.0D0
!                PHI3(2)=0.0D0
!                PHI4(2)=0.0D0
!                PHI5(2)=0.0D0
!                PHI6(2)=0.0D0
!                PHI7(2)=0.0D0
!                PHI8(2)=0.0D0
!                PHI9(2)=0.0D0
!          END SELECT
!    !
!          SELECT CASE(INDBCJ)
!          CASE(0)
!                PHI5(1)=PHI5(1)+WI(5)/(WI(5)+WJ(4)+WK(3))
!                PHI6(1)=1.0D0/(WJ(2)+1.0D0)
!                PHI8(1)=WI(3)/(WJ(5)+WI(3))
!    !
!                PHI1(2)=WJ(1)/(WJ(1)+1.0D0)
!                PHI2(2)=WJ(2)/(WJ(2)+1.0D0)
!                PHI3(2)=WJ(3)/(WK(5)+WJ(3))
!                PHI4(2)=WJ(5)/(WI(3)+WJ(5))
!                PHI5(2)=PHI5(2)+WJ(4)/(WJ(4)+WK(3)+WI(5))           
!          CASE DEFAULT
!                PHI1(2)=1.0D0
!                PHI2(2)=1.0D0
!                PHI3(2)=1.0D0
!                PHI4(2)=1.0D0
!                PHI5(2)=PHI5(2)+1.0D0          
!          END SELECT
!    !
!          SELECT CASE(INDBCK)
!          CASE(0)
!                PHI5(1)=PHI5(1)+WI(3)/(WI(3)+WJ(5)+WK(4))
!                PHI7(1)=1.0D0/(WK(1)+1.0D0)
!                PHI9(1)=WI(5)/(WK(3)+WI(5))
!    !                
!                PHI5(2)=PHI5(2)+WJ(5)/(WJ(5)+WK(4)+WI(3))
!                PHI6(2)=1.0D0/(WK(2)+1.0D0)
!                PHI8(2)=WJ(3)/(WK(5)+WJ(3))
!          CASE DEFAULT
!                CONTINUE
!          END SELECT
!!          
!          C1=2.0D0*(PHI1+PHI2+2.0D0*(PHI8+PHI9-PHI5))
!          C2=2.0D0*(2.0D0*(PHI3+PHI2+2.0D0*PHI9-PHI5-PHI4)+PHI6-PHI7)
!          C3=-8.0D0*PHI9+PHI7-3.0D0*PHI2-PHI1+4.0D0*PHI5-PHI6
!          C4=4.0D0*PHI9+PHI2-PHI7
!          C5=3.0D0*(PHI7-PHI2)-8.0D0*PHI9+4.0D0*PHI4-1.0D0
!          C6=2.0D0*(-PHI7+PHI2-2.0D0*(PHI4-PHI9)+1.0D0)
!        ELSE
! -----------------------------------   PHI_I
       PHI1(1)=0.0D0
       PHI2(1)=0.0D0
       PHI3(1)=WI(1)/(WJ(3)+WI(1))
       PHI4(1)=WI(3)/(WK(1)+WI(3))
       PHI5(1)=WI(2)/(WI(2)+WJ(1)+WK(3))+WI(3)/(WI(3)+WJ(2)+WK(1))+WI(1)/(WI(1)+WJ(3)+WK(2))
       PHI6(1)=1.0D0
       PHI7(1)=1.0D0
!       PHI1(1)=WI(1)/(WI(1)+1.0D0)
!       PHI2(1)=WI(2)/(WI(2)+1.0D0) 
!       PHI3(1)=WI(3)/(WJ(5)+WI(3))
!       PHI4(1)=WI(5)/(WK(3)+WI(5))
!       PHI5(1)=WI(4)/(WI(4)+WJ(3)+WK(5))+WI(5)/(WI(5)+WJ(4)+WK(3))+WI(3)/(WI(3)+WJ(5)+WK(4))
!       PHI6(1)=1.0D0/(WJ(2)+1.0D0)
!       PHI7(1)=1.0D0/(WK(1)+1.0D0)
! -----------------------------------   PHI_J
       PHI1(2)=0.0D0
       PHI2(2)=0.0D0
       PHI3(2)=WJ(1)/(WK(3)+WJ(1))
       PHI4(2)=WJ(3)/(WI(1)+WJ(3))
       PHI5(2)=WJ(2)/(WJ(2)+WK(1)+WI(3))+WJ(3)/(WJ(3)+WK(2)+WI(1))+WJ(1)/(WJ(1)+WK(3)+WI(2))
       PHI6(2)=1.0D0
       PHI7(2)=1.0D0
!       PHI1(2)=WJ(1)/(WJ(1)+1.0D0)
!       PHI2(2)=WJ(2)/(WJ(2)+1.0D0)
!       PHI3(2)=WJ(3)/(WK(5)+WJ(3))
!       PHI4(2)=WJ(5)/(WI(3)+WJ(5))  
!       PHI5(2)=WJ(4)/(WJ(4)+WK(3)+WI(5))+WJ(5)/(WJ(5)+WK(4)+WI(3))+WJ(3)/(WJ(3)+WK(5)+WI(4))
!       PHI6(2)=1.0D0/(WK(2)+1.0D0)
!       PHI7(2)=1.0D0/(WI(1)+1.0D0)
! -----------------------------------   197+46=243 FLOP
! -----------------------------------   C1 C2 C3 C4 C5 C6
!       CONS1=PHI2-PHI7
!       CONS2=PHI2+2.0D0*PHI4
!       CONS3=2.0D0*(PHI3-PHI5)
!       CONS4=CONS1+4.0D0*PHI4
!       CONS5=PHI1+CONS2
!       CONS6=PHI6+CONS2
!
!       C1=2.0D0*(CONS3+CONS5)
!       C2=2.0D0*(CONS6+CONS1+CONS3)
!       C3=-CONS5-CONS6-CONS4+4.0D0*PHI5
!       C4=CONS4
!       C5=-1.0D0-2.0D0*CONS1-CONS4
!       C6=2.0D0*(1.0D0+CONS1) 
!
       CONS1=-1.0D0
       CONS2=2.0D0*PHI4
       CONS3=2.0D0*(PHI3-PHI5)
       CONS4=-1.0D0+4.0D0*PHI4
       CONS5=CONS2
       CONS6=1.0D0+CONS2
!
       C1=2.0D0*(CONS3+CONS5)
       C2=2.0D0*(CONS2+CONS3)
       C3=4.0D0*(-CONS2+PHI5)
       C4=CONS4
       C5=2.0D0-4.0D0*PHI4
       C6=0.0D0 
!      ENDIF
! -----------------------------------   243+48 = 291 FLOP
!------------------------------------------  XTMP,YTMP 
       XTMP=AREA_COOR(:,1)*XY1(1)+AREA_COOR(:,2)*XY2(1)+AREA_COOR(:,3)*XY3(1)
       YTMP=AREA_COOR(:,1)*XY1(2)+AREA_COOR(:,2)*XY2(2)+AREA_COOR(:,3)*XY3(2)
!------------------------------------------ 288+NITP_IN*10 FLOP
!------------------------------------------ HIW,HIWC1,HIWC2, HJW, HJWC1,HJWC2
       HIWC1=AREA_COOR(:,2)*C1(1)+AREA_COOR(:,1)*C2(1)+C3(1)
!       HIWC2=C4(1)+AREA_COOR(:,1)*(C5(1)+AREA_COOR(:,1)*C6(1))
       HIWC2=C4(1)+AREA_COOR(:,1)*C5(1)
       HIW  =LILJ*HIWC1+AREA_COOR(:,1)*HIWC2 
!
       HJWC1=AREA_COOR(:,3)*C1(2)+AREA_COOR(:,2)*C2(2)+C3(2)
!       HJWC2=C4(2)+AREA_COOR(:,2)*(C5(2)+AREA_COOR(:,2)*C6(2))
       HJWC2=C4(2)+AREA_COOR(:,2)*C5(2)
       HJW  =LJLK*HJWC1+AREA_COOR(:,2)*HJWC2
!
       HKW  =1.0D0-HIW-HJW
!------------------------------------------ 288+NITP_IN*10+NITP_IN*24=288+NITP_IN*34 FLOP
!------------------------------------------ dLIdx,dLIdy,dLJdx,dLJdy,dLKdx,dLKdy
       dLIdx=XYC3(2)/AREA2 
       dLIdy=-XYC3(1)/AREA2
       dLJdx=-XYC2(2)/AREA2
       dLJdy=XYC2(1)/AREA2
       dLKdx=XYC1(2)/AREA2
       dLKdy=-XYC1(1)/AREA2
!------------------------------------------ 291+NITP_IN*34+6=297+NITP_IN*34 FLOP
!------------------------------------------ dHIWdx,dHIWdy,dHJWdx,dHJWdy,dLILJdx,dLILJdy,dLJLKdx,dLJLKdy
       dLILJdx=dLIdx*AREA_COOR(:,2)+dLJdx*AREA_COOR(:,1)
       dLILJdy=dLIdy*AREA_COOR(:,2)+dLJdy*AREA_COOR(:,1)
!
       dLJLKdx=dLJdx*AREA_COOR(:,3)+dLKdx*AREA_COOR(:,2)
       dLJLKdy=dLJdy*AREA_COOR(:,3)+dLKdy*AREA_COOR(:,2)
!
!       HIWC3=AREA_COOR(:,1)*(C5(1)+AREA_COOR(:,1)*(2.0D0*C6(1)))
!       HJWC3=AREA_COOR(:,2)*(C5(2)+AREA_COOR(:,2)*(2.0D0*C6(2)))
       HIWC3=AREA_COOR(:,1)*C5(1)
       HJWC3=AREA_COOR(:,2)*C5(2)
!
       dHIWdx=dLILJdx*HIWC1+dLIdx*HIWC2+LILJ*(dLJdx*C1(1)+dLIdx*C2(1))+dLIdx*HIWC3
       dHIWdy=dLILJdy*HIWC1+dLIdy*HIWC2+LILJ*(dLJdy*C1(1)+dLIdy*C2(1))+dLIdy*HIWC3        
!
       dHJWdx=dLJLKdx*HJWC1+dLJdx*HJWC2+LJLK*(dLKdx*C1(2)+dLJdx*C2(2))+dLJdx*HJWC3
       dHJWdy=dLJLKdy*HJWC1+dLJdy*HJWC2+LJLK*(dLKdy*C1(2)+dLJdy*C2(2))+dLJdy*HJWC3
!
       dHKWdx=-dHIWdx-dHJWdx
       dHKWdy=-dHIWdy-dHJWdy 
!
!       HIW=AREA_COOR(:,1)
!       HJW=AREA_COOR(:,2)
!       HKW=AREA_COOR(:,3)
!       dHIWdx=dLIdx 
!       dHIWdy=dLIdy
!       dHJWdx=dLJdx
!       dHJWdy=dLJdy
!       dHKWdx=dLKdx
!       dHKWdy=dLKdy
!------------------------------------------ 297+NITP_IN*34+NITP_IN*48=297+NITP_IN*82 FLOP
! **************************** dHImdx,dHImdy
       dHJmdx=0.0D0
       dHJmdy=0.0D0
!
       dHKmdx=0.0D0
       dHKmdy=0.0D0
       IF(INDBCI.EQ.1) THEN
          NXX=ELEPAR(3,Ele1)
          NYY=ELEPAR(4,Ele1)
!    
          CONSTX=(XTMP-X(N1))/RADI
          CONSTY=(YTMP-Y(N1))/RADI
!  
          X_BAR=NXX*CONSTX+NYY*CONSTY
          Y_BAR=-NYY*CONSTX+NXX*CONSTY
!
          SELECT CASE(NDOFI)
          CASE(1)
             dHImdx(:,1)=dHIWdx*AREA_COOR(:,1)+HIW*dLIdx
             dHImdy(:,1)=dHIWdy*AREA_COOR(:,1)+HIW*dLIdy
!
             dHJmdx(:,1)=dHIWdx*AREA_COOR(:,2)+HIW*dLJdx
             dHJmdy(:,1)=dHIWdy*AREA_COOR(:,2)+HIW*dLJdy
!
             dHKmdx(:,1)=dHIWdx*AREA_COOR(:,3)+HIW*dLKdx
             dHKmdy(:,1)=dHIWdy*AREA_COOR(:,3)+HIW*dLKdy
          CASE(2)
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
          CASE(3)
             HIWC1=HIW*X_BAR
!
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
!
             dHImdx(:,3)=dHImdx(:,2)*Y_BAR+HIWC1*(-NYY/RADI)
             dHImdy(:,3)=dHImdy(:,2)*Y_BAR+HIWC1*(NXX/RADI)
          CASE(4)
             HIWC1=HIW*X_BAR
!
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
!
             dHImdx(:,3)=dHImdx(:,2)*Y_BAR+HIWC1*(-NYY/RADI)
             dHImdy(:,3)=dHImdy(:,2)*Y_BAR+HIWC1*(NXX/RADI)
!
             dHImdx(:,4)=dHImdx(:,2)*X_BAR+HIWC1*(NXX/RADI)
             dHImdy(:,4)=dHImdy(:,2)*X_BAR+HIWC1*(NYY/RADI)             
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT         
       ELSE
          X_BAR=(XTMP-X(N1))/RADI
          Y_BAR=(YTMP-Y(N1))/RADI
          SELECT CASE(NDOFI)
          CASE(1)
             dHImdx(:,1)=dHIWdx*AREA_COOR(:,1)+HIW*dLIdx
             dHImdy(:,1)=dHIWdy*AREA_COOR(:,1)+HIW*dLIdy
!
             dHJmdx(:,1)=dHIWdx*AREA_COOR(:,2)+HIW*dLJdx
             dHJmdy(:,1)=dHIWdy*AREA_COOR(:,2)+HIW*dLJdy
!
             dHKmdx(:,1)=dHIWdx*AREA_COOR(:,3)+HIW*dLKdx
             dHKmdy(:,1)=dHIWdy*AREA_COOR(:,3)+HIW*dLKdy
          CASE(3)
             HIWC1=HIW/RADI
!
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1           
          CASE(4)
             HIWC1=HIW/RADI
!
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1
!
             dHImdx(:,4)=dHImdx(:,2)*Y_BAR
             dHImdy(:,4)=dHImdy(:,3)*X_BAR          
          CASE(6)
             HIWC1=HIW/RADI
!
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1          
!
             dHImdx(:,4)=dHImdx(:,2)*Y_BAR
             dHImdy(:,4)=dHImdy(:,3)*X_BAR
!
             dHImdx(:,5)=(dHImdx(:,2)+HIWC1)*X_BAR
             dHImdy(:,5)=dHImdy(:,2)*X_BAR
!
             dHImdx(:,6)=dHImdx(:,3)*Y_BAR
             dHImdy(:,6)=(dHImdy(:,3)+HIWC1)*Y_BAR              
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT 
       ENDIF
! **************************** dHJmdx,dHJmdy
       IF(INDBCJ.EQ.1) THEN
          NXX=ELEPAR(3,Ele2)
          NYY=ELEPAR(4,Ele2)
!
          CONSTX=(XTMP-X(N2))/RADJ
          CONSTY=(YTMP-Y(N2))/RADJ
!
          X_BAR=NXX*CONSTX+NYY*CONSTY
          Y_BAR=-NYY*CONSTX+NXX*CONSTY
!
          SELECT CASE(NDOFJ)
          CASE(1)
             dHImdx(:,1)=dHImdx(:,1)+dHJWdx*AREA_COOR(:,1)+HJW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHJWdy*AREA_COOR(:,1)+HJW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx*AREA_COOR(:,2)+HJW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy*AREA_COOR(:,2)+HJW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHJWdx*AREA_COOR(:,3)+HJW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHJWdy*AREA_COOR(:,3)+HJW*dLKdy
          CASE(2)
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
          CASE(3)
             HJWC1=HJW*X_BAR
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
!
             dHJmdx(:,3)=dHJmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADJ)
             dHJmdy(:,3)=dHJmdy(:,2)*Y_BAR+HJWC1*(NXX/RADJ)
          CASE(4)
             HJWC1=HJW*X_BAR
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
!
             dHJmdx(:,3)=dHJmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADJ)
             dHJmdy(:,3)=dHJmdy(:,2)*Y_BAR+HJWC1*(NXX/RADJ)
!
             dHJmdx(:,4)=dHJmdx(:,2)*X_BAR+HJWC1*(NXX/RADJ)
             dHJmdy(:,4)=dHJmdy(:,2)*X_BAR+HJWC1*(NYY/RADJ)              
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT         
       ELSE
          X_BAR=(XTMP-X(N2))/RADJ
          Y_BAR=(YTMP-Y(N2))/RADJ
          SELECT CASE(NDOFJ)
          CASE(1)
             dHImdx(:,1)=dHImdx(:,1)+dHJWdx*AREA_COOR(:,1)+HJW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHJWdy*AREA_COOR(:,1)+HJW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx*AREA_COOR(:,2)+HJW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy*AREA_COOR(:,2)+HJW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHJWdx*AREA_COOR(:,3)+HJW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHJWdy*AREA_COOR(:,3)+HJW*dLKdy
          CASE(3)
             HJWC1=HJW/RADJ
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1           
          CASE(4)
             HJWC1=HJW/RADJ
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1 
!
             dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
             dHJmdy(:,4)=dHJmdy(:,3)*X_BAR          
          CASE(6)
             HJWC1=HJW/RADJ
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1
!
             dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
             dHJmdy(:,4)=dHJmdy(:,3)*X_BAR
!
             dHJmdx(:,5)=(dHJmdx(:,2)+HJWC1)*X_BAR
             dHJmdy(:,5)=dHJmdy(:,2)*X_BAR
!
             dHJmdx(:,6)=dHJmdx(:,3)*Y_BAR
             dHJmdy(:,6)=(dHJmdy(:,3)+HJWC1)*Y_BAR 
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT 
       ENDIF
! **************************** dHKmdx,dHKmdy
       IF(INDBCK.EQ.1) THEN
          NXX=ELEPAR(3,Ele3)
          NYY=ELEPAR(4,Ele3)
!
          CONSTX=(XTMP-X(N3))/RADK
          CONSTY=(YTMP-Y(N3))/RADK
!
          X_BAR=NXX*CONSTX+NYY*CONSTY
          Y_BAR=-NYY*CONSTX+NXX*CONSTY
!
          SELECT CASE(NDOFK)
          CASE(1)
             dHImdx(:,1)=dHImdx(:,1)+dHKWdx*AREA_COOR(:,1)+HKW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHKWdy*AREA_COOR(:,1)+HKW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHKWdx*AREA_COOR(:,2)+HKW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHKWdy*AREA_COOR(:,2)+HKW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx*AREA_COOR(:,3)+HKW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy*AREA_COOR(:,3)+HKW*dLKdy
          CASE(2)
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
          CASE(3)
             HJWC1=HKW*X_BAR
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
!
             dHKmdx(:,3)=dHKmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADK)
             dHKmdy(:,3)=dHKmdy(:,2)*Y_BAR+HJWC1*(NXX/RADK)
          CASE(4)
             HJWC1=HKW*X_BAR
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
!
             dHKmdx(:,3)=dHKmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADK)
             dHKmdy(:,3)=dHKmdy(:,2)*Y_BAR+HJWC1*(NXX/RADK)
!
             dHKmdx(:,4)=dHKmdx(:,2)*X_BAR+HJWC1*(NXX/RADK)
             dHKmdy(:,4)=dHKmdy(:,2)*X_BAR+HJWC1*(NYY/RADK) 
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT         
       ELSE
          X_BAR=(XTMP-X(N3))/RADK
          Y_BAR=(YTMP-Y(N3))/RADK
          SELECT CASE(NDOFK)
          CASE(1)
             dHImdx(:,1)=dHImdx(:,1)+dHKWdx*AREA_COOR(:,1)+HKW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHKWdy*AREA_COOR(:,1)+HKW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHKWdx*AREA_COOR(:,2)+HKW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHKWdy*AREA_COOR(:,2)+HKW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx*AREA_COOR(:,3)+HKW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy*AREA_COOR(:,3)+HKW*dLKdy
          CASE(3)
             HJWC1=HKW/RADK
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1           
          CASE(4)
             HJWC1=HKW/RADK
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1  
!
             dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
             dHKmdy(:,4)=dHKmdy(:,3)*X_BAR          
          CASE(6)
             HJWC1=HKW/RADK
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1
!
             dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
             dHKmdy(:,4)=dHKmdy(:,3)*X_BAR
!
             dHKmdx(:,5)=(dHKmdx(:,2)+HJWC1)*X_BAR
             dHKmdy(:,5)=dHKmdy(:,2)*X_BAR
!
             dHKmdx(:,6)=dHKmdx(:,3)*Y_BAR
             dHKmdy(:,6)=(dHKmdy(:,3)+HJWC1)*Y_BAR 
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT 
       ENDIF
!------------------------------------------ 297+NITP_IN*82+NITP_IN*39=297+NITP_IN*121 FLOP
       CONST=AREA2*THIC/2.0D0
       KCS1=KCS(1)*CONST
       KCS2=KCS(2)*CONST
       KCS3=KCS(3)*CONST
!------------------------------------------ 297+NITP_IN*121+5=302+NITP_IN*121 FLOP
!------------------------------------------ STIFFII,STIFFJJ,STIFFKK
!      STIFFII
       NCOUNT=0
       DO J=1,NDOFI
          dHImdxW(:,J)=dHImdx(:,J)*GAUSS_WEI(:)
          dHImdyW(:,J)=dHImdy(:,J)*GAUSS_WEI(:)
          ADI(J)=NCOUNT
          DO K=J,NDOFI
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHImdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHImdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHImdy(:,K))  
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,1)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,1)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,1)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,1)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))           
!      STIFFJJ
       NCOUNT=0
       DO J=1,NDOFJ
          dHJmdxW(:,J)=dHJmdx(:,J)*GAUSS_WEI(:)
          dHJmdyW(:,J)=dHJmdy(:,J)*GAUSS_WEI(:)
          ADJ(J)=NCOUNT
          DO K=J,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,2)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,2)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,2)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,2)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!
!      STIFFKK
       NCOUNT=0
       DO J=1,NDOFK
          dHKmdxW(:,J)=dHKmdx(:,J)*GAUSS_WEI(:)
          dHKmdyW(:,J)=dHKmdy(:,J)*GAUSS_WEI(:)
          ADK(J)=NCOUNT
          DO K=J,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHKmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHKmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHKmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHKmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,3)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,3)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,3)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,3)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IJ
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFJ
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHJmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHJmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHJmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,4)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,4)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,4)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,4)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF IK
       NCOUNT=0
       DO J=1,NDOFI
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHImdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHImdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHImdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHImdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,5)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,5)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,5)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,5)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1)) 
!      STIFF JK
       NCOUNT=0
       DO J=1,NDOFJ
          DO K=1,NDOFK
             NCOUNT=NCOUNT+1
             DHIHJ(NCOUNT,1)=SUM(dHJmdxW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,2)=SUM(dHJmdxW(:,J)*dHKmdy(:,K))
             DHIHJ(NCOUNT,3)=SUM(dHJmdyW(:,J)*dHKmdx(:,K))
             DHIHJ(NCOUNT,4)=SUM(dHJmdyW(:,J)*dHKmdy(:,K))           
          ENDDO
       ENDDO
!
       L=NCOUNT
!
       STIFF(1:L,1,6)=(KCS1*DHIHJ(1:L,1)+KCS3*DHIHJ(1:L,4))
       STIFF(1:L,2,6)=(KCS2*DHIHJ(1:L,2)+KCS3*DHIHJ(1:L,3))
       STIFF(1:L,3,6)=(KCS2*DHIHJ(1:L,3)+KCS3*DHIHJ(1:L,2))
       STIFF(1:L,4,6)=(KCS1*DHIHJ(1:L,4)+KCS3*DHIHJ(1:L,1))    
!------------------------------------------ 302+NITP_IN*121+NITP_IN*648-78=224+NITP*769 FLOP
       CALL OFE1_QUADS_ADDBAN(Ele1,Ele2,Ele3,N1,N2,N3,NDOFI,NDOFJ,NDOFK,ADI,ADJ,ADK,STIFF)
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE OFE1_QUADS_STIFF
!
    SUBROUTINE OFE1_QUADS_ADDBAN(Ele1,Ele2,Ele3,N1,N2,N3,NDOFI,NDOFJ,NDOFK,ADI,ADJ,ADK,STIFF)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: Ele1,Ele2,Ele3,N1,N2,N3,NDOFI,NDOFJ,NDOFK
    INTEGER,DIMENSION(NPOR_MAX) :: ADI,ADJ,ADK
    REAL(8),DIMENSION(NPOR_MAX*NPOR_MAX,4,6) :: STIFF
!   IN-FUNCTION VARIABLES
    INTEGER,DIMENSION(2) :: NDOFIXY,NDOFJXY,NDOFKXY
    INTEGER,DIMENSION(2) :: IDI,IDJ,IDK
    INTEGER :: I,J,M,N
    INTEGER*8 :: LL1
    INTEGER :: KI,KJ,KK2,KK1
    LOGICAL :: LOG1,LOG2,LOG3
!
    NDOFIXY(1)=1+NDMFS(2,Ele1)
    NDOFIXY(2)=1+NDMFS(3,Ele1)
    IDI(1)=ID(1,N1)
    IDI(2)=ID(2,N1)
! 
    NDOFJXY(1)=1+NDMFS(2,Ele2)
    NDOFJXY(2)=1+NDMFS(3,Ele2)
    IDJ(1)=ID(1,N2)
    IDJ(2)=ID(2,N2)
!
    NDOFKXY(1)=1+NDMFS(2,Ele3)
    NDOFKXY(2)=1+NDMFS(3,Ele3)
    IDK(1)=ID(1,N3)
    IDK(2)=ID(2,N3)
!*******************************************************************
    LOG1=N1.LT.N2
    LOG2=N1.LT.N3
    LOG3=N2.LT.N3
!*******************************************************************
!   STIFFII
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(1)+N-NDOFIXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(2),NDOFI
       KJ=IDI(2)+M-NDOFIXY(2)
       KK2=ADI(M)
       DO N=M,NDOFI
          KK2=KK2+1
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
!
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,1)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,1)
       ENDDO
    ENDDO
!
    DO M=NDOFIXY(1),NDOFI
       KJ=IDI(1)+M-NDOFIXY(1)
       DO N=NDOFIXY(2),NDOFI
          KI=IDI(2)+N-NDOFIXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADI(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,1)
          ELSE
             KK2=ADI(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,1)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,1)
          ENDIF
       ENDDO
    ENDDO
!   STIFFJJ
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(1)+N-NDOFJXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(2),NDOFJ
       KJ=IDJ(2)+M-NDOFJXY(2)
       KK2=ADJ(M)
       DO N=M,NDOFJ
          KK2=KK2+1
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,2)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,2)
       ENDDO
    ENDDO
!
    DO M=NDOFJXY(1),NDOFJ
       KJ=IDJ(1)+M-NDOFJXY(1)
       DO N=NDOFJXY(2),NDOFJ
          KI=IDJ(2)+N-NDOFJXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADJ(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,2)
          ELSE
             KK2=ADJ(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,2)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,2)
          ENDIF
       ENDDO
    ENDDO
!   STIFFKK
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(1)+N-NDOFKXY(1)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,1,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,1,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(2),NDOFK
       KJ=IDK(2)+M-NDOFKXY(2)
       KK2=ADK(M)
       DO N=M,NDOFK
          KK2=KK2+1
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,4,3)
          KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,4,3)
       ENDDO
    ENDDO
!
    DO M=NDOFKXY(1),NDOFK
       KJ=IDK(1)+M-NDOFKXY(1)
       DO N=NDOFKXY(2),NDOFK
          KI=IDK(2)+N-NDOFKXY(2)         
          LL1=MAXA_COLSOL(KI)+KI-KJ
          IF(N.LE.M) THEN
             KK2=ADK(N)+M-N+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,3,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,3,3)
          ELSE
             KK2=ADK(M)+N-M+1
             K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,2,3)
             KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,2,3)
          ENDIF
       ENDDO
    ENDDO
!*******************************************************************
    KK1=0
    DO I=1,2
       DO J=1,2
          KK1=KK1+1
!   STIFFIJ
          IF(LOG1) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,4)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,4)                
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFJXY(J),NDOFJ
                   KI=IDJ(J)+N-NDOFJXY(J)
                   KK2=(M-1)*NDOFJ+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,4)
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,4)                
                ENDDO
             ENDDO 
          ENDIF
!   STIFFIK
          IF(LOG2) THEN
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,5) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,5)               
                ENDDO
             ENDDO 
          ELSE
             DO M=NDOFIXY(I),NDOFI
                KJ=IDI(I)+M-NDOFIXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,5) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,5)               
                ENDDO
             ENDDO
          ENDIF
!   STIFFJK
          IF(LOG3) THEN
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KI)+KI-KJ
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,6) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,6)               
                ENDDO
             ENDDO
          ELSE
             DO M=NDOFJXY(I),NDOFJ
                KJ=IDJ(I)+M-NDOFJXY(I)
                DO N=NDOFKXY(J),NDOFK
                   KI=IDK(J)+N-NDOFKXY(J)
                   KK2=(M-1)*NDOFK+N
   !
                   LL1=MAXA_COLSOL(KJ)+KJ-KI
   !                                
                   K_COLSOL(LL1)=K_COLSOL(LL1)+STIFF(KK2,KK1,6) 
                   KMATRIX(LL1)=KMATRIX(LL1)+STIFF(KK2,KK1,6)               
                ENDDO
             ENDDO
          ENDIF
!           
       ENDDO
     ENDDO
!*******************************************************************
!
    RETURN
!}}}
    END SUBROUTINE  OFE1_QUADS_ADDBAN
!
    SUBROUTINE OFE1_LINE_G(ADBC,THIC,LI)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC
    REAL(8),INTENT(IN) :: THIC
    REAL(8),INTENT(IN) :: LI(6,NITP_BC)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    LOGICAL :: EXPRN1,INDBCX,INDBCY
    INTEGER :: Ele1,Ele2,N1,N2
    REAL(8) :: LENGTH_HF,LENGTH,NX,NY,NXX,NYY
    REAL(8) :: XM,XL,YM,YL
    REAL(8) :: RADI,RADJ
    REAL(8) :: S1,S2,S3,S4
    REAL(8) :: WI,WJ
    REAL(8) :: PHII
    REAL(8),DIMENSION(2) :: FLOAD
    INTEGER,DIMENSION(2) :: FUNC
    REAL(8),DIMENSION(NITP_BC) :: HIW,HJW,XTMP,YTMP,CONSTX,CONSTY,X_BAR,Y_BAR
    REAL(8),DIMENSION(2) :: LOADD
    REAL(8),DIMENSION(NITP_BC,2) :: LOAD
    REAL(8),DIMENSION(NPOR_MAX,2):: SI,SJ
    REAL(8),DIMENSION(NITP_BC,NPOR_MAX) :: PmI,PmJ
    INTEGER :: INDM,NDOFF
    REAL(8) :: CONST1
    INTEGER :: NUM1,NUM2,JJ1,JJ2,J1
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCED ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCIND ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCPAR ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    DO I=1,NumBc
       ELe1=BCED(1,I)
       ELe2=BCED(2,I)
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
!
       IF(INDBCX) THEN
          NDMFS(4,ELe1)=1
          NDMFS(4,ELe2)=1
       ENDIF
!
       IF(INDBCY) THEN
          NDMFS(5,ELe1)=1
          NDMFS(5,ELe2)=1
       ENDIF
!
       IF(INDBCX.OR.INDBCY) THEN
          N1=NDMFS(1,ELe1)
          N2=NDMFS(1,ELe2)
!
          LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
          LENGTH=SQRT(LENGTH)
!      CALCULATE NORMAL COMPONENTS
          NX=(Y(N2)-Y(N1))/LENGTH
          NY=-(X(N2)-X(N1))/LENGTH
          ELEPAR(3,ELe1)=-NX
          ELEPAR(4,ELe1)=-NY
          ELEPAR(3,ELe2)=-NX
          ELEPAR(4,ELe2)=-NY
       ENDIF
    ENDDO
!
!   CALCULATE LOAD VECTOR
!
    DO I=1,NumBc
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
       IF(INDBCX.AND.INDBCY) CYCLE
!
       ELe1=BCED(1,I)
       ELe2=BCED(2,I)
!
       N1=NDMFS(1,ELe1)
       N2=NDMFS(1,ELe2)
!
       XM=(X(N1)+X(N2))/2.0D0
       XL=(X(N2)-X(N1))/2.0D0
       YM=(Y(N1)+Y(N2))/2.0D0
       YL=(Y(N2)-Y(N1))/2.0D0
!
       LENGTH=XL*XL+YL*YL
       LENGTH_HF=SQRT(LENGTH)
       LENGTH=LENGTH_HF*2.0D0
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
!  ------  WI -------------------------
!       S1=(/LENGTH_HF, LENGTH/)/RADI
!       S2=S1*S1
!       S3=S2*S1
!       S4=S3*S1
!!
!       S2=S2*6.0D0
!       S3=S3*8.0D0
!       S4=S4*3.0D0
!!
!       WI=1.0D0
!       WI(2:3)=WI(2:3)-S2+S3-S4
!!  ------  WJ -------------------------
!       S1=(/LENGTH_HF, LENGTH/)/RADJ
!       S2=S1*S1
!       S3=S2*S1
!       S4=S3*S1
!!
!       S2=S2*6.0D0
!       S3=S3*8.0D0
!       S4=S4*3.0D0
!!
!       WJ=1.0D0
!       WJ(2:3)=WJ(2:3)-S2+S3-S4
!***************************************
!  ------  WI -------------------------
       S1=LENGTH_HF/RADI
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WI=1.0D0-S2+S3-S4
!  ------  WJ -------------------------
       S1=LENGTH_HF/RADJ
       S2=S1*S1
       S3=S2*S1
       S4=S3*S1
!
       S2=S2*6.0D0
       S3=S3*8.0D0
       S4=S4*3.0D0
!
       WJ=1.0D0-S2+S3-S4
!-------- PHII PHIJ --------------------
!       PHII(1)=WI(2)/(WI(2)+WJ(2))
!       PHII(2)=WI(3)/(WI(3)+1.0D0)
!       PHII(2)=0.0D0
!       PHII(3)=1.0D0/(WJ(3)+1.0D0)
!       PHII(3)=1.0D0
       PHII=WI/(WI+WJ)     
!----- XTMP YTMP -----------------------
       XTMP=XM+XL*XABSC_WEIG_BC(:,1)
       YTMP=YM+YL*XABSC_WEIG_BC(:,1)
!----- HIW  HJW ------------------------
       DO J=1,NITP_BC
!          HIW(J)=LI(1,J)+DOT_PRODUCT(LI(2:4,J),PHII)
          HIW(J)=LI(1,J)+LI(2,J)*PHII+LI(4,J)
!-----   LOAD  ---------------------------
          LOADD= F_LOAD(FLOAD,FUNC,NX,NY,XTMP(J),YTMP(J),N1,N2)
          LOAD(J,1)=LOADD(1)*XABSC_WEIG_BC(J,2)
          LOAD(J,2)=LOADD(2)*XABSC_WEIG_BC(J,2)
       ENDDO
!       HIW=LI(5,:)
!       HIW=LI(5,:)
       HJW=(1.0D0-HIW)
!
       CONST1=LENGTH*THIC/2.0D0
! ------ PmI --------------
       NDOFF=NDOF(1,N1)+NDMFS(2,Ele1)
!
       PmI=1.0D0
       INDM=NDMFS(4,Ele1)
       IF(INDM.EQ.1) THEN
          NXX=ELEPAR(3,Ele1)
          NYY=ELEPAR(4,Ele1)
!
          CONSTX=(XTMP-X(N1))/RADI
          CONSTY=(YTMP-Y(N1))/RADI
          X_BAR=(NXX*CONSTX+NYY*CONSTY)
          Y_BAR=(-NYY*CONSTX+NXX*CONSTY)
          SELECT CASE(NDOFF)
          CASE(1)              
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1              
          CASE(2)
              PmI(:,2)=X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
          CASE(3)
              PmI(:,2)=X_BAR
              PmI(:,3)=PmI(:,2)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
          CASE(4)
              PmI(:,2)=X_BAR
              PmI(:,3)=PmI(:,2)*Y_BAR
              PmI(:,4)=PmI(:,2)*X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1              
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ELSE          
          X_BAR=(XTMP-X(N1))/RADI
          Y_BAR=(YTMP-Y(N1))/RADI
          SELECT CASE(NDOFF)
          CASE(1)
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
          CASE(2)
              PmI(:,2)=X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
          CASE(3)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1              
          CASE(4)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
              PmI(:,4)=PmI(:,2)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1
          CASE(6)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
              PmI(:,4)=PmI(:,2)*Y_BAR
              PmI(:,5)=PmI(:,2)*X_BAR
              PmI(:,6)=PmI(:,3)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1
!
              SI(5,1)=SUM(HIW*PmI(:,5)*LOAD(:,1))*CONST1
              SI(5,2)=SUM(HIW*PmI(:,5)*LOAD(:,2))*CONST1
!
              SI(6,1)=SUM(HIW*PmI(:,6)*LOAD(:,1))*CONST1
              SI(6,2)=SUM(HIW*PmI(:,6)*LOAD(:,2))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ENDIF
! ------- PmJ --------------------
       NDOFF=NDOF(1,N2)+NDMFS(2,Ele2)
!
       PmJ=1.0D0
       INDM=NDMFS(4,Ele2)
       IF(INDM.EQ.1) THEN
          NXX=ELEPAR(3,Ele2)
          NYY=ELEPAR(4,Ele2)
!
          CONSTX=(XTMP-X(N2))/RADJ
          CONSTY=(YTMP-Y(N2))/RADJ
          X_BAR=(NXX*CONSTX+NYY*CONSTY)
          Y_BAR=(-NYY*CONSTX+NXX*CONSTY)
          SELECT CASE(NDOFF)
          CASE(1)
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
          CASE(2)
              PmJ(:,2)=X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
          CASE(3)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=PmJ(:,2)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
          CASE(4)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=PmJ(:,2)*Y_BAR
              PmJ(:,4)=PmJ(:,2)*X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ELSE
          X_BAR=(XTMP-X(N2))/RADJ
          Y_BAR=(YTMP-Y(N2))/RADJ
          SELECT CASE(NDOFF)
          CASE(1)
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
          CASE(2)
              PmJ(:,2)=X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
          CASE(3)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
          CASE(4)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
              PmJ(:,4)=PmJ(:,2)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
          CASE(6)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
              PmJ(:,4)=PmJ(:,2)*Y_BAR
              PmJ(:,5)=PmJ(:,2)*X_BAR
              PmJ(:,6)=PmJ(:,3)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
!
              SJ(5,1)=SUM(HJW*LOAD(:,1)*PmJ(:,5))*CONST1
              SJ(5,2)=SUM(HJW*LOAD(:,2)*PmJ(:,5))*CONST1
!
              SJ(6,1)=SUM(HJW*LOAD(:,1)*PmJ(:,6))*CONST1
              SJ(6,2)=SUM(HJW*LOAD(:,2)*PmJ(:,6))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ENDIF
!   
       DO J=1,2
          IF(BCIND(J,I).EQ.1) CYCLE
!
          J1=1+J
!
          NUM1=ID(J,N1)
          NUM2=NUM1+NDOF(J,N1)-1
!
          JJ1=1+NDMFS(J1,Ele1)
          JJ2=NDOF(J,N1)+NDMFS(J1,Ele1)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SI(JJ1:JJ2,J)
!--------------------------------------------------
          NUM1=ID(J,N2)
          NUM2=NUM1+NDOF(J,N2)-1
!
          JJ1=1+NDMFS(J1,Ele2)
          JJ2=NDOF(J,N2)+NDMFS(J1,Ele2)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SJ(JJ1:JJ2,J)
       ENDDO
!
    ENDDO      
!------------------------------------------
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE OFE1_LINE_G    
!
    SUBROUTINE OFE1B_LINE_G(ADBC,THIC,HIW,HJW)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC
    REAL(8),INTENT(IN) :: THIC
    REAL(8),INTENT(IN) :: HIW(NITP_BC),HJW(NITP_BC)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: NumBc
    LOGICAL :: EXPRN1,INDBCX,INDBCY
    INTEGER :: Ele1,Ele2,N1,N2
    REAL(8) :: LENGTH_HF,LENGTH,NX,NY,NXX,NYY
    REAL(8) :: XM,XL,YM,YL
    REAL(8) :: RADI,RADJ
    REAL(8),DIMENSION(2) :: FLOAD
    INTEGER,DIMENSION(2) :: FUNC
    REAL(8),DIMENSION(NITP_BC) :: XTMP,YTMP,CONSTX,CONSTY,X_BAR,Y_BAR
    REAL(8),DIMENSION(2) :: LOADD
    REAL(8),DIMENSION(NITP_BC,2) :: LOAD
    REAL(8),DIMENSION(NPOR_MAX,2):: SI,SJ
    REAL(8),DIMENSION(NITP_BC,NPOR_MAX) :: PmI,PmJ
    INTEGER :: INDM,NDOFF
    REAL(8) :: CONST1
    INTEGER :: NUM1,NUM2,JJ1,JJ2,J1
    INTEGER :: ALLocateStatus
!------------------------------------------
!   READ
!------------------------------------------
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCED ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCIND ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY BCPAR ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    DO I=1,NumBc
       ELe1=BCED(1,I)
       ELe2=BCED(2,I)
!
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) THEN
          NDMFS(4,Ele1)=3
          NDMFS(4,Ele2)=3
          NDMFS(5,Ele1)=3
          NDMFS(5,Ele2)=3
          CYCLE
       ENDIF
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
!
       IF(INDBCX) THEN
          NDMFS(4,ELe1)=1
          NDMFS(4,ELe2)=1
       ELSE
          NDMFS(4,Ele1)=2
          NDMFS(4,Ele2)=2
       ENDIF
!
       IF(INDBCY) THEN
          NDMFS(5,ELe1)=1
          NDMFS(5,ELe2)=1
       Else
          NDMFS(5,Ele1)=2
          NDMFS(5,Ele2)=2
       ENDIF
!
       IF(INDBCX.OR.INDBCY) THEN
          N1=NDMFS(1,ELe1)
          N2=NDMFS(1,ELe2)
!
          LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
          LENGTH=SQRT(LENGTH)
!      CALCULATE NORMAL COMPONENTS
          NX=(Y(N2)-Y(N1))/LENGTH
          NY=-(X(N2)-X(N1))/LENGTH
          ELEPAR(3,ELe1)=-NX
          ELEPAR(4,ELe1)=-NY
          ELEPAR(3,ELe2)=-NX
          ELEPAR(4,ELe2)=-NY
       ENDIF
    ENDDO
!
!   CALCULATE LOAD VECTOR
!
    DO I=1,NumBc
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) CYCLE
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
       IF(INDBCX.AND.INDBCY) CYCLE
!
       ELe1=BCED(1,I)
       ELe2=BCED(2,I)
!
       N1=NDMFS(1,ELe1)
       N2=NDMFS(1,ELe2)
!
       XM=(X(N1)+X(N2))/2.0D0
       XL=(X(N2)-X(N1))/2.0D0
       YM=(Y(N1)+Y(N2))/2.0D0
       YL=(Y(N2)-Y(N1))/2.0D0
!
       LENGTH=XL*XL+YL*YL
       LENGTH_HF=SQRT(LENGTH)
       LENGTH=LENGTH_HF*2.0D0
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       RADI=ELEPAR(1,Ele1)
       RADJ=ELEPAR(1,Ele2)
!----- XTMP YTMP -----------------------
       XTMP=XM+XL*XABSC_WEIG_BC(:,1)
       YTMP=YM+YL*XABSC_WEIG_BC(:,1)
!----- HIW  HJW ------------------------
       DO J=1,NITP_BC
!-----   LOAD  ---------------------------
          LOADD= F_LOAD(FLOAD,FUNC,NX,NY,XTMP(J),YTMP(J),N1,N2)
          LOAD(J,1)=LOADD(1)*XABSC_WEIG_BC(J,2)
          LOAD(J,2)=LOADD(2)*XABSC_WEIG_BC(J,2)
       ENDDO
!
       CONST1=LENGTH*THIC/2.0D0
! ------ PmI --------------
       NDOFF=NDOF(1,N1)+NDMFS(2,Ele1)
!
       PmI=1.0D0
       INDM=NDMFS(4,Ele1)
       IF(INDM.EQ.1) THEN
          NXX=ELEPAR(3,Ele1)
          NYY=ELEPAR(4,Ele1)
!
          CONSTX=(XTMP-X(N1))/RADI
          CONSTY=(YTMP-Y(N1))/RADI
          X_BAR=(NXX*CONSTX+NYY*CONSTY)
          Y_BAR=(-NYY*CONSTX+NXX*CONSTY)
          SELECT CASE(NDOFF)
          CASE(1)              
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
          CASE(2)
              PmI(:,2)=X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
          CASE(3)
              PmI(:,2)=X_BAR
              PmI(:,3)=PmI(:,2)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
          CASE(4)
              PmI(:,2)=X_BAR
              PmI(:,3)=PmI(:,2)*Y_BAR
              PmI(:,4)=PmI(:,2)*X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1              
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ELSE          
          X_BAR=(XTMP-X(N1))/RADI
          Y_BAR=(YTMP-Y(N1))/RADI
          SELECT CASE(NDOFF)
          CASE(1)
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
          CASE(2)
              PmI(:,2)=X_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
          CASE(3)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1              
          CASE(4)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
              PmI(:,4)=PmI(:,2)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1
          CASE(6)
              PmI(:,2)=X_BAR
              PmI(:,3)=Y_BAR
              PmI(:,4)=PmI(:,2)*Y_BAR
              PmI(:,5)=PmI(:,2)*X_BAR
              PmI(:,6)=PmI(:,3)*Y_BAR
!
              SI(1,1)=SUM(HIW*LOAD(:,1))*CONST1
              SI(1,2)=SUM(HIW*LOAD(:,2))*CONST1
!
              SI(2,1)=SUM(HIW*PmI(:,2)*LOAD(:,1))*CONST1
              SI(2,2)=SUM(HIW*PmI(:,2)*LOAD(:,2))*CONST1
!
              SI(3,1)=SUM(HIW*PmI(:,3)*LOAD(:,1))*CONST1
              SI(3,2)=SUM(HIW*PmI(:,3)*LOAD(:,2))*CONST1
!
              SI(4,1)=SUM(HIW*PmI(:,4)*LOAD(:,1))*CONST1
              SI(4,2)=SUM(HIW*PmI(:,4)*LOAD(:,2))*CONST1
!
              SI(5,1)=SUM(HIW*PmI(:,5)*LOAD(:,1))*CONST1
              SI(5,2)=SUM(HIW*PmI(:,5)*LOAD(:,2))*CONST1
!
              SI(6,1)=SUM(HIW*PmI(:,6)*LOAD(:,1))*CONST1
              SI(6,2)=SUM(HIW*PmI(:,6)*LOAD(:,2))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ENDIF
! ------- PmJ --------------------
       NDOFF=NDOF(1,N2)+NDMFS(2,Ele2)
!
       PmJ=1.0D0
       INDM=NDMFS(4,Ele2)
       IF(INDM.EQ.1) THEN
          NXX=ELEPAR(3,Ele2)
          NYY=ELEPAR(4,Ele2)
!
          CONSTX=(XTMP-X(N2))/RADJ
          CONSTY=(YTMP-Y(N2))/RADJ
          X_BAR=(NXX*CONSTX+NYY*CONSTY)
          Y_BAR=(-NYY*CONSTX+NXX*CONSTY)
          SELECT CASE(NDOFF)
          CASE(1)
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
          CASE(2)
              PmJ(:,2)=X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
          CASE(3)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=PmJ(:,2)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
          CASE(4)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=PmJ(:,2)*Y_BAR
              PmJ(:,4)=PmJ(:,2)*X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ELSE
          X_BAR=(XTMP-X(N2))/RADJ
          Y_BAR=(YTMP-Y(N2))/RADJ
          SELECT CASE(NDOFF)
          CASE(1)
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
          CASE(2)
              PmJ(:,2)=X_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
          CASE(3)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
          CASE(4)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
              PmJ(:,4)=PmJ(:,2)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
          CASE(6)
              PmJ(:,2)=X_BAR
              PmJ(:,3)=Y_BAR
              PmJ(:,4)=PmJ(:,2)*Y_BAR
              PmJ(:,5)=PmJ(:,2)*X_BAR
              PmJ(:,6)=PmJ(:,3)*Y_BAR
!
              SJ(1,1)=SUM(LOAD(:,1))*CONST1-SI(1,1)
              SJ(1,2)=SUM(LOAD(:,2))*CONST1-SI(1,2)
!
              SJ(2,1)=SUM(HJW*LOAD(:,1)*PmJ(:,2))*CONST1
              SJ(2,2)=SUM(HJW*LOAD(:,2)*PmJ(:,2))*CONST1
!
              SJ(3,1)=SUM(HJW*LOAD(:,1)*PmJ(:,3))*CONST1
              SJ(3,2)=SUM(HJW*LOAD(:,2)*PmJ(:,3))*CONST1
!
              SJ(4,1)=SUM(HJW*LOAD(:,1)*PmJ(:,4))*CONST1
              SJ(4,2)=SUM(HJW*LOAD(:,2)*PmJ(:,4))*CONST1
!
              SJ(5,1)=SUM(HJW*LOAD(:,1)*PmJ(:,5))*CONST1
              SJ(5,2)=SUM(HJW*LOAD(:,2)*PmJ(:,5))*CONST1
!
              SJ(6,1)=SUM(HJW*LOAD(:,1)*PmJ(:,6))*CONST1
              SJ(6,2)=SUM(HJW*LOAD(:,2)*PmJ(:,6))*CONST1
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GENERATEM/OFE1_LINE_G!" 
             STOP
          END SELECT
       ENDIF
!   
       DO J=1,2
          IF(BCIND(J,I).EQ.1) CYCLE
!
          J1=1+J
!
          NUM1=ID(J,N1)
          NUM2=NUM1+NDOF(J,N1)-1
!
          JJ1=1+NDMFS(J1,Ele1)
          JJ2=NDOF(J,N1)+NDMFS(J1,Ele1)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SI(JJ1:JJ2,J)
!--------------------------------------------------
          NUM1=ID(J,N2)
          NUM2=NUM1+NDOF(J,N2)-1
!
          JJ1=1+NDMFS(J1,Ele2)
          JJ2=NDOF(J,N2)+NDMFS(J1,Ele2)
!
          R(NUM1:NUM2)=R(NUM1:NUM2)+SJ(JJ1:JJ2,J)
       ENDDO
!
    ENDDO      
!------------------------------------------
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE OFE1B_LINE_G  
!
    SUBROUTINE QUADS_OVLP_STIFF(NGDOF,NEleAss,NumEle,NN,THIC) 
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NGDOF,NEleAss,NumEle,NN
    REAL(8) :: THIC
!   IN-FUNCTION VARIALBLES
    INTEGER :: I,J,K,L,M,LL,KK
    INTEGER :: ELE1,ELE2,ELE3
    INTEGER :: N1,N2,N3
    REAL(8) :: AREA2
    INTEGER :: EleMin, EleMax
    REAL(8) :: DHIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: KCS1,KCS2,KCS3,CONST1
    INTEGER :: INDI,INDJ,NI,NJ,NPORI,NPORJ
    REAL(8) :: WEIGHT,XTMP,YTMP
    LOGICAL :: IND1_X,IND1_Y,IND2_X,IND2_Y
    INTEGER :: IDD,INC1,INC2,NUM
    REAL(8) :: RR(2,NPOR_MAX),FBODY_X,FBODY_Y
! -----------------------------------
    ELE1=NDINT(1,NN)
    ELE2=NDINT(2,NN)
    ELE3=NDINT(3,NN)
! -----------------------------------
    N1=NDMFS(1,ELE1)
    N2=NDMFS(1,ELE2)
    N3=NDMFS(1,ELE3)
    AREA2=X(N2)*Y(N3)-X(N3)*Y(N2)+X(N3)*Y(N1)-X(N1)*Y(N3)+X(N1)*Y(N2)-X(N2)*Y(N1)
!   CALCULATE DH_OVLP_QUADS
    CALL GAUTRIANGLE(NGDOF,N1,N2,N3)
!    
    EleMax=0
    EleMin=100000
    CALL DH_QUADS(NEleAss,NN,N1,N2,N3,AREA2,EleMin,EleMax)
!
    KCS1=KCS(1)
    KCS2=KCS(2)
    KCS3=KCS(3)
    CONST1=THIC*AREA2/2.0D0
    DO I=EleMin,EleMax
       INDI=NDMFS(2,I)
       NI=NDMFS(1,I)
       IND1_X=INDI.EQ.1
       IND1_Y=NDMFS(3,I).EQ.1
       IF(IND1_X) THEN
          NPORI=NDOF(1,NI)+1
       ELSE
          NPORI=NDOF(1,NI)
       ENDIF
!
       DO K=1,NPORI
          RR(1,K)=0.0D0
          RR(2,K)=0.0D0
       ENDDO
!   GENERATE BODY FORCE
       DO M=1,NITP_IN
          XTMP=XABSC_WEIG_IN(M,1)
          YTMP=XABSC_WEIG_IN(M,2)
          WEIGHT=XABSC_WEIG_IN(M,3)
!
          FBODY_X=0.0D0
          FBODY_Y=0.0D0
          CALL FBODY_OVLP(XTMP,YTMP,FBODY_X,FBODY_Y)
!
          DO K=1,NPORI
             RR(1,K)=RR(1,K)+H_OVLP_QUADS(K,M,I)*FBODY_X*WEIGHT
             RR(2,K)=RR(2,K)+H_OVLP_QUADS(K,M,I)*FBODY_Y*WEIGHT
          ENDDO
       ENDDO
!   CALCULATE LOAD VECTOR
       INC1=0
       IF(IND1_X) INC1=1
       INC2=0
       IF(IND1_Y) INC2=1
!
       IDD=ID(1,NI)
       IF(IDD.NE.0) THEN
          DO J=1,NDOF(1,NI)
             NUM=IDD+J-1
             KK=J+INC1
             R(NUM)=R(NUM)+RR(1,KK)*CONST1
          ENDDO
       ENDIF
!
       IDD=ID(2,NI)
       IF(IDD.NE.0) THEN
          DO J=1,NDOF(2,NI)
             NUM=IDD+J-1
             KK=J+INC2
             R(NUM)=R(NUM)+RR(2,KK)*CONST1
          ENDDO
       ENDIF
!
       DO J=I,EleMax
          INDJ=NDMFS(2,J)
          NJ=NDMFS(1,J)
          IND2_X=INDJ.EQ.1
          IND2_Y=NDMFS(3,J).EQ.1
          IF(IND2_X) THEN
             NPORJ=NDOF(1,NJ)+1
          ELSE
             NPORJ=NDOF(1,NJ)
          ENDIF
!
          KK=NPORI*NPORJ
!         INITIALIZATION
          DHIHJ=0.0D0
          STIFF_IJ=0.0D0
!
          DO M=1,NITP_IN
             WEIGHT=XABSC_WEIG_IN(M,3)
             DO K=1,NPORI
                DO L=1,NPORJ
                 LL=(K-1)*NPORJ+L
!
                 DHIHJ(LL,1)=DHIHJ(LL,1)+DH_OVLP_QUADS(K,M,I)*DH_OVLP_QUADS(L,M,J)*WEIGHT
                 DHIHJ(LL,2)=DHIHJ(LL,2)+DH_OVLP_QUADS(K,M,I)*DH_OVLP_QUADS(L+NPOR_MAX,M,J)*WEIGHT
                 DHIHJ(LL,3)=DHIHJ(LL,3)+DH_OVLP_QUADS(K+NPOR_MAX,M,I)*DH_OVLP_QUADS(L,M,J)*WEIGHT
                 DHIHJ(LL,4)=DHIHJ(LL,4)+DH_OVLP_QUADS(K+NPOR_MAX,M,I)*DH_OVLP_QUADS(L+NPOR_MAX,M,J)*WEIGHT
                ENDDO
             ENDDO
          ENDDO
!
          DO K=1,KK
             STIFF_IJ(K,1)=(KCS1*DHIHJ(K,1)+KCS3*DHIHJ(K,4))*CONST1
             STIFF_IJ(K,2)=(KCS2*DHIHJ(K,2)+KCS3*DHIHJ(K,3))*CONST1
             STIFF_IJ(K,3)=(KCS2*DHIHJ(K,3)+KCS3*DHIHJ(K,2))*CONST1
             STIFF_IJ(K,4)=(KCS1*DHIHJ(K,4)+KCS3*DHIHJ(K,1))*CONST1
          ENDDO
!
          CALL QUADS_MFS_ADDBAN_V2(I,J,NPORI,NPORJ,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))         
       ENDDO
    ENDDO
!   reset data
    DH_OVLP_QUADS(:,:,EleMIn:EleMax)=0.0D0
    H_OVLP_QUADS(:,:,EleMIn:EleMax)=0.0D0
!
    RETURN
!}}}
    END SUBROUTINE QUADS_OVLP_STIFF
!
    SUBROUTINE FBODY_OVLP(X,Y,FBODY_X,FBODY_Y)
!{{{
    IMPLICIT NONE
    REAL(8),INTENT(IN) :: X,Y
    REAL(8),INTENT(INOUT) :: FBODY_X,FBODY_Y
!   ASSISTANT VARIABLES
    REAL(8):: DUDXX,DUDXY,DUDYY,DVDXX,DVDXY,DVDYY
    REAL(8) :: EXPY,COSX,SINX,X2,Y2,X2_1,Y2_1,X22,Y22
    REAL(8) :: KCS1,KCS2,KCS3
!
    FBODY_X=0.0D0
    FBODY_Y=0.0D0
!
    KCS1=KCS(1)
    KCS2=KCS(2)
    KCS3=KCS(3)
!
    EXPY=EXP(5*Y)
    COSX=COS(5*X)
    SINX=SIN(5*X)
    X2=X*X
    Y2=Y*Y
!
    X2_1=1.0D0-X2
    Y2_1=1.0D0-Y2
!
    X22=X2_1*X2_1
    Y22=Y2_1*Y2_1
!
    DUDXX=8*X2*Y22*EXPY*COSX+(40*X2_1)*Y22*EXPY*SINX*X-(4*X2_1)*Y22*EXPY*COSX-25*X22*Y22*EXPY*COSX
    DUDXY=(16*X2_1)*Y2_1*EXPY*COSX*X*Y-(20*X2_1)*Y22*EXPY*COSX*X+20*X22*Y2_1*EXPY*SINX*Y-25*X22*Y22*EXPY*SINX
    DUDYY=8*X22*Y2*EXPY*COSX-40*X22*Y2_1*EXPY*COSX*Y-4*X22*Y2_1*EXPY*COSX+25*X22*Y22*EXPY*COSX
!
    DVDXX=8*X2*Y22*EXPY*SINX-(40*X2_1)*Y22*EXPY*COSX*X-(4*X2_1)*Y22*EXPY*SINX-25*X22*Y22*EXPY*SINX
    DVDXY=(16*X2_1)*Y2_1*EXPY*SINX*Y*X-20*X22*Y2_1*EXPY*COSX*Y-(20*X2_1)*Y22*EXPY*SINX*X+25*X22*Y22*EXPY*COSX
    DVDYY=8*X22*Y2*EXPY*SINX-40*X22*Y2_1*EXPY*SINX*Y-4*X22*Y2_1*EXPY*SINX+25*X22*Y22*EXPY*SINX
!
!    TAUXX=KCS1*DUDX+KCS2*DVDY
!    TAUYY=KCS1*DVDY+KCS2*DUDX
!    TAUXY=KCS3*(DUDY+DVDX)
!
    FBODY_X=-(KCS1*DUDXX+KCS2*DVDXY)-(KCS3*(DUDYY+DVDXY))
    FBODY_Y=-(KCS1*DVDYY+KCS2*DUDXY)-(KCS3*(DUDXY+DVDXX))
!
    RETURN
!}}}
    END SUBROUTINE FBODY_OVLP    
!
    SUBROUTINE DH_QUADS(NEleAss,NN,N1,N2,N3,AREA2,EleMin,EleMax)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEleAss,NN,N1,N2,N3
    REAL(8),INTENT(IN) :: AREA2
    INTEGER,INTENT(INOUT) :: EleMin,EleMax
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L,LL,L2
    REAL(8) :: X1,X2,X3,Y1,Y2,Y3
    REAL(8) :: XTMP,YTMP,AXYZ(3),DAXYZ_X(3),DAXYZ_Y(3)
    INTEGER :: ELE,N,IND,INDN,NDOFF
    REAL(8) :: RAD,CONSTX,CONSTY,NX,NY
    REAL(8) :: X_BAR,Y_BAR
    REAL(8) :: HH(NPOR_MAX),DHH(NPOR_MAX,2)
    REAL(8) :: WEI(NEleAss),DWEI(NEleAss,2)
    INTEGER :: ELE1,ELE2,M,INDSP
!
    X1=X(N1)
    Y1=Y(N1)
    X2=X(N2)
    Y2=Y(N2)
    X3=X(N3)
    Y3=Y(N3)
!
    DO I=1,3
       SELECT CASE(I)
       CASE(1)
           DAXYZ_X(I)=(Y2-Y3)/AREA2
           DAXYZ_Y(I)=(X3-X2)/AREA2
       CASE(2)
           DAXYZ_X(I)=(Y3-Y1)/AREA2
           DAXYZ_Y(I)=(X1-X3)/AREA2
       CASE(3)
           DAXYZ_X(I)=(Y1-Y2)/AREA2
           DAXYZ_Y(I)=(X2-X1)/AREA2
       CASE DEFAULT
           WRITE(*,*) "ERROR, /GENERATEM/OVLP_LINE_G! "
           STOP
       END SELECT
    ENDDO
!
    DO K=1,NITP_IN
       XTMP=XABSC_WEIG_IN(K,1)
       YTMP=XABSC_WEIG_IN(K,2)
!
       DO I=1,3
          AXYZ(I)=ACOOR_IN(K,I)
       ENDDO
!
       DO I=1,3
          ELE=NDINT(I,NN)
          N=NDMFS(1,ELE)
          IND=NDMFS(2,ELE)
          INDSP=NDMFS(4,ELE)
          IF ((IND.NE.0).OR.(INDSP.EQ.2)) THEN
             IF(ELE.LT.EleMin) EleMin=ELE
             IF(ELE.GT.EleMax) EleMax=ELE
             RAD=2.0D0*ELEPAR(1,ELE)
             CONSTX=(XTMP-X(N))/RAD
             CONSTY=(YTMP-Y(N))/RAD
!
             IF(IND.EQ.1) THEN
                NDOFF=NDOF(1,N)+1
             ELSE
                NDOFF=NDOF(1,N)
             ENDIF
!
             INDN=NDMFS(5,ELE)
             IF(INDN.EQ.1) THEN
                NX=ELEPAR(3,ELE)
                NY=ELEPAR(4,ELE)
                X_BAR=(NX*CONSTX+NY*CONSTY)
                Y_BAR=(-NY*CONSTX+NX*CONSTY)
!{{{
                SELECT CASE(NDOFF)
                CASE(1)
                 HH(1)=AXYZ(I)
                 DHH(1,1)=DAXYZ_X(I)
                 DHH(1,2)=DAXYZ_Y(I)
                CASE(2)
                 HH(1)=AXYZ(I)
                 DHH(1,1)=DAXYZ_X(I)
                 DHH(1,2)=DAXYZ_Y(I)
! HIX,DHIXDX,DHIXDY
                 HH(2)=HH(1)*X_BAR
                 DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                 DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
                CASE(3)
                 HH(1)=AXYZ(I)
                 DHH(1,1)=DAXYZ_X(I)
                 DHH(1,2)=DAXYZ_Y(I)
! HIX,DHIXDX,DHIXDY
                 HH(2)=HH(1)*X_BAR
                 DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                 DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
! HIXY,DHIXYDX,DHIXYDY
                 HH(3)=HH(2)*Y_BAR
                 DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                 DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
               CASE(4)
                 HH(1)=AXYZ(I)
                 DHH(1,1)=DAXYZ_X(I)
                 DHH(1,2)=DAXYZ_Y(I)
! HIX,DHIXDX,DHIXDY
                 HH(2)=HH(1)*X_BAR
                 DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                 DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
! HIXY,DHIXYDX,DHIXYDY
                 HH(3)=HH(2)*Y_BAR
                 DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                 DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
! HIXX,DHIXXDX,DHIXXDY
                 HH(4)=HH(2)*X_BAR
                 DHH(4,1)=DHH(2,1)*X_BAR+HH(2)*NX/RAD
                 DHH(4,2)=DHH(2,2)*X_BAR+HH(2)*NY/RAD
               CASE(6)
                 HH(1)=AXYZ(I)
                 DHH(1,1)=DAXYZ_X(I)
                 DHH(1,2)=DAXYZ_Y(I)
! HIX,DHIXDX,DHIXDY
                 HH(2)=HH(1)*X_BAR
                 DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                 DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
! HIXY,DHIXYDX,DHIXYDY
                 HH(3)=HH(2)*Y_BAR
                 DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                 DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
! HIXX,DHIXXDX,DHIXXDY
                 HH(4)=HH(2)*X_BAR
                 DHH(4,1)=DHH(2,1)*X_BAR+HH(2)*NX/RAD
                 DHH(4,2)=DHH(2,2)*X_BAR+HH(2)*NY/RAD
! HIXXY,DHIXXYDX,DHIXXYDY
                 HH(5)=HH(3)*X_BAR
                 DHH(5,1)=DHH(3,1)*X_BAR+HH(3)*NX/RAD
                 DHH(5,2)=DHH(3,2)*X_BAR+HH(3)*NY/RAD
! HIXYY,DHIXYYDX,DHIXYYDY
                 HH(6)=HH(3)*Y_BAR
                 DHH(6,1)=DHH(3,1)*Y_BAR+HH(3)*(-NY)/RAD
                 DHH(6,2)=DHH(3,2)*Y_BAR+HH(3)*NX/RAD
               CASE DEFAULT
                  WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                  STOP
               END SELECT
!
               DO L=1,2
                  DO J=1,NDOFF
                     LL=J+(L-1)*NPOR_MAX
                     DH_OVLP_QUADS(LL,K,ELE)=DH_OVLP_QUADS(LL,K,ELE)+DHH(J,L)
                  ENDDO
               ENDDO
!
               DO L=1,NDOFF
                  H_OVLP_QUADS(L,K,ELE)=H_OVLP_QUADS(L,K,ELE)+HH(L)
               ENDDO
!}}}
             ELSE
                IF(IND.EQ.3) THEN
                   DO J=1,3
                      ELE1=NDINT(J,NN)
!
                      HH(1)=AXYZ(I)*AXYZ(J)
                      DHH(1,1)=DAXYZ_X(I)*AXYZ(J)+AXYZ(I)*DAXYZ_X(J)
                      DHH(1,2)=DAXYZ_Y(I)*AXYZ(J)+AXYZ(I)*DAXYZ_Y(J)
!
                      DO L=1,2
                         LL=1+(L-1)*NPOR_MAX
                         DH_OVLP_QUADS(LL,K,ELE1)=DH_OVLP_QUADS(LL,K,ELE1)+DHH(1,L)
                      ENDDO
!
                      H_OVLP_QUADS(1,K,ELE1)=H_OVLP_QUADS(1,K,ELE1)+HH(1)
                   ENDDO
                ELSE
!{{{
                  SELECT CASE(NDOFF)
                  CASE(1)
                   HH(1)=AXYZ(I)
                   DHH(1,1)=DAXYZ_X(I)
                   DHH(1,2)=DAXYZ_Y(I)
                  CASE(3)
                   HH(1)=AXYZ(I)
                   DHH(1,1)=DAXYZ_X(I)
                   DHH(1,2)=DAXYZ_Y(I)
  ! HIX,DHIXDX,DHIXDY
                   HH(2)=HH(1)*CONSTX
                   DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                   DHH(2,2)=DHH(1,2)*CONSTX
  ! HIY,DHIYDX,DHIYDY
                   HH(3)=HH(1)*CONSTY
                   DHH(3,1)=DHH(1,1)*CONSTY
                   DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD             
                  CASE(4)
                   HH(1)=AXYZ(I)
                   DHH(1,1)=DAXYZ_X(I)
                   DHH(1,2)=DAXYZ_Y(I)
  ! HIX,DHIXDX,DHIXDY
                   HH(2)=HH(1)*CONSTX
                   DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                   DHH(2,2)=DHH(1,2)*CONSTX
  ! HIY,DHIYDX,DHIYDY
                   HH(3)=HH(1)*CONSTY
                   DHH(3,1)=DHH(1,1)*CONSTY
                   DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
  ! HIXY,DHIXYDX,DHIXYDY
                   HH(4)=HH(2)*CONSTY
                   DHH(4,1)=DHH(2,1)*CONSTY
                   DHH(4,2)=DHH(3,2)*CONSTX                
                 CASE(6)
                   HH(1)=AXYZ(I)
                   DHH(1,1)=DAXYZ_X(I)
                   DHH(1,2)=DAXYZ_Y(I)
  ! HIX,DHIXDX,DHIXDY
                   HH(2)=HH(1)*CONSTX
                   DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                   DHH(2,2)=DHH(1,2)*CONSTX
  ! HIY,DHIYDX,DHIYDY
                   HH(3)=HH(1)*CONSTY
                   DHH(3,1)=DHH(1,1)*CONSTY
                   DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
  ! HIXY,DHIXYDX,DHIXYDY
                   HH(4)=HH(2)*CONSTY
                   DHH(4,1)=DHH(2,1)*CONSTY
                   DHH(4,2)=DHH(3,2)*CONSTX 
  !   HIXX,DHIXXDX,DHIXXDY
                   HH(5)=HH(2)*CONSTX
                   DHH(5,1)=DHH(2,1)*CONSTX+HH(2)/RAD
                   DHH(5,2)=DHH(2,2)*CONSTX
  !   HIYY,DHIYYDX,DHIYYDY
                   HH(6)=HH(3)*CONSTY
                   DHH(6,1)=DHH(3,1)*CONSTY
                   DHH(6,2)=DHH(3,2)*CONSTY+HH(3)/RAD 
                 CASE(8)
                   HH(1)=AXYZ(I)
                   DHH(1,1)=DAXYZ_X(I)
                   DHH(1,2)=DAXYZ_Y(I)
  ! HIX,DHIXDX,DHIXDY
                   HH(2)=HH(1)*CONSTX
                   DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                   DHH(2,2)=DHH(1,2)*CONSTX
  ! HIY,DHIYDX,DHIYDY
                   HH(3)=HH(1)*CONSTY
                   DHH(3,1)=DHH(1,1)*CONSTY
                   DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
  ! HIXY,DHIXYDX,DHIXYDY
                   HH(4)=HH(2)*CONSTY
                   DHH(4,1)=DHH(2,1)*CONSTY
                   DHH(4,2)=DHH(3,2)*CONSTX 
  !   HIXX,DHIXXDX,DHIXXDY
                   HH(5)=HH(2)*CONSTX
                   DHH(5,1)=DHH(2,1)*CONSTX+HH(2)/RAD
                   DHH(5,2)=DHH(2,2)*CONSTX
  !   HIYY,DHIYYDX,DHIYYDY
                   HH(6)=HH(3)*CONSTY
                   DHH(6,1)=DHH(3,1)*CONSTY
                   DHH(6,2)=DHH(3,2)*CONSTY+HH(3)/RAD
  !   HIXXY,DHIXXYDX,DHIXXYDY
                   HH(7)=HH(4)*CONSTX
                   DHH(7,1)=DHH(4,1)*CONSTX+HH(4)/RAD
                   DHH(7,2)=DHH(4,2)*CONSTX
  !   HIXYY,DHIXYYDX,DHIXYYDY
                   HH(8)=HH(4)*CONSTY
                   DHH(8,1)=DHH(4,1)*CONSTY
                   DHH(8,2)=DHH(4,2)*CONSTY+HH(4)/RAD
                  CASE DEFAULT
                    WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                    STOP
                  END SELECT
!
                  DO L=1,2
                     DO J=1,NDOFF
                        LL=J+(L-1)*NPOR_MAX
                        DH_OVLP_QUADS(LL,K,ELE)=DH_OVLP_QUADS(LL,K,ELE)+DHH(J,L)
                     ENDDO
                  ENDDO
!
                  DO L=1,NDOFF
                     H_OVLP_QUADS(L,K,ELE)=H_OVLP_QUADS(L,K,ELE)+HH(L)
                  ENDDO
!}}}
                ENDIF
             ENDIF
          ELSE
             CALL WEI_CAL_QUADS(ELE,NEleAss,XTMP,YTMP,WEI(1),DWEI(1,1))
             DO L=1,NEleAss
                ELE1=LM(L,ELE)
                IF(ELE1.NE.0) THEN
!
                   IF(ELE1.LT.EleMin) EleMin=ELE1
                   IF(ELE1.GT.EleMax) EleMax=ELE1
!
                   M=NDMFS(1,ELE1)
                   IND=NDMFS(2,ELE1)
                   RAD=2*ELEPAR(1,ELE1)
!
                   CONSTX=(XTMP-X(M))/RAD
                   CONSTY=(YTMP-Y(M))/RAD 
!
                   IF(IND.EQ.1) THEN
                      NDOFF=NDOF(1,M)+1
                   ELSE
                      NDOFF=NDOF(1,M)
                   ENDIF
!
                   INDN=NDMFS(5,ELE1)
                   IF (INDN.EQ.1) THEN
                       NX=ELEPAR(3,ELE1)
                       NY=ELEPAR(4,ELE1)
                       X_BAR=(NX*CONSTX+NY*CONSTY)
                       Y_BAR=(-NY*CONSTX+NX*CONSTY)
!{{{
                       SELECT CASE(NDOFF)
                       CASE(1)
                          HH(1)=AXYZ(I)*WEI(L)
                          DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                          DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
                       CASE(2)
                          HH(1)=AXYZ(I)*WEI(L)
                          DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                          DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
       ! HIX,DHIXDX,DHIXDY
                          HH(2)=HH(1)*X_BAR
                          DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                          DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
                       CASE(3)
                          HH(1)=AXYZ(I)*WEI(L)
                          DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                          DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
       ! HIX,DHIXDX,DHIXDY
                          HH(2)=HH(1)*X_BAR
                          DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                          DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
       ! HIXY,DHIXYDX,DHIXYDY
                          HH(3)=HH(2)*Y_BAR
                          DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                          DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
                      CASE(4)
                          HH(1)=AXYZ(I)*WEI(L)
                          DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                          DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
       ! HIX,DHIXDX,DHIXDY
                          HH(2)=HH(1)*X_BAR
                          DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                          DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
       ! HIXY,DHIXYDX,DHIXYDY
                          HH(3)=HH(2)*Y_BAR
                          DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                          DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
       ! HIXX,DHIXXDX,DHIXXDY
                          HH(4)=HH(2)*X_BAR
                          DHH(4,1)=DHH(2,1)*X_BAR+HH(2)*NX/RAD
                          DHH(4,2)=DHH(2,2)*X_BAR+HH(2)*NY/RAD
                      CASE(6)
                          HH(1)=AXYZ(I)*WEI(L)
                          DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                          DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
         ! HIX,DHIXDX,DHIXDY
                          HH(2)=HH(1)*X_BAR
                          DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                          DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
         ! HIXY,DHIXYDX,DHIXYDY
                          HH(3)=HH(2)*Y_BAR
                          DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                          DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
         ! HIXX,DHIXXDX,DHIXXDY
                          HH(4)=HH(2)*X_BAR
                          DHH(4,1)=DHH(2,1)*X_BAR+HH(2)*NX/RAD
                          DHH(4,2)=DHH(2,2)*X_BAR+HH(2)*NY/RAD
         ! HIXXY,DHIXXYDX,DHIXXYDY
                          HH(5)=HH(3)*X_BAR
                          DHH(5,1)=DHH(3,1)*X_BAR+HH(3)*NX/RAD
                          DHH(5,2)=DHH(3,2)*X_BAR+HH(3)*NY/RAD
         ! HIXYY,DHIXYYDX,DHIXYYDY
                          HH(6)=HH(3)*Y_BAR
                          DHH(6,1)=DHH(3,1)*Y_BAR+HH(3)*(-NY)/RAD
                          DHH(6,2)=DHH(3,2)*Y_BAR+HH(3)*NX/RAD
                       CASE DEFAULT
                         WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                         STOP
                       END SELECT
!}}}  
!
                       DO L2=1,2
                          DO J=1,NDOFF
                             LL=J+(L2-1)*NPOR_MAX
                             DH_OVLP_QUADS(LL,K,ELE1)=DH_OVLP_QUADS(LL,K,ELE1)+DHH(J,L2)
                          ENDDO
                       ENDDO
!
                       DO L2=1,NDOFF
                          H_OVLP_QUADS(L2,K,ELE1)=H_OVLP_QUADS(L2,K,ELE1)+HH(L2)
                       ENDDO
                   ELSE
!{{{
                      IF(IND.EQ.3) THEN 
                         DO J=1,3
                            ELE2=NDINT(J,NN)
   !
                            HH(1)=AXYZ(I)*AXYZ(J)*WEI(L)
                            DHH(1,1)=(DAXYZ_X(I)*AXYZ(J)+AXYZ(I)*DAXYZ_X(J))*WEI(L)+AXYZ(I)*AXYZ(J)*DWEI(L,1)
                            DHH(1,2)=(DAXYZ_Y(I)*AXYZ(J)+AXYZ(I)*DAXYZ_Y(J))*WEI(L)+AXYZ(I)*AXYZ(J)*DWEI(L,2)
   !
                            DO L2=1,2
                               LL=1+(L2-1)*NPOR_MAX
                               DH_OVLP_QUADS(LL,K,ELE2)=DH_OVLP_QUADS(LL,K,ELE2)+DHH(1,L2)
                            ENDDO
!
                            H_OVLP_QUADS(1,K,ELE2)=H_OVLP_QUADS(1,K,ELE2)+HH(1)
                         ENDDO
                      ELSE
!
                         SELECT CASE(NDOFF)
                         CASE(1)
                            HH(1)=AXYZ(I)*WEI(L)
                            DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                            DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
                         CASE(3)
                            HH(1)=AXYZ(I)*WEI(L)
                            DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                            DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
      ! HIX,DHIXDX,DHIXDY
                            HH(2)=HH(1)*CONSTX
                            DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                            DHH(2,2)=DHH(1,2)*CONSTX
      ! HIY,DHIYDX,DHIYDY
                            HH(3)=HH(1)*CONSTY
                            DHH(3,1)=DHH(1,1)*CONSTY
                            DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD  
                         CASE(4)
                            HH(1)=AXYZ(I)*WEI(L)
                            DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                            DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
      ! HIX,DHIXDX,DHIXDY
                            HH(2)=HH(1)*CONSTX
                            DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                            DHH(2,2)=DHH(1,2)*CONSTX
      ! HIY,DHIYDX,DHIYDY
                            HH(3)=HH(1)*CONSTY
                            DHH(3,1)=DHH(1,1)*CONSTY
                            DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
      ! HIXY,DHIXYDX,DHIXYDY
                            HH(4)=HH(2)*CONSTY
                            DHH(4,1)=DHH(2,1)*CONSTY
                            DHH(4,2)=DHH(3,2)*CONSTX                
                        CASE(6)
                            HH(1)=AXYZ(I)*WEI(L)
                            DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                            DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
      ! HIX,DHIXDX,DHIXDY
                            HH(2)=HH(1)*CONSTX
                            DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                            DHH(2,2)=DHH(1,2)*CONSTX
      ! HIY,DHIYDX,DHIYDY
                            HH(3)=HH(1)*CONSTY
                            DHH(3,1)=DHH(1,1)*CONSTY
                            DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
      ! HIXY,DHIXYDX,DHIXYDY
                            HH(4)=HH(2)*CONSTY
                            DHH(4,1)=DHH(2,1)*CONSTY
                            DHH(4,2)=DHH(3,2)*CONSTX 
      ! HIXX,DHIXXDX,DHIXXDY
                            HH(5)=HH(2)*CONSTX
                            DHH(5,1)=DHH(2,1)*CONSTX+HH(2)/RAD
                            DHH(5,2)=DHH(2,2)*CONSTX
      ! HIYY,DHIYYDX,DHIYYDY
                            HH(6)=HH(3)*CONSTY
                            DHH(6,1)=DHH(3,1)*CONSTY
                            DHH(6,2)=DHH(3,2)*CONSTY+HH(3)/RAD 
                         CASE(8)
                            HH(1)=AXYZ(I)*WEI(L)
                            DHH(1,1)=DAXYZ_X(I)*WEI(L)+AXYZ(I)*DWEI(L,1)
                            DHH(1,2)=DAXYZ_Y(I)*WEI(L)+AXYZ(I)*DWEI(L,2)
          ! HIX,DHIXDX,DHIXDY
                           HH(2)=HH(1)*CONSTX
                           DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                           DHH(2,2)=DHH(1,2)*CONSTX
          ! HIY,DHIYDX,DHIYDY
                           HH(3)=HH(1)*CONSTY
                           DHH(3,1)=DHH(1,1)*CONSTY
                           DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD 
          ! HIXY,DHIXYDX,DHIXYDY
                           HH(4)=HH(2)*CONSTY
                           DHH(4,1)=DHH(2,1)*CONSTY
                           DHH(4,2)=DHH(3,2)*CONSTX 
          !   HIXX,DHIXXDX,DHIXXDY
                           HH(5)=HH(2)*CONSTX
                           DHH(5,1)=DHH(2,1)*CONSTX+HH(2)/RAD
                           DHH(5,2)=DHH(2,2)*CONSTX
          !   HIYY,DHIYYDX,DHIYYDY
                           HH(6)=HH(3)*CONSTY
                           DHH(6,1)=DHH(3,1)*CONSTY
                           DHH(6,2)=DHH(3,2)*CONSTY+HH(3)/RAD
          !   HIXXY,DHIXXYDX,DHIXXYDY
                           HH(7)=HH(4)*CONSTX
                           DHH(7,1)=DHH(4,1)*CONSTX+HH(4)/RAD
                           DHH(7,2)=DHH(4,2)*CONSTX
          !   HIXYY,DHIXYYDX,DHIXYYDY
                           HH(8)=HH(4)*CONSTY
                           DHH(8,1)=DHH(4,1)*CONSTY
                           DHH(8,2)=DHH(4,2)*CONSTY+HH(4)/RAD
                         CASE DEFAULT
                             WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                             STOP
                         END SELECT
!
                         DO L2=1,2
                            DO J=1,NDOFF
                               LL=J+(L2-1)*NPOR_MAX
                               DH_OVLP_QUADS(LL,K,ELE1)=DH_OVLP_QUADS(LL,K,ELE1)+DHH(J,L2)
                            ENDDO
                         ENDDO
!
                         DO L2=1,NDOFF
                            H_OVLP_QUADS(L2,K,ELE1)=H_OVLP_QUADS(L2,K,ELE1)+HH(L2)
                         ENDDO
                      ENDIF
!}}}
                   ENDIF
!
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE DH_QUADS 
!
    SUBROUTINE WEI_CAL_QUADS(ELE,NEleAss,XTMP,YTMP,WEI,DWEI)
!{{{
    IMPLICIT NONE
    INTEGER :: ELE,NEleAss
    REAL(8) :: XTMP,YTMP
    REAL(8) :: WEI(NEleAss),DWEI(NEleAss,2)
!   ASSISTANT VARIABLES
    INTEGER :: I
    INTEGER :: ELE1,N
    REAL(8) :: RAD
    REAL(8) :: S,S1,SS(4)
    REAL(8) :: WEI_SUM,DWEI_SUM(2)
    REAL(8) :: AAA,BBB
!   INITIALIZATION
    WEI=0.0D0
    DWEI=0.D0
!    
    WEI_SUM=0.0D0
    DWEI_SUM=0.0D0
!
    RAD=2.0D0*ELEPAR(1,ELE)
    DO I=1,NEleAss
       ELE1=LM(I,ELE)
       IF(ELE1.NE.0) THEN
          N=NDMFS(1,ELE1)
          S1=(XTMP-X(N))*(XTMP-X(N))+(YTMP-Y(N))*(YTMP-Y(N))
          S1=SQRT(S1)
          S=S1/RAD
          IF ((1.0D0-S).LT.1E-12) THEN
              WRITE(*,*) "ERROR, RADIUS IS TOO SMALL FOR ELEMENT=",ELE
              STOP
          ENDIF
          SS(1)=S*S
          SS(2)=SS(1)*S
          SS(3)=SS(2)*S
          SS(4)=-12.0D0*S+24.0D0*SS(1)-12.0D0*SS(2)
!
          IF (S.LT.1E-15) THEN
             WEI(I)=1.0D0
             DWEI(I,1)=0.0D0
             DWEI(I,2)=0.0D0
          ELSE
             WEI(I)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
             DWEI(I,1)=SS(4)*(XTMP-X(N))/S1/RAD
             DWEI(I,2)=SS(4)*(YTMP-Y(N))/S1/RAD
!
          ENDIF
             WEI_SUM=WEI_SUM+WEI(I)
             DWEI_SUM(1)=DWEI_SUM(1)+DWEI(I,1)
             DWEI_SUM(2)=DWEI_SUM(2)+DWEI(I,2)
        ENDIF
    ENDDO
!
    DO I=1,NEleAss
       WEI(I)=WEI(I)/WEI_SUM
       DWEI(I,1)=DWEI(I,1)/WEI_SUM-WEI(I)/WEI_SUM*DWEI_SUM(1)
       DWEI(I,2)=DWEI(I,2)/WEI_SUM-WEI(I)/WEI_SUM*DWEI_SUM(2)
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE WEI_CAL_QUADS
!
    SUBROUTINE OVLP_LINE_G(ADBC,NEleAss,THIC)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC,NEleAss
    REAL(8),INTENT(INOUT) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L,JJ,KK
    INTEGER :: NumBc
    INTEGER :: ALLocateStatus
    LOGICAL :: EXPRN1,INDBCX,INDBCY
    INTEGER :: L1,L2,N1,N2,N
    REAL(8) :: LENGTH,NX,NY
    REAL(8) :: XM,XL,YM,YL,XTMP,YTMP,WEIGHT
    REAL(8) :: FLOAD(2)
    INTEGER :: FUNC(2) 
    INTEGER :: NUM
    REAL(8) :: LOAD(2)
    INTEGER :: K1,K2
    INTEGER :: EleMin, EleMax
    INTEGER :: IND_X,IND_Y,INC_X,INC_Y,IDD_X,IDD_Y
    REAL(8) :: S(NPOR_MAX,2)
!   READ 
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    DO I=1,NumBc
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) CYCLE
!
       L1=BCED(1,I)
       L2=BCED(2,I)
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
!
       IF(INDBCX) THEN
          NDMFS(5,L1)=1
          NDMFS(5,L2)=1
       ENDIF
!
       IF(INDBCY) THEN
          NDMFS(6,L1)=1
          NDMFS(6,L2)=1
       ENDIF
!
       IF(INDBCX.OR.INDBCY) THEN
!
          N1=NDMFS(1,L1)
          N2=NDMFS(1,L2)
!
          LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
          LENGTH=SQRT(LENGTH)
!      CALCULATE NORMAL COMPONENTS
          NX=(Y(N2)-Y(N1))/LENGTH
          NY=-(X(N2)-X(N1))/LENGTH
          ELEPAR(3,L1)=-NX
          ELEPAR(4,L1)=-NY
          ELEPAR(3,L2)=-NX
          ELEPAR(4,L2)=-NY
       ENDIF
    ENDDO
!
!   CALCULATE LOAD VECTOR
!
    DO I=1,NumBc
       EleMin=100000
       EleMax=0
!
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) CYCLE
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
       IF(INDBCX.AND.INDBCY) CYCLE
!
       L1=BCED(1,I)
       L2=BCED(2,I)
!
       N1=NDMFS(1,L1)
       N2=NDMFS(1,L2)
!
       XM=(X(N1)+X(N2))/2.0D0
       XL=(X(N1)-X(N2))/2.0D0
       YM=(Y(N1)+Y(N2))/2.0D0
       YL=(Y(N1)-Y(N2))/2.0D0
!
       LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
       LENGTH=SQRT(LENGTH)
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       CALL H_QUADS(3,I,EleMin,EleMax,LENGTH,XM,XL,YM,YL)
!
       DO J=EleMin,EleMax
          N=NDMFS(1,J)
          DO K=1,NPOR_MAX
             S(K,1)=0.0D0
             S(K,2)=0.0D0
          ENDDO
!
          DO K=1,NITP_BC
             XTMP = XM+XL*XABSC_WEIG_BC(K,1) 
             YTMP = YM+YL*XABSC_WEIG_BC(K,1)
             LOAD=F_LOAD(FLOAD,FUNC,NX,NY,XTMP,YTMP,N1,N2)
             DO L=1,NPOR_MAX
                S(L,1)=S(L,1)+H_OVLP_QUADS(L,K,J)*LOAD(1)
                S(L,2)=S(L,2)+H_OVLP_QUADS(L,K,J)*LOAD(2)
             ENDDO
          ENDDO
!
          DO K=1,NPOR_MAX
             S(K,1)=S(K,1)*LENGTH/2.0D0*THIC
             S(K,2)=S(K,2)*LENGTH/2.0D0*THIC
          ENDDO
!
          IND_X=NDMFS(2,J)
          IND_Y=NDMFS(3,J)
          IF (IND_X.EQ.1) THEN
             INC_X=1
          ELSE
             INC_X=0
          ENDIF
!
          IF (IND_Y.EQ.1) THEN
             INC_Y=1
          ELSE
             INC_Y=0
          ENDIF
!
          IDD_X=ID(1,N)
          IDD_Y=ID(2,N)
          IF(IDD_X.NE.0) THEN
             DO K=1,NDOF(1,N)
                JJ=K+INC_X
                NUM=IDD_X+K-1
                R(NUM)=R(NUM)+S(JJ,1)
             ENDDO
          ENDIF
!
          IF(IDD_Y.NE.0) THEN
             DO K=1,NDOF(2,N)
                NUM=IDD_Y+K-1
                JJ=K+INC_Y
                R(NUM)=R(NUM)+S(JJ,2)
             ENDDO
          ENDIF
       ENDDO
!
       DO J=EleMin,EleMax
          DO K=1,NITP_BC
             DO L=1,NPOR_MAX
                H_OVLP_QUADS(L,K,J)=0.0D0
             ENDDO
          ENDDO
       ENDDO
!
    ENDDO      
!
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE OVLP_LINE_G    
!
    SUBROUTINE H_QUADS(NEleAss,NN,EleMin,EleMax,LENGTH,XM,XL,YM,YL)
!{{{
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NEleAss, NN
    INTEGER, INTENT(INOUT) :: EleMin,  EleMax
    REAL(8), INTENT(IN) :: LENGTH
    REAL(8), INTENT(IN) :: XM,XL,YM,YL
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    REAL(8) :: XTMP,YTMP,WEIGHT
    INTEGER :: ELE,N
    REAL(8) :: S,RAD,CONSTX,CONSTY,NX,NY,X_BAR,Y_BAR
    INTEGER :: IND,INDM,NDOFF
    REAL(8) :: H(8)
!
    DO I=1,NITP_BC
       XTMP = XM+XL*XABSC_WEIG_BC(I,1) 
       YTMP = YM+YL*XABSC_WEIG_BC(I,1)
       WEIGHT = XABSC_WEIG_BC(I,2)
       DO J=1,2
          ELE=BCED(J,NN)
!
          IF(ELE.LT.EleMin) EleMin=ELE
          IF(ELE.GT.EleMax) EleMax=ELE
!
          N=NDMFS(1,ELE)
          S=(XTMP-X(N))*(XTMP-X(N))+(YTMP-Y(N))*(YTMP-Y(N))
          S=SQRT(S)
!
          IND=NDMFS(2,ELE)
          IF(IND.EQ.1) THEN
             NDOFF=NDOF(1,N)+1
          ELSE
             NDOFF=NDOF(1,N)
          ENDIF
!
          RAD=2.0D0*ELEPAR(1,ELE)
          CONSTX=(XTMP-X(N))/RAD
          CONSTY=(YTMP-Y(N))/RAD
!
          INDM=NDMFS(5,ELE)
          IF(INDM.EQ.1) THEN
             NX=ELEPAR(3,ELE)
             NY=ELEPAR(4,ELE)
             X_BAR=(NX*CONSTX+NY*CONSTY)
             Y_BAR=(-NY*CONSTX+NX*CONSTY)
             SELECT CASE(NDOFF)
             CASE(1)
                  H(1)=1.0D0-S/LENGTH
!
                  H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
             CASE(2)
                  H(1)=1.0D0-S/LENGTH
                  H(2)=H(1)*X_BAR
!
                  H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                  H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
             CASE(3)
                  H(1)=1.0D0-S/LENGTH
                  H(2)=H(1)*X_BAR
                  H(3)=H(2)*Y_BAR
!
                  H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                  H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                  H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
             CASE(4)
                  H(1)=1.0D0-S/LENGTH
                  H(2)=H(1)*X_BAR
                  H(3)=H(2)*Y_BAR
                  H(4)=H(2)*X_BAR
!
                  H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                  H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                  H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
                  H_OVLP_QUADS(4,I,ELE)=H(4)*WEIGHT
             CASE(6)
                  H(1)=1.0D0-S/LENGTH
                  H(2)=H(1)*X_BAR
                  H(3)=H(2)*Y_BAR
                  H(4)=H(2)*X_BAR
                  H(5)=H(3)*X_BAR
                  H(6)=H(3)*Y_BAR
!
                  H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                  H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                  H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
                  H_OVLP_QUADS(4,I,ELE)=H(4)*WEIGHT
                  H_OVLP_QUADS(5,I,ELE)=H(5)*WEIGHT
                  H_OVLP_QUADS(6,I,ELE)=H(6)*WEIGHT
             CASE DEFAULT
                  WRITE(*,*) "ERROR, /GENERATEM/H_QUADS!"
                  STOP
             END SELECT
          ELSE
!
            SELECT CASE(NDOFF) 
            CASE(1)
                 H(1)=(1.0D0-S/LENGTH)
!
                 H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
            CASE(3)
                 H(1)=(1.0D0-S/LENGTH)
                 H(2)=H(1)*CONSTX
                 H(3)=H(1)*CONSTY
!
                 H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                 H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                 H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
            CASE(4)
                 H(1)=(1.0D0-S/LENGTH)
                 H(2)=H(1)*CONSTX
                 H(3)=H(1)*CONSTY
                 H(4)=H(2)*CONSTY
!
                 H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                 H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                 H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
                 H_OVLP_QUADS(4,I,ELE)=H(4)*WEIGHT
            CASE(6)
                 H(1)=(1.0D0-S/LENGTH)
                 H(2)=H(1)*CONSTX
                 H(3)=H(1)*CONSTY
                 H(4)=H(2)*CONSTY
                 H(5)=H(2)*CONSTX
                 H(6)=H(3)*CONSTY
!
                 H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                 H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                 H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
                 H_OVLP_QUADS(4,I,ELE)=H(4)*WEIGHT
                 H_OVLP_QUADS(5,I,ELE)=H(5)*WEIGHT
                 H_OVLP_QUADS(6,I,ELE)=H(6)*WEIGHT
           CASE(8)
                 H(1)=(1.0D0-S/LENGTH)
                 H(2)=H(1)*CONSTX
                 H(3)=H(1)*CONSTY
                 H(4)=H(2)*CONSTY
                 H(5)=H(2)*CONSTX
                 H(6)=H(3)*CONSTY
                 H(7)=H(4)*CONSTX
                 H(8)=H(4)*CONSTY
!
                 H_OVLP_QUADS(1,I,ELE)=H(1)*WEIGHT
                 H_OVLP_QUADS(2,I,ELE)=H(2)*WEIGHT
                 H_OVLP_QUADS(3,I,ELE)=H(3)*WEIGHT
                 H_OVLP_QUADS(4,I,ELE)=H(4)*WEIGHT
                 H_OVLP_QUADS(5,I,ELE)=H(5)*WEIGHT
                 H_OVLP_QUADS(6,I,ELE)=H(6)*WEIGHT
                 H_OVLP_QUADS(7,I,ELE)=H(7)*WEIGHT
                 H_OVLP_QUADS(8,I,ELE)=H(8)*WEIGHT
            CASE DEFAULT
                 WRITE(*,*) "ERROR, /GENERATEM/H_QUADS!"
                 STOP
            END SELECT
!
          ENDIF
       ENDDO       
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE H_QUADS
!
    SUBROUTINE WEI_CAL_QUADS2(ELE,NEleAss,XTMP,YTMP,WEI)
!{{{
      IMPLICIT NONE
      INTEGER :: ELE,NEleAss
      REAL(8) :: XTMP,YTMP
      REAL(8) :: WEI(NEleAss)
!   ASSISTANT VARIABLES
      INTEGER :: I
      INTEGER :: ELE1,N
      REAL(8) :: RAD
      REAL(8) :: S,S1,SS(3)
      REAL(8) :: WEI_SUM
!
      DO I=1,NEleAss
         WEI(I)=0.0D0
      ENDDO
      WEI(1)=1.0D0
!    
      WEI_SUM=0.0D0
!
      RAD=2.0D0*ELEPAR(1,ELE)
      DO I=1,NEleAss
         ELE1=LM(I,ELE)
         IF(ELE1.NE.0) THEN
!
            N=NDMFS(1,ELE1)
            S1=(XTMP-X(N))*(XTMP-X(N))+(YTMP-Y(N))*(YTMP-Y(N))
            S1=SQRT(S1)
            S=S1/RAD
            IF ((1.0D0-S).LT.1E-12) THEN
                WRITE(*,*) "ERROR, RADIUS IS TOO SMALL FOR ELEMENT=",ELE
                STOP
            ENDIF
!
            SS(1)=S*S
            SS(2)=SS(1)*S
            SS(3)=SS(2)*S
!
            WEI(I)=1.0D0-6.0D0*SS(1)+8.0D0*SS(2)-3.0D0*SS(3)
!
            WEI_SUM=WEI_SUM+WEI(I)
         ENDIF
      ENDDO
!
      DO I=1,NEleAss
         WEI(I)=WEI(I)/WEI_SUM
      ENDDO
!
      RETURN
!}}}
    END SUBROUTINE WEI_CAL_QUADS2
!
    SUBROUTINE LeastSQ(IND,NN,MM,ELE1,Rmatrix,Qmatrix,Wmatrix)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IND,NN,MM,ELE1
    REAL(8),DIMENSION(NN,NN) :: Rmatrix
    REAL(8),DIMENSION(NN,MM) :: Qmatrix
    REAL(8),DIMENSION(MM)    :: Wmatrix
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: N,N1,ELE
    REAL(8) :: RAD1,RAD,SS,S(3)
    REAL(8),DIMENSION(MM,NN) :: Amatrix
    REAL(8),DIMENSION(MM,NN) :: Vmatrix
    REAL(8) :: XNORM2,XNORM,VNORM,X1
!
!   CONSTRUCT Q MATRIX
    DO I=1,MM
       DO J=1,NN
          Qmatrix(J,I)=0D0
       ENDDO
    ENDDO
!
    DO I=1,NN
       Qmatrix(I,I)=1.0D0
    ENDDO
!   Reset R matrix
    DO I=1,NN
       DO J=1,NN
          Rmatrix(J,I)=0.0D0
       ENDDO
    ENDDO
    IF(IND.NE.0) RETURN
!   Construct weighted matrix and A matrix
    N1=NDMFS(1,ELE1)
    RAD1=ELEPAR(1,ELE1)
    DO I=1,MM
       ELE=LM(I,ELE1)
       N=NDMFS(1,ELE)
       SS=ABS(X(N)-X(N1))/RAD1
       IF (SS.LT.1E-12) THEN
           Wmatrix(I)=1.0D0
       ELSE
           IF ((1.0-SS).GT.1E-12) THEN
               S(1)=SS*SS
               S(2)=S(1)*SS
               S(3)=S(2)*SS
               Wmatrix(I)=1-6*S(1)+8*S(2)-3*S(3)
           ENDIF
       ENDIF
!
       SELECT CASE(NN)
       CASE(2)
           Amatrix(I,1)=1.0D0*Wmatrix(I)
           Amatrix(I,2)=(X(N)-X(N1))/RAD1*Wmatrix(I)
       CASE(3)
           Amatrix(I,1)=1.0D0*Wmatrix(I)
           Amatrix(I,2)=(X(N)-X(N1))/RAD1
           Amatrix(I,3)=Amatrix(I,2)*Amatrix(I,2)
!
           Amatrix(I,2)=Amatrix(I,2)*Wmatrix(I)
           Amatrix(I,3)=Amatrix(I,3)*Wmatrix(I)
       CASE DEFAULT
           WRITE(*,*) "ERROR, /GENERATEM/LeastSQ!"
           STOP
       END SELECT
    ENDDO
!    WRITE(*,*) "*********Amatrix"
!    DO I=1,MM
!       WRITE(*,*) Amatrix(I,:)
!    ENDDO
!
!   Householder QR Factorization
!
    DO K=1,NN
       XNORM2=0.0D0
       DO I=K,MM
          XNORM2=XNORM2+Amatrix(I,K)*Amatrix(I,K)
       ENDDO
       X1=Amatrix(K,K)
       XNORM=SQRT(XNORM2)
       VNORM=2*XNORM2+2*SIGN(1.0D0,X1)*XNORM*X1
       VNORM=SQRT(VNORM)
!
       Vmatrix(K,K)=(SIGN(1.0D0,X1)*XNORM+X1)/VNORM
       DO I=K+1,MM
          Vmatrix(I,K)=Amatrix(I,K)/VNORM
       ENDDO
!   Calculate the R matrix
       DO I=K,NN
          XNORM2=0.0D0
          DO J=K,MM
             XNORM2=XNORM2+Vmatrix(J,K)*Amatrix(J,I)
          ENDDO
          DO J=K,MM
             Amatrix(J,I)=Amatrix(J,I)-2*Vmatrix(J,K)*XNORM2
          ENDDO
       ENDDO
!
       DO I=1,K
          Rmatrix(I,K)=Amatrix(I,K)
       ENDDO
    ENDDO
!   Calculate Q matrix
    DO K=NN,1,-1        
       DO I=K,NN
          XNORM2=0.0D0
          IF (I.EQ.K) THEN
             XNORM2=Vmatrix(K,K)
          ELSE
             DO J=K,MM
                XNORM2=XNORM2+Vmatrix(J,K)*Qmatrix(I,J)
             ENDDO
          ENDIF
          DO J=K,MM
             Qmatrix(I,J)=Qmatrix(I,J)-2*Vmatrix(J,K)*XNORM2         
          ENDDO
       ENDDO
    ENDDO
!
!    WRITE(*,*) "*********Rmatrix"
!    DO I=1,NN
!       WRITE(*,*) Rmatrix(I,:)
!    ENDDO
!
!    WRITE(*,*) "*********Qmatrix"
!    DO I=1,NN
!       WRITE(*,*) Qmatrix(I,:)
!    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE LeastSQ 
!
    SUBROUTINE TRUSS_OVLP_STIFF(NITG,NEleAss,NN,E) 
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NITG,NEleAss,NN
    REAL(8),INTENT(IN) :: E
!   IN-FUNCTION VARIALBLES
    INTEGER :: I,J,K,L,M
    INTEGER :: KK,LL
    INTEGER :: II,JJ,NII,NJJ,IND1,IND2
    INTEGER :: ELE1,ELE2,N1,N2
    REAL(8) :: XM,XL
    LOGICAL :: EXPN
!
    ELE1=NDINT(1,NN)
    ELE2=NDINT(2,NN)
    N1=NDMFS(1,ELE1)
    N2=NDMFS(1,ELE2)
!
    XM=(X(N1)+X(N2))/2.0
    XL=ABS(X(N1)-X(N2))/2.0
!   GENERATE DHI_OVLP_TRUSS
    SELECT CASE(NPAR(1))
    CASE(6)
         CALL TRUSS_OVLP_DHI(NITG,NEleAss,NN) 
    CASE(7)
         CALL TRUSS_OVLP2_DHI(NITG,NEleAss,NN)
    CASE(8)
         CONTINUE
    CASE DEFAULT
         WRITE(*,*) "*** ERROR ***, /GENERATEM/TRUSS_OVLP_STIFF!"
         Stop
    END SELECT            
!   BODY FORCE
!    CALL TRUSS_OVLP_BODY(NITG,NEleAss,NN)
!   NODE 1 & NODE 1
    DO I=1,2
       ELE1=NDINT(I,NN)
       DO J=I,2
          ELE2=NDINT(J,NN)
          IF(I.EQ.J) THEN
             EXPN=.FALSE.
             DO K=1,NEleAss 
                II=LM(K,ELE1)
                IF(II.NE.0) THEN 
                   KK=(I-1)*NEleAss+K
                   NII=NDMFS(1,II)
                   IND1=NDMFS(2,II)
                   DO L=K,NEleAss
                      JJ=LM(L,ELE2)
                      IF (JJ.NE.0) THEN
                         LL=(J-1)*NEleAss+L
                         NJJ=NDMFS(1,JJ)
                         IND2=NDMFS(2,JJ)
                         CALL TRUSS_OVLP_STIFF_IJ(KK,LL,NITG,E,XM,XL)
                         CALL TRUSS_OVLP_ADDBAN(IND1,IND2,NII,NJJ,STIFF_IJ(1,1),EXPN)
                      ENDIF
                   ENDDO 
                ENDIF
             ENDDO 
          ELSE
             DO K=1,NEleAss 
                II=LM(K,ELE1)
                IF(II.NE.0) THEN
                   KK=(I-1)*NEleAss+K
                   NII=NDMFS(1,II)
                   IND1=NDMFS(2,II)
                   DO L=1,NEleAss
                      JJ=LM(L,ELE2)
                      IF (JJ.NE.0) THEN
                         LL=(J-1)*NEleAss+L
                         NJJ=NDMFS(1,JJ)
                         IND2=NDMFS(2,JJ)
                         CALL TRUSS_OVLP_STIFF_IJ(KK,LL,NITG,E,XM,XL)
                         IF (II.EQ.JJ) THEN
                             EXPN=.TRUE.
                         ELSE
                             EXPN=.FALSE.
                         ENDIF
                         CALL TRUSS_OVLP_ADDBAN(IND1,IND2,NII,NJJ,STIFF_IJ(1,1),EXPN)
                      ENDIF
                   ENDDO 
                ENDIF
             ENDDO 
          ENDIF
       ENDDO
    ENDDO    
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_STIFF
!
    SUBROUTINE TRUSS_OVLP_STIFF3(NITG,NEleAss,NN,E) 
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NITG,NEleAss,NN
    REAL(8),INTENT(IN) :: E
!   IN-FUNCTION VARIALBLES
    INTEGER :: I,J,K
    INTEGER :: LL
    INTEGER :: ELE1,ELE2,N1,N2
    REAL(8) :: XM,XL,LENGTH,XTMP
    INTEGER :: IDDMin,IDDMax
    REAL(8) :: STIFF,AREA,RLOAD
    REAL(8) :: WEIGHT
!
    ELE1=NDINT(1,NN)
    ELE2=NDINT(2,NN)
    N1=NDMFS(1,ELE1)
    N2=NDMFS(1,ELE2)
!
    XM=(X(N1)+X(N2))/2.0D0
    XL=ABS(X(N1)-X(N2))/2.0D0
    LENGTH=ABS(X(N1)-X(N2))
    
!   GENERATE TRUSS_OVLP_H
    DO I=1,NITG
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)          
       CALL H_TRUSS(I,NEleAss,NN,XTMP,LENGTH)
    ENDDO 
!    
    IDDMin=NDINT(3,NN)
    IDDMax=NDINT(4,NN)
!    
    DO I=IDDMin,IDDMax
!   BODY FORCE
!       RLOAD=0.0D0
!       DO K=1,NITG
!          AREA=1.0D0
!          XTMP=XM+XL*XABSC_WEIG_IN(K,1)
!          WEIGHT=XABSC_WEIG_IN(K,2)
!          RLOAD=RLOAD+H_OVLP_TRUSS(1,K,I)*XTMP*WEIGHT*AREA
!       ENDDO
!       R(I)=R(I)+RLOAD*XL
!
       DO J=I,IDDMax
          STIFF=0.0D0
          DO K=1,NITG
             XTMP=XM+XL*XABSC_WEIG_IN(K,1)
             AREA=1.0D0+XTMP
             WEIGHT=XABSC_WEIG_IN(K,2) 
             STIFF=STIFF+H_OVLP_TRUSS(2,K,I)*H_OVLP_TRUSS(2,K,J)*E*AREA*WEIGHT
          ENDDO
          STIFF=STIFF*XL
!
          LL=MAXA_COLSOL(J)+(J-I)
          K_COLSOL(LL)=K_COLSOL(LL)+STIFF
          KMATRIX(LL)=KMATRIX(LL)+STIFF        
       ENDDO
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_STIFF3
!
    SUBROUTINE H_TRUSS(II,NEleAss,NN,XTMP,LENGTH)
!{{{
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: II,NEleAss,NN
    REAL(8), INTENT(IN) :: XTMP,LENGTH
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: ELE,N,IND,INC,N1,IND1
    REAL(8) :: RAD,X_BAR,XX1,RAD1
    REAL(8),DIMENSION(NPOR_MAX+1) :: HH,DHH
    INTEGER :: M
!   Least square fitting
    REAL(8) :: HI, DHI
    REAL(8),DIMENSION(3) :: POLY_B,DPOLY_B
    REAL(8),DIMENSION(NEleAss) :: POLY_DISP,DPOLY_DISP
    INTEGER :: ELE1,IDD
    INTEGER :: NDOFF
! 
    DO K=1,2
       ELE=NDINT(K,NN)
       N=NDMFS(1,ELE)
       IND=NDMFS(2,ELE)
       RAD=2*ELEPAR(1,ELE)
       NDOFF=NDOF(1,N)
       IDD=ID(1,N)
       IF(IND.NE.0) THEN
          IF(IDD.NE.0) THEN
            IF(IND.EQ.1) THEN
                INC=1
            ELSE
                INC=0
            ENDIF
!
            X_BAR=(XTMP-X(N))/RAD
            SELECT CASE(NDOFF)
            CASE(1)
              XX1=XTMP-X(N)
              HH(1)=(1.0-ABS(XX1)/LENGTH)
              HH(2)=HH(1)*X_BAR
!
              DHH(1)=-SIGN(1.0D0,XX1)/LENGTH
              DHH(2)=DHH(1)*X_BAR+HH(1)/RAD
!
              H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1+INC)
              H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1+INC)
            CASE(2)
              XX1=XTMP-X(N)
              HH(1)=(1.0-ABS(XX1)/LENGTH)
              HH(2)=HH(1)*X_BAR
              HH(3)=HH(2)*X_BAR
!
              DHH(1)=-SIGN(1.0D0,XX1)/LENGTH
              DHH(2)=DHH(1)*X_BAR+HH(1)/RAD
              DHH(3)=DHH(2)*X_BAR+HH(2)/RAD
!
              H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1+INC)
              H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2+INC)
              H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1+INC)
              H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2+INC)
            CASE(3)
              XX1=XTMP-X(N)
              HH(1)=(1.0-ABS(XX1)/LENGTH)
              HH(2)=HH(1)*X_BAR
              HH(3)=HH(2)*X_BAR
!
              DHH(1)=-SIGN(1.0D0,XX1)/LENGTH
              DHH(2)=DHH(1)*X_BAR+HH(1)/RAD
              DHH(3)=DHH(2)*X_BAR+HH(2)/RAD
!
              H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1)
              H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2)
              H_OVLP_TRUSS(1,II,IDD+2)=H_OVLP_TRUSS(1,II,IDD+2)+HH(3)
!
              H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1)
              H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2)
              H_OVLP_TRUSS(2,II,IDD+2)=H_OVLP_TRUSS(2,II,IDD+2)+DHH(3)
            CASE DEFAULT
              WRITE(*,*) "ERROR, NDOFF.GT.3!"
              STOP
            END SELECT
          ENDIF
       ELSE
          M=NumLM(ELE)
          RAD=ELEPAR(1,ELE)
!
          XX1=XTMP-X(N)
          X_BAR=XX1/RAD
!
          HI=(1.0-ABS(XX1)/LENGTH)
          DHI=-SIGN(1.0D0,XX1)/LENGTH
!
          POLY_B(1)=HI
          POLY_B(2)=HI*X_BAR
          POLY_B(3)=POLY_B(2)*X_BAR
!
          DPOLY_B(1)=DHI
          DPOLY_B(2)=DHI*X_BAR+HI/RAD
          DPOLY_B(3)=DPOLY_B(2)*X_BAR+POLY_B(2)/RAD
!       (R-1) T * 
!       Forward substitution
          DO I=1,3
             POLY_B(I)=POLY_B(I)/Rmatrix_all(I,I,ELE)
             DPOLY_B(I)=DPOLY_B(I)/Rmatrix_all(I,I,ELE)
             DO J=I+1,3
                POLY_B(J)=POLY_B(J)-Rmatrix_all(I,J,ELE)*POLY_B(I)
                DPOLY_B(J)=DPOLY_B(J)-Rmatrix_all(I,J,ELE)*DPOLY_B(I)
             ENDDO         
          ENDDO
!       W * Q *
          DO I=1,M
             POLY_DISP(I)=0.0D0
             DPOLY_DISP(I)=0.0D0
             DO J=1,3
                POLY_DISP(I)=POLY_DISP(I)+Qmatrix_all(J,I,ELE)*POLY_B(J)
                DPOLY_DISP(I)=DPOLY_DISP(I)+Qmatrix_all(J,I,ELE)*DPOLY_B(J)
             ENDDO
             POLY_DISP(I)=POLY_DISP(I)*Wmatrix_all(I,ELE)
             DPOLY_DISP(I)=DPOLY_DISP(I)*Wmatrix_all(I,ELE)
          ENDDO
!   ****************************************
!            X_BAR=(XTMP-X(N))/RAD
!            SELECT CASE(NDOFF-1)
!            CASE(0) 
!              CONTINUE
!            CASE(1)
!              XX1=XTMP-X(N)
!              HH(1)=(1.0-ABS(XX1)/LENGTH)
!              HH(2)=HH(1)*X_BAR*X_BAR*X_BAR
!!
!              DHH(1)=-SIGN(1.0D0,XX1)/LENGTH
!              DHH(2)=DHH(1)*X_BAR*X_BAR*X_BAR+3*HH(1)*X_BAR*X_BAR/RAD
!!
!              H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2)
!              H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2)
!            CASE(2)
!              XX1=XTMP-X(N)
!              HH(1)=(1.0-ABS(XX1)/LENGTH)
!              HH(2)=HH(1)*X_BAR*X_BAR*X_BAR
!              HH(3)=HH(2)*X_BAR
!!
!              DHH(1)=-SIGN(1.0D0,XX1)/LENGTH
!              DHH(2)=DHH(1)*X_BAR*X_BAR*X_BAR+3*HH(1)*X_BAR*X_BAR/RAD
!              DHH(3)=DHH(2)*X_BAR+HH(2)/RAD
!!
!              H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2)
!              H_OVLP_TRUSS(1,II,IDD+2)=H_OVLP_TRUSS(1,II,IDD+2)+HH(3)
!              H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2)
!              H_OVLP_TRUSS(2,II,IDD+2)=H_OVLP_TRUSS(2,II,IDD+2)+DHH(3)
!            CASE DEFAULT
!              WRITE(*,*) "ERROR, NDOFF.GT.3!"
!              STOP
!            END SELECT
!   ****************************************          
          DO I=1,M
             ELE1=LM(I,ELE)
             N1=NDMFS(1,ELE1)
             IND1=NDMFS(2,ELE1)
             IDD=ID(1,N1)
!
             IF(IND1.EQ.1) THEN
                INC=1
                NDOFF=0
             ELSEIF(IND1.EQ.2) THEN
                INC=0
                NDOFF=1
             ELSE
                INC=0
                NDOFF=NDOF(1,N1)
             ENDIF
!
             RAD1=2*ELEPAR(1,ELE1)
             X_BAR=(XTMP-X(N1))/RAD1
!          IND1.EQ.0
             SELECT CASE(NDOFF)
             CASE(0)
                CONTINUE
             CASE(1)
                HH(1)=POLY_DISP(I)
!                HH(2)=HH(1)*X_BAR

                DHH(1)=DPOLY_DISP(I)
!                DHH(2)=DHH(1)*X_BAR+HH(1)/RAD1

                H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1+INC)
                H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1+INC)  
             CASE(2)
                HH(1)=POLY_DISP(I)
                HH(2)=HH(1)*X_BAR
!                HH(3)=HH(2)*X_BAR*X_BAR

                DHH(1)=DPOLY_DISP(I)
                DHH(2)=DHH(1)*X_BAR+HH(1)/RAD1
!                DHH(3)=DHH(2)*X_BAR+HH(2)/RAD1

                H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1+INC)
                H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2+INC)

                H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1+INC)  
                H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2+INC)             
             CASE(3)
                HH(1)=POLY_DISP(I)
                HH(2)=HH(1)*X_BAR
                HH(3)=HH(2)*X_BAR

                DHH(1)=DPOLY_DISP(I)
                DHH(2)=DHH(1)*X_BAR+HH(1)/RAD1
                DHH(3)=DHH(2)*X_BAR+HH(2)/RAD1

                H_OVLP_TRUSS(1,II,IDD)=H_OVLP_TRUSS(1,II,IDD)+HH(1)
                H_OVLP_TRUSS(1,II,IDD+1)=H_OVLP_TRUSS(1,II,IDD+1)+HH(2)
                H_OVLP_TRUSS(1,II,IDD+2)=H_OVLP_TRUSS(1,II,IDD+2)+HH(3)

                H_OVLP_TRUSS(2,II,IDD)=H_OVLP_TRUSS(2,II,IDD)+DHH(1)  
                H_OVLP_TRUSS(2,II,IDD+1)=H_OVLP_TRUSS(2,II,IDD+1)+DHH(2) 
                H_OVLP_TRUSS(2,II,IDD+2)=H_OVLP_TRUSS(2,II,IDD+2)+DHH(3)                 
             CASE DEFAULT
                WRITE(*,*) "ERROR, /GENERATEM/H_TRUSS!"
                STOP
             END SELECT
          ENDDO
       ENDIF
    ENDDO   
!
    RETURN
!}}}
    END SUBROUTINE H_TRUSS 
!
    SUBROUTINE TRUSS_OVLP_BODY(NITG,NEleAss,NN)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NITG,NEleAss,NN
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: ELE,ELE1,I2
    INTEGER :: N1,N2
    INTEGER :: N,IND,IDD,IDD2
    INTEGER :: AREA
    REAL(8) :: WEIGHT
    REAL(8) :: FORCE,RLOAD
    REAL(8) :: XL,XM,XTMP
    INTEGER :: NDOFF,INC
!
    FORCE=1.0
!
    ELE=NDINT(1,NN)
    ELE1=NDINT(2,NN)
    N1=NDMFS(1,ELE)
    N2=NDMFS(1,ELE1)
    XM=(X(N1)+X(N2))/2.0
    XL=ABS(X(N1)-X(N2))/2.0
    AREA=INTPAR(1,NN)
!   
    DO I=1,2
       ELE=NDINT(I,NN)
       DO J=1,NEleAss
          I2=(I-1)*NEleAss+J
          ELE1=LM(J,ELE)
          IF(ELE1.NE.0) THEN
             N=NDMFS(1,ELE1)
             IND=NDMFS(2,ELE1)
             IF(IND.EQ.1) THEN
                NDOFF=NDOF(1,N)+1
                INC=2
             ELSE
                NDOFF=NDOF(1,N)
                INC=1
             ENDIF
!
             DO K=INC,NDOFF
                IDD2=ID(1,N)+K-INC
                RLOAD=0.0
                DO L=1,NITG
                   XTMP=XM+XL*XABSC_WEIG_IN(L,1)
                   WEIGHT=XABSC_WEIG_IN(L,2)
                   RLOAD=RLOAD+HI_OVLP_TRUSS(L,K,I2)*WEIGHT*(1.0d0+XTMP)
                ENDDO
                RLOAD=RLOAD*XL*AREA
                R(IDD2)=R(IDD2)+RLOAD
             ENDDO
          ENDIF      
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_BODY
!
    SUBROUTINE TRUSS_OVLP_STIFF_IJ(NEL1,NEL2,NITG,E,XM,XL)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NEL1,NEL2,NITG
    REAL(8),INTENT(IN) :: E,XM,XL
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J,K
    INTEGER :: L1,L2
    REAL(8) :: WEIGHT
    REAL(8) :: XTMP
    REAL(8) :: AREA
!
    L1=NPOR_MAX*NPOR_MAX
!
    DO I=1,L1
       STIFF_IJ(I,1)=0.0D0
    ENDDO
!
    DO I=1,NITG
      XTMP=XM+XL*XABSC_WEIG_IN(I,1)
      AREA=1.0D0+XTMP
      WEIGHT=XABSC_WEIG_IN(I,2)
      DO J=1,NPOR_MAX
         DO K=1,NPOR_MAX
          L2=(J-1)*NPOR_MAX+K
          STIFF_IJ(L2,1)=STIFF_IJ(L2,1)+DHI_OVLP_TRUSS(I,J,NEL1)*DHI_OVLP_TRUSS(I,K,NEL2)*WEIGHT*AREA
         ENDDO
      ENDDO
    ENDDO
!
    DO I=1,L1
       STIFF_IJ(I,1)=STIFF_IJ(I,1)*XL*E
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_STIFF_IJ
!
    SUBROUTINE TRUSS_OVLP_ADDBAN(INDN,INDM,N,M,S,EXPN)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: INDN,INDM,N,M
    REAL(8) :: S(NPOR_MAX*NPOR_MAX)
    LOGICAL :: EXPN
!   ASSISTANT VARIABLES
    INTEGER :: K,INDK,L,INDL
    INTEGER :: NPOR1,NPOR2
    INTEGER :: KI,KJ,KK2,LL1
    REAL(8) :: CONST
!
    NPOR1=NDOF(1,N)
    NPOR2=NDOF(1,M)
!
    K=1
    INDK=0
    IF (INDN.EQ.1) THEN
        INDK=1
    ENDIF
    DO WHILE(K.LE.NPOR1)
       L=1
       INDL=0
       IF (INDM.EQ.1) THEN
           INDL=1
       ENDIF
!
       DO WHILE(L.LE.NPOR2)
          KI=ID(1,N)+K-1
          KJ=ID(1,M)+L-1
          KK2=(K+INDK-1)*NPOR_MAX+L+INDL
          IF (EXPN.AND.((K+INDK).EQ.(L+INDL))) THEN
              CONST=2.0
          ELSE
              CONST=1.0
          ENDIF
          IF (EXPN.OR.((N.NE.M).OR.(KJ.GE.KI))) THEN
             IF(KJ.LE.KI) THEN
                LL1=MAXA_COLSOL(KI)+ABS(KJ-KI)
             ELSE
                LL1=MAXA_COLSOL(KJ)+ABS(KJ-KI)
             ENDIF
             K_COLSOL(LL1)=K_COLSOL(LL1)+S(KK2)*CONST
             KMATRIX(LL1)=KMATRIX(LL1)+S(KK2)*CONST
          ENDIF
          L=L+1
       ENDDO
       K=K+1
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_ADDBAN
!
    SUBROUTINE TRUSS_OVLP_DHI(NITG,NEleAss,NN)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NITG,NEleAss,NN
!   IN-FUNCTION VARIABLES
    INTEGER :: ELE1,ELE2,N1,N2
    REAL(8) :: XM,XL,XTMP,DIV
    INTEGER :: ELE,N
    REAL(8) :: RAD
    INTEGER :: I,J
!
    ELE1=NDINT(1,NN)
    ELE2=NDINT(2,NN)
    N1=NDMFS(1,ELE1)
    N2=NDMFS(1,ELE2)
!
    XM=(X(N1)+X(N2))/2.0
    XL=ABS(XM-X(N1))
! 
    DO I=1,NITG
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)          
       CALL DHI_TRUSS(I,XTMP,X(N1),X(N2),ELE1,ELE2,NEleAss)
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP_DHI
!
    SUBROUTINE TRUSS_OVLP2_DHI(NITG,NEleAss,NN)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NITG,NEleAss,NN
!   IN-FUNCTION VARIABLES
    INTEGER :: ELE1,ELE2,N1,N2
    REAL(8) :: XM,XL,XTMP,DIV
    INTEGER :: ELE,N
    REAL(8) :: RAD
    INTEGER :: I,J
    REAL(8) :: WEI(NEleAss,2),DWEI(NEleAss,2)
!
    ELE1=NDINT(1,NN)
    ELE2=NDINT(2,NN)
    N1=NDMFS(1,ELE1)
    N2=NDMFS(1,ELE2)
!
    XM=(X(N1)+X(N2))/2.0
    XL=ABS(XM-X(N1))
! 
    DO I=1,NITG
       XTMP=XM+XL*XABSC_WEIG_IN(I,1)          
       CALL WEI_CAL(NN,NEleAss,XTMP,WEI(1,1),DWEI(1,1))
       CALL DHI_TRUSS2(I,XTMP,X(N1),X(N2),ELE1,ELE2,NEleAss,WEI(1,1),DWEI(1,1))
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE TRUSS_OVLP2_DHI
!
    SUBROUTINE WEI_CAL(NN,NEleAss,XTMP,WEI,DWEI)
!{{{
    IMPLICIT NONE
    INTEGER :: NN,NEleAss
    REAL(8) :: XTMP
    REAL(8) :: WEI(NEleAss,2),DWEI(NEleAss,2)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,N
    INTEGER :: ELE1,ELE2
    INTEGER :: IND
    REAL(8) :: WEIGHT
    REAL(8) :: RAD
    REAL(8) :: S,SS(4)
!
    DO I=1,NEleAss
       WEI(I,1)=0.0D0
       WEI(I,2)=0.0D0
       DWEI(I,1)=0.0D0
       DWEI(I,2)=0.0D0
    ENDDO
!
    DO I=1,2
       ELE1=NDINT(I,NN)
       IND=NDMFS(2,ELE1)
       RAD=ELEPAR(1,ELE1)
       IF(IND.EQ.0) THEN
          DO J=1,NEleAss
             ELE2=LM(J,ELE1)
             IF(ELE2.NE.0) THEN
                N=NDMFS(1,ELE2)
                WEIGHT=(XTMP-X(N))/2.0D0/RAD
                S=ABS(WEIGHT)
                IF ((1-S).LT.1E-12) CYCLE
                SS(1)=S*S
                SS(2)=SS(1)*S
                SS(3)=SS(2)*S
                SS(4)=-12*S+24*SS(1)-12*SS(2)
!  WEIGHT 2
                WEI(J,I)=1-6*SS(1)+8*SS(2)-3*SS(3)
                DWEI(J,I)=SS(4)*SIGN(1.0D0,WEIGHT)/2.0D0/RAD
!  WEIGHT 1
!                WEIGHT=(XTMP-X(N))/2.0D0/RAD
!                IF ((1-S).LT.1E-12) CYCLE
!                WEI(J,I)=1.0-S
!                DWEI(J,I)=-SIGN(1.0D0,WEIGHT)/2.0D0/RAD
             ENDIF
          ENDDO
       ELSE
          DO J=1,NEleAss
             ELE2=LM(J,ELE1)
             IF(ELE2.EQ.ELE1) WEI(J,I)=1.0D0
          ENDDO
       ENDIF
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE WEI_CAL 
!
    SUBROUTINE DHI_TRUSS(II,XTMP,XX1,XX2,ELE1,ELE2,NEleAss)
!{{{
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: II
    REAL(8), INTENT(IN) :: XTMP,XX1,XX2
    INTEGER, INTENT(IN) :: ELE1,ELE2,NEleAss
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J
    INTEGER :: ELE,N,IND,NDOFF
    REAL(8) :: RAD,LENGTH
    REAL(8) :: WEIGHT
    INTEGER :: I2
    REAL(8) :: X_BAR
!
    LENGTH=ABS(XX1-XX2)
    DO I=1,NEleAss
       ELE=LM(I,ELE1)
       IF(ELE.NE.0) THEN
          WEIGHT=WEIND(I,ELE1)
          N=NDMFS(1,ELE)
          IND=NDMFS(2,ELE)
          RAD=ELEPAR(1,ELE)
!
          IF(IND.EQ.1) THEN
             NDOFF=NDOF(1,N)+1
          ELSE
             NDOFF=NDOF(1,N)
          ENDIF
!
          X_BAR=(XTMP-X(N))/RAD
!
          SELECT CASE(NDOFF)
          CASE(1)
             HI_OVLP_TRUSS(II,1,I)=(1.0-(XTMP-XX1)/LENGTH)*WEIGHT
!
             DHI_OVLP_TRUSS(II,1,I)=-WEIGHT/LENGTH
          CASE(2)
             HI_OVLP_TRUSS(II,1,I)=1.0-(XTMP-XX1)/LENGTH
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,1,I)*X_BAR

             DHI_OVLP_TRUSS(II,1,I)=-1.0/LENGTH
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,1,I)*X_BAR+HI_OVLP_TRUSS(II,1,I)/RAD

             HI_OVLP_TRUSS(II,1,I)=HI_OVLP_TRUSS(II,1,I)*WEIGHT
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,2,I)*WEIGHT
             DHI_OVLP_TRUSS(II,1,I)=DHI_OVLP_TRUSS(II,1,I)*WEIGHT
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,2,I)*WEIGHT
          CASE(3)
             HI_OVLP_TRUSS(II,1,I)=1.0-(XTMP-XX1)/LENGTH
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,1,I)*X_BAR
             HI_OVLP_TRUSS(II,3,I)=HI_OVLP_TRUSS(II,2,I)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I)=-1.0/LENGTH
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,1,I)*X_BAR+HI_OVLP_TRUSS(II,1,I)/RAD
             DHI_OVLP_TRUSS(II,3,I)=DHI_OVLP_TRUSS(II,2,I)*X_BAR+HI_OVLP_TRUSS(II,2,I)/RAD
!
             HI_OVLP_TRUSS(II,1,I)=HI_OVLP_TRUSS(II,1,I)*WEIGHT
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,2,I)*WEIGHT
             HI_OVLP_TRUSS(II,3,I)=HI_OVLP_TRUSS(II,3,I)*WEIGHT 
             DHI_OVLP_TRUSS(II,1,I)=DHI_OVLP_TRUSS(II,1,I)*WEIGHT
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,2,I)*WEIGHT
             DHI_OVLP_TRUSS(II,3,I)=DHI_OVLP_TRUSS(II,3,I)*WEIGHT          
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GNERATEM/DHI_TRUSS!"
             STOP
          END SELECT
       ENDIF
!
       ELE=LM(I,ELE2)
       IF(ELE.NE.0) THEN
          WEIGHT=WEIND(I,ELE2)
          N=NDMFS(1,ELE)
          IND=NDMFS(2,ELE)
          RAD=ELEPAR(1,ELE)
          IF(IND.EQ.1) THEN
             NDOFF=NDOF(1,N)+1
          ELSE
             NDOFF=NDOF(1,N)
          ENDIF
!
          X_BAR=(XTMP-X(N))/RAD
!
          I2=I+NEleAss
          SELECT CASE(NDOFF)
          CASE(1)
             HI_OVLP_TRUSS(II,1,I2)=(1.0+(XTMP-XX2)/LENGTH)*WEIGHT
             DHI_OVLP_TRUSS(II,1,I2)=WEIGHT/LENGTH
          CASE(2)
             HI_OVLP_TRUSS(II,1,I2)=1.0+(XTMP-XX2)/LENGTH
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,1,I2)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I2)=1.0/LENGTH
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,1,I2)*X_BAR+HI_OVLP_TRUSS(II,1,I2)/RAD
!
             HI_OVLP_TRUSS(II,1,I2)=HI_OVLP_TRUSS(II,1,I2)*WEIGHT
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,2,I2)*WEIGHT
             DHI_OVLP_TRUSS(II,1,I2)=DHI_OVLP_TRUSS(II,1,I2)*WEIGHT
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,2,I2)*WEIGHT
          CASE(3)
             HI_OVLP_TRUSS(II,1,I2)=1.0+(XTMP-XX2)/LENGTH
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,1,I2)*X_BAR
             HI_OVLP_TRUSS(II,3,I2)=HI_OVLP_TRUSS(II,2,I2)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I2)=1.0/LENGTH
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,1,I2)*X_BAR+HI_OVLP_TRUSS(II,1,I2)/RAD
             DHI_OVLP_TRUSS(II,3,I2)=DHI_OVLP_TRUSS(II,2,I2)*X_BAR+HI_OVLP_TRUSS(II,2,I2)/RAD
!
             HI_OVLP_TRUSS(II,1,I2)=HI_OVLP_TRUSS(II,1,I2)*WEIGHT
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,2,I2)*WEIGHT
             HI_OVLP_TRUSS(II,3,I2)=HI_OVLP_TRUSS(II,3,I2)*WEIGHT
             DHI_OVLP_TRUSS(II,1,I2)=DHI_OVLP_TRUSS(II,1,I2)*WEIGHT
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,2,I2)*WEIGHT
             DHI_OVLP_TRUSS(II,3,I2)=DHI_OVLP_TRUSS(II,3,I2)*WEIGHT
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GNERATEM/DHI_TRUSS!"
             STOP
          END SELECT
       ENDIF
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE DHI_TRUSS
!
    SUBROUTINE DHI_TRUSS2(II,XTMP,XX1,XX2,ELE1,ELE2,NEleAss,WEI,DWEI)
!{{{
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: II
    REAL(8), INTENT(IN) :: XTMP,XX1,XX2
    REAL(8) :: WEI(NEleAss,2),DWEI(NEleAss,2)
    INTEGER, INTENT(IN) :: ELE1,ELE2,NEleAss
!   IN-FUNCTION VARIABLES
    INTEGER :: I,J
    INTEGER :: ELE,N,IND,NDOFF
    REAL(8) :: RAD,LENGTH
    INTEGER :: I2
    REAL(8) :: X_BAR
    REAL(8) :: WEI_SUM(2),DWEI_SUM(2)
    REAL(8) :: H
!
    DO I=1,2
       WEI_SUM(I)=0.0D0
       DWEI_SUM(I)=0.0D0
       DO J=1,NEleAss
          WEI_SUM(I)=WEI_SUM(I)+WEI(J,I)
          DWEI_SUM(I)=DWEI_SUM(I)+DWEI(J,I)
       ENDDO
    ENDDO
!
    DO I=1,2
       DO J=1,NEleAss
          WEI(J,I)=WEI(J,I)/WEI_SUM(I)
          DWEI(J,I)=DWEI(J,I)/WEI_SUM(I)-WEI(J,I)/WEI_SUM(I)*DWEI_SUM(I)
       ENDDO
    ENDDO
!
    LENGTH=ABS(XX1-XX2)
    DO I=1,NEleAss
       ELE=LM(I,ELE1)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IND=NDMFS(2,ELE)
          RAD=2*ELEPAR(1,ELE)
!
          IF(IND.EQ.1) THEN
             NDOFF=NDOF(1,N)+1
          ELSE
             NDOFF=NDOF(1,N)
          ENDIF
!
          X_BAR=(XTMP-X(N))/RAD
!
          SELECT CASE(NDOFF)
          CASE(1)
             H=(1.0-(XTMP-XX1)/LENGTH)
             HI_OVLP_TRUSS(II,1,I)=H*WEI(I,1)
!
             DHI_OVLP_TRUSS(II,1,I)=-1.0D0/LENGTH*WEI(I,1)+H*DWEI(I,1)
          CASE(2)
             H=1.0-(XTMP-XX1)/LENGTH
             HI_OVLP_TRUSS(II,1,I)=H*WEI(I,1)
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,1,I)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I)=-1.0D0/LENGTH*WEI(I,1)+H*DWEI(I,1)
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,1,I)*X_BAR+HI_OVLP_TRUSS(II,1,I)/RAD
          CASE(3)
             H=1.0-(XTMP-XX1)/LENGTH
             HI_OVLP_TRUSS(II,1,I)=H*WEI(I,1)
             HI_OVLP_TRUSS(II,2,I)=HI_OVLP_TRUSS(II,1,I)*X_BAR
             HI_OVLP_TRUSS(II,3,I)=HI_OVLP_TRUSS(II,2,I)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I)=-1.0D0/LENGTH*WEI(I,1)+H*DWEI(I,1)
             DHI_OVLP_TRUSS(II,2,I)=DHI_OVLP_TRUSS(II,1,I)*X_BAR+HI_OVLP_TRUSS(II,1,I)/RAD
             DHI_OVLP_TRUSS(II,3,I)=DHI_OVLP_TRUSS(II,2,I)*X_BAR+HI_OVLP_TRUSS(II,2,I)/RAD         
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GNERATEM/DHI_TRUSS!"
             STOP
          END SELECT
       ENDIF
!
       ELE=LM(I,ELE2)
       IF(ELE.NE.0) THEN
          N=NDMFS(1,ELE)
          IND=NDMFS(2,ELE)
          RAD=2*ELEPAR(1,ELE)
          IF(IND.EQ.1) THEN
             NDOFF=NDOF(1,N)+1
          ELSE
             NDOFF=NDOF(1,N)
          ENDIF
!
          X_BAR=(XTMP-X(N))/RAD
!
          I2=I+NEleAss
          SELECT CASE(NDOFF)
          CASE(1)
             H=1.0+(XTMP-XX2)/LENGTH
             HI_OVLP_TRUSS(II,1,I2)=H*WEI(I,2)
             DHI_OVLP_TRUSS(II,1,I2)=1.0D0/LENGTH*WEI(I,2)+H*DWEI(I,2)
          CASE(2)
             H=1.0+(XTMP-XX2)/LENGTH
             HI_OVLP_TRUSS(II,1,I2)=H*WEI(I,2)
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,1,I2)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I2)=1.0D0/LENGTH*WEI(I,2)+H*DWEI(I,2)
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,1,I2)*X_BAR+HI_OVLP_TRUSS(II,1,I2)/RAD
          CASE(3)
             H=1.0D0+(XTMP-XX2)/LENGTH
             HI_OVLP_TRUSS(II,1,I2)=H*WEI(I,2)
             HI_OVLP_TRUSS(II,2,I2)=HI_OVLP_TRUSS(II,1,I2)*X_BAR
             HI_OVLP_TRUSS(II,3,I2)=HI_OVLP_TRUSS(II,2,I2)*X_BAR
!
             DHI_OVLP_TRUSS(II,1,I2)=1.0D0/LENGTH*WEI(I,2)+H*DWEI(I,2)
             DHI_OVLP_TRUSS(II,2,I2)=DHI_OVLP_TRUSS(II,1,I2)*X_BAR+HI_OVLP_TRUSS(II,1,I2)/RAD
             DHI_OVLP_TRUSS(II,3,I2)=DHI_OVLP_TRUSS(II,2,I2)*X_BAR+HI_OVLP_TRUSS(II,2,I2)/RAD
          CASE DEFAULT
             WRITE(*,*) "ERROR, /GNERATEM/DHI_TRUSS!"
             STOP
          END SELECT
       ENDIF
    ENDDO 
!
    RETURN
!}}}
    END SUBROUTINE DHI_TRUSS2
!
    SUBROUTINE MFS_LINE_G(ADBC,NEleAss,THIC)
!{{{
!--------------------------------------------------
!   GENERATE LOAD VECTOR
!      LOAD      VECTOR
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADBC,NEleAss
    REAL(8),INTENT(INOUT) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,JJ,KK,JJJ
    INTEGER :: NumBc
    INTEGER :: ALLocateStatus
!  
    INTEGER :: L1,L2,N1,N2,L,N
    REAL(8) :: XM,XL,YM,YL,XX,YY
    REAL(8) :: FLOAD(2)
    INTEGER :: FUNC(2)
    INTEGER :: INDBC,IND
    LOGICAL :: INDBCX,INDBCY,EXPRN1
    INTEGER :: NUM
    INTEGER :: NPOR,NPORMAX 
    REAL(8) :: DIV
    REAL(8) :: S(NPOR_MAX,2),HI(NPOR_MAX),LENGTH
    REAL(8) :: S1,S2,S3,DHIDX,DHIDY
    REAL(8) :: NX,NY
    REAL(8) :: LOAD(2)
!   READ 
    NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
    ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
    ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
    READ(ADBC) ((BCED(J,I),J=1,2),I=1,NumBc),((BCIND(J,I),J=1,4),I=1,NumBc),((BCPAR(J,I),J=1,2),I=1,NumBc)
!
    IF(PENALTY.EQ.0) THEN   
      DO I=1,NumBc
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) CYCLE
!
       L1=BCED(1,I)
       L2=BCED(2,I)
!
       N1=NDMFS(1,L1)
       N2=NDMFS(1,L2)
!
       LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
       LENGTH=SQRT(LENGTH)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       NDMFS(7,L1)=BCIND(1,I)
       NDMFS(7,L2)=BCIND(1,I)
!
       NDMFS(8,L1)=BCIND(2,I)
       NDMFS(8,L2)=BCIND(2,I)
!
       ELEPAR(5,L1)=-NX
       ELEPAR(6,L1)=-NY
       ELEPAR(5,L2)=-NX
       ELEPAR(6,L2)=-NY
      ENDDO
    ENDIF
!
    DO I=1,NumBc
!      CALCULATE LOAD VECTOR
       L1=BCED(1,I)
       L2=BCED(2,I)
!
       NDMFS(5,L1)=L2
       NDMFS(4,L2)=L1
!
       EXPRN1=BCIND(1,I).EQ.3
       IF(EXPRN1) CYCLE
!
       N1=NDMFS(1,L1)
       N2=NDMFS(1,L2)
!
       LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
       LENGTH=SQRT(LENGTH)
!
       XM=0.5*(X(N1)+X(N2))
       XL=0.5*(X(N1)-X(N2))
       YM=0.5*(Y(N1)+Y(N2))
       YL=0.5*(Y(N1)-Y(N2))
!
       FLOAD(1)=BCPAR(1,I)
       FLOAD(2)=BCPAR(2,I)
       FUNC(1)=BCIND(3,I)
       FUNC(2)=BCIND(4,I)
!
       NX=(Y(N2)-Y(N1))/LENGTH
       NY=-(X(N2)-X(N1))/LENGTH
!
       INDBCX=BCIND(1,I).EQ.1
       INDBCY=BCIND(2,I).EQ.1
       EXPRN1=INDBCX.OR.INDBCY 
!       
       DO K=1,2
          L=BCED(K,I)
          N=NDMFS(1,L)
          IF(NDMFS(2,L).EQ.1) THEN
            NPORMAX=NDOF(1,N)+1
          ELSE
            NPORMAX=NDOF(1,N)
          ENDIF
!
          DO J=1,2
             DO JJ=1,NPORMAX
                S(JJ,J)=0.0D0
             ENDDO
          ENDDO
!
          DO JJ=1,NITP_BC
             XX = XM+XL*XABSC_WEIG_BC(JJ,1)
             YY = YM+YL*XABSC_WEIG_BC(JJ,1)
             DIV= XABSC_WEIG_BC(JJ,2)
!
             HI = QUADS_MFS_HIHJ(NEleAss,L,NPORMAX,XX,YY)
             LOAD=F_LOAD(FLOAD,FUNC,NX,NY,XX,YY,N1,N2)
             DO KK=1,NPORMAX
                S(KK,1)= S(KK,1) + HI(KK)*DIV*LOAD(1)
                S(KK,2)= S(KK,2) + HI(KK)*DIV*LOAD(2)
             ENDDO
          ENDDO
!
          DO J=1,2
             INDBC=BCIND(J,I)
             SELECT CASE(INDBC)
             CASE(2)
                 IND=NDMFS(J+1,L)
                 DO JJ=1,NDOF(J,N)
                    NUM=ID(J,N)+JJ-1
                    IF(IND.EQ.1) THEN
                       KK=JJ+1
                    ELSE
                       KK=JJ
                    ENDIF
                    R(NUM)=R(NUM)+S(KK,J)*LENGTH/2.*THIC
                 ENDDO
             CASE(1) 
                 CONTINUE
             CASE DEFAULT
                 WRITE(*,*) "ERROR, IND.GT.1 IN SUBROUTINE MFS_LINE_G!"
                 STOP
             END SELECT
          ENDDO       
       ENDDO    
!******************************
! IMPOSITION OF DIRICHLET Bcs
!******************************
       IF(PENALTY.EQ.1) THEN
         IF(EXPRN1) THEN
          DO J=1,2
             DO K=J,2
                L1=BCED(J,I)
                L2=BCED(K,I)
!******************************
!  PENALTY METHOD
!******************************
                CALL MFS_LINE_STIFF_PENALTY(L1,L2,I,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
             ENDDO
          ENDDO
         ENDIF
       ENDIF
    ENDDO
!
    DEALLOCATE(BCED)
    DEALLOCATE(BCIND)
    DEALLOCATE(BCPAR) 
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_G
!
    SUBROUTINE MFS_LINE_STIFF_WEAK(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: LL
    INTEGER :: N1,N2
!   REAL(8) :: S(2,2)
    INTEGER :: NPOR1,NPOR2
!
    LOGICAL :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: HI(3,NPOR_MAX),HJ(3,NPOR_MAX)
    REAL(8) :: XTMP,YTMP
    REAL(8) :: WEIGHT
    REAL(8) :: CONST1
    REAL(8) :: KCS1,KCS2,KCS3
!
!    INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II).EQ.1
    IND1_Y  =  NDMFS(3,II).EQ.1
    IND2_X  =  NDMFS(2,JJ).EQ.1
    IND2_Y  =  NDMFS(3,JJ).EQ.1
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N1)+1
    ELSE
        NPOR1=NDOF(1,N1)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,N2)+1
    ELSE
        NPOR2=NDOF(1,N2)
    ENDIF
!
!    DO I=1,NPOR_MAX*NPOR_MAX
!       STIFF_IJ(1,I)=0.0D0
!       STIFF_IJ(2,I)=0.0D0
!       STIFF_IJ(3,I)=0.0D0
!       STIFF_IJ(4,I)=0.0D0
!    ENDDO
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJ(I,1)=0.0D0
       HIHJ(I,2)=0.0D0
       HIHJ(I,3)=0.0D0
       HIHJ(I,4)=0.0D0
    ENDDO
!
!    DO I=1,3
!       DO J=1,NPOR_MAX
!          HI(J,I)=0.0D0
!          HJ(J,I)=0.0D0
!       ENDDO
!    ENDDO
!
      DO I=1,NITP_IN
         XTMP = XM+XL*XABSC_WEIG_IN(I,1)
         YTMP = YM+YL*XABSC_WEIG_IN(I,1)
!
         CALL QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XTMP,YTMP,NPOR1,NPOR2,HI(1,1),HJ(1,1))
!
         WEIGHT=XABSC_WEIG_IN(I,2)
         DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL,1)=HIHJ(LL,1)+HI(1,J)*HJ(2,K)*WEIGHT
               HIHJ(LL,2)=HIHJ(LL,2)+HI(1,J)*HJ(3,K)*WEIGHT
               HIHJ(LL,3)=HIHJ(LL,3)+HI(2,J)*HJ(1,K)*WEIGHT
               HIHJ(LL,4)=HIHJ(LL,4)+HI(3,J)*HJ(1,K)*WEIGHT
            ENDDO
         ENDDO
      ENDDO
!
      CONST1=LENGTH/2.0*THIC
      KCS1=KCS(1)
      KCS2=KCS(2)
      KCS3=KCS(3)
      DO I=1,L
         STIFF_IJ(I,1)=(HIHJ(I,1)+HIHJ(I,3))*NX*KCS1+(HIHJ(I,2)+HIHJ(I,4))*NY*KCS3
         STIFF_IJ(I,2)=(HIHJ(I,2)*NX+HIHJ(I,3)*NY)*KCS2+(HIHJ(I,1)*NY+HIHJ(I,4)*NX)*KCS3
         STIFF_IJ(I,3)=(HIHJ(I,1)*NY+HIHJ(I,4)*NX)*KCS2+(HIHJ(I,2)*NX+HIHJ(I,3)*NY)*KCS3
         STIFF_IJ(I,4)=(HIHJ(I,2)+HIHJ(I,4))*NY*KCS1+(HIHJ(I,1)+HIHJ(I,3))*NX*KCS3
!
         STIFF_IJ(I,1)=STIFF_IJ(I,1)*CONST1
         STIFF_IJ(I,2)=STIFF_IJ(I,2)*CONST1
         STIFF_IJ(I,3)=STIFF_IJ(I,3)*CONST1
         STIFF_IJ(I,4)=STIFF_IJ(I,4)*CONST1
      ENDDO  
!
      CALL MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF_WEAK
!
    SUBROUTINE MFS_LINE_STIFF_PENALTY(II,JJ,KK,NEleAss,NX,NY,XM,XL,YM,YL,LENGTH,THIC)
!{{{
!--------------------------------------------------
!   GENERATE STIFFNESS MATRIX RESULT FROM BCs
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER :: II,JJ,KK,NEleAss
    REAL(8) :: NX,NY,XM,XL,YM,YL,LENGTH
    REAL(8) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: LL
    INTEGER :: N1,N2
!   REAL(8) :: S(2,2)
    INTEGER :: NPOR1,NPOR2
!
    LOGICAL :: IND1_X,IND1_Y,IND2_X,IND2_Y 
    REAL(8) :: HIHJ(NPOR_MAX*NPOR_MAX)
    REAL(8) :: HI(NPOR_MAX),HJ(NPOR_MAX)
    REAL(8) :: XTMP,YTMP
    REAL(8) :: WEIGHT
    REAL(8) :: CONST1
    REAL(8) :: KCS1,KCS2,KCS3
!
!    INTEGER :: KI,KJ
!
    N1=NDMFS(1,II)
    N2=NDMFS(1,JJ)
!
    IND1_X  =  NDMFS(2,II).EQ.1
    IND1_Y  =  NDMFS(3,II).EQ.1
    IND2_X  =  NDMFS(2,JJ).EQ.1
    IND2_Y  =  NDMFS(3,JJ).EQ.1
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N1)+1
    ELSE
        NPOR1=NDOF(1,N1)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,N2)+1
    ELSE
        NPOR2=NDOF(1,N2)
    ENDIF
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJ(I)=0.0D0
    ENDDO
!
      DO I=1,NITP_IN
         XTMP = XM+XL*XABSC_WEIG_IN(I,1)
         YTMP = YM+YL*XABSC_WEIG_IN(I,1)
         WEIGHT=XABSC_WEIG_IN(I,2)
!
         HI = QUADS_MFS_HIHJ(NEleAss,II,NPOR1,XTMP,YTMP)
         HJ = QUADS_MFS_HIHJ(NEleAss,JJ,NPOR2,XTMP,YTMP)
         DO J=1,NPOR1
            DO K=1,NPOR2
               LL=(J-1)*NPOR2+K
               HIHJ(LL)=HIHJ(LL)+HI(J)*HJ(K)*WEIGHT
            ENDDO
         ENDDO
      ENDDO
!
      CONST1=LENGTH/2.0*THIC*BETA
      DO I=1,L
         STIFF_IJ(I,1)=HIHJ(I)*CONST1
         STIFF_IJ(I,2)=0.0D0
         STIFF_IJ(I,3)=0.0D0
         STIFF_IJ(I,4)=HIHJ(I)*CONST1
      ENDDO  
!
      CALL MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!
    RETURN
!}}}
    END SUBROUTINE MFS_LINE_STIFF_PENALTY
!
    SUBROUTINE MFS_LINE_ADDBAN(II,JJ,KK,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,KK,NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L
    INTEGER :: JMIN,JMAX
    INTEGER :: N,M
    LOGICAL :: INDBCX,INDBCY
    INTEGER :: INDK,INDL,KK1,KK2
    INTEGER :: KI,KJ
    INTEGER :: LL
!
    INDBCX=(BCIND(1,KK).EQ.2)
    INDBCY=(BCIND(2,KK).EQ.2)
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       IF(INDBCX) THEN
          JMIN=2
       ELSE
          JMIN=1
       ENDIF
       IF(INDBCY) THEN
          JMAX=1
       ELSE
          JMAX=2
       ENDIF
!
       DO J=JMIN,JMAX
          K=1
          INDK=1
          IF((IND1_X).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y).AND.(J.EQ.2)) THEN
                 L=2
                 INDL=2
             ENDIF 
             DO WHILE(L.LE.NPOR2)
                KI=ID(I,N)+K-INDK
                KJ=ID(J,M)+L-INDL
                KK1=(I-1)*2+J
                KK2=(K-1)*NPOR2+L
                IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                   IF(KJ.LE.KI) THEN
                      LL=MAXA_COLSOL(KI)+ABS(KJ-KI)
                   ELSE
                      LL=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                   ENDIF
                   K_COLSOL(LL)=K_COLSOL(LL)-S(KK2,KK1)
                ENDIF
                L=L+1
             ENDDO
             K=K+1
          ENDDO
       ENDDO
    ENDDO
!
    RETURN    
!}}}
    END SUBROUTINE MFS_LINE_ADDBAN
!
    SUBROUTINE QUADS_MFS_STIFF_V2(II,JJ,NEleAss,THIC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE STIFFNESS MATRIX OF TRUSS_MFS                         .
! .          A(NWK)       :  GLOBAL STIFFNESS MATRIX                  .
! .          AMAXA(NEQ+1) :  ADDRESS MATRIX                           .
! .          N,M          :  NODE NUMBER                              .
! .          E            :  YOUNG'S MODULUS                          .
! .          PR           :  POisson ratio                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ
    INTEGER            :: NEleAss
    REAL(8)            :: THIC
!   ASSISTANT VARIABLES
    INTEGER            :: I,J,K,L,K1
    INTEGER            :: N,M
    INTEGER            :: NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    LOGICAL            :: EXPN1
    INTEGER            :: IND1,IND2
!
    REAL(8)            :: XM,YM,XL0,YL0,XL,YL,LENGTHN,LENGTHM,XTMP,WEIGHT
    REAL(8)            :: HIHJH(NPOR_MAX*NPOR_MAX,4),HIHJG(NPOR_MAX*NPOR_MAX,4)
    REAL(8)            :: KCS1,KCS2,KCS3,CONST1
    REAL(8)            :: RAD1J,RAD2J
!
    N       =  NDMFS(1,II)
    M       =  NDMFS(1,JJ)
!
    IND1    = NDMFS(2,II) 
    IND2    = NDMFS(2,JJ)
!
    IND1_X  = (IND1.EQ.1)
    IND1_Y  = (NDMFS(3,II).EQ.1)
    IND2_X  = (IND2.EQ.1)
    IND2_Y  = (NDMFS(3,JJ).EQ.1)
!
    EXPN1=((IND1.NE.0).AND.(IND2.NE.0))
!
    IF (IND1_X) THEN
        NPOR1=NDOF(1,N)+1
    ELSE
        NPOR1=NDOF(1,N)
    ENDIF
!
    IF (IND2_X) THEN
        NPOR2=NDOF(1,M)+1
    ELSE
        NPOR2=NDOF(1,M)
    ENDIF
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJH(I,1)=0.0D0
       HIHJH(I,2)=0.0D0
       HIHJH(I,3)=0.0D0
       HIHJH(I,4)=0.0D0
    ENDDO
!
    RAD1J=ELEPAR(3,JJ)
    RAD1J=1.0/RAD1J
    RAD2J=ELEPAR(4,JJ)
    RAD2J=1.0/RAD2J
!
    XL0=ELEPAR(3,II)
    XL =XL0/NIDIV_IN
    YL0=ELEPAR(4,II)
    YL =YL0/NIDIV_IN
!
    DO K1=1,NIDIV_IN
       XM = X(N)-XL0+XL*(2*K1-1)
       DO I=1,NITP_IN
          XTMP=XM+XL*XABSC_WEIG_IN(I,1)
          WEIGHT=XABSC_WEIG_IN(I,2)
!
          CALL QUADS_MFS_GAUSSG_V2(N,M,II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YL0,YL,HIHJG(1,1),EXPN1,L,RAD1J,RAD2J)
          DO J=1,L
             HIHJH(J,1)=HIHJH(J,1)+HIHJG(J,1)*WEIGHT
             HIHJH(J,2)=HIHJH(J,2)+HIHJG(J,2)*WEIGHT
             HIHJH(J,3)=HIHJH(J,3)+HIHJG(J,3)*WEIGHT
             HIHJH(J,4)=HIHJH(J,4)+HIHJG(J,4)*WEIGHT
          ENDDO
       ENDDO
    ENDDO
!
    KCS1=KCS(1)
    KCS2=KCS(2)
    KCS3=KCS(3)
    CONST1=THIC*XL
    DO I=1,L
       STIFF_IJ(I,1)=(KCS1*HIHJH(I,1)+KCS3*HIHJH(I,4))*CONST1
       STIFF_IJ(I,2)=(KCS2*HIHJH(I,2)+KCS3*HIHJH(I,3))*CONST1
       STIFF_IJ(I,3)=(KCS2*HIHJH(I,3)+KCS3*HIHJH(I,2))*CONST1
       STIFF_IJ(I,4)=(KCS1*HIHJH(I,4)+KCS3*HIHJH(I,1))*CONST1
    ENDDO
!
    CALL QUADS_MFS_ADDBAN_V2(II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,STIFF_IJ(1,1))
!   
    RETURN
!!}}}
    END SUBROUTINE QUADS_MFS_STIFF_V2
!
    SUBROUTINE COUP1_STIFF(II,JJ,NN,THIC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE STIFFNESS MATRIX OF COUPLING                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,NN
    REAL(8),INTENT(IN) :: THIC
!   ASSISTANT VARIABLES
    INTEGER :: N1,N2
    INTEGER :: NPOR1,NPOR2
    REAL(8) :: RM,RTMP
    INTEGER :: I,J,K,L
    REAL(8) :: WEIGHT
    REAL(8) :: HIHJH(NPOR_MAX*NPOR_MAX,4),HIHJG(NPOR_MAX*NPOR_MAX,4)
    REAL(8) :: KCS1,KCS2,KCS3,CONST1
!
    N1=NDMFS(II,NN)
    N2=NDMFS(JJ,NN)
!
    NPOR1=NDOF(1,N1)
    NPOR2=NDOF(1,N2)
!
    L=NPOR1*NPOR2
    DO I=1,L
       HIHJH(I,1)=0.0D0
       HIHJH(I,2)=0.0D0
       HIHJH(I,3)=0.0D0
       HIHJH(I,4)=0.0D0
    ENDDO
!
    DO K=1,NIDIV_IN
       RM = -1.0 + (2*K-1.0)/NIDIV_IN
       DO I=1,NITP_IN
          RTMP=RM+XABSC_WEIG_IN(I,1)/NIDIV_IN
          WEIGHT=XABSC_WEIG_IN(I,2)
          CALL COUP1_GAUSSG(II,JJ,NN,NPOR1,NPOR2,RTMP,L,HIHJG(1,1))
          DO J=1,L
             HIHJH(J,1)=HIHJH(J,1)+HIHJG(J,1)*WEIGHT
             HIHJH(J,2)=HIHJH(J,2)+HIHJG(J,2)*WEIGHT
             HIHJH(J,3)=HIHJH(J,3)+HIHJG(J,3)*WEIGHT
             HIHJH(J,4)=HIHJH(J,4)+HIHJG(J,4)*WEIGHT
          ENDDO
       ENDDO
    ENDDO
!
    KCS1=KCS(1)
    KCS2=KCS(2)
    KCS3=KCS(3)
    CONST1=THIC/NIDIV_IN
    DO I=1,L
       STIFF_IJ(I,1)=(KCS1*HIHJH(I,1)+KCS3*HIHJH(I,4))*CONST1
       STIFF_IJ(I,2)=(KCS2*HIHJH(I,2)+KCS3*HIHJH(I,3))*CONST1
       STIFF_IJ(I,3)=(KCS2*HIHJH(I,3)+KCS3*HIHJH(I,2))*CONST1
       STIFF_IJ(I,4)=(KCS1*HIHJH(I,4)+KCS3*HIHJH(I,1))*CONST1
    ENDDO
!
    CALL COUP1_ADDBAN(II,JJ,NN,NPOR1,NPOR2,STIFF_IJ(1,1))
!   
    RETURN
!}}}
    END SUBROUTINE COUP1_STIFF
!
    SUBROUTINE QUADS_MFS_GAUSSG_V2(N,M,II,JJ,NPOR1,NPOR2,NEleAss,XTMP,YL0,YL,HIHJG,EXPN1,LL,RAD1J,RAD2J)
!{{{
    IMPLICIT NONE
    INTEGER :: N,M,II,JJ,NPOR1,NPOR2,NEleAss
    REAL(8) :: XTMP,YL0,YL
    REAL(8) :: HIHJG(NPOR_MAX*NPOR_MAX,4)
    LOGICAL :: EXPN1
    INTEGER :: LL
    REAL(8) :: RAD1J,RAD2J
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L,K1
    INTEGER :: L1,L2,N1,N2
    REAL(8) :: YM,YTMP
    REAL(8) :: HI(3,NPOR1),HJ(3,NPOR2)
    REAL(8) :: X1,X2,Y1,Y2,XX1,YY1,WEIGHT
    LOGICAL :: REST1,REST2,REST3,LOGIC1,LOGIC2
    LOGICAL :: EXPRN1
    REAL(8) :: SX,SY
!
    DO I=1,4
       DO J=1,LL
        HIHJG(J,I)=0.0D0
       ENDDO
    ENDDO
!
    XX1=XTMP-X(N)
    SX=1.0-ABS(XTMP-X(M))*RAD1J
    EXPRN1=SX.LT.1E-12
    IF(EXPRN1) RETURN
    DO K1=1,NIDIV_IN
        YM=Y(N)-YL0+YL*(2*K1-1)
        DO I=1,NITP_IN     
           YTMP=YM+YL*XABSC_WEIG_IN(I,1)
           YY1=YTMP-Y(N)
!
           SY=1.0-ABS(YTMP-Y(M))*RAD2J
           EXPRN1=SY.LT.1E-12
           IF(EXPRN1) CYCLE               
!
           WEIGHT=XABSC_WEIG_IN(I,2)
           IF (EXPN1) THEN
              L1=NDMFS(4,II)
              L2=NDMFS(5,II)
              N1=NDMFS(1,L1)
              N2=NDMFS(1,L2)
! 
              X2=X(N2)-X(N)
              X1=X(N1)-X(N)
              Y2=Y(N2)-Y(N)
              Y1=Y(N1)-Y(N)
!  
              REST1=areClockwise(X1,X2,Y1,Y2)
              REST2=areClockwise(X2,XX1,Y2,YY1)
              REST3=areClockwise(XX1,X1,YY1,Y1)
              LOGIC1=(REST1.AND.(REST2.OR.REST3))
              LOGIC2=(REST2.AND.REST3)
              IF (LOGIC1.OR.LOGIC2) CYCLE
            ENDIF
            CALL QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XTMP,YTMP,NPOR1,NPOR2,HI(1,1),HJ(1,1))
!
            DO J=1,NPOR1
               DO K=1,NPOR2
                  L=(J-1)*NPOR2+K
                  HIHJG(L,1)=HIHJG(L,1)+HI(2,J)*HJ(2,K)*WEIGHT
                  HIHJG(L,2)=HIHJG(L,2)+HI(2,J)*HJ(3,K)*WEIGHT
                  HIHJG(L,3)=HIHJG(L,3)+HI(3,J)*HJ(2,K)*WEIGHT
                  HIHJG(L,4)=HIHJG(L,4)+HI(3,J)*HJ(3,K)*WEIGHT
               ENDDO             
            ENDDO  
!
        ENDDO
    ENDDO
!
    DO I=1,LL
       HIHJG(I,1)=HIHJG(I,1)*YL
       HIHJG(I,2)=HIHJG(I,2)*YL
       HIHJG(I,3)=HIHJG(I,3)*YL
       HIHJG(I,4)=HIHJG(I,4)*YL
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_GAUSSG_V2
!
    SUBROUTINE COUP1_GAUSSG(II,JJ,NN,NPOR1,NPOR2,RTMP,LL,HIHJG)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,NN,NPOR1,NPOR2,LL
    REAL(8),INTENT(IN) :: RTMP
    REAL(8) :: HIHJG(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L,K1
    REAL(8) :: DET_J
    REAL(8) :: SM,STMP,WEIGHT
    REAL(8) :: HI(3,NPOR1),HJ(3,NPOR2)
!
    DO I=1,4
       DO J=1,LL
        HIHJG(J,I)=0.0D0
       ENDDO
    ENDDO
!
    DO K1 = 1,NIDIV_IN
       SM =-1+(2*K1-1.0)/NIDIV_IN
       DO I=1,NITP_IN     
          STMP=SM+XABSC_WEIG_IN(I,1)/NIDIV_IN
          WEIGHT=XABSC_WEIG_IN(I,2)
          CALL COUP1_HIHJ(II,JJ,NN,RTMP,STMP,NPOR1,NPOR2,HI(1,1),HJ(1,1),DET_J)
          DO J=1,NPOR1
             DO K=1,NPOR2
                L=(J-1)*NPOR2+K
                HIHJG(L,1)=HIHJG(L,1)+HI(2,J)*HJ(2,K)*WEIGHT*DET_J
                HIHJG(L,2)=HIHJG(L,2)+HI(2,J)*HJ(3,K)*WEIGHT*DET_J
                HIHJG(L,3)=HIHJG(L,3)+HI(3,J)*HJ(2,K)*WEIGHT*DET_J
                HIHJG(L,4)=HIHJG(L,4)+HI(3,J)*HJ(3,K)*WEIGHT*DET_J
             ENDDO             
          ENDDO  
       ENDDO
    ENDDO
!
    DO I=1,LL
       HIHJG(I,1)=HIHJG(I,1)/NIDIV_IN
       HIHJG(I,2)=HIHJG(I,2)/NIDIV_IN
       HIHJG(I,3)=HIHJG(I,3)/NIDIV_IN
       HIHJG(I,4)=HIHJG(I,4)/NIDIV_IN
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COUP1_GAUSSG
!
    SUBROUTINE QUADS_MFS_ADDBAN_V2(II,JJ,NPOR1,NPOR2,IND1_X,IND1_Y,IND2_X,IND2_Y,S)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,NPOR1,NPOR2
    LOGICAL            :: IND1_X,IND1_Y,IND2_X,IND2_Y
    REAL(8)            :: S(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: N,M
!
    INTEGER :: I,J,K,L
    INTEGER :: KI,KJ
    INTEGER :: LL,LL1,KK1,KK2
    INTEGER :: INDK,INDL
!
    N         =  NDMFS(1,II)
    M         =  NDMFS(1,JJ)
!
    DO I=1,2
       DO J=1,2
          K=1
          INDK=1
          IF((IND1_X).AND.(I.EQ.1)) THEN
             K=2
             INDK=2
          ENDIF
          IF((IND1_Y).AND.(I.EQ.2)) THEN
             K=2
             INDK=2
          ENDIF          
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             IF((IND2_X).AND.(J.EQ.1)) THEN
                 L=2
                 INDL=2
             ENDIF
             IF((IND2_Y).AND.(J.EQ.2)) THEN
                 L=2
                 INDL=2
             ENDIF 
             DO WHILE(L.LE.NPOR2)
                KI=ID(I,N)+K-INDK
                KJ=ID(J,M)+L-INDL
                KK1=(I-1)*2+J
                KK2=(K-1)*NPOR2+L
                IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                   IF(KJ.LE.KI) THEN
!                     LL=AMAXA(KI)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KI)+ABS(KJ-KI)
                   ELSE
!                     LL=AMAXA(KJ)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                   ENDIF
!                  A(LL)=A(LL)+S(KK1,KK2)
                   K_COLSOL(LL1)=K_COLSOL(LL1)+S(KK2,KK1)
                   KMATRIX(LL1)=KMATRIX(LL1)+S(KK2,KK1)
                ENDIF
                L=L+1
             ENDDO
             K=K+1
          ENDDO
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE QUADS_MFS_ADDBAN_V2
!
    SUBROUTINE COUP1_ADDBAN(II,JJ,NN,NPOR1,NPOR2,S)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: II,JJ,NN,NPOR1,NPOR2
    REAL(8)            :: S(NPOR_MAX*NPOR_MAX,4)
!   ASSISTANT VARIABLES
    INTEGER :: N,M
    INTEGER :: I,J,K,L
    INTEGER :: KI,KJ
    INTEGER :: LL,LL1,KK1,KK2
    INTEGER :: INDK,INDL
!
    N         =  NDMFS(II,NN)
    M         =  NDMFS(JJ,NN)
    DO I=1,2
       DO J=1,2
          K=1
          INDK=1         
          DO WHILE(K.LE.NPOR1)
             L=1
             INDL=1
             DO WHILE(L.LE.NPOR2)
                KI=ID(I,N)+K-INDK
                KJ=ID(J,M)+L-INDL
                KK1=(I-1)*2+J
                KK2=(K-1)*NPOR2+L
                IF ((N.NE.M).OR.(KJ.GE.KI)) THEN
                   IF(KJ.LE.KI) THEN
!                     LL=AMAXA(KI)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KI)+ABS(KJ-KI)
                   ELSE
!                     LL=AMAXA(KJ)+ABS(KJ-KI)
                      LL1=MAXA_COLSOL(KJ)+ABS(KJ-KI)
                   ENDIF
!                  A(LL)=A(LL)+S(KK1,KK2)
                   K_COLSOL(LL1)=K_COLSOL(LL1)+S(KK2,KK1)
                   KMATRIX(LL1)=KMATRIX(LL1)+S(KK2,KK1)
                ENDIF
                L=L+1
             ENDDO
             K=K+1
          ENDDO
       ENDDO
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE COUP1_ADDBAN
!
!**************************************************   version 2
!
    SUBROUTINE ADDBAN (S,NND,LMM,ND,IND)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble upper triangular element stiffness into             .
! .   compacted global stiffness                                      .
! .                                                                   .
! .         A = GLOBAL STIFFNESS                                      .
! .         S = ELEMENT STIFFNESS                                     .
! .         ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS              .
! .                                                                   .
! .                   S(1)        S(2)        S(3)        . . .       .
! .         S   =                 S(ND+1)     S(ND+2)     . . .       .
! .                                           S(2*ND)     . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! .                   A(1)        A(3)        A(6)        . . .       .
! .         A   =                 A(2)        A(5)        . . .       .
! .                                           A(4)        . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ND,NND
      REAL(8)    :: S(NND)
      INTEGER    :: LMM(ND)
      INTEGER    :: NDI, I,II, MI, KS, J, JJ, IJ, KK, KSS
      INTEGER    :: IND
!
      NDI=0
      DO I=1,ND
         II=LMM(I)
         IF (II .GT. 0) THEN
            MI=MAXA_COLSOL(II)
            KS=I
            DO J=1,ND
               JJ=LMM(J)
               IF (JJ .GT. 0) THEN
                  IJ=II - JJ
                  IF (IJ .GE. 0) THEN
                     KK=MI + IJ
                     KSS=KS
                     IF (J.GE.I) KSS=J + NDI
                     IF(IND.EQ.1) THEN
                        K_COLSOL(KK)=K_COLSOL(KK) + S(KSS)
                        KMATRIX(KK)=KMATRIX(KK)+S(KSS)
                     ELSE
                        M_COLSOL(KK)=M_COLSOL(KK) + S(KSS)
                     ENDIF
                  END IF
               END IF
               KS=KS + ND - J
            END DO
         END IF
         NDI=NDI + ND - I
      END DO
!
      RETURN
!}}}
    END SUBROUTINE ADDBAN
!
    SUBROUTINE OUTPUTKM(ADK,ADM)
!{{{
!--------------------------------------------------
!   OUTPUT STIFFNESS AND MASS MATRIX
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ADK,ADM
    INTEGER :: I,J,K,L
    REAL(8) :: KI(NEQ)
!--- ASSISTANT VARIABLES
    INTEGER :: II,JJ
!
    KI=0.0D0
!
!    DO I=1,NEQ
!       DO L=I,NEQ
!          JJ=MHT(L)+1
!          J=L-I+1
!          IF(JJ.GE.J) THEN
!            II=MAXA_COLSOL(L)+J-1
!            KI(L)=K_COLSOL(II)
!          ENDIF
!       ENDDO
!       WRITE(ADK,1000) (KI(K),K=1,NEQ)
!       KI=0.0D0
!    ENDDO
!
    DO I=1,NEQ
       JJ=MHT(I)+1
       DO J=1,JJ
          II=MAXA_COLSOL(I)+J-1
          KI(J)=K_COLSOL(II)
       ENDDO
       WRITE(ADK,1000) (KI(K),K=1,JJ)
    ENDDO    
!
    RETURN
!
1000 FORMAT(*(ES20.13,1X))
!}}}
    END SUBROUTINE OUTPUTKM
!      
    SUBROUTINE ADDBAN_ZLB(ND,S,LMM,IND)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To assemble upper triangular element stiffness into             .
! .   compacted global stiffness                                      .
! .                                                                   .
! .         A = GLOBAL STIFFNESS                                      .
! .         S = ELEMENT STIFFNESS                                     .
! .         ND = DEGREES OF FREEDOM IN ELEMENT STIFFNESS              .
! .                                                                   .
! .                   S(1)        S(2)        S(3)        . . .       .
! .         S   =                 S(ND+1)     S(ND+2)     . . .       .
! .                                           S(2*ND)     . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! .                   A(1)        A(3)        A(6)        . . .       .
! .         A   =                 A(2)        A(5)        . . .       .
! .                                           A(4)        . . .       .
! .                                                       . . .       .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: ND
   REAL(8),DIMENSION(ND,ND),INTENT(IN) :: S
   INTEGER,DIMENSION(ND),INTENT(IN) :: LMM
   INTEGER :: IND
!
   INTEGER :: I,J,L,M,K
!
   DO I=1,ND
      L=LMM(I)
      IF(L.GT.0) THEN
        DO J=1,ND
           M=LMM(J)
           IF(M.GT.0) THEN
              K=M-L
              IF(K.GE.0) THEN
                K=MAXA_COLSOL(M)+K
                IF(IND.EQ.1) THEN
                   K_COLSOL(K)=K_COLSOL(K)+S(I,J)
                   KMATRIX(K)=KMATRIX(K)+S(I,J)
                ELSE
                   M_COLSOL(K)=M_COLSOL(K)+S(I,J)
                ENDIF                
              END IF
           ENDIF
        ENDDO
      ENDIF
   ENDDO
!
   RETURN
!}}}
    END SUBROUTINE ADDBAN_ZLB
!
end module GENERATEM
!{{{
