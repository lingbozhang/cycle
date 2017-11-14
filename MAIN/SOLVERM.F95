module SOLVERM
    USE WAREHOUSE
    implicit none
    PUBLIC  :: COLHT
!    PUBLIC  :: COLHT_COUP
    PUBLIC  :: COLHT_MFS
    PUBLIC  :: COLSOL_PRE
    PUBLIC  :: COLSOL
    PUBLIC  :: COL_ALLOCATE
!    PRIVATE :: BUILDKM
!
    CONTAINS
!**************************************************
!
!   SUBROUTINES FOR COLSOL
!
!**************************************************
    SUBROUTINE COLSOL_PRE()
!{{{
!--------------------------------------------------                                 
!     SOLVER: COLSOL
!     OUTPUT:
!       MAXA_COLSOL     --   ADDRESSES OF DIAGONAL ELEMENTS
!       NWK_COLSOL
!       K_COLSOL
!       M_COLSOL
!--------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I
!
      DO I=1,NEQ
         DISP(I)=R(I)
      ENDDO
!
!      CALL BUILDKM(KE,CMAXA,CMAXA_COLSOL,CK_COLSOL)
!
      RETURN
!}}}
    END SUBROUTINE COLSOL_PRE
!
    SUBROUTINE COLHT(ND,LMM)
!{{{
!--------------------------------------------------                                 
!    To calculate column heights                                     
!--------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ND
      INTEGER,INTENT(IN) :: LMM(ND)
      INTEGER :: I,LS, II, ME
!
      LS=100000
      DO I=1,ND
         IF (LMM(I) .NE. 0) THEN
            IF (LMM(I)-LS .LT. 0) LS=LMM(I)
         END IF
      END DO
!
      DO I=1,ND
         II=LMM(I)
         IF (II.NE.0) THEN
            ME=II - LS
            IF (ME.GT.MHT(II)) MHT(II)=ME
         END IF
      END DO
! 
      RETURN
!}}}
    END SUBROUTINE COLHT
!
!    SUBROUTINE COLHT_COUP(ND,LMM)
!!{{{
!!--------------------------------------------------                                 
!!    To calculate column heights                                     
!!--------------------------------------------------
!      IMPLICIT NONE
!      INTEGER,INTENT(IN) :: ND
!      INTEGER,INTENT(IN) :: LMM(ND)
!      INTEGER :: I,LS, II, ME
!      INTEGER :: J,RANG,KK
!!
!      LS=100000
!      DO I=1,ND
!         IF (LMM(I) .NE. 0) THEN
!            IF (LMM(I)-LS .LT. 0) LS=LMM(I)
!         END IF
!      END DO
!!
!      DO I=1,4
!         II=LMM(2*I-1)
!         IF(II.NE.0) THEN
!            RANG=2*(LMM(2*I)-II)
!            DO J=1,RANG
!               KK=II+J-1
!               ME=KK-LS
!               IF (ME.GT.MHT(KK)) MHT(KK)=ME
!            ENDDO
!         ENDIF
!      ENDDO
!!
!      RETURN
!!}}}
!    END SUBROUTINE COLHT_COUP
!
    SUBROUTINE COLHT_MFS(N,M)
!{{{
!--------------------------------------------------                                 
!    To calculate column heights                                     
!--------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: N,M
      INTEGER :: I,J,LS, II, ME
!      
      LS=100000
      DO I=1,3
         II=ID(I,N)
         IF (II.NE.0) THEN
            IF(II-LS.LT.0) LS=II
         ENDIF
!
         II=ID(I,M)
         IF (II.NE.0) THEN
            IF(II-LS.LT.0) LS=II
         ENDIF
      END DO
!
      DO I=1,3
         II=ID(I,N)
         IF(II.NE.0) THEN
            DO J=1,NDOF(I,N)
               II=ID(I,N)+J-1
               ME=II-LS
               IF (ME.GT.MHT(II)) MHT(II)=ME
            ENDDO
         ENDIF
!
         II=ID(I,M)
         IF(II.NE.0) THEN
            DO J=1,NDOF(I,M)
               II=ID(I,M)+J-1
               ME=II-LS
               IF(ME.GT.MHT(II)) MHT(II)=ME
            ENDDO
         ENDIF
      ENDDO
! 
      RETURN
!}}}
    END SUBROUTINE COLHT_MFS
!
!    SUBROUTINE BUILDKM(KE,CMAXA,CMAXA_COLSOL,CK_COLSOL)
!!{{{
!!--------------------------------------------------
!!   CONSTRUCT STIFFNESS AND MASS MATRIX FOR COLSOL
!!--------------------------------------------------
!    IMPLICIT NONE
!    REAL(8),DIMENSION(NWK),INTENT(IN) :: KE
!!
!    INTEGER,DIMENSION(NEQ+1),INTENT(IN)           :: CMAXA,CMAXA_COLSOL
!    REAL(8),DIMENSION(NWK_COLSOL),INTENT(INOUT)   :: CK_COLSOL
!!   ASSISTANT VARIABLES
!    INTEGER :: I,J
!    INTEGER :: Num1,Num2
!!
!    DO I=1,NEQ 
!       Num1=CMAXA_COLSOL(I)
!       Num2=CMAXA(I)
!       DO J=1,MHT(I)+1
!          CK_COLSOL(Num1)=KE(Num2)
!          Num1=Num1+1
!          Num2=Num2+1
!       END DO
!    ENDDO
!!
!    RETURN
!!}}}
!    END SUBROUTINE BUILDKM
!
    SUBROUTINE COLSOL(A,V,MAXAA,NN,NWKK,NNM,KKK)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To solve finite element static equilibrium equations in         .
! .   core, using compacted storage and column reduction scheme       .
! .                                                                   .
! .  - - Input variables - -                                          .
! .        A(NWK)    = Stiffness matrix stored in compacted form      .
! .        V(NN)     = Right-hand-side load vector                    .
! .        MAXA(NNM) = Vector containing addresses of diagonal        .
! .                    elements of stiffness matrix in a              .
! .        NN        = Number of equations                            .
! .        NWK       = Number of elements below skyline of matrix     .
! .        NNM       = NN + 1                                         .
! .        KKK       = Input flag                                     .
! .            EQ. 1   Triangularization of stiffness matrix          .
! .            EQ. 2   Reduction and back-substitution of load vector .
! .                                                                   .
! .  - - OUTPUT - -                                                   .
! .        A(NWK)    = D and L - Factors of stiffness matrix          .
! .        V(NN)     = Displacement vector                            .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
      IMPLICIT NONE
      INTEGER(8) :: MAXAA(1),NWKK
      INTEGER    :: NN,NNM,KKK
      REAL(8)    :: A(NWKK),V(1),C,B
      INTEGER    :: N,K,KN,KL,KU,KH,IC,KLT,KI,J,ND,KK,L
      INTEGER    :: II
      INTEGER    :: MIN0
!
      II=1
!
!     Perform L*D*L(T) factorization of stiffness matrix
!
      IF (KKK-2) 40,150,150

   40 DO N=1,NN
         KN=MAXAA(N)
         KL=KN + 1
         KU=MAXAA(N+1) - 1
         KH=KU - KL
         
         IF (KH) 110,90,50
         
   50    K=N - KH
         IC=0
         KLT=KU
         DO J=1,KH
            IC=IC + 1
            KLT=KLT - 1
            KI=MAXAA(K)
            ND=MAXAA(K+1) - KI - 1
            IF (ND .GT. 0) THEN
               KK=MIN0(IC,ND)
               C=0.
               DO L=1,KK
                  C=C + A(KI+L)*A(KLT+L)
               END DO
               A(KLT)=A(KLT) - C
            END IF
            K=K + 1
         END DO
         
   90    K=N
         B=0.
         DO KK=KL,KU
            K=K - 1
            KI=MAXAA(K)
            C=A(KK)/A(KI)
            B=B + C*A(KK)
            A(KK)=C
         END DO
         A(KN)=A(KN) - B
         
  110    IF ((A(KN) .LE. 0).AND.(II.EQ.1)) THEN
            WRITE (*,2000) N,A(KN)
            II=2
            CONTINUE
         END IF
      END DO   
      RETURN
!
!     REDUCE RIGHT-HAND-SIDE LOAD VECTOR
!
  150 DO N=1,NN
         KL=MAXAA(N) + 1
         KU=MAXAA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            C=0.
            DO KK=KL,KU
               K=K - 1
               C=C + A(KK)*V(K)
            END DO
            V(N)=V(N) - C
         END IF
      END DO
!
!     BACK-SUBSTITUTE
!
      DO N=1,NN
         K=MAXAA(N)
         V(N)=V(N)/A(K)
      END DO

      IF (NN.EQ.1) RETURN

      N=NN
      DO L=2,NN
         KL=MAXAA(N) + 1
         KU=MAXAA(N+1) - 1
         IF (KU-KL .GE. 0) THEN
            K=N
            DO KK=KL,KU
               K=K - 1
               V(K)=V(K) - A(KK)*V(N)
            END DO
         END IF
         N=N - 1
      END DO
!
      RETURN
!
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,&
                ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,&
                ' PIVOT = ',E20.12 )
!}}}
    END SUBROUTINE COLSOL
!
    SUBROUTINE COL_ALLOCATE()           
!**************************************************
!   WAREHOUSE FOR SUBROUTINE SOLVER
!**************************************************
! Build DISP(NEQ), MAXA_COLSOL, K_COLSOL
!{{{
    IMPLICIT NONE
!
    INTEGER :: ALLocateStatus
    INTEGER :: I,NN
!   ASSISTANT VARIABLES
    INTEGER :: J,K
    INTEGER :: Num1
!
!---- ALLOCATE VARIABLES
      ALLOCATE (DISP(NEQ+NCQ),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
      DISP=0.0D0
!
      ALLOCATE (MAXA_COLSOL(NEQ+NCQ+1),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
      ALLOCATE (MAXA_ND(NUMNP+1),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!
!     Clear array maxa
!
      MAXA_COLSOL=0
      MAXA_ND=0
!
      MAXA_COLSOL(1)=1
      MAXA_COLSOL(2)=2
      IF ((NEQ+NCQ).GE.2) THEN
         DO I=2,NEQ+NCQ
           MAXA_COLSOL(I+1)=MAXA_COLSOL(I) + MHT(I) + 1
         END DO
      ENDIF
      NWK_COLSOL=MAXA_COLSOL(NEQ+NCQ+1) - MAXA_COLSOL(1)
!
      MAXA_ND(1)=1
      MAXA_ND(2)=2
      IF(NUMNP.GE.2) THEN
         DO I=2,NUMNP
            MAXA_ND(I+1)=MAXA_ND(I)+MHT_ND(I)+1
         ENDDO
      ENDIF
      NWK_ND=MAXA_ND(NUMNP+1)-MAXA_ND(1)
!
      WRITE(*,*) "NUNMBER OF EQUATIONS      ", NEQ 
      WRITE(*,*) "STORAGE ENTRIES IN KMATRIX    ", NWK_COLSOL
      WRITE(*,*) "AVERAGE BANDWIDTH is (including zeros) ",NWK_COLSOL/NEQ
!---- ALLOCATE VARIABLES
      ALLOCATE (K_COLSOL(NWK_COLSOL),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
      ALLOCATE (BAND_FLAG(NWK_ND),STAT=ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
      ALLOCATE (KMATRIX(NWK_COLSOL),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!      ALLOCATE (M_COLSOL(NWK_COLSOL),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
!      ALLOCATE (C_COLSOL(NWK_COLSOL),STAT = ALLocateStatus)
!      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
      ALLOCATE (HESTRAIN(NEQ),STAT = ALLocateStatus)
      IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"      
! ----- INITIALIZATION -------
      HESTRAIN=0.0D0
      K_COLSOL=0.0D0
      BAND_FLAG=.FALSE.
      KMATRIX=0.0D0
!      M_COLSOL=0.0D0
!      C_COLSOL=0.0D0
! ----------------------------
    RETURN
!}}} 
    END SUBROUTINE COL_ALLOCATE
!**************************************************
!
!   SUBROUTINES FOR COLSOL
!
!**************************************************
end module SOLVERM
