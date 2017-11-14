! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Operation module  for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
       MODULE OPERATIONM
!--------------------------------------------------
!       USER DEFINED OPERATIONS
!--------------------------------------------------
        USE WAREHOUSE
        IMPLICIT NONE
! SUBROUTINES
        PUBLIC SHEPARD_2D
        PUBLIC SHEPARD_COUP1
        PUBLIC SHEPARD_TRUSS
        PUBLIC QUADS_MFS_HIHJ_V2
        PUBLIC COUP1_HIHJ
        PUBLIC OUTPUT_2DMFS
        PUBLIC OUTPUT_COUP1
        PUBLIC OUTPUT_2DMFS_SHEPARD
        PUBLIC OUTPUT_COUP1_SHEPARD
        PUBLIC G_KCS
        PUBLIC HALF_BANDWIDTH
! FUNCTIONS
        PUBLIC MATRIXB_VECTOR
        PUBLIC QUADS_MFS_HIHJ 
        PUBLIC QUADS_OVLP_HIHJ
! * * * * * * * * * * * * * * * * * * 
        PUBLIC TRUSS_MFS_HIHJ
        PUBLIC TRUSS_MFS_HIHJ_V2
        PUBLIC areClockwise
!
        CONTAINS
!
        FUNCTION MATRIXB_VECTOR(A,NMAXA,NMHT,V,NWKK,NDOFF) RESULT(H)
!{{{
          IMPLICIT NONE
          INTEGER*8 :: NWKK
          INTEGER :: NDOFF
          INTEGER*8 :: NMAXA(NDOFF+1)
          INTEGER :: NMHT(NDOFF)
          REAL(8) :: A(NWKK),V(NDOFF)
          REAL(8) :: H(NDOFF)
!         ASSISTANT VARIABLES
          INTEGER :: I,J
          INTEGER :: JJ
!
          DO I=1,NDOFF
             H(I)=0.0D0
          ENDDO
!
          DO I=1,NDOFF
             JJ=NMAXA(I)
             H(I)=H(I)+A(JJ)*V(I)
             DO J=1,NMHT(I)
                JJ=JJ+1
                H(I)=H(I)+A(JJ)*V(I-J)
                H(I-J)=H(I-J)+A(JJ)*V(I)
             ENDDO
          ENDDO
!}}}
        END FUNCTION MATRIXB_VECTOR
!
        FUNCTION areClockwise(X1,X2,Y1,Y2) RESULT(REST)
!{{{
        IMPLICIT NONE
        REAL(8) :: X1,X2,Y1,Y2
        LOGICAL :: REST
!
        REST= -X1*Y2+Y1*X2>=-1e-12  
!}}}
        END FUNCTION areClockwise
!
      SUBROUTINE SHEPARD_2D(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND_BLEND)
!{{{
      IMPLICIT NONE
      INTEGER :: NEleAss, II,JJ
      REAL(8) :: XX,YY
      REAL(8) :: W(5),DWDX(5),DWDY(5)
      INTEGER :: IND_BLEND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: M
      INTEGER :: INDSP
      REAL(8) :: RAD1,RAD2
      REAL(8) :: S,SS(4),SX,SY,S1
      REAL(8) :: CONSTX,CONSTY
      REAL(8) :: CONSTXD,CONSTYD
      REAL(8) :: WW,DWDXX,DWDYY
      LOGICAL :: EXPRN1
      REAL(8) :: W3,DWDX3,DWDY3
      INTEGER :: INDSPI,INDSPJ
!     
      DO I=1,5
         W(I)=0.0D0
         DWDX(I)=0.0D0
         DWDY(I)=0.0D0
      ENDDO
!
      W3=0.0D0
      DWDX3=0.0D0
      DWDY3=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(6,J)
             RAD1=ELEPAR(3,J)
             RAD1=1/RAD1
!
             IF(INDSP.GE.3) THEN
                RAD2=ELEPAR(4,J)
                RAD2=1.0/RAD2
                CONSTX=(XX-X(M))
!
                SX=1.0-ABS(CONSTX)*RAD1
                EXPRN1=SX.LT.1E-12
                IF(EXPRN1) CYCLE
                CONSTY=(YY-Y(M))
                SY=1.0-ABS(CONSTY)*RAD2
                EXPRN1=SY.LT.1E-12
                IF(EXPRN1) CYCLE
             ELSE
                CONSTX=(XX-X(M))
                CONSTY=(YY-Y(M))
!
                S1=SQRT(CONSTX*CONSTX+CONSTY*CONSTY)
                S=S1*RAD1
                EXPRN1=(1.0-S).LT.1E-12
                IF(EXPRN1) CYCLE
                IF(S.LT.1E-12) THEN
                   CONSTXD=0.0D0
                   CONSTYD=0.0D0
                ELSE
                   CONSTXD=CONSTX/S1*RAD1
                   CONSTYD=CONSTY/S1*RAD1
                ENDIF
             ENDIF
!
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  WW=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDXX=SS(4)*CONSTXD
                  DWDYY=SS(4)*CONSTYD
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
!
                     WW=2/3.0-4*SS(2)+4*SS(3)
                     DWDXX=SS(4)*CONSTXD
                     DWDYY=SS(4)*CONSTYD
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
!
                     WW=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDXX=SS(4)*CONSTXD
                     DWDYY=SS(4)*CONSTYD
                  ENDIF
                CASE(3)
                     WW=SX*SY
                     DWDXX=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDYY=-SX*RAD2*SIGN(1.0D0,CONSTY)
                CASE(4)
                     IND_BLEND=1
!
                     WW=SX*SY
                     DWDXX=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDYY=-SX*RAD2*SIGN(1.0D0,CONSTY)
!
                     W3=W3+WW
                     DWDX3=DWDX3+DWDXX
                     DWDY3=DWDY3+DWDYY
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT                                 
!
             IF(II.EQ.J) THEN
                W(1)=WW
                DWDX(1)=DWDXX
                DWDY(1)=DWDYY
             ENDIF
!
             IF(JJ.EQ.J) THEN
                W(2)=WW
                DWDX(2)=DWDXX
                DWDY(2)=DWDYY
             ENDIF
!
             W(3)=W(3)+WW
             DWDX(3)=DWDX(3)+DWDXX
             DWDY(3)=DWDY(3)+DWDYY
         ENDIF
      ENDDO
!
      IF(IND_BLEND.EQ.0) RETURN
!
      INDSPI=NDMFS(6,II)
      INDSPJ=NDMFS(6,JJ)
!
      IF(INDSPI.EQ.4) THEN
         W(4)=W(1)*W3
         DWDX(4)=DWDX(1)*W3+W(1)*DWDX3
         DWDY(4)=DWDY(1)*W3+W(1)*DWDY3
      ENDIF
      IF(INDSPI.EQ.3) THEN
         W(4)=W(1)+W3*W(1)
         DWDX(4)=DWDX(1)+DWDX3*W(1)+W3*DWDX(1)
         DWDY(4)=DWDY(1)+DWDY3*W(1)+W3*DWDY(1)
      ENDIF
!
      SX=1.0-ABS(CONSTX)*RAD1
      SY=1.0-ABS(CONSTY)*RAD2
      IF(INDSPJ.EQ.4) THEN
         W(5)=W(2)*W3
         DWDX(5)=DWDX(2)*W3+W(2)*DWDX3
         DWDY(5)=DWDY(2)*W3+W(2)*DWDY3
      ENDIF
      IF(INDSPJ.EQ.3) THEN   
         W(5)=W(2)+W3*W(2)
         DWDX(5)=DWDX(2)+DWDX3*W(2)+W3*DWDX(2)
         DWDY(5)=DWDY(2)+DWDY3*W(2)+W3*DWDY(2)
      ENDIF
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD_2D
!
      SUBROUTINE SHEPARD_COUP1(II,JJ,NN,RR,SS,XX,YY,W,DWDX,DWDY,IND_BLEND,DET_J)
!{{{
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: II,JJ,NN
      REAL(8),INTENT(IN) :: RR,SS
      REAL(8) :: W(5),DWDX(5),DWDY(5)
      INTEGER :: IND_BLEND
      REAL(8),INTENT(OUT) :: XX,YY,DET_J
!     ASSISTANT VARIABLES
      INTEGER :: I
      INTEGER :: IND,N,N1,N2,N3,N4
      REAL(8) :: RR1,RR2,SS1,SS2
      REAL(8) :: DXDR,DYDR,DXDS,DYDS
      REAL(8) :: RAD1,S,S1,SSS(4)
      REAL(8) :: CONSTX,CONSTY,CONSTXD,CONSTYD
      LOGICAL :: EXPRN1
      REAL(8) :: WW,DWDXX,DWDYY
      REAL(8) :: W3,DWDX3,DWDY3
!     
      DO I=1,5
         W(I)=0.0D0
         DWDX(I)=0.0D0
         DWDY(I)=0.0D0
      ENDDO
!
      N1=NDMFS(1,NN)
      N2=NDMFS(2,NN)
      N3=NDMFS(3,NN)
      N4=NDMFS(4,NN)
!     GLOBAL COORDINATES
      RR1=1+RR
      RR2=1-RR
      SS1=1+SS
      SS2=1-SS
!
      XX=RR1*SS1*X(N1)/4.0
      XX=XX+RR2*SS1*X(N2)/4.0
      XX=XX+RR2*SS2*X(N3)/4.0
      XX=XX+RR1*SS2*X(N4)/4.0
!
      YY=RR1*SS1*Y(N1)/4.0
      YY=YY+RR2*SS1*Y(N2)/4.0
      YY=YY+RR2*SS2*Y(N3)/4.0
      YY=YY+RR1*SS2*Y(N4)/4.0
!
      DXDR=(SS1*X(N1)-SS1*X(N2)-SS2*X(N3)+SS2*X(N4))/4.0
      DXDS=(RR1*X(N1)+RR2*X(N2)-RR2*X(N3)-RR1*X(N4))/4.0
      DYDR=(SS1*Y(N1)-SS1*Y(N2)-SS2*Y(N3)+SS2*Y(N4))/4.0
      DYDS=(RR1*Y(N1)+RR2*Y(N2)-RR2*Y(N3)-RR1*Y(N4))/4.0
!
      DET_J=DXDR*DYDS-DYDR*DXDS

      W3=0.0D0
      DWDX3=0.0D0
      DWDY3=0.0D0
!
      N1=NDMFS(II,NN)
      N2=NDMFS(JJ,NN)
!
      DO I=1,4
         N=NDMFS(I,NN)
         IND=NDMFS(I+4,NN)
         IF(IND.EQ.1) THEN
            RAD1=ELEPAR(I,NN)
            RAD1=1.0/RAD1
!
            CONSTX=(XX-X(N))
            CONSTY=(YY-Y(N))
!
            S1=SQRT(CONSTX*CONSTX+CONSTY*CONSTY)
            S=S1*RAD1
            EXPRN1=(1.0-S).LT.1E-12
            IF(EXPRN1) CYCLE
            IF(S.LT.1E-12) THEN
              CONSTXD=0.0D0
              CONSTYD=0.0D0
            ELSE
              CONSTXD=CONSTX/S1*RAD1
              CONSTYD=CONSTY/S1*RAD1
            ENDIF
         ENDIF
!
         SELECT CASE(IND)
         CASE(1)
             SSS(1)=S*S
             SSS(2)=SSS(1)*S
             SSS(3)=SSS(2)*S
             SSS(4)=-12*S+24*SSS(1)-12*SSS(2)
!
             WW=1-6*SSS(1)+8*SSS(2)-3*SSS(3)
             DWDXX=SSS(4)*CONSTXD
             DWDYY=SSS(4)*CONSTYD
         CASE(3)
             SELECT CASE(I)
             CASE(1)
                WW=RR1*SS1/4.0
                DWDXX=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
                DWDYY=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
             CASE(2)
                WW=RR2*SS1/4.0
                DWDXX=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
                DWDYY=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
             CASE(3)
                WW=RR2*SS2/4.0
                DWDXX=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
                DWDYY=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
             CASE(4)
                WW=RR1*SS2/4.0
                DWDXX=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
                DWDYY=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
             CASE DEFAULT
                WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
                STOP
             END SELECT
         CASE(4)
             IND_BLEND=1
             SELECT CASE(I)
             CASE(1)
                WW=RR1*SS1/4.0
                DWDXX=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
                DWDYY=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
             CASE(2)
                WW=RR2*SS1/4.0
                DWDXX=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
                DWDYY=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
             CASE(3)
                WW=RR2*SS2/4.0
                DWDXX=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
                DWDYY=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
             CASE(4)
                WW=RR1*SS2/4.0
                DWDXX=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
                DWDYY=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
             CASE DEFAULT
                WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
                STOP
             END SELECT
!
             W3=W3+WW
             DWDX3=DWDX3+DWDXX
             DWDY3=DWDY3+DWDYY
         CASE DEFAULT
             WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
             WRITE(*,*) "IND",IND
             STOP
         END SELECT                                 
!
         IF(N1.EQ.N) THEN
            W(1)=WW
            DWDX(1)=DWDXX
            DWDY(1)=DWDYY
         ENDIF
         IF(N2.EQ.N) THEN
            W(2)=WW
            DWDX(2)=DWDXX
            DWDY(2)=DWDYY
         ENDIF
!
         W(3)=W(3)+WW
         DWDX(3)=DWDX(3)+DWDXX
         DWDY(3)=DWDY(3)+DWDYY
      ENDDO
!
      IF(IND_BLEND.EQ.0) RETURN
!
      IND=NDMFS(II+4,NN)
      SELECT CASE(IND)
      CASE(4)
         W(4)=W(1)*W3
         DWDX(4)=DWDX(1)*W3+W(1)*DWDX3
         DWDY(4)=DWDY(1)*W3+W(1)*DWDY3
      CASE(3)
         W(4)=W(1)+W3*W(1)
         DWDX(4)=DWDX(1)+DWDX3*W(1)+W3*DWDX(1)
         DWDY(4)=DWDY(1)+DWDY3*W(1)+W3*DWDY(1)
      CASE(1)
         SELECT CASE(II)
         CASE(1)
            WW=RR1*SS1/4.0
            DWDXX=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
            DWDYY=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
         CASE(2)
            WW=RR2*SS1/4.0
            DWDXX=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
            DWDYY=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
         CASE(3)
            WW=RR2*SS2/4.0
            DWDXX=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
            DWDYY=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
         CASE(4)
            WW=RR1*SS2/4.0
            DWDXX=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
            DWDYY=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
         CASE DEFAULT
            WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
            STOP
         END SELECT
         W(4)=W(1)+W3*WW
         DWDX(4)=DWDX(1)+DWDX3*WW+W3*DWDXX
         DWDY(4)=DWDY(1)+DWDY3*WW+W3*DWDYY
      CASE DEFAULT
         WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
         STOP
      END SELECT
!
      IND=NDMFS(JJ+4,NN)
      SELECT CASE(IND)
      CASE(4)
         W(5)=W(2)*W3
         DWDX(5)=DWDX(2)*W3+W(2)*DWDX3
         DWDY(5)=DWDY(2)*W3+W(2)*DWDY3
      CASE(3)
         W(5)=W(2)+W3*W(2)
         DWDX(5)=DWDX(2)+DWDX3*W(2)+W3*DWDX(2)
         DWDY(5)=DWDY(2)+DWDY3*W(2)+W3*DWDY(2)
      CASE(1)
         SELECT CASE(JJ)
         CASE(1)
            WW=RR1*SS1/4.0
            DWDXX=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
            DWDYY=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
         CASE(2)
            WW=RR2*SS1/4.0
            DWDXX=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
            DWDYY=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
         CASE(3)
            WW=RR2*SS2/4.0
            DWDXX=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
            DWDYY=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
         CASE(4)
            WW=RR1*SS2/4.0
            DWDXX=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
            DWDYY=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
         CASE DEFAULT
            WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
            STOP
         END SELECT
         W(5)=W(2)+W3*WW
         DWDX(5)=DWDX(2)+DWDX3*WW+W3*DWDXX
         DWDY(5)=DWDY(2)+DWDY3*WW+W3*DWDYY
      CASE DEFAULT
         WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
         STOP
      END SELECT
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD_COUP1
!
      SUBROUTINE SHEPARD_TRUSS(NEleAss,II,XX,W,DWDX,IND)
!{{{
!---------------------------------------------------------------
!     W(1):  WEIGHTING FUNCTION CORRESPONDING TO NODE II
!     DWDX(1): ...
!     W(2):  SUM OF WEIGHTING FUNCTIONS
!     DWDX(2): ...
!     W(3):  SUM OF WEIGHTING FUNCTIONS FOR BLENDING ELEMENT
!     DWDX(3): ...
!---------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: NEleAss, II
      REAL(8) :: XX
      REAL(8) :: W(3),DWDX(3)
      INTEGER :: IND
!     ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: M
      INTEGER :: INDSP
      REAL(8) :: Radius
      REAL(8) :: S,SS(4)
      REAL(8) :: CONSTM1,CONSTMX
      REAL(8) :: W2,DWDX2,W3,DWDX3
!
      W(2)=0.0D0
      DWDX(2)=0.0D0
      W3=0.0D0
      DWDX3=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(3,J)
             Radius=ELEPAR(2,J)
!
             CONSTM1=1/Radius
             CONSTMX=(XX-X(M))
!
             S=ABS(CONSTMX)*CONSTM1
!
             IF(S.LE.1) THEN
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  W2=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
                     W2=2/3.0-4*SS(2)+4*SS(3)
                     DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX)
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
                     W2=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDX2=SS(4)*CONSTM1*SIGN(1.0D0,CONSTMX) 
                  ENDIF
                CASE(3)
                  W2=1-S
                  DWDX2=-CONSTM1*SIGN(1.0D0,CONSTMX)
                CASE(4)
                  W2=1-S
                  DWDX2=-CONSTM1*SIGN(1.0D0,CONSTMX)
!
                  W3=W3+W2
                  DWDX3=DWDX3+DWDX2
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT
!
                IF(J.EQ.II) THEN
                   W(1)=W2
                   DWDX(1)=DWDX2
                ENDIF
!
                W(2)=W(2)+W2 
                DWDX(2)=DWDX(2)+DWDX2
             ENDIF
         ENDIF
      ENDDO
!
      INDSP=NDMFS(3,II)
      Radius=ELEPAR(2,II)
!
      IF(INDSP.EQ.4) THEN
         W(3)=W(1)*W3
         DWDX(3)=DWDX(1)*W3+W(1)*DWDX3
      ELSE
         M=NDMFS(1,II)
         CONSTM1=1/Radius
         CONSTMX=(XX-X(M))
         S=ABS(CONSTMX)*CONSTM1
!
         W(3)=W(1)+W3*(1-S)
         DWDX(3)=DWDX(1)+DWDX3*(1-S)+W3*(-CONSTM1*SIGN(1.0D0,CONSTMX))
      ENDIF
!
      RETURN
!}}}
      END SUBROUTINE SHEPARD_TRUSS
!
      SUBROUTINE QUADS_MFS_HIHJ_V2(NEleAss,II,JJ,XX,YY,NPORMAX1,NPORMAX2,HI,HJ) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,JJ,NPORMAX1,NPORMAX2
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8)            :: HI(3,NPORMAX1),HJ(3,NPORMAX2)
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: W(5),DWDX(5),DWDY(5)
        REAL(8) :: RAD_IX,RAD_IY,RAD_JX,RAD_JY
        REAL(8) :: CONSTW2,CONSTW2_X,CONSTW2_Y
        REAL(8) :: CONSTX_I,CONSTX_J,CONSTY_I,CONSTY_J
        INTEGER :: IND_BLEND
        INTEGER :: INDSPI,INDSPJ
        INTEGER :: INDN,INDM
        REAL(8) :: NX,NY
        REAL(8) :: X_BAR,Y_BAR
!
        IND_BLEND=0
!
        RAD_IX=ELEPAR(3,II)
        RAD_IX=1.0/RAD_IX
!
        RAD_IY=ELEPAR(4,II)
        RAD_IY=1.0/RAD_IY
!
        RAD_JX=ELEPAR(3,JJ)
        RAD_JX=1.0/RAD_JX
!
        RAD_JY=ELEPAR(4,JJ)
        RAD_JY=1.0/RAD_JY
!
        N=NDMFS(1,II)
        M=NDMFS(1,JJ)
!
        CALL SHEPARD_2D(NEleAss,II,JJ,XX,YY,W,DWDX,DWDY,IND_BLEND)
!
        CONSTW2=1/W(3)
        CONSTW2_X=DWDX(3)*CONSTW2*CONSTW2
        CONSTW2_Y=DWDY(3)*CONSTW2*CONSTW2
!
        CONSTX_I=(XX-X(N))*RAD_IX
        CONSTY_I=(YY-Y(N))*RAD_IY
        CONSTX_J=(XX-X(M))*RAD_JX
        CONSTY_J=(YY-Y(M))*RAD_JY
!
        INDN=NDMFS(7,II)
        INDM=NDMFS(7,JJ)
        IF((INDN.EQ.1).AND.(PENALTY.EQ.0)) THEN
          NX=ELEPAR(5,II)
          NY=ELEPAR(6,II)
          X_BAR=(NX*CONSTX_I+NY*CONSTY_I)
          Y_BAR=(-NY*CONSTX_I+NX*CONSTY_I)
          SELECT CASE(NPORMAX1)
!{{{
          CASE(1)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
          CASE(2)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
          CASE(3)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HI(1,3)=HI(1,2)*Y_BAR
           HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
           HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
         CASE(4)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HI(1,2)=HI(1,1)*X_BAR
           HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
           HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HI(1,3)=HI(1,2)*Y_BAR
           HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
           HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
!   HIXX,DHIXXDX,DHIXXDY
           HI(1,4)=HI(1,2)*X_BAR
           HI(2,4)=HI(2,2)*X_BAR+HI(1,2)*RAD_IX*NX
           HI(3,4)=HI(3,2)*X_BAR+HI(1,2)*RAD_IY*NY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ELSE
          SELECT CASE(NPORMAX1)
!{{{
          CASE(1)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
          CASE(3)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
          CASE(4)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
          CASE(6)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
!   HIXX,DHIXXDX,DHIXXDY
                    HI(1,5)=HI(1,2)*CONSTX_I
                    HI(2,5)=HI(2,2)*CONSTX_I+HI(1,2)*RAD_IX
                    HI(3,5)=HI(3,2)*CONSTX_I
!   HIYY,DHIYYDX,DHIYYDY
                    HI(1,6)=HI(1,3)*CONSTY_I
                    HI(2,6)=HI(2,3)*CONSTY_I
                    HI(3,6)=HI(3,3)*CONSTY_I+HI(1,3)*RAD_IY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ENDIF
!
        IF((INDM.EQ.1).AND.(PENALTY.EQ.0)) THEN
          NX=ELEPAR(5,JJ)
          NY=ELEPAR(6,JJ)
          X_BAR=(NX*CONSTX_J+NY*CONSTY_J)
          Y_BAR=(-NY*CONSTX_J+NX*CONSTY_J)
          SELECT CASE(NPORMAX2)
!{{{
          CASE(1)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
          CASE(2)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
          CASE(3)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HJ(1,3)=HJ(1,2)*Y_BAR
           HJ(2,3)=HJ(2,2)*Y_BAR+HJ(1,2)*(-NY)*RAD_JX
           HJ(3,3)=HJ(3,2)*Y_BAR+HJ(1,2)*NX*RAD_JY
         CASE(4)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
           HJ(1,2)=HJ(1,1)*X_BAR
           HJ(2,2)=HJ(2,1)*X_BAR+HJ(1,1)*RAD_JX*NX
           HJ(3,2)=HJ(3,1)*X_BAR+HJ(1,1)*RAD_JY*NY
!   HIXY,DHIXYDX,DHIXYDY
           HJ(1,3)=HJ(1,2)*Y_BAR
           HJ(2,3)=HJ(2,2)*Y_BAR+HJ(1,2)*(-NY)*RAD_JX
           HJ(3,3)=HJ(3,2)*Y_BAR+HJ(1,2)*NX*RAD_JY
!   HIXX,DHIXXDX,DHIXXDY
           HJ(1,4)=HJ(1,2)*X_BAR
           HJ(2,4)=HJ(2,2)*X_BAR+HJ(1,2)*RAD_JX*NX
           HJ(3,4)=HJ(3,2)*X_BAR+HJ(1,2)*RAD_JY*NY
          CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
          END SELECT
        ELSE
         SELECT CASE(NPORMAX2)
!{{{
        CASE(1)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
        CASE(3)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
        CASE(4)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                 HJ(1,4)=HJ(1,2)*CONSTY_J
                 HJ(2,4)=HJ(2,2)*CONSTY_J
                 HJ(3,4)=HJ(3,3)*CONSTX_J
        CASE(6)
           HJ(1,1)=W(2)*CONSTW2
           HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
           HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HJ(1,2)=HJ(1,1)*CONSTX_J
              HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
              HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
              HJ(1,3)=HJ(1,1)*CONSTY_J
              HJ(2,3)=HJ(2,1)*CONSTY_J
              HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                 HJ(1,4)=HJ(1,2)*CONSTY_J
                 HJ(2,4)=HJ(2,2)*CONSTY_J
                 HJ(3,4)=HJ(3,3)*CONSTX_J
!   HIXX,DHIXXDX,DHIXXDY
                    HJ(1,5)=HJ(1,2)*CONSTX_J
                    HJ(2,5)=HJ(2,2)*CONSTX_J+HJ(1,2)*RAD_JX
                    HJ(3,5)=HJ(3,2)*CONSTX_J
!   HIYY,DHIYYDX,DHIYYDY
                    HJ(1,6)=HJ(1,3)*CONSTY_J
                    HJ(2,6)=HJ(2,3)*CONSTY_J
                    HJ(3,6)=HJ(3,3)*CONSTY_J+HJ(1,3)*RAD_JY
        CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            STOP
!}}}
         END SELECT
        ENDIF
!
        IF(IND_BLEND.EQ.0) RETURN
!
           INDSPI=NDMFS(6,II)
           INDSPJ=NDMFS(6,JJ)
           IF(INDSPI.GE.3) THEN
              HI(1,1)=W(4)*CONSTW2
              HI(2,1)=DWDX(4)*CONSTW2-W(4)*CONSTW2_X
              HI(3,1)=DWDY(4)*CONSTW2-W(4)*CONSTW2_Y
           ENDIF
!
           IF(INDSPJ.GE.3) THEN
              HJ(1,1)=W(5)*CONSTW2
              HJ(2,1)=DWDX(5)*CONSTW2-W(5)*CONSTW2_X
              HJ(3,1)=DWDY(5)*CONSTW2-W(5)*CONSTW2_Y
           ENDIF
!
        RETURN
!}}}
      END SUBROUTINE QUADS_MFS_HIHJ_V2
!
      SUBROUTINE COUP1_HIHJ(II,JJ,NN,RR,SS,NPORMAX1,NPORMAX2,HI,HJ,DET_J) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: II,JJ,NN,NPORMAX1,NPORMAX2
        REAL(8),INTENT(IN)  :: RR,SS
        REAL(8)             :: HI(3,NPORMAX1),HJ(3,NPORMAX2)
        REAL(8),INTENT(OUT) :: DET_J
!   ASSISTANT VARIABLES
        INTEGER :: IND_BLEND
        REAL(8) :: RAD_IX,RAD_IY,RAD_JX,RAD_JY
        REAL(8) :: XX,YY
        INTEGER :: IND
        INTEGER :: N,M
        REAL(8) :: W(5),DWDX(5),DWDY(5)
        REAL(8) :: CONSTW2,CONSTW2_X,CONSTW2_Y
        REAL(8) :: CONSTX_I,CONSTX_J,CONSTY_I,CONSTY_J
!
        IND_BLEND=0
!
        RAD_IX=ELEPAR(II,NN)
        RAD_IX=1.0/RAD_IX
        RAD_IY=RAD_IX
!
        RAD_JX=ELEPAR(JJ,NN)
        RAD_JX=1.0/RAD_JX
        RAD_JY=RAD_JX
!
        N=NDMFS(II,NN)
        M=NDMFS(JJ,NN)
!
        XX=0.0D0
        YY=0.0D0
        DET_J=0.0D0
        CALL SHEPARD_COUP1(II,JJ,NN,RR,SS,XX,YY,W,DWDX,DWDY,IND_BLEND,DET_J)
!
        CONSTW2=1/W(3)
        CONSTW2_X=DWDX(3)*CONSTW2*CONSTW2
        CONSTW2_Y=DWDY(3)*CONSTW2*CONSTW2
!
        CONSTX_I=(XX-X(N))*RAD_IX
        CONSTY_I=(YY-Y(N))*RAD_IY
        CONSTX_J=(XX-X(M))*RAD_JX
        CONSTY_J=(YY-Y(M))*RAD_JY
!
         SELECT CASE(NPORMAX1)
!{{{
         CASE(1)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
         CASE(3)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
         CASE(4)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
         CASE(6)
           HI(1,1)=W(1)*CONSTW2
           HI(2,1)=DWDX(1)*CONSTW2-W(1)*CONSTW2_X
           HI(3,1)=DWDY(1)*CONSTW2-W(1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
              HI(1,2)=HI(1,1)*CONSTX_I
              HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
              HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
              HI(1,3)=HI(1,1)*CONSTY_I
              HI(2,3)=HI(2,1)*CONSTY_I
              HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                 HI(1,4)=HI(1,2)*CONSTY_I
                 HI(2,4)=HI(2,2)*CONSTY_I
                 HI(3,4)=HI(3,3)*CONSTX_I
!   HIXX,DHIXXDX,DHIXXDY
                    HI(1,5)=HI(1,2)*CONSTX_I
                    HI(2,5)=HI(2,2)*CONSTX_I+HI(1,2)*RAD_IX
                    HI(3,5)=HI(3,2)*CONSTX_I
!   HIYY,DHIYYDX,DHIYYDY
                    HI(1,6)=HI(1,3)*CONSTY_I
                    HI(2,6)=HI(2,3)*CONSTY_I
                    HI(3,6)=HI(3,3)*CONSTY_I+HI(1,3)*RAD_IY
         CASE DEFAULT
            WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
            WRITE(*,*) "NPORMAX1",NPORMAX1
            STOP
!}}}
         END SELECT
!
         SELECT CASE(NPORMAX2)
!{{{
         CASE(1)
            HJ(1,1)=W(2)*CONSTW2
            HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
            HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
         CASE(3)
            HJ(1,1)=W(2)*CONSTW2
            HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
            HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HJ(1,2)=HJ(1,1)*CONSTX_J
               HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
               HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
               HJ(1,3)=HJ(1,1)*CONSTY_J
               HJ(2,3)=HJ(2,1)*CONSTY_J
               HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
         CASE(4)
            HJ(1,1)=W(2)*CONSTW2
            HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
            HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HJ(1,2)=HJ(1,1)*CONSTX_J
               HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
               HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
               HJ(1,3)=HJ(1,1)*CONSTY_J
               HJ(2,3)=HJ(2,1)*CONSTY_J
               HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                  HJ(1,4)=HJ(1,2)*CONSTY_J
                  HJ(2,4)=HJ(2,2)*CONSTY_J
                  HJ(3,4)=HJ(3,3)*CONSTX_J
         CASE(6)
            HJ(1,1)=W(2)*CONSTW2
            HJ(2,1)=DWDX(2)*CONSTW2-W(2)*CONSTW2_X
            HJ(3,1)=DWDY(2)*CONSTW2-W(2)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HJ(1,2)=HJ(1,1)*CONSTX_J
               HJ(2,2)=HJ(2,1)*CONSTX_J+HJ(1,1)*RAD_JX
               HJ(3,2)=HJ(3,1)*CONSTX_J
!   HIY,DHIYDX,DHIYDY
               HJ(1,3)=HJ(1,1)*CONSTY_J
               HJ(2,3)=HJ(2,1)*CONSTY_J
               HJ(3,3)=HJ(3,1)*CONSTY_J+HJ(1,1)*RAD_JY
!   HIXY,DHIXYDX,DHIXYDY
                  HJ(1,4)=HJ(1,2)*CONSTY_J
                  HJ(2,4)=HJ(2,2)*CONSTY_J
                  HJ(3,4)=HJ(3,3)*CONSTX_J
!   HIXX,DHIXXDX,DHIXXDY
                     HJ(1,5)=HJ(1,2)*CONSTX_J
                     HJ(2,5)=HJ(2,2)*CONSTX_J+HJ(1,2)*RAD_JX
                     HJ(3,5)=HJ(3,2)*CONSTX_J
!   HIYY,DHIYYDX,DHIYYDY
                     HJ(1,6)=HJ(1,3)*CONSTY_J
                     HJ(2,6)=HJ(2,3)*CONSTY_J
                     HJ(3,6)=HJ(3,3)*CONSTY_J+HJ(1,3)*RAD_JY
         CASE DEFAULT
             WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
             WRITE(*,*) "NPORMAX2",NPORMAX2
             STOP
!}}}
         END SELECT
!
         IF(IND_BLEND.EQ.0) RETURN
!
         HI(1,1)=W(4)*CONSTW2
         HI(2,1)=DWDX(4)*CONSTW2-W(4)*CONSTW2_X
         HI(3,1)=DWDY(4)*CONSTW2-W(4)*CONSTW2_Y
!
         HJ(1,1)=W(5)*CONSTW2
         HJ(2,1)=DWDX(5)*CONSTW2-W(5)*CONSTW2_X
         HJ(3,1)=DWDY(5)*CONSTW2-W(5)*CONSTW2_Y
!
         RETURN
!}}}
      END SUBROUTINE COUP1_HIHJ
!
      FUNCTION QUADS_MFS_HIHJ(NEleAss,II,NPORMAX,XX,YY) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8) :: H(NPORMAX)
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        INTEGER :: INDSP
        REAL(8) :: W(3)
        REAL(8) :: Radius1,Radius2,RAD1,RAD2
        REAL(8) :: S,SS(3),SX,SY
        REAL(8) :: CONSTX,CONSTY
        REAL(8) :: CONST1
        REAL(8) :: WW,W3
        LOGICAL :: EXPRN1
        INTEGER :: IND_BLEND
        INTEGER :: INDN
        REAL(8) :: NX,NY,X_BAR,Y_BAR
!
        DO I=1,NPORMAX
           H(I)=0.0D0
        ENDDO
!
        IND_BLEND=0
!
        Radius1=ELEPAR(3,II)
        Radius2=ELEPAR(4,II)
!
        N=NDMFS(1,II)
        CONSTX=(XX-X(N))/Radius1
        CONSTY=(YY-Y(N))/Radius2
!
        W(1)=0.0D0
        W(2)=0.0D0
        W3=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(6,J)
             RAD1=ELEPAR(3,J)
             RAD1=1.0/RAD1
             IF(INDSP.GE.3) THEN
               RAD2=ELEPAR(4,J)
               RAD2=1.0/RAD2
!
               SX=1.0-ABS(XX-X(M))*RAD1
               EXPRN1=SX.LT.1E-12
               IF(EXPRN1) CYCLE
               SY=1.0-ABS(YY-Y(M))*RAD2
               EXPRN1=SY.LT.1E-12
               IF(EXPRN1) CYCLE
             ELSE
               S=(XX-X(M))*(XX-X(M))+(YY-Y(M))*(YY-Y(M))
               S=SQRT(S)*RAD1
               EXPRN1=(1-S).LT.1E-12
               IF(EXPRN1) CYCLE
             ENDIF
!
               SELECT CASE(INDSP)
               CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  WW=1-6*SS(1)+8*SS(2)-3*SS(3)
               CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     WW=2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     WW=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
               CASE(3)
                  WW=SX*SY
               CASE(4)
                  IND_BLEND=1
                  WW=SX*SY                  
!
                  W3=W3+WW                  
               CASE DEFAULT
                  WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                  WRITE(*,*) INDSP
                  STOP
               END SELECT 
!
             IF(J.EQ.II) W(1)=WW
             W(2)=W(2)+WW                
           ENDIF
        ENDDO
!
        CONST1=1/W(2)
        INDN=NDMFS(7,II)
        IF((INDN.EQ.1).AND.(PENALTY.EQ.0)) THEN
          NX=ELEPAR(5,II)
          NY=ELEPAR(6,II)
          X_BAR=(NX*CONSTX+NY*CONSTY)
          Y_BAR=(-NY*CONSTX+NX*CONSTY)
          SELECT CASE(NPORMAX) 
          CASE(1)
           H(1)=W(1)*CONST1
          CASE(2)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
          CASE(3)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
           H(3)=H(2)*Y_BAR
          CASE(4)
           H(1)=W(1)*CONST1
           H(2)=H(1)*X_BAR
           H(3)=H(2)*Y_BAR
           H(4)=H(2)*X_BAR
          CASE DEFAULT
           WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
           STOP
          END SELECT
        ELSE
!
          SELECT CASE(NPORMAX) 
          CASE(1)
           H(1)=W(1)*CONST1
          CASE(3)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
          CASE(4)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
           H(4)=H(2)*CONSTY
          CASE(6)
           H(1)=W(1)*CONST1
           H(2)=H(1)*CONSTX
           H(3)=H(1)*CONSTY
           H(4)=H(2)*CONSTY
           H(5)=H(2)*CONSTX
           H(6)=H(3)*CONSTY
          CASE DEFAULT
           WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
           STOP
          END SELECT
        ENDIF
!
        IF(IND_BLEND.EQ.0) RETURN
!
        INDSP=NDMFS(6,II)
        IF(INDSP.EQ.4) W(3)=W(1)*W3
        IF(INDSP.EQ.3) THEN
           SX=1.0-ABS(CONSTX)
           SY=1.0-ABS(CONSTY)
!
           W(3)=W(1)+W3*SX*SY
        ENDIF
!
        IF(INDSP.GE.3) THEN
          H(1)=W(3)*CONST1 
        ENDIF  
!
        RETURN
!}}}
      END FUNCTION QUADS_MFS_HIHJ
!    
      FUNCTION QUADS_OVLP_HIHJ(II,NDOFF,XX,YY,LENGTH) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: II,NDOFF
        REAL(8),INTENT(IN) :: XX,YY,LENGTH
        REAL(8) :: H(NPOR_MAX)
!   ASSISTANT VARIABLES
        INTEGER :: I,N
        REAL(8) :: RAD1,RAD2,CONSTX,CONSTY
        REAL(8) :: S
        INTEGER :: IND
        REAL(8) :: NX,NY,X_BAR,Y_BAR
!
        DO I=1,NPOR_MAX
           H(I)=0.0D0
        ENDDO
!
        RAD1=2.0D0*ELEPAR(1,II)
        RAD2=2.0D0*ELEPAR(2,II)
!
        N=NDMFS(1,II)
        CONSTX=(XX-X(N))/RAD1
        CONSTY=(YY-Y(N))/RAD2
        S=(XX-X(N))*(XX-X(N))+(YY-Y(N))*(YY-Y(N))
        S=SQRT(S)
!
        IND=NDMFS(5,II)
        IF(IND.EQ.1) THEN
          NX=ELEPAR(3,II)
          NY=ELEPAR(4,II)
          X_BAR=(NX*CONSTX+NY*CONSTY)
          Y_BAR=(-NY*CONSTX+NX*CONSTY)
          SELECT CASE(NDOFF) 
          CASE(1)
               H(1)=1.0D0-S/LENGTH
          CASE(2)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*X_BAR
          CASE(3)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*X_BAR
               H(3)=H(2)*Y_BAR
          CASE(4)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*X_BAR
               H(3)=H(2)*Y_BAR
               H(4)=H(2)*X_BAR
          CASE DEFAULT
               WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
               STOP
          END SELECT
        ELSE
          SELECT CASE(NDOFF) 
          CASE(1)
               H(1)=1.0D0-S/LENGTH
          CASE(3)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*CONSTX
               H(3)=H(1)*CONSTY
          CASE(4)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*CONSTX
               H(3)=H(1)*CONSTY
               H(4)=H(2)*CONSTY
          CASE(6)
               H(1)=1.0D0-S/LENGTH
               H(2)=H(1)*CONSTX
               H(3)=H(1)*CONSTY
               H(4)=H(2)*CONSTY
               H(5)=H(2)*CONSTX
               H(6)=H(3)*CONSTY
          CASE DEFAULT
               WRITE(*,*) "ERROR, NPOR SHOULD LT 4!"
               STOP
          END SELECT
        ENDIF
!
        RETURN
!}}}
      END FUNCTION QUADS_OVLP_HIHJ
!
      FUNCTION TRUSS_MFS_HIHJ(NEleAss,II,NPORMAX,XX,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NEleAss
      INTEGER, INTENT(IN) :: II,NPORMAX
      INTEGER, INTENT(IN) :: IND
      REAL(8), INTENT(IN) :: XX
      REAL(8)             :: H(NPOR_MAX)
!   ASSISTANT VARIABLES
      INTEGER :: I,J
      INTEGER :: N,M
      INTEGER :: INDSP
      REAL(8) :: W(3)
      REAL(8) :: Radius1,Radius
      REAL(8) :: S,SS(3)
      REAL(8) :: CONSTX,CONST2
      REAL(8) :: W2
      REAL(8) :: W3
!
      Radius1=ELEPAR(2,II)
      N=NDMFS(1,II)
      CONSTX=(XX-X(N))/Radius1
!
        W(2)=0.0D0
        W3=0.0D0
        DO I=1,NEleAss
           J=LM(I,II)
           IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(3,J)
             Radius=ELEPAR(2,J)
             S=ABS(XX-X(M))/Radius
             IF(S.LE.1) THEN
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
!
                  W2=1-6*SS(1)+8*SS(2)-3*SS(3)
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     W2=2/3.0-4*SS(2)+4*SS(3)
                  ELSE
                     W2=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                  ENDIF
                CASE(3)
                  W2=1-S
                CASE(4)
                  W2=1-S
!
                  W3=W3+W2
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT 
                IF(J.EQ.II) W(1)=W2
                W(2)=W(2)+W2
             ENDIF
           ENDIF
        ENDDO
!
        INDSP=NDMFS(3,II)
        S=ABS(CONSTX)
        IF(INDSP.EQ.4) THEN
           W(3)=W(1)*W3
        ELSE
           W(3)=W(1)+W3*(1-S)
        ENDIF
!
        CONST2=1/W(2)
        SELECT CASE(NPORMAX)
        CASE(1)
           H(1)=W(3)*CONST2
        CASE(2)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(1)=W(3)*CONST2
        CASE(3)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(3)=H(2)*CONSTX
           H(1)=W(3)*CONST2
        CASE(4)
           H(1)=W(1)*CONST2
           H(2)=H(1)*CONSTX
           H(3)=H(2)*CONSTX
           H(4)=H(3)*CONSTX
           H(1)=W(3)*CONST2
        CASE DEFAULT
             WRITE(*,*) "NPORMAX.GT.4 IN TRUSS_MFS_HIHJ"
             STOP
        END SELECT     
!
      RETURN
!}}}
      END FUNCTION TRUSS_MFS_HIHJ
!
      FUNCTION TRUSS_MFS_HIHJ_V2(NEleAss,II,NPORMAX,XX,IND) RESULT(H)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,II,NPORMAX
        REAL(8),INTENT(IN) :: XX
        REAL(8)            :: H(NPORMAX)
        INTEGER :: IND
!   ASSISTANT VARIABLES
        INTEGER :: I,J
        INTEGER :: N,M
        REAL(8) :: HH1,HH2,HH3
        REAL(8) :: YY
        REAL(8) :: W(3),DWDX(3)
        REAL(8) :: Radius1
        REAL(8) :: CONST1,CONST2,CONSTX
!
        Radius1=ELEPAR(2,II)
        CONST1=1/Radius1
        N=NDMFS(1,II)
!
        CALL SHEPARD_TRUSS(NEleAss,II,XX,W,DWDX,IND)
        CONST2=1/W(2)
!
        CONSTX=(XX-X(N))*CONST1
!
        SELECT CASE(NPORMAX)
        CASE(1)
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(2)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           H(2)=H(1)*CONSTX+HH1*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(3)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           HH2=HH1*CONSTX
           H(2)=H(1)*CONSTX+HH1*CONST1
           H(3)=H(2)*CONSTX+HH2*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE(4)
           HH1=W(1)*CONST2
           H(1)=DWDX(1)*CONST2-W(1)*DWDX(2)*CONST2*CONST2           
           HH2=HH1*CONSTX
           H(2)=H(1)*CONSTX+HH1*CONST1
           HH3=HH2*CONSTX
           H(3)=H(2)*CONSTX+HH2*CONST1
           H(4)=H(3)*CONSTX+HH3*CONST1
           H(1)=DWDX(3)*CONST2-W(3)*DWDX(2)*CONST2*CONST2
        CASE DEFAULT
             WRITE(*,*) "NPORMAX.GT.4 IN QUADS_MFS_HIHJ2"
             STOP
        END SELECT
!
        RETURN
!}}}
      END FUNCTION TRUSS_MFS_HIHJ_V2
!
      SUBROUTINE OUTPUT_2DMFS_SHEPARD(NEleAss,II,XX,YY,IND_BLEND,W,DWDX,DWDY)
!{{{
      IMPLICIT NONE
      INTEGER :: NEleAss,II
      REAL(8) :: XX,YY
      INTEGER :: IND_BLEND
      REAL(8) :: W(NEleAss,2),DWDX(NEleAss,2),DWDY(NEleAss,2)
!     ASSISTANT VARIABLES
      INTEGER :: I,J,M
      INTEGER :: INDSP
      REAL(8) :: RAD1,RAD2
      REAL(8) :: S,SS(4),SX,SY,S1
      REAL(8) :: CONSTX,CONSTY
      REAL(8) :: CONSTXD,CONSTYD
      REAL(8) :: W3,DWDX3,DWDY3
      LOGICAL :: EXPRN1
!
      DO J=1,2
         DO I=1,NEleAss
            W(I,J)=0.0D0
            DWDX(I,J)=0.0D0
            DWDY(I,J)=0.0D0
         ENDDO
      ENDDO
!
      W3=0.0D0
      DWDX3=0.0D0
      DWDY3=0.0D0
      DO I=1,NEleAss
         J=LM(I,II)
         IF (J.NE.0) THEN
             M=NDMFS(1,J)
             INDSP=NDMFS(6,J)
             RAD1=ELEPAR(3,J)
             RAD1=1/RAD1
!
             IF(INDSP.GE.3) THEN
                RAD2=ELEPAR(4,J)
                RAD2=1.0/RAD2
                CONSTX=(XX-X(M))
!
                SX=1.0-ABS(CONSTX)*RAD1
                EXPRN1=SX.LT.1E-12
                IF(EXPRN1) CYCLE
                CONSTY=(YY-Y(M))
                SY=1.0-ABS(CONSTY)*RAD2
                EXPRN1=SY.LT.1E-12
                IF(EXPRN1) CYCLE
             ELSE
                CONSTX=(XX-X(M))
                CONSTY=(YY-Y(M))
!
                S1=SQRT(CONSTX*CONSTX+CONSTY*CONSTY)
                S=S1*RAD1
                EXPRN1=(1.0-S).LT.1E-12
                IF(EXPRN1) CYCLE
                IF(S.LT.1E-12) THEN
                   CONSTXD=0.0D0
                   CONSTYD=0.0D0
                ELSE
                   CONSTXD=CONSTX/S1*RAD1
                   CONSTYD=CONSTY/S1*RAD1
                ENDIF
             ENDIF
!
                SELECT CASE(INDSP)
                CASE(1)
                  SS(1)=S*S
                  SS(2)=SS(1)*S
                  SS(3)=SS(2)*S
                  SS(4)=-12*S+24*SS(1)-12*SS(2)
!
                  W(I,1)=1-6*SS(1)+8*SS(2)-3*SS(3)
                  DWDX(I,1)=SS(4)*CONSTXD
                  DWDY(I,1)=SS(4)*CONSTYD
                CASE(2)
                  SS(1)=S
                  SS(2)=S*S
                  SS(3)=SS(2)*S
                  IF(S.LE.0.5) THEN
                     SS(4)=-8*SS(1)+12*SS(2)
!
                     W(I,1)=2/3.0-4*SS(2)+4*SS(3)
                     DWDX(I,1)=SS(4)*CONSTXD
                     DWDY(I,1)=SS(4)*CONSTYD
                  ELSE
                     SS(4)=-4+8*SS(1)-4*SS(2)
!
                     W(I,1)=4/3.0-4*SS(1)+4*SS(2)-4/3.0*SS(3)
                     DWDX(I,1)=SS(4)*CONSTXD
                     DWDY(I,1)=SS(4)*CONSTYD
                  ENDIF
                CASE(3)
                     W(I,1)=SX*SY
                     DWDX(I,1)=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDY(I,1)=-SX*RAD2*SIGN(1.0D0,CONSTY)
                CASE(4)
                     IND_BLEND=1
!
                     W(I,1)=SX*SY
                     DWDX(I,1)=-SY*RAD1*SIGN(1.0D0,CONSTX)
                     DWDY(I,1)=-SX*RAD2*SIGN(1.0D0,CONSTY)
!
                     W3=W3+W(I,1)
                     DWDX3=DWDX3+DWDX(I,1)
                     DWDY3=DWDY3+DWDY(I,1)
                CASE DEFAULT
                   WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
                   STOP
                END SELECT                                 
         ENDIF
      ENDDO
!
      IF(IND_BLEND.EQ.0) RETURN
!
      DO I=1,NEleAss         
         J=LM(I,II)                
         IF(J.NE.0) THEN
            IF(W(I,1).LT.1E-12) CYCLE
            INDSP=NDMFS(6,J)
            IF(INDSP.EQ.4) THEN
               W(I,2)=W(I,1)*W3
               DWDX(I,2)=DWDX(I,1)*W3+W(I,1)*DWDX3
               DWDY(I,2)=DWDY(I,1)*W3+W(I,1)*DWDY3            
            ENDIF
            IF(INDSP.EQ.3) THEN
               M=NDMFS(1,J)
               RAD1=ELEPAR(3,J)
               RAD1=1.0/RAD1  
               RAD2=ELEPAR(4,J)
               RAD2=1.0/RAD2
               CONSTX=(XX-X(M))
               CONSTY=(YY-Y(M))
               SX=1.0-ABS(CONSTX)*RAD1
               SY=1.0-ABS(CONSTY)*RAD2
!
               W(I,2)=W(I,1)+W3*SX*SY
               DWDX(I,2)=DWDX(I,1)+DWDX3*SX*SY-W3*SY*RAD1*SIGN(1.0D0,CONSTX)
               DWDY(I,2)=DWDY(I,1)+DWDY3*SX*SY-W3*SX*RAD2*SIGN(1.0D0,CONSTY)
            ENDIF
         ENDIF         
      ENDDO
!
      RETURN
!}}}
      END SUBROUTINE OUTPUT_2DMFS_SHEPARD    
!
      SUBROUTINE OUTPUT_2DMFS(NEleAss,KK,XX,YY,D,STRESS) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: NEleAss,KK
        REAL(8),INTENT(IN) :: XX,YY 
        REAL(8)            :: D(2),STRESS(4)
!   ASSISTANT VARIABLES
        INTEGER :: I,J,K
        INTEGER :: N,II
        INTEGER :: NPORMAX
        REAL(8) :: W(NEleAss,2),DWDX(NEleAss,2),DWDY(NEleAss,2)
        REAL(8) :: HI(3,6)
        REAL(8) :: RAD_IX,RAD_IY
        REAL(8) :: CONSTW2,CONSTW2_X,CONSTW2_Y
        REAL(8) :: CONSTX_I,CONSTY_I
        INTEGER :: IND_BLEND
        INTEGER :: INDSPI
        INTEGER :: IND,INDN,INC
        REAL(8) :: NX,NY
        REAL(8) :: X_BAR,Y_BAR
        REAL(8) :: WW,DWDXX,DWDYY
        INTEGER :: NUM,NPOR,NDOFF
        REAL(8) :: HU(4)
!
        IND_BLEND=0
!
        CALL OUTPUT_2DMFS_SHEPARD(NEleAss,KK,XX,YY,IND_BLEND,W,DWDX,DWDY)
!
        WW=0.0D0
        DWDXX=0.0D0
        DWDYY=0.0D0
        DO I=1,NEleAss
           WW=WW+W(I,1)
           DWDXX=DWDXX+DWDX(I,1)
           DWDYY=DWDYY+DWDY(I,1)
        ENDDO
!
        DO I=1,4
           HU(I)=0.0D0
        ENDDO
!
        DO K=1,NEleAss
           II=LM(K,KK)
           IF(II.NE.0) THEN
              N=NDMFS(1,II)
              DO I=1,6
                DO J=1,3
                 HI(J,I)=0.0D0
                ENDDO
              ENDDO
!
!              IF (W(K,1).LT.1E-12) CYCLE
!
              IND=NDMFS(2,II)
              IF(IND.EQ.1) THEN
                 NPORMAX=NDOF(1,N)+1
              ELSE
                 NPORMAX=NDOF(1,N)
              ENDIF
!
              RAD_IX=ELEPAR(3,II)
              RAD_IX=1.0/RAD_IX
!
              RAD_IY=ELEPAR(4,II)
              RAD_IY=1.0/RAD_IY
!
              CONSTW2=1/WW
              CONSTW2_X=DWDXX*CONSTW2*CONSTW2
              CONSTW2_Y=DWDYY*CONSTW2*CONSTW2
!
              CONSTX_I=(XX-X(N))*RAD_IX
              CONSTY_I=(YY-Y(N))*RAD_IY
!
              INDN=NDMFS(7,II)
              IF((INDN.EQ.1).AND.(PENALTY.EQ.0)) THEN
                NX=ELEPAR(5,II)
                NY=ELEPAR(6,II)
                X_BAR=(NX*CONSTX_I+NY*CONSTY_I)
                Y_BAR=(-NY*CONSTX_I+NX*CONSTY_I)
                SELECT CASE(NPORMAX)
!{{{
                CASE(1)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
                CASE(2)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                 HI(1,2)=HI(1,1)*X_BAR
                 HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
                 HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
                CASE(3)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                 HI(1,2)=HI(1,1)*X_BAR
                 HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
                 HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
      !   HIXY,DHIXYDX,DHIXYDY
                 HI(1,3)=HI(1,2)*Y_BAR
                 HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
                 HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
               CASE(4)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                 HI(1,2)=HI(1,1)*X_BAR
                 HI(2,2)=HI(2,1)*X_BAR+HI(1,1)*RAD_IX*NX
                 HI(3,2)=HI(3,1)*X_BAR+HI(1,1)*RAD_IY*NY
      !   HIXY,DHIXYDX,DHIXYDY
                 HI(1,3)=HI(1,2)*Y_BAR
                 HI(2,3)=HI(2,2)*Y_BAR+HI(1,2)*(-NY)*RAD_IX
                 HI(3,3)=HI(3,2)*Y_BAR+HI(1,2)*NX*RAD_IY
      !   HIXX,DHIXXDX,DHIXXDY
                 HI(1,4)=HI(1,2)*X_BAR
                 HI(2,4)=HI(2,2)*X_BAR+HI(1,2)*RAD_IX*NX
                 HI(3,4)=HI(3,2)*X_BAR+HI(1,2)*RAD_IY*NY
                CASE DEFAULT
                  WRITE(*,*) "ERROR, NPORMAX DOES NOT EXIST!"
                  STOP
!}}}
                END SELECT
              ELSE
                SELECT CASE(NPORMAX)
      !{{{
                CASE(1)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
                CASE(3)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                    HI(1,2)=HI(1,1)*CONSTX_I
                    HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
                    HI(3,2)=HI(3,1)*CONSTX_I
      !   HIY,DHIYDX,DHIYDY
                    HI(1,3)=HI(1,1)*CONSTY_I
                    HI(2,3)=HI(2,1)*CONSTY_I
                    HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
                CASE(4)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                    HI(1,2)=HI(1,1)*CONSTX_I
                    HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
                    HI(3,2)=HI(3,1)*CONSTX_I
      !   HIY,DHIYDX,DHIYDY
                    HI(1,3)=HI(1,1)*CONSTY_I
                    HI(2,3)=HI(2,1)*CONSTY_I
                    HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
      !   HIXY,DHIXYDX,DHIXYDY
                       HI(1,4)=HI(1,2)*CONSTY_I
                       HI(2,4)=HI(2,2)*CONSTY_I
                       HI(3,4)=HI(3,3)*CONSTX_I
                CASE(6)
                 HI(1,1)=W(K,1)*CONSTW2
                 HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
                 HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
      !   HIX,DHIXDX,DHIXDY
                    HI(1,2)=HI(1,1)*CONSTX_I
                    HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
                    HI(3,2)=HI(3,1)*CONSTX_I
      !   HIY,DHIYDX,DHIYDY
                    HI(1,3)=HI(1,1)*CONSTY_I
                    HI(2,3)=HI(2,1)*CONSTY_I
                    HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
      !   HIXY,DHIXYDX,DHIXYDY
                       HI(1,4)=HI(1,2)*CONSTY_I
                       HI(2,4)=HI(2,2)*CONSTY_I
                       HI(3,4)=HI(3,3)*CONSTX_I
      !   HIXX,DHIXXDX,DHIXXDY
                          HI(1,5)=HI(1,2)*CONSTX_I
                          HI(2,5)=HI(2,2)*CONSTX_I+HI(1,2)*RAD_IX
                          HI(3,5)=HI(3,2)*CONSTX_I
      !   HIYY,DHIYYDX,DHIYYDY
                          HI(1,6)=HI(1,3)*CONSTY_I
                          HI(2,6)=HI(2,3)*CONSTY_I
                          HI(3,6)=HI(3,3)*CONSTY_I+HI(1,3)*RAD_IY
                CASE DEFAULT
                  WRITE(*,*) "ERROR, NPORMAX DOES NOT EXIST!"
                  STOP
      !}}}
                END SELECT
              ENDIF
!
              IF(IND_BLEND.NE.0) THEN
                 INDSPI=NDMFS(6,II)
                 IF(INDSPI.GE.3) THEN
                    HI(1,1)=W(K,2)*CONSTW2
                    HI(2,1)=DWDX(K,2)*CONSTW2-W(K,2)*CONSTW2_X
                    HI(3,1)=DWDY(K,2)*CONSTW2-W(K,2)*CONSTW2_Y
                 ENDIF
              ENDIF
!
              DO I=1,2
                 IND=NDMFS(I+1,II)
                 IF(IND.NE.1) THEN
                    INC=0
                 ELSE
                    INC=1
                 ENDIF
                 NDOFF=NDOF(I,N)
                 IF(IND.EQ.1) THEN
                    D(I)=D(I)+RESS(I,N)*HI(1,1)
                    HU(I)=HU(I)+HI(I+1,1)*RESS(I,N)
                    HU(I+2)=HU(I+2)+HI((4-I),1)*RESS(I,N)
                 ENDIF
                 IF(NDOFF.GT.0) THEN
                   DO J=1,NDOFF                           
                      NUM=ID(I,N)+J-1
                      NPOR=J+INC
                      D(I)=D(I)+DISP(NUM)*HI(1,NPOR) 
                      HU(I)=HU(I)+HI(I+1,NPOR)*DISP(NUM)
                      HU(I+2)=HU(I+2)+HI((4-I),NPOR)*DISP(NUM)
                   ENDDO
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        STRESS(1)=KCS(1)*HU(1)+ KCS(2)*HU(2)
        STRESS(2)=KCS(1)*HU(2)+ KCS(2)*HU(1)
        STRESS(3)=KCS(3)*HU(4)+ KCS(3)*HU(3)
        STRESS(4)=STRESS(1)*STRESS(1)+STRESS(2)*STRESS(2)-STRESS(1)*STRESS(2)
        STRESS(4)=STRESS(4)+3*STRESS(3)*STRESS(3)
        STRESS(4)=SQRT(STRESS(4))
!
        RETURN
!}}}
      END SUBROUTINE OUTPUT_2DMFS
!
      SUBROUTINE OUTPUT_COUP1(RR,SS,NN,XX,YY,D,STRESS) 
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .   CALCULATE VALUE OF HI,HJ,DHI,DHJ                                .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        IMPLICIT NONE
        REAL(8),INTENT(IN)  :: RR,SS
        INTEGER,INTENT(IN)  :: NN
        REAL(8)             :: D(2),STRESS(4)
        REAL(8),INTENT(OUT) :: XX,YY
!   ASSISTANT VARIABLES
        INTEGER :: I,J,K
        INTEGER :: N
        INTEGER :: NPORMAX
        REAL(8) :: W(4,2),DWDX(4,2),DWDY(4,2)
        REAL(8) :: HI(3,6)
        REAL(8) :: RAD_IX,RAD_IY
        REAL(8) :: CONSTW2,CONSTW2_X,CONSTW2_Y
        REAL(8) :: CONSTX_I,CONSTY_I
        INTEGER :: IND_BLEND
        REAL(8) :: WW,DWDXX,DWDYY
        INTEGER :: NUM,NPOR,NDOFF
        REAL(8) :: HU(4)
!
        IND_BLEND=0
!
        XX=0.0D0
        YY=0.0D0
        CALL OUTPUT_COUP1_SHEPARD(RR,SS,NN,XX,YY,IND_BLEND,W,DWDX,DWDY)
!
        WW=0.0D0
        DWDXX=0.0D0
        DWDYY=0.0D0
        DO I=1,4
           WW=WW+W(I,1)
           DWDXX=DWDXX+DWDX(I,1)
           DWDYY=DWDYY+DWDY(I,1)
        ENDDO
!
        DO I=1,4
           HU(I)=0.0D0
        ENDDO
!
        DO K=1,4
           N=NDMFS(K,NN)
           DO I=1,6
              DO J=1,3
                 HI(J,I)=0.0D0
              ENDDO
           ENDDO
!
           NPORMAX=NDOF(1,N)
!
           RAD_IX=ELEPAR(K,NN)
           RAD_IX=1.0/RAD_IX
!
           RAD_IY=RAD_IX
!
           CONSTW2=1/WW
           CONSTW2_X=DWDXX*CONSTW2*CONSTW2
           CONSTW2_Y=DWDYY*CONSTW2*CONSTW2
!
           CONSTX_I=(XX-X(N))*RAD_IX
           CONSTY_I=(YY-Y(N))*RAD_IY
!
           SELECT CASE(NPORMAX)
!{{{
           CASE(1)
            HI(1,1)=W(K,1)*CONSTW2
            HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
            HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
           CASE(3)
            HI(1,1)=W(K,1)*CONSTW2
            HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
            HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HI(1,2)=HI(1,1)*CONSTX_I
               HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
               HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
                  HI(1,3)=HI(1,1)*CONSTY_I
                  HI(2,3)=HI(2,1)*CONSTY_I
                  HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
           CASE(4)
            HI(1,1)=W(K,1)*CONSTW2
            HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
            HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HI(1,2)=HI(1,1)*CONSTX_I
               HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
               HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
               HI(1,3)=HI(1,1)*CONSTY_I
               HI(2,3)=HI(2,1)*CONSTY_I
               HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                  HI(1,4)=HI(1,2)*CONSTY_I
                  HI(2,4)=HI(2,2)*CONSTY_I
                  HI(3,4)=HI(3,3)*CONSTX_I
           CASE(6)
            HI(1,1)=W(K,1)*CONSTW2
            HI(2,1)=DWDX(K,1)*CONSTW2-W(K,1)*CONSTW2_X
            HI(3,1)=DWDY(K,1)*CONSTW2-W(K,1)*CONSTW2_Y
!   HIX,DHIXDX,DHIXDY
               HI(1,2)=HI(1,1)*CONSTX_I
               HI(2,2)=HI(2,1)*CONSTX_I+HI(1,1)*RAD_IX
               HI(3,2)=HI(3,1)*CONSTX_I
!   HIY,DHIYDX,DHIYDY
               HI(1,3)=HI(1,1)*CONSTY_I
               HI(2,3)=HI(2,1)*CONSTY_I
               HI(3,3)=HI(3,1)*CONSTY_I+HI(1,1)*RAD_IY
!   HIXY,DHIXYDX,DHIXYDY
                  HI(1,4)=HI(1,2)*CONSTY_I
                  HI(2,4)=HI(2,2)*CONSTY_I
                  HI(3,4)=HI(3,3)*CONSTX_I
!   HIXX,DHIXXDX,DHIXXDY
                     HI(1,5)=HI(1,2)*CONSTX_I
                     HI(2,5)=HI(2,2)*CONSTX_I+HI(1,2)*RAD_IX
                     HI(3,5)=HI(3,2)*CONSTX_I
!   HIYY,DHIYYDX,DHIYYDY
                     HI(1,6)=HI(1,3)*CONSTY_I
                     HI(2,6)=HI(2,3)*CONSTY_I
                     HI(3,6)=HI(3,3)*CONSTY_I+HI(1,3)*RAD_IY
           CASE DEFAULT
             WRITE(*,*) "ERROR, NPORMAX DOES NOT EXIST!"
             STOP
!}}}
           END SELECT
!
           IF(IND_BLEND.NE.0) THEN
              HI(1,1)=W(K,2)*CONSTW2
              HI(2,1)=DWDX(K,2)*CONSTW2-W(K,2)*CONSTW2_X
              HI(3,1)=DWDY(K,2)*CONSTW2-W(K,2)*CONSTW2_Y
           ENDIF
!
           DO I=1,2
              NDOFF=NDOF(I,N)
              DO J=1,NDOFF                           
                 NUM=ID(I,N)+J-1
                 NPOR=J
                 D(I)=D(I)+DISP(NUM)*HI(1,NPOR) 
                 HU(I)=HU(I)+HI(I+1,NPOR)*DISP(NUM)
                 HU(I+2)=HU(I+2)+HI((4-I),NPOR)*DISP(NUM)
              ENDDO
           ENDDO
        ENDDO
        STRESS(1)=KCS(1)*HU(1)+ KCS(2)*HU(2)
        STRESS(2)=KCS(1)*HU(2)+ KCS(2)*HU(1)
        STRESS(3)=KCS(3)*HU(4)+ KCS(3)*HU(3)
        STRESS(4)=(STRESS(1)-STRESS(2))*(STRESS(1)-STRESS(2))
        STRESS(4)=STRESS(4)+6*STRESS(3)*STRESS(3)
        STRESS(4)=1.0/2*SQRT(STRESS(4))
!
        RETURN
!}}}
      END SUBROUTINE OUTPUT_COUP1
!
      SUBROUTINE OUTPUT_COUP1_SHEPARD(RR,SS,NN,XX,YY,IND_BLEND,W,DWDX,DWDY)
!{{{
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NN
      REAL(8),INTENT(IN) :: RR,SS
      INTEGER :: IND_BLEND
      REAL(8) :: W(4,2),DWDX(4,2),DWDY(4,2)
      REAL(8),INTENT(OUT) :: XX,YY
!     ASSISTANT VARIABLES
      INTEGER :: N
      REAL(8) :: CONSTX,CONSTY
      REAL(8) :: CONSTXD,CONSTYD
      REAL(8) :: W3,DWDX3,DWDY3
      REAL(8) :: WW,DWDXX,DWDYY
      LOGICAL :: EXPRN1
      REAL(8) :: S1,S,RAD1,SSS(4)
      INTEGER :: IND
      INTEGER :: I,J,N1,N2,N3,N4
      REAL(8) :: RR1,RR2,SS1,SS2 
      REAL(8) :: DXDR,DYDR,DXDS,DYDS
      REAL(8) :: DET_J
!
      DO J=1,2
         DO I=1,4
            W(I,J)=0.0D0
            DWDX(I,J)=0.0D0
            DWDY(I,J)=0.0D0
         ENDDO
      ENDDO
!
      N1=NDMFS(1,NN)
      N2=NDMFS(2,NN)
      N3=NDMFS(3,NN)
      N4=NDMFS(4,NN)
!     GLOBAL COORDINATES
      RR1=1+RR
      RR2=1-RR
      SS1=1+SS
      SS2=1-SS
!
      XX=RR1*SS1*X(N1)/4.0
      XX=XX+RR2*SS1*X(N2)/4.0
      XX=XX+RR2*SS2*X(N3)/4.0
      XX=XX+RR1*SS2*X(N4)/4.0
!
      YY=RR1*SS1*Y(N1)/4.0
      YY=YY+RR2*SS1*Y(N2)/4.0
      YY=YY+RR2*SS2*Y(N3)/4.0
      YY=YY+RR1*SS2*Y(N4)/4.0
!
      DXDR=(SS1*X(N1)-SS1*X(N2)-SS2*X(N3)+SS2*X(N4))/4.0
      DXDS=(RR1*X(N1)+RR2*X(N2)-RR2*X(N3)-RR1*X(N4))/4.0
      DYDR=(SS1*Y(N1)-SS1*Y(N2)-SS2*Y(N3)+SS2*Y(N4))/4.0
      DYDS=(RR1*Y(N1)+RR2*Y(N2)-RR2*Y(N3)-RR1*Y(N4))/4.0
!
      DET_J=DXDR*DYDS-DYDR*DXDS
!
      W3=0.0D0
      DWDX3=0.0D0
      DWDY3=0.0D0
      DO I=1,4
         N=NDMFS(I,NN)
         IND=NDMFS(I+4,NN)
         IF(IND.EQ.1) THEN
            RAD1=ELEPAR(I,NN)
            RAD1=1.0/RAD1
!
            CONSTX=(XX-X(N))
            CONSTY=(YY-Y(N))
!
            S1=SQRT(CONSTX*CONSTX+CONSTY*CONSTY)
            S=S1*RAD1
            EXPRN1=(1.0-S).LT.1E-12
            IF(EXPRN1) CYCLE
            IF(S.LT.1E-12) THEN
              CONSTXD=0.0D0
              CONSTYD=0.0D0
            ELSE
              CONSTXD=CONSTX/S1*RAD1
              CONSTYD=CONSTY/S1*RAD1
            ENDIF
         ENDIF
!
         SELECT CASE(IND)
         CASE(1)
             SSS(1)=S*S
             SSS(2)=SSS(1)*S
             SSS(3)=SSS(2)*S
             SSS(4)=-12*S+24*SSS(1)-12*SSS(2)
!
             W(I,1)=1-6*SSS(1)+8*SSS(2)-3*SSS(3)
             DWDX(I,1)=SSS(4)*CONSTXD
             DWDY(I,1)=SSS(4)*CONSTYD
         CASE(3)
             SELECT CASE(I)
             CASE(1)
                W(I,1)=RR1*SS1/4.0
                DWDX(I,1)=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
                DWDY(I,1)=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
             CASE(2)
                W(I,1)=RR2*SS1/4.0
                DWDX(I,1)=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
                DWDY(I,1)=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
             CASE(3)
                W(I,1)=RR2*SS2/4.0
                DWDX(I,1)=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
                DWDY(I,1)=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
             CASE(4)
                W(I,1)=RR1*SS2/4.0
                DWDX(I,1)=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
                DWDY(I,1)=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
             CASE DEFAULT
                WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
                STOP
             END SELECT
         CASE(4)
             IND_BLEND=1
             SELECT CASE(I)
             CASE(1)
                W(I,1)=RR1*SS1/4.0
                DWDX(I,1)=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
                DWDY(I,1)=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
             CASE(2)
                W(I,1)=RR2*SS1/4.0
                DWDX(I,1)=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
                DWDY(I,1)=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
             CASE(3)
                W(I,1)=RR2*SS2/4.0
                DWDX(I,1)=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
                DWDY(I,1)=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
             CASE(4)
                W(I,1)=RR1*SS2/4.0
                DWDX(I,1)=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
                DWDY(I,1)=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
             CASE DEFAULT
                WRITE(*,*) "ERROR, OPERATIONM/SHEPARD_COUP1 !"
                STOP
             END SELECT
!
             W3=W3+W(I,1)
             DWDX3=DWDX3+DWDX(I,1)
             DWDY3=DWDY3+DWDY(I,1)
         CASE DEFAULT
             WRITE(*,*) "ERROR, SHEPARD FUNCITON DOES NOT EXIST!"
             WRITE(*,*) "IND",IND
             STOP
         END SELECT                                 
      ENDDO
!
      IF(IND_BLEND.EQ.0) RETURN
!
      DO I=1,4
         IND=NDMFS(I+4,NN)
         SELECT CASE(IND)
         CASE(4)
            W(I,2)=W(I,1)*W3
            DWDX(I,2)=DWDX(I,1)*W3+W(I,1)*DWDX3
            DWDY(I,2)=DWDY(I,1)*W3+W(I,1)*DWDY3
         CASE(3)
            W(I,2)=W(I,1)+W3*W(I,1)
            DWDX(I,2)=DWDX(I,1)+DWDX3*W(I,1)+W3*DWDX(I,1)
            DWDY(I,2)=DWDY(I,1)+DWDY3*W(I,1)+W3*DWDY(I,1)
         CASE(1)
            SELECT CASE(I)
            CASE(1)
               WW=RR1*SS1/4.0
               DWDXX=(SS1*DYDS-RR1*DYDR)/4.0/DET_J
               DWDYY=(-SS1*DXDS+RR1*DXDR)/4.0/DET_J
            CASE(2)
               WW=RR2*SS1/4.0
               DWDXX=(-SS1*DYDS-RR2*DYDR)/4.0/DET_J
               DWDYY=(SS1*DXDS+RR2*DXDR)/4.0/DET_J
            CASE(3)
               WW=RR2*SS2/4.0
               DWDXX=(-SS2*DYDS+RR2*DYDR)/4.0/DET_J
               DWDYY=(SS2*DXDS-RR2*DXDR)/4.0/DET_J
            CASE(4)
               WW=RR1*SS2/4.0
               DWDXX=(SS2*DYDS+RR1*DYDR)/4.0/DET_J
               DWDYY=(-SS2*DXDS-RR1*DXDR)/4.0/DET_J
            CASE DEFAULT
               WRITE(*,*) "ERROR, OPERATIONM/OUTPUT_COUP1_SHEPARD !"
               STOP
            END SELECT
            W(I,2)=W(I,1)+W3*WW
            DWDX(I,2)=DWDX(I,1)+DWDX3*WW+W3*DWDXX
            DWDY(I,2)=DWDY(I,1)+DWDY3*WW+W3*DWDYY
         CASE DEFAULT
            WRITE(*,*) "ERROR, OPERATIONM/OUTPUT_COUP1_SHEPARD !"
            STOP
         END SELECT
      ENDDO
!
      RETURN
!}}}
      END SUBROUTINE OUTPUT_COUP1_SHEPARD    
!
      SUBROUTINE G_KCS(DIMEN,ITYPE,E,PR)
!{{{
!--------------------------------------------------
!   GENERATE CONSTITUTIVE MATRIX
!--------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: DIMEN,ITYPE
      REAL(8),INTENT(IN)  :: E,PR
      REAL(8) :: KK
!
      INTEGER :: I
      INTEGER :: ALLocateStatus
!
      SELECT CASE(DIMEN)
      CASE(2)
         ALLOCATE (KCS(3),STAT = ALLocateStatus)
         IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY KCS ***" 
!
         DO I=1,3
            KCS(I)=0.
         ENDDO
!
         SELECT CASE(ITYPE)
         CASE(1)
            KCS(1)=E*(1-PR)/(1+PR)/(1-2*PR)
            KCS(2)=E*PR/(1+PR)/(1-2*PR)
            KCS(3)=E/2/(1+PR)
         CASE(2)
            KCS(1)=E/(1-PR*PR)
            KCS(2)=E*PR/(1-PR*PR)
            KCS(3)=E/2/(1+PR)
         CASE DEFAULT
            WRITE(*,*) "ERROR, ITYPE.GT.2!"
            STOP
         END SELECT
      CASE(3)
         ALLOCATE (K3CS(6,6),STAT = ALLocateStatus)
         IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY KCS ***" 

         KK=E*(1.0D0-PR)/(1.0D0+PR)/(1.0D0-2.0D0*PR)
         K3CS=0.0D0
         K3CS(1,1)=1.0D0
         K3CS(1,2)=PR/(1.0D0-PR)
         K3CS(1,3)=PR/(1.0D0-PR)
         K3CS(2,1)=PR/(1.0D0-PR)
         K3CS(2,2)=1.0D0
         K3CS(2,3)=PR/(1.0D0-PR)
         K3CS(3,1)=PR/(1.0D0-PR)
         K3CS(3,2)=PR/(1.0D0-PR)
         K3CS(3,3)=1.0D0
         K3CS(4,4)=(1.0D0-2.0D0*PR)/2.0D0/(1.0D0-PR)
         K3CS(5,5)=(1.0D0-2.0D0*PR)/2.0D0/(1.0D0-PR)
         K3CS(6,6)=(1.0D0-2.0D0*PR)/2.0D0/(1.0D0-PR)
!
         K3CS=K3CS*KK
      CASE DEFAULT
           WRITE(*,*) "ERROR, DIMEN.GT.2!"
           STOP
      END SELECT
!
      RETURN
!}}}
      END SUBROUTINE G_KCS
!
      FUNCTION F_LOAD(FLOAD,FUNC,NX,NY,XX,YY,N1,N2) RESULT(LOAD)
!{{{
      IMPLICIT NONE
      REAL(8) :: LOAD(2)
      REAL(8),INTENT(IN)  :: FLOAD(2),NX,NY,XX,YY
      INTEGER,INTENT(IN)  :: FUNC(2), N1,N2
!     ASSISTANT VARIABLES
      INTEGER :: I
      REAL(8) :: FX,FY
      REAL(8) :: X1,Y1,X2,Y2,MAG,MAG1,MAG2
      REAL(8) :: LENGTH,LENGTH2
      REAL(8) :: THETA1,THETA2
!
      FX=-FLOAD(2)*NX+FLOAD(1)*NY
      FY=-FLOAD(2)*NY-FLOAD(1)*NX
!
      SELECT CASE(FUNC(1))
      CASE(1)
         MAG=1 
      CASE(2)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         MAG1=1.0D0-2.0D0*ACOS(ABS(Y1)/3.0D0)/Pi
         MAG2=1.0D0-2.0D0*ACOS(ABS(Y2)/3.0D0)/Pi
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(3)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         MAG1=1.0D0-2.0D0*ACOS(ABS(X1-27.0D0)/3.0D0)/Pi
         MAG2=1.0D0-2.0D0*ACOS(ABS(X2-27.0D0)/3.0D0)/Pi
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(4)  !LINEAR
         MAG=2.0D0/13.0D0*YY-1.0D0        
      CASE(5)  !QUADRATIC
         MAG=-4.0D0/169.0D0*YY*(YY-13.0D0)
      CASE(6)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         THETA1=ASIN(ABS(Y1)/4.0D0)
         THETA2=ASIN(ABS(Y2)/4.0D0)
         IF(THETA1.LE.Pi/4.0D0) THEN
            MAG1=4.0D0*THETA1/Pi
         ELSE
            MAG1=2.0D0-4.0D0*THETA1/Pi
         ENDIF
!
         IF(THETA2.LE.Pi/4.0D0) THEN
            MAG2=4.0D0*THETA2/Pi
         ELSE
            MAG2=2.0D0-4.0D0*THETA2/Pi
         ENDIF
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(7)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         THETA1=ASIN(ABS(X1)/4.0D0)
         THETA2=ASIN(ABS(X2)/4.0D0)
         IF(THETA1.LE.Pi/4.0D0) THEN
            MAG1=4.0D0*THETA1/Pi
         ELSE
            MAG1=2.0D0-4.0D0*THETA1/Pi
         ENDIF
!
         IF(THETA2.LE.Pi/4.0D0) THEN
            MAG2=4.0D0*THETA2/Pi
         ELSE
            MAG2=2.0D0-4.0D0*THETA2/Pi
         ENDIF
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE DEFAULT
        WRITE(*,*) "***STOP***, /OPERATIONM/LOAD!"
        STOP
      END SELECT
      LOAD(1)=FX*MAG
!
      SELECT CASE(FUNC(2))
      CASE(1)
         MAG=1 
      CASE(2)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         MAG1=1.0D0-2.0D0*ACOS(ABS(Y1)/3.0D0)/Pi
         MAG2=1.0D0-2.0D0*ACOS(ABS(Y2)/3.0D0)/Pi
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(3)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         MAG1=1.0D0-2.0D0*ACOS(ABS(X1-27.0D0)/3.0D0)/Pi
         MAG2=1.0D0-2.0D0*ACOS(ABS(X2-27.0D0)/3.0D0)/Pi
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(4)   ! LINEAR
         MAG=2.0D0/13.0D0*YY-1.0D0 
      CASE(5) ! QUADRATIC
         MAG=-4.0D0/169.0D0*YY*(YY-13.0D0)
      CASE(6)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         THETA1=ASIN(ABS(Y1)/4.0D0)
         THETA2=ASIN(ABS(Y2)/4.0D0)
         IF(THETA1.LE.Pi/4.0D0) THEN
            MAG1=4.0D0*THETA1/Pi
         ELSE
            MAG1=2.0D0-4.0D0*THETA1/Pi
         ENDIF
!
         IF(THETA2.LE.Pi/4.0D0) THEN
            MAG2=4.0D0*THETA2/Pi
         ELSE
            MAG2=2.0D0-4.0D0*THETA2/Pi
         ENDIF
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE(7)
         X1=X(N1)
         Y1=Y(N1)
         X2=X(N2)
         Y2=Y(N2)
         LENGTH=SQRT((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))
         LENGTH2=SQRT((XX-X2)*(XX-X2)+(YY-Y2)*(YY-Y2))
!
         THETA1=ASIN(ABS(X1)/4.0D0)
         THETA2=ASIN(ABS(X2)/4.0D0)
         IF(THETA1.LE.Pi/4.0D0) THEN
            MAG1=4.0D0*THETA1/Pi
         ELSE
            MAG1=2.0D0-4.0D0*THETA1/Pi
         ENDIF
!
         IF(THETA2.LE.Pi/4.0D0) THEN
            MAG2=4.0D0*THETA2/Pi
         ELSE
            MAG2=2.0D0-4.0D0*THETA2/Pi
         ENDIF
         MAG=MAG2+(MAG1-MAG2)*LENGTH2/LENGTH
      CASE DEFAULT
        WRITE(*,*) "***STOP***, /OPERATIONM/LOAD!"
        STOP
      END SELECT
      LOAD(2)=FY*MAG
!
      RETURN
!}}}
      END FUNCTION F_LOAD
!
      SUBROUTINE HALF_BANDWIDTH(A,NEQQ,NWKK)
!{{{
!------------------------------------------
!     CALCULATE THE HALF-BANDWIDTH WITH 
!     ZERO ENTRIES ELIMINATED
!------------------------------------------
      IMPLICIT NONE
      INTEGER(8) :: NWKK
      INTEGER :: NEQQ
      REAL(8) :: A(NWKK)
!     ASSISTANT VARIABLES
      INTEGER(8) :: I
      REAL(8) :: NWKK_R
!
      NWKK_R=0
      DO I=1,NWKK
         IF(ABS(A(I)).NE.0.0d0) THEN
           NWKK_R=NWKK_R+1.0
         ENDIF
      ENDDO

      HALF_BAND=NWKK_R/NEQQ-1
      RETURN
!}}}
      END SUBROUTINE HALF_BANDWIDTH
!
      END MODULE OPERATIONM




