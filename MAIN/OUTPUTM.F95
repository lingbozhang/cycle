! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Output moduel for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .
module OUTPUTM
    USE BUILDINS
    USE INTERPRETER
    USE WAREHOUSE
    USE GENERATEM
    USE OPERATIONM
    IMPLICIT NONE
    PUBLIC  :: OUTPUT
    PUBLIC  :: OUTPUTD
    PRIVATE :: OUTDISPD1
    PRIVATE :: OUTDISPD2
    PRIVATE :: OUTDISP_STRESS
    PRIVATE :: OUTDISP_FEM
    PRIVATE :: OUTDISP_2D_MFS
    PRIVATE :: OUTDISP_TRUSS_MFS
    PRIVATE :: OUTPUT_2DMFS_MATLAB
    PRIVATE :: OUTPUT_COUP1_MATLAB
    PRIVATE :: TRUSS_OUT
    PRIVATE :: OUTPUT_TRUSS_OVLP
    PRIVATE :: OUTPUT_TRUSS_OVLP2
    PRIVATE :: WEI_CAL_OUT
    PRIVATE :: OUTPUT_2D_OVLP
    PRIVATE :: OUTPUT_2D_OFE1
    PRIVATE :: OUTPUT_COUPLE_OFE1
    PRIVATE :: OUTPUT_QUADS
    CONTAINS
!   
    SUBROUTINE OUTPUT(IIN,IOUT,IDISP,ISTRE,ADELE,ADBC)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT AND STRESS
!----------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,ISTRE,ADELE,ADBC
!   OUTPUT DISPLACEMENT AND STRESS
!
!   Calculate strain energy 
!
    CALL OUTDISP_STRESS(IIN,IOUT,IDISP,ISTRE,ADELE,ADBC)    
!
    RETURN
!
 1050 FORMAT (' TOTAL STRAIN ENERGY =',5X,ES18.6)
!}}}
    END SUBROUTINE OUTPUT
!
    SUBROUTINE OUTPUTD(IND,IDISP,NODE,IDIRN)
!{{{
!----------------------------------------
!   OUTPUT DATA FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IND,IDISP
    INTEGER, INTENT(IN),OPTIONAL :: NODE,IDIRN
!   OUTPUT DISPLACEMENT FOR PLOTTING
    SELECT CASE(IND)
    CASE(1)
      CALL OUTDISPD1(IDISP,NODE,IDIRN)
    CASE(2)
      CALL OUTDISPD2(IDISP)
    CASE DEFAULT
      WRITE(*,*) "ERROR, WRONG INPUT IN OUTPUTD"
      STOP
    END SELECT
    RETURN
!}}}
    END SUBROUTINE OUTPUTD
!
    SUBROUTINE OUTDISPD1(IDISP,NODE,IDIRN)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDISP
    INTEGER, INTENT(IN) :: NODE,IDIRN
!
    INTEGER :: M
    REAL(8) :: DD,VV,AA
!
    M=ID(IDIRN,NODE)
    IF(M.EQ.0) THEN
       DD = 0
       VV = 0
       AA = 0
    ELSE
       DD = U(M)
       VV = V(M)
       AA = ACC(M)
    ENDIF
!
    WRITE(IDISP,2000) TIME,DD,VV,AA
    RETURN
!
 2000 FORMAT(F10.5,5X,3F13.6)
!}}}
    END SUBROUTINE OUTDISPD1
!
    SUBROUTINE OUTDISPD2(IDISP)
!{{{
!----------------------------------------
!   OUTPUT DISPLACEMENT FOR PLOTTING
!----------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDISP
!
    INTEGER :: I
    REAL(8) :: DD,VV,AA
!
    WRITE(IDISP,1000) TIME
    DO I=1,NEQ
       DD=U(I)
       VV=V(I)
       AA=ACC(I)
       WRITE(IDISP,2000) I,DD,VV,AA
    END DO
!
    RETURN
!
 1000 FORMAT(//,"TIME=",F10.6)
 2000 FORMAT(2X,'NEQ',1X,I5,5X,ES13.6,2X,ES13.6,2X,ES13.6)
!}}}
    END SUBROUTINE OUTDISPD2
!
    SUBROUTINE OUTDISP_STRESS(IIN,IOUT,IDISP,ISTRE,ADELE,ADBC)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,ISTRE,ADELE,ADBC 
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
!   TRUSS
    INTEGER :: EleType,NumEle,INDSP
    INTEGER :: BcType,NumBc
    INTEGER :: L1,L2
    INTEGER :: N1,N2
    REAL(8) :: LENGTH,NX,NY
    LOGICAL :: INDBCX,INDBCY,EXPRN1
!   QUADS
!   MFS
    INTEGER :: NEleAss
    INTEGER :: ITYPE,Mtype
    REAL(8) :: E,PR
!   OVLP
    INTEGER :: NumInt,NumOFE,NumCPL
    INTEGER :: NDIV, NOP
!   
    INTEGER :: ALLocateStatus
!
    CALL L_OPEN("DISP_M")
    CALL L_OPEN("STRESS_M")
!
    REWIND(ADELE)
    REWIND(ADBC)
    DO I=1,NEGROUP
       READ(ADELE) (NPAR(J),J=1,12)
       EleType=NPAR(1)
       SELECT CASE(EleType)
       CASE(1)
!--------------------------------------------------
!      TRUSS
!--------------------------------------------------
!{{{
           NumEle = NPAR(3)
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (XYZ(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
           DO J=1,NumEle
             DO K=1,6
                XYZ(K,J)=0
                LM(K,J)=0
             END DO
             ELEPAR(1,J)=0.
           END DO
!--------------------------------------------------
           READ(ADELE) ((XYZ(K,J),K=1,6),J=1,NumEle),((LM(K,J),K=1,6),J=1,NumEle),(ELEPAR(1,J),J=1,NPAR(3))
!--------------------------------------------------
!}}}
           CALL OUTDISP_FEM(IDISP)
           CALL TRUSS_OUT(ISTRE,ADELE,I)
!--------------------------------------------------
           DEALLOCATE(XYZ)
           DEALLOCATE(LM)
           DEALLOCATE(ELEPAR)
       CASE(2)
!--------------------------------------------------
!      QUADS
!--------------------------------------------------
!{{{
!      READ
!--------------------------------------------------
           NumEle = NPAR(4)
           ITYPE  = NPAR(2)
           Mtype  = NPAR(3)
           NDIV=4
!
           E  = MATP(1,Mtype)
           PR = MATP(2,Mtype)
           CALL G_KCS(2,ITYPE,E,PR)           
!--------------------------------------------------
!      ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (XYZ(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(1,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"             
!--------------------------------------------------
           DO J=1,NumEle
              DO K=1,8
                 XYZ(K,J)=0
                 LM(K,J)=0
              END DO
              ELEPAR(1,J)=0.
           END DO  
!--------------------------------------------------
           READ(ADELE) ((XYZ(K,J),K=1,8),J=1,NumEle),((LM(K,J),K=1,8),J=1,NumEle),(ELEPAR(1,K),K=1,NumEle)
!}}}
           CALL OUTPUT_QUADS(NDIV,NumEle)
!--------------------------------------------------
           DEALLOCATE(XYZ)
           DEALLOCATE(LM)
           DEALLOCATE(ELEPAR)
! **************************
           DEALLOCATE(KCS)           
       CASE(3)
!--------------------------------------------------
!      MFS_BAR
!--------------------------------------------------
!{{{
           NEleAss= NPAR(4)
           NumEle = NPAR(3)  
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (NDMFS(3,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(2,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              NDMFS(1,K)=0
              NDMFS(2,K)=0
              NDMFS(3,K)=0
              ELEPAR(1,K)=0.
              ELEPAR(2,K)=0.
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,3),J=1,NumEle),((ELEPAR(J,K),J=1,2),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
           CALL OUTDISP_TRUSS_MFS(IDISP,NEleAss,NumEle)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM)
       CASE(4)
!--------------------------------------------------
!      MFS_QUADS
!--------------------------------------------------
!{{{
           ITYPE = NPAR(2)
           Mtype = NPAR(3)
           E =MATP(1,Mtype)
           PR=MATP(2,Mtype)
           NEleAss= NPAR(5)
           NumEle = NPAR(4)
           CALL G_KCS(2,ITYPE,E,PR)             
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(6,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (LM(NEleAss,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!--------------------------------------------------
           DO K=1,NumEle
              DO J=1,8
                 NDMFS(J,K)=0
              ENDDO
!
              DO J=1,6
                 ELEPAR(J,K)=0.
              ENDDO
!
              DO J=1,NEleAss
                 LM(J,K)=0
              ENDDO
           END DO
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,8),J=1,NumEle),((ELEPAR(J,K),J=1,6),K=1,NumEle),((LM(J,K),J=1,NEleAss),K=1,NumEle)
!}}}
!{{{
           READ(ADBC) (NPAR(J),J=1,7)
           BcType=NPAR(1)
           SELECT CASE(BcType)
           CASE(1)
             NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
             ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
             ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
             ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
             IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
             DO K=1,NumBc
                DO J=1,4
                  BCIND(J,K)=0
                ENDDO
                DO J=1,2
                  BCED(J,K)=0
                  BCPAR(J,K)=0.
                ENDDO
             ENDDO
!
             READ(ADBC) ((BCED(J,K),J=1,2),K=1,NumBc),((BCIND(J,K),J=1,4),K=1,NumBc),((BCPAR(J,K),J=1,2),K=1,NumBc)
!
             IF(PENALTY.EQ.0) THEN   
               DO K=1,NumBc
                EXPRN1=BCIND(1,K).EQ.3
                IF(EXPRN1) CYCLE
!
                L1=BCED(1,K)
                L2=BCED(2,K)
!
                NDMFS(5,L1)=L2
                NDMFS(4,L2)=L1
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
                NDMFS(7,L1)=BCIND(1,K)
                NDMFS(7,L2)=BCIND(1,K)
!
                NDMFS(8,L1)=BCIND(2,K)
                NDMFS(8,L2)=BCIND(2,K)
!
                ELEPAR(5,L1)=-NX
                ELEPAR(6,L1)=-NY
                ELEPAR(5,L2)=-NX
                ELEPAR(6,L2)=-NY
               ENDDO
             ELSE
               DO K=1,NumBc
                  L1=BCED(1,K)
                  L2=BCED(2,K)
!
                  NDMFS(5,L1)=L2
                  NDMFS(4,L2)=L1
               ENDDO
             ENDIF
           CASE DEFAULT
              WRITE(*,*), "ERROR, WRONG INPUT IN GENERATE BOUNDARY PROCEDURE"
              STOP
           END SELECT
           DEALLOCATE(BCED)
           DEALLOCATE(BCIND)
           DEALLOCATE(BCPAR)
!}}}
           CALL OUTPUT_2DMFS_MATLAB(IIN,IOUT,IDISP,NEleAss,NumEle)
!          CALL OUTDISP_2D_MFS(IDISP,NEleAss,NumEle)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM)
           DEALLOCATE(KCS)
       CASE(5)
! COUPLING ELEMENT
           ITYPE = NPAR(2)
           Mtype = NPAR(3)
           E =MATP(1,Mtype)
           PR=MATP(2,Mtype)
           NumEle = NPAR(4)
           CALL G_KCS(2,ITYPE,E,PR)
!--------------------------------------------------
!   ALLOCATE VARIABLES
!--------------------------------------------------
           ALLOCATE (NDMFS(8,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***" 
!--------------------------------------------------
           DO K=1,NumEle
              DO J=1,8
                 NDMFS(J,K)=0
              ENDDO
              DO J=1,4
                 ELEPAR(J,K)=0.
              ENDDO
           END DO
!--------------------------------------------------
           READ(ADELE) ((NDMFS(K,J),K=1,8),J=1,NumEle),((ELEPAR(J,K),J=1,4),K=1,NumEle)  
           CALL OUTPUT_COUP1_MATLAB(IIN,IOUT,IDISP,NumEle)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(KCS) 
       CASE(6)
! New Overlapping ELEMENT1
!{{{
           NumInt = NPAR(3)
           NumEle = NPAR(4)
           NEleAss= NPAR(5)
!
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
           DO J=1,NumInt
              NDINT(1,J)=0
              NDINT(2,J)=0
              INTPAR(1,J)=0.0D0
           END DO
       !
           DO J=1,NumEle
              DO K=1,NEleAss
                 LM(K,J)=0
                 WEIND(K,J)=0.0D0
              END DO
              ELEPAR(1,J)=0.0D0
              NDMFS(1,J)=0
              NDMFS(2,J)=0
           END DO
       !   READ ELEMENT DATA
           READ(ADELE) ((NDMFS(J,K),J=1,2),K=1,NPAR(4)),(ELEPAR(1,J),J=1,NPAR(4)),((LM(J,K),J=1,NPAR(5)),K=1,NPAR(4)),&
           ((WEIND(J,K),J=1,NPAR(5)),K=1,NPAR(4)),((NDINT(J,K),J=1,2),K=1,NPAR(3)),(INTPAR(1,J),J=1,NPAR(3))  
!
           CALL OUTPUT_TRUSS_OVLP(IIN,IOUT,IDISP,NumEle,NEleAss)          
!
           DEALLOCATE(NDINT)
           DEALLOCATE(INTPAR)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM) 
           DEALLOCATE(WEIND)
!}}}
       CASE(7)
! New Overlapping ELEMENT2
!{{{
           NumInt = NPAR(3)
           NumEle = NPAR(4)
           NEleAss= NPAR(5)
!
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
           DO J=1,NumInt
              NDINT(1,J)=0
              NDINT(2,J)=0
              INTPAR(1,J)=0.0D0
           END DO
       !
           DO J=1,NumEle
              DO K=1,NEleAss
                 LM(K,J)=0
              END DO
              ELEPAR(1,J)=0.0D0
              NDMFS(1,J)=0
              NDMFS(2,J)=0
           END DO
       !   READ ELEMENT DATA
           READ(ADELE) ((NDMFS(J,K),J=1,2),K=1,NPAR(4)),(ELEPAR(1,J),J=1,NPAR(4)),((LM(J,K),J=1,NPAR(5)),K=1,NPAR(4)),&
           ((NDINT(J,K),J=1,2),K=1,NPAR(3)),(INTPAR(1,J),J=1,NPAR(3))  
!
           CALL OUTPUT_TRUSS_OVLP2(IIN,IOUT,IDISP,NumEle,NEleAss)          
!
           DEALLOCATE(NDINT)
           DEALLOCATE(INTPAR)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM) 
!}}}
       CASE(9)
!{{{
           ITYPE   = NPAR(7)
           Mtype   = NPAR(2)
           E  = MATP(1,Mtype)
           PR = MATP(2,Mtype)
           NEleAss = NPAR(5)
           NumEle  = NPAR(4)
           NDIV = 10
           NumInt  = NPAR(3)
           CALL G_KCS(2,ITYPE,E,PR)           
!
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
!
           DO J=1,NumInt
              DO K=1,3
                 NDINT(K,J)=0
              ENDDO
              INTPAR(1,J)=0.0D0
           END DO
       !
           DO J=1,NumEle
              NumLM(J)=0
              ELEPAR(1,J)=0.0D0
              ELEPAR(2,J)=0.0D0
              ELEPAR(3,J)=0.0D0
              ELEPAR(4,J)=0.0D0
       !
              DO K=1,NEleAss
                 LM(K,J)=0
              END DO
       !
              DO K=1,6
                 NDMFS(K,J)=0
              ENDDO
           END DO
!   READ ELEMENT DATA
           READ(ADELE) ((NDMFS(J,K),J=1,6),K=1,NPAR(4)),((ELEPAR(J,K),J=1,4),K=1,NPAR(4)),((LM(J,K),J=1,NPAR(5)),K=1,NPAR(4)),&
           ((NDINT(J,K),J=1,3),K=1,NPAR(3)),(INTPAR(1,J),J=1,NPAR(3)),(NumLM(J),J=1,NPAR(4))
!}}}
!{{{
           READ(ADBC) (NPAR(J),J=1,7)
           NumBc  = NPAR(2)
       !   ALLOCATE VARIABLES
           ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
           ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
       !
           READ(ADBC) ((BCED(J,K),J=1,2),K=1,NumBc),((BCIND(J,K),J=1,4),K=1,NumBc),((BCPAR(J,K),J=1,2),K=1,NumBc)
       !
           DO J=1,NumBc
              EXPRN1=BCIND(1,J).EQ.3
              IF(EXPRN1) CYCLE
       !
              L1=BCED(1,J)
              L2=BCED(2,J)
       !
              INDBCX=BCIND(1,J).EQ.1
              INDBCY=BCIND(2,J).EQ.1
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
           DEALLOCATE(BCED)
           DEALLOCATE(BCIND)
           DEALLOCATE(BCPAR)
!}}}
           CALL OUTPUT_2D_OVLP(NDIV,NumInt,NumEle,NEleAss)
!
           DEALLOCATE(NDINT)
           DEALLOCATE(INTPAR)
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(LM) 
           DEALLOCATE(NumLM)
! **************************
           DEALLOCATE(KCS)
       CASE(12)
!{{{
           ITYPE   = NPAR(6)
           Mtype   = NPAR(2)
           E  = MATP(1,Mtype)
           PR = MATP(2,Mtype)
           NumEle  = NPAR(4)
           NumOFE  = NPAR(3)
           NUMCPL  = NPAR(8)
! ---------------------------------------
           NDIV = 10
           NOP  = 66
           CALL G_KCS(2,ITYPE,E,PR)           
! ----- Allocate memory -------
           ALLOCATE (NDMFS(5,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDMFS***"
           ALLOCATE (ELEPAR(4,NumEle),STAT = ALLocateStatus)
           IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ELEPAR***"        
!   READ ELEMENT DATA
!   INITIALIZATION
           READ(ADELE) ((NDMFS(J,K),J=1,5),K=1,NumEle),((ELEPAR(J,K),J=1,4),K=1,NumEle)
!
           IF(OFE1BCONOFF) THEN
              READ(ADBC) (NPAR(J),J=1,7)
              NumBc  = NPAR(2)
!   ALLOCATE VARIABLES
              ALLOCATE (BCED(2,NumBc),STAT = ALLocateStatus)
              IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
              ALLOCATE (BCIND(4,NumBc),STAT = ALLocateStatus)
              IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"
              ALLOCATE (BCPAR(2,NumBc),STAT = ALLocateStatus)
              IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY ***"   
!
              READ(ADBC) ((BCED(J,K),J=1,2),K=1,NumBc),((BCIND(J,K),J=1,4),K=1,NumBc),((BCPAR(J,K),J=1,2),K=1,NumBc)
!
              DO J=1,NumBc
                 L1=BCED(1,J)
                 L2=BCED(2,J)
!
                 EXPRN1=BCIND(1,J).EQ.3
                 IF(EXPRN1) THEN
                    NDMFS(4,L1)=3
                    NDMFS(4,L2)=3
                    NDMFS(5,L1)=3
                    NDMFS(5,L2)=3
                    CYCLE
                 ENDIF
!                     
                 INDBCX=BCIND(1,J).EQ.1
                 INDBCY=BCIND(2,J).EQ.1
!           
                 IF(INDBCX) THEN
                    NDMFS(4,L1)=1
                    NDMFS(4,L2)=1
                 ENDIF
!           
                 IF(INDBCY) THEN
                    NDMFS(5,L1)=1
                    NDMFS(5,L2)=1
                 ENDIF
!           
                 IF(INDBCX.OR.INDBCY) THEN
                    N1=NDMFS(1,L1)
                    N2=NDMFS(1,L2)
!           
                    LENGTH=(X(N1)-X(N2))*(X(N1)-X(N2))+(Y(N1)-Y(N2))*(Y(N1)-Y(N2))
                    LENGTH=SQRT(LENGTH)
!-----------   CALCULATE NORMAL COMPONENTS
                    NX=(Y(N2)-Y(N1))/LENGTH
                    NY=-(X(N2)-X(N1))/LENGTH
                    ELEPAR(3,L1)=-NX
                    ELEPAR(4,L1)=-NY
                    ELEPAR(3,L2)=-NX
                    ELEPAR(4,L2)=-NY
                 ENDIF
              ENDDO
!
              DEALLOCATE(BCED)
              DEALLOCATE(BCIND)
              DEALLOCATE(BCPAR) 
!
           ENDIF           
!}}}
           IF (NumOFE.GT.0) THEN
              ALLOCATE (NDINT(3,NumOFE),STAT = ALLocateStatus)
              IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
              READ(ADELE) ((NDINT(J,K),J=1,3),K=1,NumOFE)
!
              CALL OUTPUT_2D_OFE1(NDIV,NOP,NumOFE)
!
              DEALLOCATE(NDINT)
           ENDIF
!
           IF(NumCPL.GT.0) THEN
              ALLOCATE (CPLINT(4,NumCPL),STAT = ALLocateStatus)
              IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY NDINT***"
!
              READ(ADELE) ((CPLINT(J,K),J=1,4),K=1,NumCPL)
!   Calculate contributions from coupling regions       
              CALL OUTPUT_COUPLE_OFE1(NDIV,NUMCPL)
!**************************
              DEALLOCATE(CPLINT)
           ENDIF   
! **************************
           DEALLOCATE(NDMFS)
           DEALLOCATE(ELEPAR)
           DEALLOCATE(KCS)
       CASE DEFAULT
           WRITE(*,*) "ERROR, WRONG INPUT IN OUTDISP SUBROUTINE!"
           STOP
       END SELECT
   ENDDO
!
   CLOSE(IFILE("DISP_M"))
   CLOSE(IFILE("STRESS_M"))
!}}}
    END SUBROUTINE OUTDISP_STRESS
!
    SUBROUTINE OUTPUT_COUPLE_OFE1(NDIV,NUMCPL)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NDIV,NUMCPL
!   ASSISTANT VARIABLES
    INTEGER :: ALLocateStatus
    INTEGER :: I,J,K,L,NCOUNT
    INTEGER :: NS,NE
    REAL(8),DIMENSION(NDIV*NDIV) :: HI,HJ,HK,HL,HC1
    REAL(8),DIMENSION(NDIV*NDIV) :: dHIdx,dHIdy,dHJdx,dHJdy,dHKdx,dHKdy,dHLdx,dHLdy
    REAL(8),DIMENSION(NDIV*NDIV) :: dHIdxC,dHIdyC,dHJdxC,dHJdyC,dHKdxC,dHKdyC
    REAL(8),DIMENSION(NDIV*NDIV) :: XTMP,YTMP,X_BAR,Y_BAR
    REAL(8),DIMENSION(NDIV*NDIV) :: XTMPC,YTMPC
    REAL(8),DIMENSION(NDIV*NDIV) :: GAUSS_WEI
!
    REAL(8),DIMENSION(NDIV*NDIV,NPOR_MAX) :: dHImdx,dHImdy,dHJmdx,dHJmdy,dHKmdx,dHKmdy,dHLmdx,dHLmdy
    REAL(8),DIMENSION(NDIV*NDIV,NPOR_MAX) :: HIm,HJm,HKm,HLm
!
    INTEGER :: Ele1,Ele2,Ele3,Ele4
    INTEGER :: N1,N2,N3,N4
    INTEGER :: NDOFI,NDOFJ,NDOFK,NDOFL
    REAL(8) :: RADI,RADJ,RADK,RADL
    REAL(8) :: LENGTH,LENGTH2,XM,YM
!
    REAL(8),DIMENSION(NDIV*NDIV,2) :: D
    REAL(8),DIMENSION(NDIV*NDIV,4) :: HU,STRESS
    INTEGER :: IDDIX,IDDIY,IDDJX,IDDJY,IDDKX,IDDKY,IDDLX,IDDLY
    INTEGER :: NUM,NDOFF
!
!    WRITE(IFILE("DISP_M"),1000)
!    WRITE(IFILE("DISP_M"),2000)
!    WRITE(IFILE("STRESS_M"),1000)
!    WRITE(IFILE("STRESS_M"),2010) 
!
    ALLOCATE (XABSC_WEIG_IN(NDIV,2),STAT = ALLocateStatus)
    IF (ALLocateStatus /=0) STOP "*** NOT ENOUGH MEMORY XABSC_WEIG_BC ***"   
!    
    CALL GAULEG(NDIV,XABSC_WEIG_IN(1,1),XABSC_WEIG_IN(1,2))
!
    DO I=1,NDIV
       NS=(I-1)*NDIV+1
       NE=I*NDIV
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
          HJm=0.0D0
       dHJmdx=0.0D0
       dHJmdy=0.0D0
!
          HKm=0.0D0
       dHKmdx=0.0D0
       dHKmdy=0.0D0
!
          HLm=0.0D0
       dHLmdx=0.0D0
       dHLmdy=0.0D0     
!------------------------------------------
!     HIm,dHImdx,dHImdy
!------------------------------------------
      X_BAR=(XTMP-X(N1))/RADI
      Y_BAR=(YTMP-Y(N1))/RADI
      SELECT CASE(NDOFI)
      CASE(1)
         HIm(:,1)=HI*HI
         HJm(:,1)=HI*HJ
         HKm(:,1)=HI*HK
         HLm(:,1)=HI*HL
!
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
         HIm(:,1)=HI
!
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
         HIm(:,1)=HI
!
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
         HIm(:,1)=HI
!
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
!     HJ,dHJmdx,dHJmdy
!---------------------------------------
      X_BAR=(XTMP-X(N2))/RADJ
      Y_BAR=(YTMP-Y(N2))/RADJ
      SELECT CASE(NDOFJ)
      CASE(1)
         HIm(:,1)=HIm(:,1)+HJ*HI
         HJm(:,1)=HJm(:,1)+HJ*HJ
         HKm(:,1)=HKm(:,1)+HJ*HK
         HLm(:,1)=HLm(:,1)+HJ*HL
!
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
         HJm(:,1)=HJm(:,1)+HJ
!
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
         HJm(:,1)=HJm(:,1)+HJ
!
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
         HJm(:,1)=HJm(:,1)+HJ
!
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
!     HKm,dHKmdx,dHKmdy
!---------------------------------------
      X_BAR=(XTMP-X(N3))/RADK
      Y_BAR=(YTMP-Y(N3))/RADK
      SELECT CASE(NDOFK)
      CASE(1)
         HIm(:,1)=HIm(:,1)+HK*HI
         HJm(:,1)=HJm(:,1)+HK*HJ
         HKm(:,1)=HKm(:,1)+HK*HK
         HLm(:,1)=HLm(:,1)+HK*HL
!
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
         HKm(:,1)=HKm(:,1)+HK
!
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
         HKm(:,1)=HKm(:,1)+HK
!
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
         HKm(:,1)=HKm(:,1)+HK
!
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
!     HL,dHLmdx,dHLmdy
!---------------------------------------
      X_BAR=(XTMP-X(N4))/RADL
      Y_BAR=(YTMP-Y(N4))/RADL
      SELECT CASE(NDOFL)
      CASE(1)
         HIm(:,1)=HIm(:,1)+HL*HI
         HJm(:,1)=HJm(:,1)+HL*HJ
         HKm(:,1)=HKm(:,1)+HL*HK
         HLm(:,1)=HLm(:,1)+HL*HL
!
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
         HLm(:,1)=HLm(:,1)+HL
!
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
         HLm(:,1)=HLm(:,1)+HL
!
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
         HLm(:,1)=HLm(:,1)+HL
!
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
!
      IDDIX=ID(1,N1)
      IDDIY=ID(2,N1)
      IDDJX=ID(1,N2)
      IDDJY=ID(2,N2)
      IDDKX=ID(1,N3)
      IDDKY=ID(2,N3)
      IDDLX=ID(1,N4)
      IDDLY=ID(2,N4)
! ------ INITIALIZATION ------
      D=0.0D0
      HU=0.0D0
! ------ N1 ---------------
      IF(IDDIX.NE.0) THEN
         DO J=1,NDOF(1,N1)                   
            NDOFF=J+NDMFS(2,Ele1)
            NUM=IDDIX+J-1
            D(:,1) =D(:,1) +DISP(NUM)*HIm(:,NDOFF)
            HU(:,1)=HU(:,1)+DISP(NUM)*dHImdx(:,NDOFF)
            HU(:,3)=HU(:,3)+DISP(NUM)*dHImdy(:,NDOFF)
         ENDDO
      ENDIF
!
      IF(IDDIY.NE.0) THEN
         DO J=1,NDOF(2,N1)                   
            NDOFF=J+NDMFS(3,Ele1)
            NUM=IDDIY+J-1
            D(:,2) =D(:,2) +DISP(NUM)*HIm(:,NDOFF)
            HU(:,2)=HU(:,2)+DISP(NUM)*dHImdy(:,NDOFF)
            HU(:,4)=HU(:,4)+DISP(NUM)*dHImdx(:,NDOFF)
         ENDDO
      ENDIF
! ------ N2 ---------------       
      IF(IDDJX.NE.0) THEN
         DO J=1,NDOF(1,N2)                   
            NDOFF=J+NDMFS(2,Ele2)
            NUM=IDDJX+J-1
            D(:,1) =D(:,1) +DISP(NUM)*HJm(:,NDOFF)
            HU(:,1)=HU(:,1)+DISP(NUM)*dHJmdx(:,NDOFF)
            HU(:,3)=HU(:,3)+DISP(NUM)*dHJmdy(:,NDOFF)
         ENDDO
      ENDIF
!
      IF(IDDJY.NE.0) THEN
         DO J=1,NDOF(2,N2)                   
            NDOFF=J+NDMFS(3,Ele2)
            NUM=IDDJY+J-1
            D(:,2) =D(:,2) +DISP(NUM)*HJm(:,NDOFF)
            HU(:,2)=HU(:,2)+DISP(NUM)*dHJmdy(:,NDOFF)
            HU(:,4)=HU(:,4)+DISP(NUM)*dHJmdx(:,NDOFF)
         ENDDO
      ENDIF
! ------ N3 ---------------
      IF(IDDKX.NE.0) THEN
         DO J=1,NDOF(1,N3)                   
            NDOFF=J+NDMFS(2,Ele3)
            NUM=IDDKX+J-1
            D(:,1) =D(:,1) +DISP(NUM)*HKm(:,NDOFF)
            HU(:,1)=HU(:,1)+DISP(NUM)*dHKmdx(:,NDOFF)
            HU(:,3)=HU(:,3)+DISP(NUM)*dHKmdy(:,NDOFF)
         ENDDO
      ENDIF
!
      IF(IDDKY.NE.0) THEN
         DO J=1,NDOF(2,N3)                   
            NDOFF=J+NDMFS(3,Ele3)
            NUM=IDDKY+J-1
            D(:,2) =D(:,2) +DISP(NUM)*HKm(:,NDOFF)
            HU(:,2)=HU(:,2)+DISP(NUM)*dHKmdy(:,NDOFF)
            HU(:,4)=HU(:,4)+DISP(NUM)*dHKmdx(:,NDOFF)
         ENDDO
      ENDIF
! ------ N4 ---------------
      IF(IDDLX.NE.0) THEN
         DO J=1,NDOF(1,N4)                   
            NDOFF=J+NDMFS(2,Ele4)
            NUM=IDDLX+J-1
            D(:,1) =D(:,1) +DISP(NUM)*HLm(:,NDOFF)
            HU(:,1)=HU(:,1)+DISP(NUM)*dHLmdx(:,NDOFF)
            HU(:,3)=HU(:,3)+DISP(NUM)*dHLmdy(:,NDOFF)
         ENDDO
      ENDIF
!
      IF(IDDLY.NE.0) THEN
         DO J=1,NDOF(2,N4)                   
            NDOFF=J+NDMFS(3,Ele4)
            NUM=IDDLY+J-1
            D(:,2) =D(:,2) +DISP(NUM)*HLm(:,NDOFF)
            HU(:,2)=HU(:,2)+DISP(NUM)*dHLmdy(:,NDOFF)
            HU(:,4)=HU(:,4)+DISP(NUM)*dHLmdx(:,NDOFF)
         ENDDO
      ENDIF
! ------ CALCULATE STRESS ------
      STRESS(:,1)=KCS(1)*HU(:,1)+KCS(2)*HU(:,2)
      STRESS(:,2)=KCS(1)*HU(:,2)+KCS(2)*HU(:,1)
      STRESS(:,3)=KCS(3)*HU(:,4)+KCS(3)*HU(:,3)
      STRESS(:,4)=STRESS(:,1)*STRESS(:,1)+STRESS(:,2)*STRESS(:,2)-STRESS(:,1)*STRESS(:,2)
      STRESS(:,4)=STRESS(:,4)+3.0D0*STRESS(:,3)*STRESS(:,3)
      STRESS(:,4)=SQRT(STRESS(:,4))
! ------ OUTPUT DISPLACEMENT AND STRESS DATA -----
      DO J=1,NDIV*NDIV
         WRITE(IFILE("DISP_M"),2020) I,XTMP(J),YTMP(J),D(J,1),D(J,2)
         WRITE(IFILE("STRESS_M"),2030) I,XTMP(J),YTMP(J),STRESS(J,1),STRESS(J,2),STRESS(J,3),STRESS(J,4) 
      ENDDO
    ENDDO
!
    DEALLOCATE(XABSC_WEIG_IN) 
! **************************
    RETURN 
!
 1000 FORMAT (' R E S U L T S   O F   T H E   O V E R L A P P I N G   E L E M E N T S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_COUPLE_OFE1
!
    SUBROUTINE OUTPUT_2D_OFE1(NDIV,NOP,NumOFE)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NDIV,NOP,NumOFE
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: LOOP1,LOOP2,NCOUNT
    INTEGER :: Ele1,Ele2,Ele3,N1,N2,N3,NDOFI,NDOFJ,NDOFK
    INTEGER :: INDBCI,INDBCJ,INDBCK
!
    REAL(8),DIMENSION(2) :: XY1,XY2,XY3,XYC1,XYC2,XYC3,XYCM1,XYCM2,XYCM3
    REAL(8)              :: RADI,RADJ,RADK
    REAL(8)              :: L1,L2,L3,L1_HF,L2_HF,L3_HF,LM1,LM2,LM3,AREA2
    REAL(8),DIMENSION(3) :: S1,S2,S3,S4,WI,WJ,WK
    REAL(8),DIMENSION(2) :: PHI1,PHI2,PHI3,PHI4,PHI5,PHI6,PHI7,PHI8,PHI9,CONS1,CONS2,CONS3,CONS4,CONS5,CONS6
    REAL(8),DIMENSION(2) :: C1,C2,C3,C4,C5,C6
    REAL(8),DIMENSION(NOP)   :: XTMP,YTMP,X_BAR,Y_BAR,CONSTX,CONSTY
    REAL(8),DIMENSION(NOP,3) :: AREA_COOR
    REAL(8),DIMENSION(NOP)   :: LILJ,LJLK
!
    REAL(8),DIMENSION(NOP) :: HIW,dHIWdx,dHIWdy,HIWC1,HIWC2,HIWC3
    REAL(8),DIMENSION(NOP) :: HJW,dHJWdx,dHJWdy,HJWC1,HJWC2,HJWC3
    REAL(8),DIMENSION(NOP) :: HKW,dHKWdx,dHKWdy
    REAL(8),DIMENSION(NOP) :: dLILJdx,dLILJdy,dLJLKdx,dLJLKdy
    REAL(8)                :: dLIdx,dLIdy,dLJdx,dLJdy,dLKdx,dLKdy
!
    REAL(8),DIMENSION(NOP,NPOR_MAX) :: HIm,HJm,HKm
    REAL(8),DIMENSION(NOP,NPOR_MAX) :: dHImdx,dHImdy
    REAL(8),DIMENSION(NOP,NPOR_MAX) :: dHJmdx,dHJmdy
    REAL(8),DIMENSION(NOP,NPOR_MAX) :: dHKmdx,dHKmdy
!
    REAL(8) :: NXX,NYY
!
    REAL(8),DIMENSION(NOP,2) :: D
    REAL(8),DIMENSION(NOP,4) :: HU,STRESS
    INTEGER                  :: IDDIX,IDDIY,IDDJX,IDDJY,IDDKX,IDDKY
    INTEGER                  :: NUM,NDOFF
!
    WRITE(IFILE("DISP_M"),1000)
    WRITE(IFILE("DISP_M"),2000)
    WRITE(IFILE("STRESS_M"),1000)
    WRITE(IFILE("STRESS_M"),2010) 
!
    NCOUNT=0
    LOOP1=NDIV+1
    DO I=1,LOOP1
       LOOP2=NDIV+2-I
       DO J=1,LOOP2
          NCOUNT=NCOUNT+1
          AREA_COOR(NCOUNT,1)=(I-1.0D0)/NDIV
          AREA_COOR(NCOUNT,2)=(J-1.0D0)/NDIV
          AREA_COOR(NCOUNT,3)=1.0D0-AREA_COOR(NCOUNT,1)-AREA_COOR(NCOUNT,2)
       ENDDO
    ENDDO
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
!      HIWC2=C4(1)+AREA_COOR(:,1)*(C5(1)+AREA_COOR(:,1)*C6(1))
       HIWC2=C4(1)+AREA_COOR(:,1)*C5(1)
       HIW  =LILJ*HIWC1+AREA_COOR(:,1)*HIWC2 
!
       HJWC1=AREA_COOR(:,3)*C1(2)+AREA_COOR(:,2)*C2(2)+C3(2)
!      HJWC2=C4(2)+AREA_COOR(:,2)*(C5(2)+AREA_COOR(:,2)*C6(2))
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
!      HIWC3=AREA_COOR(:,1)*(C5(1)+AREA_COOR(:,1)*(2.0D0*C6(1)))
!      HJWC3=AREA_COOR(:,2)*(C5(2)+AREA_COOR(:,2)*(2.0D0*C6(2)))
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
!       DO J=1,NOP
!          HIW(J)=AREA_COOR(J,1)
!          HJW(J)=AREA_COOR(J,2)
!          HKW(J)=AREA_COOR(J,3)
!          dHIWdx(J)=dLIdx
!          dHIWdy(J)=dLIdy
!          dHJWdx(J)=dLJdx
!          dHJWdy(J)=dLJdy
!          dHKWdx(J)=dLKdx
!          dHKWdy(J)=dLKdy
!       ENDDO
!
          HJm=0.0D0
       dHJmdx=0.0D0
       dHJmdy=0.0D0
!
          HKm=0.0D0
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
             HIm(:,1)=HIW*AREA_COOR(:,1)
             HJm(:,1)=HIW*AREA_COOR(:,2)
             HKm(:,1)=HIW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHIWdx*AREA_COOR(:,1)+HIW*dLIdx
             dHImdy(:,1)=dHIWdy*AREA_COOR(:,1)+HIW*dLIdy
!
             dHJmdx(:,1)=dHIWdx*AREA_COOR(:,2)+HIW*dLJdx
             dHJmdy(:,1)=dHIWdy*AREA_COOR(:,2)+HIW*dLJdy
!
             dHKmdx(:,1)=dHIWdx*AREA_COOR(:,3)+HIW*dLKdx
             dHKmdy(:,1)=dHIWdy*AREA_COOR(:,3)+HIW*dLKdy
          CASE(2)
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
          CASE(3)
             HIWC1=HIW*X_BAR
!
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
!
             HIm(:,3)=HIm(:,2)*Y_BAR
             dHImdx(:,3)=dHImdx(:,2)*Y_BAR+HIWC1*(-NYY/RADI)
             dHImdy(:,3)=dHImdy(:,2)*Y_BAR+HIWC1*(NXX/RADI)
          CASE(4)
             HIWC1=HIW*X_BAR
!
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIW*(NXX/RADI)
             dHImdy(:,2)=dHIWdy*X_BAR+HIW*(NYY/RADI)
!
             HIm(:,3)=HIm(:,2)*Y_BAR
             dHImdx(:,3)=dHImdx(:,2)*Y_BAR+HIWC1*(-NYY/RADI)
             dHImdy(:,3)=dHImdy(:,2)*Y_BAR+HIWC1*(NXX/RADI)
!
             HIm(:,4)=HIm(:,2)*X_BAR
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
             HIm(:,1)=HIW*AREA_COOR(:,1)
             HJm(:,1)=HIW*AREA_COOR(:,2)
             HKm(:,1)=HIW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHIWdx*AREA_COOR(:,1)+HIW*dLIdx
             dHImdy(:,1)=dHIWdy*AREA_COOR(:,1)+HIW*dLIdy
!
             dHJmdx(:,1)=dHIWdx*AREA_COOR(:,2)+HIW*dLJdx
             dHJmdy(:,1)=dHIWdy*AREA_COOR(:,2)+HIW*dLJdy
!
             dHKmdx(:,1)=dHIWdx*AREA_COOR(:,3)+HIW*dLKdx
             dHKmdy(:,1)=dHIWdy*AREA_COOR(:,3)+HIW*dLKdy
!             HIm(:,1)=HIW
!             dHImdx(:,1)=dHIWdx
!             dHImdy(:,1)=dHIWdy
          CASE(3)
             HIWC1=HIW/RADI
!
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             HIm(:,3)=HIW*Y_BAR
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1           
          CASE(4)
             HIWC1=HIW/RADI
!
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             HIm(:,3)=HIW*Y_BAR
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1
!
             HIm(:,4)=HIm(:,2)*Y_BAR
             dHImdx(:,4)=dHImdx(:,2)*Y_BAR
             dHImdy(:,4)=dHImdy(:,3)*X_BAR          
          CASE(6)
             HIWC1=HIW/RADI
!
             HIm(:,1)=HIW
             dHImdx(:,1)=dHIWdx
             dHImdy(:,1)=dHIWdy
!
             HIm(:,2)=HIW*X_BAR
             dHImdx(:,2)=dHIWdx*X_BAR+HIWC1
             dHImdy(:,2)=dHIWdy*X_BAR
!
             HIm(:,3)=HIW*Y_BAR
             dHImdx(:,3)=dHIWdx*Y_BAR
             dHImdy(:,3)=dHIWdy*Y_BAR+HIWC1      
!
             HIm(:,4)=HIm(:,2)*Y_BAR
             dHImdx(:,4)=dHImdx(:,2)*Y_BAR
             dHImdy(:,4)=dHImdy(:,3)*X_BAR
!
             HIm(:,5)=HIm(:,2)*X_BAR
             dHImdx(:,5)=(dHImdx(:,2)+HIWC1)*X_BAR
             dHImdy(:,5)=dHImdy(:,2)*X_BAR
!
             HIm(:,6)=HIm(:,3)*Y_BAR
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
             HIm(:,1)=HIm(:,1)+HJW*AREA_COOR(:,1)
             HJm(:,1)=HJm(:,1)+HJW*AREA_COOR(:,2)
             HKm(:,1)=HKm(:,1)+HJW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHImdx(:,1)+dHJWdx*AREA_COOR(:,1)+HJW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHJWdy*AREA_COOR(:,1)+HJW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx*AREA_COOR(:,2)+HJW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy*AREA_COOR(:,2)+HJW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHJWdx*AREA_COOR(:,3)+HJW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHJWdy*AREA_COOR(:,3)+HJW*dLKdy
          CASE(2)
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
          CASE(3)
             HJWC1=HJW*X_BAR
!
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
!
             HJm(:,3)=HJm(:,2)*Y_BAR
             dHJmdx(:,3)=dHJmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADJ)
             dHJmdy(:,3)=dHJmdy(:,2)*Y_BAR+HJWC1*(NXX/RADJ)
          CASE(4)
             HJWC1=HJW*X_BAR
!
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJW*(NXX/RADJ)
             dHJmdy(:,2)=dHJWdy*X_BAR+HJW*(NYY/RADJ)
!
             HJm(:,3)=HJm(:,2)*Y_BAR
             dHJmdx(:,3)=dHJmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADJ)
             dHJmdy(:,3)=dHJmdy(:,2)*Y_BAR+HJWC1*(NXX/RADJ)
!
             HJm(:,4)=HJm(:,2)*X_BAR
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
             HIm(:,1)=HIm(:,1)+HJW*AREA_COOR(:,1)
             HJm(:,1)=HJm(:,1)+HJW*AREA_COOR(:,2)
             HKm(:,1)=HKm(:,1)+HJW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHImdx(:,1)+dHJWdx*AREA_COOR(:,1)+HJW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHJWdy*AREA_COOR(:,1)+HJW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx*AREA_COOR(:,2)+HJW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy*AREA_COOR(:,2)+HJW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHJWdx*AREA_COOR(:,3)+HJW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHJWdy*AREA_COOR(:,3)+HJW*dLKdy
!             HJm(:,1)=HJW
!             dHJmdx(:,1)=dHJWdx
!             dHJmdy(:,1)=dHJWdy
          CASE(3)
             HJWC1=HJW/RADJ
!
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             HJm(:,3)=HJW*Y_BAR
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1           
          CASE(4)
             HJWC1=HJW/RADJ
!
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             HJm(:,3)=HJW*Y_BAR
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1 
!
             HJm(:,4)=HJm(:,2)*Y_BAR
             dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
             dHJmdy(:,4)=dHJmdy(:,3)*X_BAR          
          CASE(6)
             HJWC1=HJW/RADJ
!
             HJm(:,1)=HJm(:,1)+HJW
             dHJmdx(:,1)=dHJmdx(:,1)+dHJWdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHJWdy
!
             HJm(:,2)=HJW*X_BAR
             dHJmdx(:,2)=dHJWdx*X_BAR+HJWC1
             dHJmdy(:,2)=dHJWdy*X_BAR
!
             HJm(:,3)=HJW*Y_BAR
             dHJmdx(:,3)=dHJWdx*Y_BAR
             dHJmdy(:,3)=dHJWdy*Y_BAR+HJWC1       
!
             HJm(:,4)=HJm(:,2)*Y_BAR
             dHJmdx(:,4)=dHJmdx(:,2)*Y_BAR
             dHJmdy(:,4)=dHJmdy(:,3)*X_BAR
!
             HJm(:,5)=HJm(:,2)*X_BAR
             dHJmdx(:,5)=(dHJmdx(:,2)+HJWC1)*X_BAR
             dHJmdy(:,5)=dHJmdy(:,2)*X_BAR
!
             HJm(:,6)=HJm(:,3)*Y_BAR
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
          SELECT CASE(NDOFK)
          CASE(1)
             HIm(:,1)=HIm(:,1)+HKW*AREA_COOR(:,1)
             HJm(:,1)=HJm(:,1)+HKW*AREA_COOR(:,2)
             HKm(:,1)=HKm(:,1)+HKW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHImdx(:,1)+dHKWdx*AREA_COOR(:,1)+HKW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHKWdy*AREA_COOR(:,1)+HKW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHKWdx*AREA_COOR(:,2)+HKW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHKWdy*AREA_COOR(:,2)+HKW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx*AREA_COOR(:,3)+HKW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy*AREA_COOR(:,3)+HKW*dLKdy
          CASE(2)
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
          CASE(3)
             HJWC1=HKW*X_BAR
!
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
!
             HKm(:,3)=HKm(:,2)*Y_BAR
             dHKmdx(:,3)=dHKmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADK)
             dHKmdy(:,3)=dHKmdy(:,2)*Y_BAR+HJWC1*(NXX/RADK)
          CASE(4)
             HJWC1=HKW*X_BAR
!
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HKW*(NXX/RADK)
             dHKmdy(:,2)=dHKWdy*X_BAR+HKW*(NYY/RADK)
!
             HKm(:,3)=HKm(:,2)*Y_BAR
             dHKmdx(:,3)=dHKmdx(:,2)*Y_BAR+HJWC1*(-NYY/RADK)
             dHKmdy(:,3)=dHKmdy(:,2)*Y_BAR+HJWC1*(NXX/RADK)
!
             HKm(:,4)=HKm(:,2)*X_BAR
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
             HIm(:,1)=HIm(:,1)+HKW*AREA_COOR(:,1)
             HJm(:,1)=HJm(:,1)+HKW*AREA_COOR(:,2)
             HKm(:,1)=HKm(:,1)+HKW*AREA_COOR(:,3)
!
             dHImdx(:,1)=dHImdx(:,1)+dHKWdx*AREA_COOR(:,1)+HKW*dLIdx
             dHImdy(:,1)=dHImdy(:,1)+dHKWdy*AREA_COOR(:,1)+HKW*dLIdy
!
             dHJmdx(:,1)=dHJmdx(:,1)+dHKWdx*AREA_COOR(:,2)+HKW*dLJdx
             dHJmdy(:,1)=dHJmdy(:,1)+dHKWdy*AREA_COOR(:,2)+HKW*dLJdy
!
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx*AREA_COOR(:,3)+HKW*dLKdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy*AREA_COOR(:,3)+HKW*dLKdy
!             HKm(:,1)=HKW
!             dHKmdx(:,1)=dHKWdx
!             dHKmdy(:,1)=dHKWdy
          CASE(3)
             HJWC1=HKW/RADK
!
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             HKm(:,3)=HKW*Y_BAR
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1           
          CASE(4)
             HJWC1=HKW/RADK
!
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!       
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             HKm(:,3)=HKW*Y_BAR
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1  
!
             HKm(:,4)=HKm(:,2)*Y_BAR
             dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
             dHKmdy(:,4)=dHKmdy(:,3)*X_BAR          
          CASE(6)
             HJWC1=HKW/RADK
!
             HKm(:,1)=HKm(:,1)+HKW
             dHKmdx(:,1)=dHKmdx(:,1)+dHKWdx
             dHKmdy(:,1)=dHKmdy(:,1)+dHKWdy
!
             HKm(:,2)=HKW*X_BAR
             dHKmdx(:,2)=dHKWdx*X_BAR+HJWC1
             dHKmdy(:,2)=dHKWdy*X_BAR
!
             HKm(:,3)=HKW*Y_BAR
             dHKmdx(:,3)=dHKWdx*Y_BAR
             dHKmdy(:,3)=dHKWdy*Y_BAR+HJWC1
!
             HKm(:,4)=HKm(:,2)*Y_BAR
             dHKmdx(:,4)=dHKmdx(:,2)*Y_BAR
             dHKmdy(:,4)=dHKmdy(:,3)*X_BAR
!
             HKm(:,5)=HKm(:,2)*X_BAR
             dHKmdx(:,5)=(dHKmdx(:,2)+HJWC1)*X_BAR
             dHKmdy(:,5)=dHKmdy(:,2)*X_BAR
!
             HKm(:,6)=HKm(:,3)*Y_BAR
             dHKmdx(:,6)=dHKmdx(:,3)*Y_BAR
             dHKmdy(:,6)=(dHKmdy(:,3)+HJWC1)*Y_BAR 
          CASE DEFAULT
             WRITE(*,*) "ERROR, GENERATEM/OFE1_QUADS_STIFF!"
             STOP
          END SELECT 
       ENDIF
!
       IDDIX=ID(1,N1)
       IDDIY=ID(2,N1)
       IDDJX=ID(1,N2)
       IDDJY=ID(2,N2)
       IDDKX=ID(1,N3)
       IDDKY=ID(2,N3)
! ------ INITIALIZATION ------
       D=0.0D0
       HU=0.0D0
! ------ N1 ---------------
       IF(IDDIX.NE.0) THEN
          DO J=1,NDOF(1,N1)                   
             NDOFF=J+NDMFS(2,Ele1)
             NUM=IDDIX+J-1
             D(:,1) =D(:,1) +DISP(NUM)*HIm(:,NDOFF)
             HU(:,1)=HU(:,1)+DISP(NUM)*dHImdx(:,NDOFF)
             HU(:,3)=HU(:,3)+DISP(NUM)*dHImdy(:,NDOFF)
          ENDDO
       ENDIF
!
       IF(IDDIY.NE.0) THEN
          DO J=1,NDOF(2,N1)                   
             NDOFF=J+NDMFS(3,Ele1)
             NUM=IDDIY+J-1
             D(:,2) =D(:,2) +DISP(NUM)*HIm(:,NDOFF)
             HU(:,2)=HU(:,2)+DISP(NUM)*dHImdy(:,NDOFF)
             HU(:,4)=HU(:,4)+DISP(NUM)*dHImdx(:,NDOFF)
          ENDDO
       ENDIF
! ------ N2 ---------------       
       IF(IDDJX.NE.0) THEN
          DO J=1,NDOF(1,N2)                   
             NDOFF=J+NDMFS(2,Ele2)
             NUM=IDDJX+J-1
             D(:,1) =D(:,1) +DISP(NUM)*HJm(:,NDOFF)
             HU(:,1)=HU(:,1)+DISP(NUM)*dHJmdx(:,NDOFF)
             HU(:,3)=HU(:,3)+DISP(NUM)*dHJmdy(:,NDOFF)
          ENDDO
       ENDIF
!
       IF(IDDJY.NE.0) THEN
          DO J=1,NDOF(2,N2)                   
             NDOFF=J+NDMFS(3,Ele2)
             NUM=IDDJY+J-1
             D(:,2) =D(:,2) +DISP(NUM)*HJm(:,NDOFF)
             HU(:,2)=HU(:,2)+DISP(NUM)*dHJmdy(:,NDOFF)
             HU(:,4)=HU(:,4)+DISP(NUM)*dHJmdx(:,NDOFF)
          ENDDO
       ENDIF
! ------ N3 ---------------
       IF(IDDKX.NE.0) THEN
          DO J=1,NDOF(1,N3)                   
             NDOFF=J+NDMFS(2,Ele3)
             NUM=IDDKX+J-1
             D(:,1) =D(:,1) +DISP(NUM)*HKm(:,NDOFF)
             HU(:,1)=HU(:,1)+DISP(NUM)*dHKmdx(:,NDOFF)
             HU(:,3)=HU(:,3)+DISP(NUM)*dHKmdy(:,NDOFF)
          ENDDO
       ENDIF
!
       IF(IDDKY.NE.0) THEN
          DO J=1,NDOF(2,N3)                   
             NDOFF=J+NDMFS(3,Ele3)
             NUM=IDDKY+J-1
             D(:,2) =D(:,2) +DISP(NUM)*HKm(:,NDOFF)
             HU(:,2)=HU(:,2)+DISP(NUM)*dHKmdy(:,NDOFF)
             HU(:,4)=HU(:,4)+DISP(NUM)*dHKmdx(:,NDOFF)
          ENDDO
       ENDIF
! ------ CALCULATE STRESS ------
       STRESS(:,1)=KCS(1)*HU(:,1)+KCS(2)*HU(:,2)
       STRESS(:,2)=KCS(1)*HU(:,2)+KCS(2)*HU(:,1)
       STRESS(:,3)=KCS(3)*HU(:,4)+KCS(3)*HU(:,3)
       STRESS(:,4)=STRESS(:,1)*STRESS(:,1)+STRESS(:,2)*STRESS(:,2)-STRESS(:,1)*STRESS(:,2)
       STRESS(:,4)=STRESS(:,4)+3.0D0*STRESS(:,3)*STRESS(:,3)
       STRESS(:,4)=SQRT(STRESS(:,4))
! ------ OUTPUT DISPLACEMENT AND STRESS DATA -----
       DO J=1,NOP
          WRITE(IFILE("DISP_M"),2020) I,XTMP(J),YTMP(J),D(J,1),D(J,2)
          WRITE(IFILE("STRESS_M"),2030) I,XTMP(J),YTMP(J),STRESS(J,1),STRESS(J,2),STRESS(J,3),STRESS(J,4) 
       ENDDO
!
    ENDDO
! **************************
    RETURN 
!
 1000 FORMAT (' R E S U L T S   O F   T H E   O V E R L A P P I N G   E L E M E N T S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_2D_OFE1
!
    SUBROUTINE OUTPUT_2D_OVLP(NDIV,NumInt,NumEle,NEleAss)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NDIV,NumInt,NumEle,NEleAss
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K,L,I1,I2
    INTEGER :: ELE1,ELE2,ELE3,N1,N2,N3,ELEE,NN,LL
    REAL(8) :: X1,Y1,X2,Y2,X3,Y3,AREA2
    INTEGER :: LOOP1,LOOP2
    REAL(8) :: DAXYZ(2,3),AXYZ,DAXYZ_X,DAXYZ_Y,XI(3)
    REAL(8) :: XTMP,YTMP
    REAL(8) :: HH(6),DHH(6,2)
    REAL(8) :: D(2),HU(4),STRESS(4),WEI(NEleAss),DWEI(NEleAss,2)
    INTEGER :: N,M,ELE,IND,INDN,INDSP
    REAL(8) :: X_BAR,Y_BAR,CONSTX,CONSTY,RAD,NX,NY
    INTEGER :: NDOFF,NUM
    INTEGER :: IDD_X,IDD_Y,IND_X,IND_Y,INC
!
    WRITE(IFILE("DISP_M"),1000)
    WRITE(IFILE("DISP_M"),2000)
    WRITE(IFILE("STRESS_M"),1000)
    WRITE(IFILE("STRESS_M"),2010) 
!
    DO I=1,NumInt
       ELE1=NDINT(1,I)
       ELE2=NDINT(2,I)
       ELE3=NDINT(3,I)
! -----------------------------------
       N1=NDMFS(1,ELE1)
       N2=NDMFS(1,ELE2)
       N3=NDMFS(1,ELE3)
!
       X1=X(N1)
       Y1=Y(N1)
       X2=X(N2)
       Y2=Y(N2)
       X3=X(N3)
       Y3=Y(N3)
!
       AREA2=X2*Y3-X3*Y2+X3*Y1-X1*Y3+X1*Y2-X2*Y1
!
       DAXYZ(1,1)=(Y2-Y3)/AREA2
       DAXYZ(2,1)=(X3-X2)/AREA2
       DAXYZ(1,2)=(Y3-Y1)/AREA2
       DAXYZ(2,2)=(X1-X3)/AREA2
       DAXYZ(1,3)=(Y1-Y2)/AREA2
       DAXYZ(2,3)=(X2-X1)/AREA2
!
       LOOP1=NDIV+1
       DO I1=1,LOOP1
          XI(1)=(I1-1.0D0)/NDIV
          LOOP2=NDIV+2-I1
          DO I2=1,LOOP2
             XI(2)=(I2-1.0D0)/NDIV
             XI(3)=1.0D0-XI(1)-XI(2)
!
             D(1)=0.0D0
             D(2)=0.0D0
!
             DO J=1,4
                HU(J)=0.0D0
             ENDDO
!
             XTMP=XI(1)*X1+XI(2)*X2+XI(3)*X3
             YTMP=XI(1)*Y1+XI(2)*Y2+XI(3)*Y3
             DO J=1,3
                AXYZ=XI(J)
                DAXYZ_X=DAXYZ(1,J)
                DAXYZ_Y=DAXYZ(2,J)
!
                ELE=NDINT(J,I)
                N=NDMFS(1,ELE)
                IND=NDMFS(2,ELE)
                INDSP=NDMFS(4,ELE)
                IF ((IND.NE.0).OR.(INDSP.EQ.2)) THEN
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
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
                      CASE(2)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
! HIX,DHIXDX,DHIXDY
                       HH(2)=HH(1)*X_BAR
                       DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                       DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
                      CASE(3)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
! HIX,DHIXDX,DHIXDY
                       HH(2)=HH(1)*X_BAR
                       DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                       DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
! HIXY,DHIXYDX,DHIXYDY
                       HH(3)=HH(2)*Y_BAR
                       DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                       DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
                     CASE(4)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
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
                     CASE DEFAULT
                        WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                        STOP
                     END SELECT
!}}}
                   ELSE
!{{{
                    IF(IND.EQ.3) THEN
                      DO L=1,3
                         ELEE=NDINT(L,I)
                         NN=NDMFS(1,ELEE)

                         HH(1)=AXYZ*XI(L)
                         DHH(1,1)=DAXYZ_X*XI(L)+AXYZ*DAXYZ(1,L)
                         DHH(1,2)=DAXYZ_Y*XI(L)+AXYZ*DAXYZ(2,L)
!
                         IND_X=NDMFS(2,ELEE)
                         IND_Y=NDMFS(3,ELEE)
!
                         IF(IND_X.EQ.1) THEN
                            INC=1
                         ELSE
                            INC=0
                         ENDIF
!
                         IF(INC.EQ.0) THEN                  
                             NUM=ID(1,NN)
                             D(1)=D(1)+DISP(NUM)*HH(1)
                             HU(1)=HU(1)+DISP(NUM)*DHH(1,1)
                             HU(3)=HU(3)+DISP(NUM)*DHH(1,2)
                         ENDIF
!
                         IF(IND_Y.EQ.1) THEN
                            INC=1
                         ELSE
                            INC=0
                         ENDIF
!
                         IF(INC.EQ.0) THEN                  
                             NUM=ID(2,NN)
                             D(2)=D(2)+DISP(NUM)*HH(1)
                             HU(2)=HU(2)+DISP(NUM)*DHH(1,2)
                             HU(4)=HU(4)+DISP(NUM)*DHH(1,1)
                         ENDIF
                      ENDDO
                    ELSE
                      SELECT CASE(NDOFF)
                      CASE(1)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
                      CASE(3)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
! HIX,DHIXDX,DHIXDY
                       HH(2)=HH(1)*CONSTX
                       DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                       DHH(2,2)=DHH(1,2)*CONSTX
! HIY,DHIYDX,DHIYDY
                       HH(3)=HH(1)*CONSTY
                       DHH(3,1)=DHH(1,1)*CONSTY
                       DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD             
                      CASE(4)
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
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
                       HH(1)=AXYZ
                       DHH(1,1)=DAXYZ_X
                       DHH(1,2)=DAXYZ_Y
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
                      CASE DEFAULT
                        WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                        STOP
                      END SELECT
!}}}
                    ENDIF
                   ENDIF
!
                   IF(IND.NE.3) THEN
                      IDD_X=ID(1,N)
                      IDD_Y=ID(2,N)
   !
                      IND_X=NDMFS(2,ELE)
                      IND_Y=NDMFS(3,ELE)
                      IF(IND_X.EQ.1) THEN
                         INC=1
                      ELSE
                         INC=0
                      ENDIF
                      IF(IDD_X.NE.0) THEN
                         DO K=1,NDOF(1,N)                   
                            NDOFF=K+INC
                            NUM=IDD_X+K-1
                            D(1)=D(1)+DISP(NUM)*HH(NDOFF)
                            HU(1)=HU(1)+DISP(NUM)*DHH(NDOFF,1)
                            HU(3)=HU(3)+DISP(NUM)*DHH(NDOFF,2)
                         ENDDO
                      ENDIF
         !
                      IF(IND_Y.EQ.1) THEN
                         INC=1
                      ELSE
                         INC=0
                      ENDIF
                      IF(IDD_Y.NE.0) THEN
                         DO K=1,NDOF(2,N)                   
                            NDOFF=K+INC
                            NUM=IDD_Y+K-1
                            D(2)=D(2)+DISP(NUM)*HH(NDOFF)
                            HU(2)=HU(2)+DISP(NUM)*DHH(NDOFF,2)
                            HU(4)=HU(4)+DISP(NUM)*DHH(NDOFF,1)
                         ENDDO
                      ENDIF
                   ENDIF
                ELSE
                   CALL WEI_CAL_QUADS(ELE,NEleAss,XTMP,YTMP,WEI(1),DWEI(1,1))
                   DO L=1,NEleAss
                      ELE1=LM(L,ELE)
                      IF(ELE1.NE.0) THEN
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
                                HH(1)=AXYZ*WEI(L)
                                DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
                             CASE(2)
                                HH(1)=AXYZ*WEI(L)
                                DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
! HIX,DHIXDX,DHIXDY
                                HH(2)=HH(1)*X_BAR
                                DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                                DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
                             CASE(3)
                                HH(1)=AXYZ*WEI(L)
                                DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
! HIX,DHIXDX,DHIXDY
                                HH(2)=HH(1)*X_BAR
                                DHH(2,1)=DHH(1,1)*X_BAR+HH(1)*NX/RAD
                                DHH(2,2)=DHH(1,2)*X_BAR+HH(1)*NY/RAD
! HIXY,DHIXYDX,DHIXYDY
                                HH(3)=HH(2)*Y_BAR
                                DHH(3,1)=DHH(2,1)*Y_BAR+HH(2)*(-NY)/RAD
                                DHH(3,2)=DHH(2,2)*Y_BAR+HH(2)*NX/RAD  
                            CASE(4)
                                HH(1)=AXYZ*WEI(L)
                                DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
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
                             CASE DEFAULT
                               WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                               STOP
                             END SELECT
!}}}       
                         ELSE
!{{{
                            IF(IND.EQ.3) THEN
                              DO LL=1,3
                                 ELEE=NDINT(LL,I)
                                 NN=NDMFS(1,ELEE)

                                 HH(1)=AXYZ*XI(LL)*WEI(L)
                                 DHH(1,1)=(DAXYZ_X*XI(LL)+AXYZ*DAXYZ(1,LL))*WEI(L)+AXYZ*XI(LL)*DWEI(L,1)
                                 DHH(1,2)=(DAXYZ_Y*XI(LL)+AXYZ*DAXYZ(2,LL))*WEI(L)+AXYZ*XI(LL)*DWEI(L,2)
!
                                 IND_X=NDMFS(2,ELEE)
                                 IND_Y=NDMFS(3,ELEE)
                                 IF(IND_X.EQ.1) THEN
                                    INC=1
                                 ELSE
                                    INC=0
                                 ENDIF
                                 IF(INC.EQ.0) THEN                 
                                    NUM=ID(1,NN)
                                    D(1)=D(1)+DISP(NUM)*HH(1)
                                    HU(1)=HU(1)+DISP(NUM)*DHH(1,1)
                                    HU(3)=HU(3)+DISP(NUM)*DHH(1,2)
                                 ENDIF
        !
                                 IF(IND_Y.EQ.1) THEN
                                    INC=1
                                 ELSE
                                    INC=0
                                 ENDIF
                                 IF(INC.EQ.0) THEN                  
                                    NUM=ID(2,NN)
                                    D(2)=D(2)+DISP(NUM)*HH(1)
                                    HU(2)=HU(2)+DISP(NUM)*DHH(1,2)
                                    HU(4)=HU(4)+DISP(NUM)*DHH(1,1)
                                 ENDIF
                              ENDDO
                            ELSE
                              SELECT CASE(NDOFF)
                              CASE(1)
                                 HH(1)=AXYZ*WEI(L)
                                 DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                 DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
                              CASE(3)
                                 HH(1)=AXYZ*WEI(L)
                                 DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                 DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
  ! HIX,DHIXDX,DHIXDY
                                 HH(2)=HH(1)*CONSTX
                                 DHH(2,1)=DHH(1,1)*CONSTX+HH(1)/RAD
                                 DHH(2,2)=DHH(1,2)*CONSTX
  ! HIY,DHIYDX,DHIYDY
                                 HH(3)=HH(1)*CONSTY
                                 DHH(3,1)=DHH(1,1)*CONSTY
                                 DHH(3,2)=DHH(1,2)*CONSTY+HH(1)/RAD  
                              CASE(4)
                                 HH(1)=AXYZ*WEI(L)
                                 DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                 DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
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
                                 HH(1)=AXYZ*WEI(L)
                                 DHH(1,1)=DAXYZ_X*WEI(L)+AXYZ*DWEI(L,1)
                                 DHH(1,2)=DAXYZ_Y*WEI(L)+AXYZ*DWEI(L,2)
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
                              CASE DEFAULT
                                  WRITE(*,*) "ERROR, /GENERATEM/DH_QUADS!"
                                  STOP
                              END SELECT
!}}}
                            ENDIF
                         ENDIF

                         IF(IND.NE.3) THEN
                            IDD_X=ID(1,M)
                            IDD_Y=ID(2,M)
         !
                            IND_X=NDMFS(2,ELE1)
                            IND_Y=NDMFS(3,ELE1)
                            IF(IND_X.EQ.1) THEN
                               INC=1
                            ELSE
                               INC=0
                            ENDIF
                            IF(IDD_X.NE.0) THEN
                               DO K=1,NDOF(1,M)                   
                                  NDOFF=K+INC
                                  NUM=IDD_X+K-1
                                  D(1)=D(1)+DISP(NUM)*HH(NDOFF)
                                  HU(1)=HU(1)+DISP(NUM)*DHH(NDOFF,1)
                                  HU(3)=HU(3)+DISP(NUM)*DHH(NDOFF,2)
                               ENDDO
                            ENDIF
               !
                            IF(IND_Y.EQ.1) THEN
                               INC=1
                            ELSE
                               INC=0
                            ENDIF
                            IF(IDD_Y.NE.0) THEN
                               DO K=1,NDOF(2,M)                   
                                  NDOFF=K+INC
                                  NUM=IDD_Y+K-1
                                  D(2)=D(2)+DISP(NUM)*HH(NDOFF)
                                  HU(2)=HU(2)+DISP(NUM)*DHH(NDOFF,2)
                                  HU(4)=HU(4)+DISP(NUM)*DHH(NDOFF,1)
                               ENDDO
                            ENDIF
!
                         ENDIF
                      ENDIF
                    ENDDO
!
                ENDIF
!
             ENDDO
!
             IF  ((abs(XTMP-100.0D0).LT.1E-12) .AND. (abs(YTMP-10.0D0).LT.1E-12)) THEN
                 WRITE(*,*) XTMP,YTMP,D(2)
             ENDIF
!
             WRITE(IFILE("DISP_M"),2020) I,XTMP,YTMP,D(1),D(2)
             STRESS(1)=KCS(1)*HU(1)+ KCS(2)*HU(2)
             STRESS(2)=KCS(1)*HU(2)+ KCS(2)*HU(1)
             STRESS(3)=KCS(3)*HU(4)+ KCS(3)*HU(3)
             STRESS(4)=STRESS(1)*STRESS(1)+STRESS(2)*STRESS(2)-STRESS(1)*STRESS(2)
             STRESS(4)=STRESS(4)+3.0D0*STRESS(3)*STRESS(3)
             STRESS(4)=SQRT(STRESS(4))
             WRITE(IFILE("STRESS_M"),2030) I,XTMP,YTMP,STRESS(1),STRESS(2),STRESS(3),STRESS(4)   
          ENDDO
       ENDDO
!
    ENDDO
! **************************
    RETURN 
!
 1000 FORMAT (' R E S U L T S   O F   T H E   O V E R L A P P I N G   E L E M E N T S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_2D_OVLP
!
    SUBROUTINE OUTPUT_QUADS(NDIV,NumEle)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NDIV,NumEle
!   ASSISTANT VARIABLES
    INTEGER :: NDIV1
    INTEGER :: NN,I,J,K,KK,KK1,NUM
    REAL(8) :: RI,RI1,RI2,SI,SI1,SI2,DRS
!
    REAL(8) :: D(2),HU(4),STRESS(4)
    REAL(8) :: HH(1,4),DHH(2,4)
    REAL(8) :: DXDR,DXDS,DYDR,DYDS,DETJ
    REAL(8) :: XTMP,YTMP
!
!    WRITE(IFILE("DISP_M"),1000)
!    WRITE(IFILE("DISP_M"),2000)
!    WRITE(IFILE("STRESS_M"),1000)
!    WRITE(IFILE("STRESS_M"),2010) 
!
    NDIV1=NDIV+1
    DRS=2.0D0/NDIV
!
    DO NN=1,NumEle
       DO I=1,NDIV1
          RI=-1.0D0+(I-1)*DRS
          DO J=1,NDIV1
             SI=-1.0D0+(J-1)*DRS
!            RESET DATA
             D(1)=0.0D0
             D(2)=0.0D0
!
             DO K=1,4
                HU(K)=0.0D0
             ENDDO
! 
             RI1=1.0D0+RI
             RI2=1.0D0-RI
             SI1=1.0D0+SI
             SI2=1.0D0-SI
!
             DXDR=(SI1*(XYZ(1,NN)-XYZ(3,NN))-SI2*(XYZ(5,NN)-XYZ(7,NN)))/4.0D0
             DXDS=(RI1*(XYZ(1,NN)-XYZ(7,NN))+RI2*(XYZ(3,NN)-XYZ(5,NN)))/4.0D0
             DYDR=(SI1*(XYZ(2,NN)-XYZ(4,NN))-SI2*(XYZ(6,NN)-XYZ(8,NN)))/4.0D0
             DYDS=(RI1*(XYZ(2,NN)-XYZ(8,NN))+RI2*(XYZ(4,NN)-XYZ(6,NN)))/4.0D0
!
             DETJ=DXDR*DYDS-DYDR*DXDS
!
             HH(1,1)=RI1*SI1/4.0D0
             DHH(1,1)=(SI1*DYDS-RI1*DYDR)/4.0D0/DETJ
             DHH(2,1)=(-SI1*DXDS+RI1*DXDR)/4.0D0/DETJ
!
             HH(1,2)=RI2*SI1/4.0D0
             DHH(1,2)=(-SI1*DYDS-RI2*DYDR)/4.0/DETJ
             DHH(2,2)=(SI1*DXDS+RI2*DXDR)/4.0/DETJ
!             
             HH(1,3)=RI2*SI2/4.0D0
             DHH(1,3)=(-SI2*DYDS+RI2*DYDR)/4.0/DETJ
             DHH(2,3)=(SI2*DXDS-RI2*DXDR)/4.0/DETJ
!
             HH(1,4)=RI1*SI2/4.0D0
             DHH(1,4)=(SI2*DYDS+RI1*DYDR)/4.0/DETJ
             DHH(2,4)=(-SI2*DXDS-RI1*DXDR)/4.0/DETJ
!
             XTMP=0.0D0
             YTMP=0.0D0
             DO K=1,4
                KK=2*K-1
                KK1=KK+1
                XTMP=XTMP+HH(1,K)*XYZ(KK,NN)
                YTMP=YTMP+HH(1,K)*XYZ(KK1,NN)
!
                NUM=LM(KK,NN) 
                IF(NUM.NE.0) THEN  
                   D(1)=D(1)+DISP(NUM)*HH(1,K)
                   HU(1)=HU(1)+DISP(NUM)*DHH(1,K)
                   HU(3)=HU(3)+DISP(NUM)*DHH(2,K)                   
                ENDIF
                NUM=LM(KK1,NN)
                IF(NUM.NE.0) THEN  
                   D(2)=D(2)+DISP(NUM)*HH(1,K)
                   HU(2)=HU(2)+DISP(NUM)*DHH(2,K)
                   HU(4)=HU(4)+DISP(NUM)*DHH(1,K)
                ENDIF
             ENDDO
             WRITE(IFILE("DISP_M"),2020) NN,XTMP,YTMP,D(1),D(2)
             STRESS(1)=KCS(1)*HU(1)+ KCS(2)*HU(2)
             STRESS(2)=KCS(1)*HU(2)+ KCS(2)*HU(1)
             STRESS(3)=KCS(3)*HU(4)+ KCS(3)*HU(3)
             STRESS(4)=STRESS(1)*STRESS(1)+STRESS(2)*STRESS(2)-STRESS(1)*STRESS(2)
             STRESS(4)=STRESS(4)+3.0D0*STRESS(3)*STRESS(3)
             STRESS(4)=SQRT(STRESS(4))
             WRITE(IFILE("STRESS_M"),2030) NN,XTMP,YTMP,STRESS(1),STRESS(2),STRESS(3),STRESS(4)            
          ENDDO
       ENDDO
    ENDDO
! **************************
    RETURN 
!
 1000 FORMAT (' R E S U L T S   O F   T H E   Q U A D S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_QUADS
!
    SUBROUTINE OUTPUT_TRUSS_OVLP(IIN,IOUT,IDISP,NumEle,NEleAss)
!{{{
    IMPLICIT NONE 
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,NumEle,NEleAss
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: II,KK,IDD,IND,INCK
    INTEGER :: N,ELE,M
    INTEGER :: IC
    REAL(8) :: D(3)
    REAL(8) :: WEIGHT,RAD
    REAL(8) :: DISPP
    REAL(8) :: X_BAR
!
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      DO I=1,3
         D(I)=0.0D0
      END DO
!
      DO II=1,NumEle
         N=NDMFS(1,II)
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         D(1)=0.0D0
!
         KK=ID(1,N)
         IF (KK.NE.0) THEN
             DO J=1,NEleAss
                ELE=LM(J,II)
                IF(ELE.NE.0) THEN
                   M=NDMFS(1,ELE)
                   WEIGHT=WEIND(J,II)
                   IND=NDMFS(2,ELE)
                   RAD=ELEPAR(1,ELE)
                   IF(IND.EQ.1) THEN
                      INCK=1
                   ELSE
                      INCK=0
                   ENDIF
                   DO K=1,NDOF(1,M)
                      IDD=ID(1,M)+K-1
                      SELECT CASE(K+INCK)
                      CASE(1)
                        DISPP=DISP(IDD)
                      CASE(2)
                        X_BAR=(X(N)-X(M))/RAD
                        DISPP=X_BAR*DISP(IDD)
                      CASE(3)
                        X_BAR=(X(N)-X(M))/RAD
                        DISPP=X_BAR*X_BAR*DISP(IDD)
                      CASE DEFAULT
                           WRITE(*,*) "ERROR, /OUTPUTM/OUTPUT_TRUSS_OVLP"
                           STOP
                      END SELECT
                      D(1)=D(1)+WEIGHT*DISPP
                   ENDDO
                ENDIF
             ENDDO
         ENDIF
!
         WRITE (IDISP,2010) N,D
      END DO
!
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,ES14.7) 
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)
!}}}
    END SUBROUTINE OUTPUT_TRUSS_OVLP
!
    SUBROUTINE OUTPUT_TRUSS_OVLP2(IIN,IOUT,IDISP,NumEle,NEleAss)
!{{{
    IMPLICIT NONE 
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,NumEle,NEleAss
!   ASSISTANT VARIABLES
    INTEGER :: I,J,K
    INTEGER :: II,KK,IDD,IND,INCK
    INTEGER :: N,ELE,M
    INTEGER :: IC
    REAL(8) :: D(3)
    REAL(8) :: WEIGHT,RAD
    REAL(8) :: DISPP
    REAL(8) :: X_BAR
    REAL(8) :: WEI(NEleAss)
!
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      DO I=1,3
         D(I)=0.0D0
      END DO
!
      DO II=1,NumEle
         N=NDMFS(1,II)
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         D(1)=0.0D0
!
         KK=ID(1,N)
         IF (KK.NE.0) THEN
             CALL WEI_CAL_OUT(II,NEleAss,X(N),WEI(1))
             DO J=1,NEleAss
                ELE=LM(J,II)
                IF(ELE.NE.0) THEN
                   M=NDMFS(1,ELE)
                   WEIGHT=WEI(J)
                   IND=NDMFS(2,ELE)
                   RAD=ELEPAR(1,ELE)
                   IF(IND.EQ.1) THEN
                      INCK=1
                   ELSE
                      INCK=0
                   ENDIF
                   DO K=1,NDOF(1,M)
                      IDD=ID(1,M)+K-1
                      SELECT CASE(K+INCK)
                      CASE(1)
                        DISPP=DISP(IDD)
                      CASE(2)
                        X_BAR=(X(N)-X(M))/RAD
                        DISPP=X_BAR*DISP(IDD)
                      CASE(3)
                        X_BAR=(X(N)-X(M))/RAD
                        DISPP=X_BAR*X_BAR*DISP(IDD)
                      CASE DEFAULT
                           WRITE(*,*) "ERROR, /OUTPUTM/OUTPUT_TRUSS_OVLP"
                           STOP
                      END SELECT
                      D(1)=D(1)+WEIGHT*DISPP
                   ENDDO
                ENDIF
             ENDDO
         ENDIF
!
         WRITE (IDISP,2010) N,D
      END DO
!
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,ES14.7) 
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)
!}}}
    END SUBROUTINE OUTPUT_TRUSS_OVLP2
!
    SUBROUTINE WEI_CAL_OUT(II,NEleAss,XTMP,WEI)
!{{{
    IMPLICIT NONE
    INTEGER :: II,NEleAss
    REAL(8) :: XTMP
    REAL(8) :: WEI(NEleAss)
!   ASSISTANT VARIABLES
    INTEGER :: I,J,N
    INTEGER :: ELE1,ELE2
    INTEGER :: IND
    REAL(8) :: WEIGHT
    REAL(8) :: WEI_SUM
    REAL(8) :: S,SS(4)
    REAL(8) :: RAD
!
    DO I=1,NEleAss
       WEI(I)=0.0D0
    ENDDO
!
    IND=NDMFS(2,II)
    RAD=ELEPAR(1,II)
    IF(IND.EQ.0) THEN
       DO J=1,NEleAss
          ELE1=LM(J,II)
          IF(ELE1.NE.0) THEN
             N=NDMFS(1,ELE1)
             WEIGHT=(XTMP-X(N))/2.0D0/RAD
             S=ABS(WEIGHT)
             IF ((1-S).LT.1E-12) CYCLE
             SS(1)=S*S
             SS(2)=SS(1)*S
             SS(3)=SS(2)*S
             SS(4)=-12*S+24*SS(1)-12*SS(2)
!
             WEI(J)=1-6*SS(1)+8*SS(2)-3*SS(3)
          ENDIF
       ENDDO
    ELSE
       DO J=1,NEleAss
          ELE1=LM(J,II)
          IF(ELE1.EQ.II) WEI(J)=1.0D0
       ENDDO
    ENDIF
!
    WEI_SUM=0.0D0
    DO I=1,NEleAss
       WEI_SUM=WEI_SUM+WEI(I)
    ENDDO
!
    DO I=1,NEleAss
       WEI(I)=WEI(I)/WEI_SUM
    ENDDO
!
    RETURN
!}}}
    END SUBROUTINE WEI_CAL_OUT
!
    SUBROUTINE OUTDISP_FEM(IDISP)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: IDISP
!
      REAL(8) :: D(3)
      INTEGER :: IC,II,I,KK,IL
      REAL(8) :: H(NEQ)
      LOGICAL :: EXPN1,EXPN2
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NUMNP
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
         DO I=1,3
            KK=ID(I,II)
            EXPN2=KK.NE.0
            IL=I
            EXPN1=NDOF(I,II).EQ.0
            IF (EXPN1.AND.EXPN2) THEN
                D(IL)=DISP(KK)
            ENDIF
         END DO
!
         WRITE (IDISP,2010) II,D
      END DO
!
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,ES14.7) 
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)
!}}}
    END SUBROUTINE OUTDISP_FEM
!
    SUBROUTINE OUTDISP_2D_MFS(IDISP,NEleAss,NumEle)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP
    INTEGER :: NEleAss,NumEle
!
    REAL(8) :: D(3)
    REAL(8) :: H(NEQ)
    REAL(8) :: HI(NPOR_MAX)
    INTEGER :: IC,II,I,KK,IL
    INTEGER :: I1
    INTEGER :: J,JJ,K,KMAX
    INTEGER :: N,M,NUM
    INTEGER :: INDD,NPOR
    REAL(8) :: RAD1,RAD2,SX,SY,S
    LOGICAL :: EXPRN1
    INTEGER :: INC
    INTEGER :: INDSP
!
!
!     Print displacements
!
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NumEle
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
         DO I=1,2
            I1=I+1
            N=NDMFS(1,II)
            INDSP=NDMFS(6,II)
            KK=NDMFS(I1,II)
            SELECT CASE(KK)
            CASE(1)
               D(I)=RESS(I,N)
            CASE(2)
               NUM=ID(I,N)
               D(I)=DISP(NUM)
            CASE(3)
               NUM=ID(I,N)
               D(I)=DISP(NUM)
            CASE(0)
               DO J=1,NEleAss
                  JJ=LM(J,II)
                  IF(JJ.EQ.0) CYCLE
                  M=NDMFS(1,JJ)
                  INDSP=NDMFS(6,JJ)
!
                  RAD1=ELEPAR(3,JJ)
                  RAD1=1.0/RAD1
                  IF(INDSP.EQ.3) THEN
                     RAD2=ELEPAR(4,JJ)
                     RAD2=1.0/RAD2
!
                     SX=1.0-ABS(X(N)-X(M))*RAD1
                     EXPRN1=SX.LT.1E-8
                     IF(EXPRN1) CYCLE
                     SY=1.0-ABS(Y(N)-Y(M))*RAD2
                     EXPRN1=SY.LT.1E-8
                     IF(EXPRN1) CYCLE
                  ELSE
                     S=(X(N)-X(M))*(X(N)-X(M))+(Y(N)-Y(M))*(Y(N)-Y(M))
                     S=SQRT(S)*RAD1
                     EXPRN1=S.GT.1.0
                     IF(EXPRN1) CYCLE
                  ENDIF
!
                  INDD=NDMFS(I+1,JJ)
                  IF(INDD.NE.1) THEN 
                     INC=0
                  ELSE
                     INC=1
                  ENDIF
                  KMAX=NDOF(I,M)
                  HI=QUADS_MFS_HIHJ(NEleAss,JJ,KMAX+INC,X(N),Y(N))                      
                  IF(INDD.EQ.1) D(I)=D(I)+RESS(I,M)*HI(1)
                  DO K=1,KMAX                             
                     NUM=ID(I,M)+K-1
                     NPOR=K+INC
                     D(I)=D(I)+DISP(NUM)*HI(NPOR) 
                  ENDDO
               ENDDO
            CASE DEFAULT 
               WRITE(*,*) "ERROR, IN OUTDISP_MFS!"
               STOP
            END SELECT
         ENDDO 
!
         WRITE (IDISP,2010) N,D
      END DO    
!
    RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)

!}}}
    END SUBROUTINE OUTDISP_2D_MFS
!
    SUBROUTINE OUTPUT_2DMFS_MATLAB(IIN,IOUT,IDISP,NEleAss,NumEle)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,NEleAss,NumEle
!
    REAL(8) :: LENGTHX,LENGTHY
    REAL(8) :: DIRN(3)
    REAL(8) :: XX,YY
!
    INTEGER :: N,INDSP,IND
    INTEGER :: L1,L2,N1,N2
    REAL(8) :: RAD1,RAD2,SX,SY,S
    REAL(8) :: X1,X2,Y1,Y2,XX1,YY1
    LOGICAL :: EXPRN1
    INTEGER :: KMAX,INC,NPOR,NUM
    REAL(8) :: HI(NPOR_MAX)
!
    REAL(8) :: D(2),STRESS(4),STRESSN(4)
    LOGICAL :: REST1,REST2,REST3,LOGIC1,LOGIC2,EXPRN2
!
    INTEGER :: NCOUNT,NCOUNT2,NCOUNT3
    INTEGER :: I,J,K,I1,J1
!
    DIRN(1)=1.0D0
    DIRN(2)=-1.0D0
    DIRN(3)=0.0D0 
!
    WRITE(IFILE("DISP_M"),1000)
    WRITE(IFILE("DISP_M"),2000)
    WRITE(IFILE("STRESS_M"),1000)
    WRITE(IFILE("STRESS_M"),2010)
!
    NCOUNT=0
    NCOUNT2=0
    DO K=1,NumEle
       N=NDMFS(1,K)
       IND=NDMFS(2,K)
       INDSP=NDMFS(6,K)
       RAD1=ELEPAR(3,K)
       LENGTHX=RAD1
       RAD1=1.0/RAD1
       RAD2=ELEPAR(4,K)
       LENGTHY=RAD2
       RAD2=1.0/RAD2
       NCOUNT3=0
!
       DO I=1,4
          STRESSN(I)=0.0D0
       ENDDO
!
       DO I=1,3
          XX=X(N)+LENGTHX*0.2*DIRN(I)
          DO J=1,3
             YY=Y(N)+LENGTHY*0.2*DIRN(J)
!
             D(1)=0.0D0
             D(2)=0.0D0
!
             DO I1=1,4
                STRESS(I1)=0.0D0
             ENDDO
!
             EXPRN2=INDSP.GE.3
             IF(EXPRN2) THEN 
                SX=1.0-ABS(XX-X(N))*RAD1
                EXPRN1=SX.LT.1E-8
                IF(EXPRN1) CYCLE
                SY=1.0-ABS(YY-Y(N))*RAD2
                EXPRN1=SY.LT.1E-8
                IF(EXPRN1) CYCLE
             ELSE
                S=(XX-X(N))*(XX-X(N))+(YY-Y(N))*(YY-Y(N))
                S=1-SQRT(S)*RAD1
                EXPRN1=S.LT.1E-8
                IF(EXPRN1) CYCLE
             ENDIF
!   
             XX1=XX-X(N)
             YY1=YY-Y(N)
             REST1=ABS(XX1).GT.1E-12
             REST2=ABS(YY1).GT.1E-12
             IF ((IND.NE.0).AND.(REST1.OR.REST2)) THEN
                L1=NDMFS(4,K)
                L2=NDMFS(5,K)
                N1=NDMFS(1,L1)
                N2=NDMFS(1,L2)
!
                X2=X(N2)-X(N)
                X1=X(N1)-X(N)
                Y2=Y(N2)-Y(N)
                Y1=Y(N1)-Y(N)
!  
                REST1=areClockwise(X1,X2,Y1,Y2)
                REST2=areClockwise(XX1,X2,YY1,Y2)
                REST3=areClockwise(X1,XX1,Y1,YY1)
                LOGIC1=(REST1.AND.(.NOT.(REST2.AND.REST3)))
                LOGIC2=((.NOT.REST3).AND.(.NOT.REST2))
                IF (LOGIC1.OR.LOGIC2) CYCLE
             ENDIF
!
             CALL OUTPUT_2DMFS(NEleAss,K,XX,YY,D,STRESS) 
             NCOUNT=NCOUNT+1
             WRITE(IFILE("DISP_M"),2020) NCOUNT,XX,YY,D(1),D(2)
             IF((EXPRN2).AND.((I.EQ.3).OR.(J.EQ.3))) CYCLE
             NCOUNT2=NCOUNT2+1
             WRITE(IFILE("STRESS_M"),2030) NCOUNT2,XX,YY,STRESS(1),STRESS(2),STRESS(3),STRESS(4)
             IF(EXPRN2) THEN
                NCOUNT3=NCOUNT3+1
                STRESSN(1)=STRESSN(1)+STRESS(1)
                STRESSN(2)=STRESSN(2)+STRESS(2)
                STRESSN(3)=STRESSN(3)+STRESS(3)
                STRESSN(4)=STRESSN(4)+STRESS(4)
             ENDIF
         ENDDO
       ENDDO
       IF(EXPRN2) THEN
          NCOUNT2=NCOUNT2+1
          WRITE(IFILE("STRESS_M"),2030) NCOUNT2,X(N),Y(N),&
                  STRESSN(1)/NCOUNT3,STRESSN(2)/NCOUNT3,STRESSN(3)/NCOUNT3,STRESSN(4)/NCOUNT3
       ENDIF
    ENDDO
!
 1000 FORMAT (' R E S U L T S   O F   T H E   O V E R L A P P I N G   E L E M E N T S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_2DMFS_MATLAB
!
    SUBROUTINE OUTPUT_COUP1_MATLAB(IIN,IOUT,IDISP,NumEle)
!{{{
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IIN,IOUT,IDISP,NumEle
! ASSISTANT VARIABLE
    REAL(8) :: R(3),S(3),RR,SS
    INTEGER :: NCOUNT
    REAL(8) :: D(2),STRESS(4)
    INTEGER :: N,I,J,K
    REAL(8) :: XX,YY
!
    WRITE(IFILE("DISP_M"),1000)
    WRITE(IFILE("DISP_M"),2000)
    WRITE(IFILE("STRESS_M"),1000)
    WRITE(IFILE("STRESS_M"),2010)
!
    R(1)=0.7746
    R(2)=0
    R(3)=-0.7746
    S(1)=0.7746
    S(2)=0
    S(3)=-0.7746    
!
    NCOUNT=0
    DO N=1,NumEle
       DO I=1,3
          RR=R(I)
          DO J=1,3
             SS=S(J)
!
             DO K=1,4
                STRESS(K)=0.0D0
             ENDDO
!
             D(1)=0.0D0
             D(2)=0.0D0
             CALL OUTPUT_COUP1(RR,SS,N,XX,YY,D,STRESS)
             NCOUNT=NCOUNT+1
             WRITE(IFILE("DISP_M"),2020) NCOUNT,XX,YY,D(1),D(2)
             WRITE(IFILE("STRESS_M"),2030) NCOUNT,XX,YY,STRESS(1),STRESS(2),STRESS(3),STRESS(4)
          ENDDO
       ENDDO
    ENDDO
!
 1000 FORMAT (' R E S U L T S   O F   T H E   C O U P L I N G   E L E M E N T S')
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE    X-DISPLACEMENT   Y-DISPLACEMENT')
 2010 FORMAT (' S T R E S S ',//,'  NUMBER ',10X,&
              'X-COORDINATE    Y-COORDINATE      STRESS_XX       STRESS_YY       STRESS_XY     STRESS_EFFECTIVE')
 2020 FORMAT (1X,I5,8X,4ES18.6)
 2030 FORMAT (1X,I5,8X,6ES18.6)
!!}}}
    END SUBROUTINE OUTPUT_COUP1_MATLAB
!
    SUBROUTINE OUTDISP_TRUSS_MFS(IDISP,NEleAss,NumEle)
!{{{
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To print displacements                                          .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: IDISP
    INTEGER :: NEleAss,NumEle
!
    REAL(8) :: D(3)
    REAL(8) :: H(NEQ)
    REAL(8) :: HI(NPOR_MAX)
    INTEGER :: IC,II,I,KK,IL
    INTEGER :: I1
    INTEGER :: J,JJ,K,KMAX
    INTEGER :: N,M,NUM
    INTEGER :: INDD,NPOR
    REAL(8) :: S,Radius
    INTEGER :: INC
    INTEGER :: INDSP
!
!     Print displacements
!
      WRITE (IDISP,1200) 
      WRITE (IDISP,1000) TIME
      WRITE (IDISP,2000)
      IC=4
!
      DO II=1,NumEle
         IC=IC + 1
         IF (IC.GE.56) THEN
            WRITE (IDISP,1000) TIME
            WRITE (IDISP,2000)
            IC=4
         END IF
!
         DO I=1,3
            D(I)=0.
         END DO
!
            N=NDMFS(1,II)
            KK=NDMFS(2,II)
            SELECT CASE(KK)
            CASE(1)
               D(1)=RESS(1,N)
            CASE(2)
               NUM=ID(1,N)
               D(1)=DISP(NUM)
            CASE(3)
               NUM=ID(1,N)
               D(1)=DISP(NUM)
            CASE(0)
               DO J=1,NEleAss
                  JJ=LM(J,II)
                  IF(JJ.EQ.0) CYCLE
                     M=NDMFS(1,JJ)
                     Radius=ELEPAR(2,JJ)
                     S=ABS(X(N)-X(M))
                     S=S/Radius
                     IF(S.LE.1) THEN
                        INDD=NDMFS(2,JJ)
                        IF(INDD.NE.1) THEN 
                           INC=0
                        ELSE
                           INC=1
                        ENDIF
                        KMAX=NDOF(1,M)
                        INDSP=NDMFS(3,JJ)
                        HI=TRUSS_MFS_HIHJ(NEleAss,JJ,KMAX+INC,X(N),INDSP)       
                        IF(INDD.EQ.1) D(1)=D(1)+RESS(1,M)*HI(1)
                        DO K=1,KMAX                             
                             NUM=ID(1,M)+K-1
                             NPOR=K+INC
                             D(1)=D(1)+DISP(NUM)*HI(K) 
                        ENDDO
                     ENDIF
               ENDDO
            CASE DEFAULT 
               WRITE(*,*) "ERROR, IN OUTDISP_MFS!"
               STOP
            END SELECT
!
         WRITE (IDISP,2010) N,D
      END DO    
!
    RETURN
!
 1200 FORMAT (/,'TRUSS(METHOD OF FINITE SPHERES)')
 1000 FORMAT (/,' TIME =',5X,F13.3)
 1050 FORMAT (' TOTAL STRAIN ENERGY =',5X,ES18.6)
 2000 FORMAT (' D I S P L A C E M E N T S',//,'  NODE ',10X,&
              'X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT')
 2010 FORMAT (1X,I3,8X,3ES18.6)

!}}}
    END SUBROUTINE OUTDISP_TRUSS_MFS
!
    SUBROUTINE TRUSS_OUT(ISTRE,ADELE,NG)
!{{{
!--------------------------------------------------
!   CALCULATE AND OUTPUT STRESS OF TRUSS ELEMENT
!--------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ISTRE,ADELE,NG
!
    INTEGER :: NumEle,Mtype
!   ASSISTANT VARIABLES
    REAL(8) :: E,DEN,AREA
    INTEGER :: ALLocateStatus
!   VARIABLES FOR STRESS CALCULATION
    INTEGER :: IPRINT
    INTEGER :: N,L,I,J
    REAL(8) :: XL2,D(3),ST(6),STR,P

!   READ 
    NumEle = NPAR(3)
    Mtype  = NPAR(2)
!   SET MATERIAL PARAMETERS
    E    = MATP(1,Mtype)
    DEN  = MATP(3,Mtype)
!
!   CALCULATE AND OUTPUT STRESS
!   
      IPRINT=0
      DO N=1,NumEle
         AREA=ELEPAR(1,N)
         IPRINT=IPRINT + 1
         IF (IPRINT.GT.50) IPRINT=1
         IF (IPRINT.EQ.1) THEN
            WRITE (ISTRE,1000) TIME
            WRITE (ISTRE,2060) NG
         END IF
         XL2=0.
         DO L=1,3
            D(L) = XYZ(L,N) - XYZ(L+3,N)
            XL2=XL2 + D(L)*D(L)
         END DO
         DO L=1,3
            ST(L)=(D(L)/XL2)*E
            ST(L+3)=-ST(L)
         END DO
         STR=0.0         
         DO L=1,3
            I=LM(L,N)
            IF (I.GT.0) STR=STR + ST(L)*DISP(I)
            J=LM(L+3,N)
            IF (J.GT.0) STR=STR + ST(L+3)*DISP(J)
         END DO
         P=STR*AREA
         WRITE (ISTRE,2070) N,P,STR
      END DO 
      RETURN
!
 1000 FORMAT (/,' TIME =',5X,F13.3)
 2060 FORMAT (' S T R E S S  C A L C U L A T I O N S  F O R  ',&
              'E L E M E N T  G R O U P',I4,//,&
              '  ELEMENT',13X,'FORCE',12X,'STRESS',/,&
              '  NUMBER',/)
 2070 FORMAT (1X,I5,11X,E13.6,4X,E13.6)
!}}}
    END SUBROUTINE TRUSS_OUT
!
end module OUTPUTM
