! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       
! .                            CYCLE  V0.3                                
! .                                                                       
! .    Gaussian integration module for CYCLE 
! .
! .    Copyright(c) 2016-2056 Lingbo Zhang      
! .    Author: Lingbo Zhang 
! .    06/03/2016 at Massachusetts Institute of Technology
! .    EMAIL: lingboz2015@gmail.com
! .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

module GAUSSM
   USE WAREHOUSE
   IMPLICIT NONE
   PUBLIC gauleg
   PUBLIC GAUTRIANGLE
   PUBLIC GAUTRIANGLE_OFE1
!   PUBLIC NINTG_TETRA
!   PRIVATE SET_VALUE
!
   INTEGER,PARAMETER :: DBP = SELECTED_REAL_KIND(15,307)
   REAL(DBP) :: NEWV
   REAL(DBP) :: EPS,M_PI
   PARAMETER (EPS=3.0D-15)         ! EPS IS THE RELATIVE PRECISION
   PARAMETER (M_PI=3.141592654D0)  ! PI VALUE
!
   CONTAINS
!
   SUBROUTINE  GAULEG(ngp, xabsc, weig)
!{{{
!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!********************************************************************************
      implicit none
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
!
      m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

      do i = 1, m       ! Loop over the desired roots */

         z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root
!*          we enter the main loop of refinement by NEWTON'S method   */
100      p1 = 1.0d0
         p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

         do j = 1, ngp
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
         enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
         pp = ngp*(z*p1-p2)/(z*z-1.0d0)
         z1 = z
         z = z1 - p1/pp             ! Newton's Method  */

         if (dabs(z-z1) .gt. EPS) GOTO  100

         xabsc(i) =  - z                        ! Roots will be bewteen -1.0 & 1.0 */
         xabsc(ngp+1-i) =  + z                  ! and symmetric about the origin  */
         weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
         weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop
!
   RETURN
!}}}
   End subroutine GAULEG
!********************************************************************************
!*     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
   recursive function qgauss(func, x1, x2, ngp) RESULT(intgrl)
!{{{
      implicit none
      REAL(dbp)  intgrl, x1, x2, func
      REAL(dbp)  xm, xl
      INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)
!
      call gauleg(ngp, xabsc, weig)
!
      intgrl = 0.0d0
      xm = 0.5 * (x2 + x1)
      xl = 0.5 * (x2 - x1)
      do j = 1, ngp
        intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) )
      END do

      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
!}}}
   END function qgauss
!
   SUBROUTINE GAUTRIANGLE(NGDOF,N1,N2,N3)
!{{{
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: NGDOF,N1,N2,N3
!  ASSISTANT VARIABLES
   INTEGER :: I
   REAL(8) :: X1,Y1,X2,Y2,X3,Y3
   REAL(8) :: XG(3,79),WGT(79)
!
   X1=X(N1)
   Y1=Y(N1)
   X2=X(N2)
   Y2=Y(N2)
   X3=X(N3)
   Y3=Y(N3) 
!
   SELECT CASE(NGDOF)
   CASE(1)
        WGT(1)  = 1.0D0
        XG(1,1) = 0.166666666666667D0
        XG(2,1) = 0.666666666666667D0
        XG(3,1) = 0.166666666666667D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3

   CASE(2)
!{{{
        WGT(1)  = 0.333333333333333D0
        XG(1,1) = 0.666666666666667D0
        XG(2,1) = 0.166666666666667D0
        XG(3,1) = 0.166666666666667D0
!
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(2,3)=WGT(1)
        XABSC_WEIG_IN(3,3)=WGT(1)
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(2,1)=XG(2,1)
        ACOOR_IN(2,2)=XG(1,1)
        ACOOR_IN(2,3)=XG(3,1)
        XABSC_WEIG_IN(2,1)=XG(2,1)*X1+XG(1,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(2,2)=XG(2,1)*Y1+XG(1,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(3,1)=XG(2,1)
        ACOOR_IN(3,2)=XG(3,1)
        ACOOR_IN(3,3)=XG(1,1)
        XABSC_WEIG_IN(3,1)=XG(2,1)*X1+XG(3,1)*X2+XG(1,1)*X3
        XABSC_WEIG_IN(3,2)=XG(2,1)*Y1+XG(3,1)*Y2+XG(1,1)*Y3
!}}}
   CASE(3)
!{{{
!         WGT(1)  =-0.562500000000000D0
!         XG(1,1) = 0.333333333333333D0
!         XG(2,1) = 0.333333333333333D0
!         XG(3,1) = 0.333333333333333D0
!
!         WGT(2)  = 0.520833333333333D0
!         XG(1,2) = 0.600000000000000D0
!         XG(2,2) = 0.200000000000000D0
!         XG(3,2) = 0.200000000000000D0
!!
!         XABSC_WEIG_IN(1,3)=WGT(1)
!         ACOOR_IN(1,1)=XG(1,1)
!         ACOOR_IN(1,2)=XG(2,1)
!         ACOOR_IN(1,3)=XG(3,1)
!         XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
!         XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!!
!         DO I=2,4
!            XABSC_WEIG_IN(I,3)=WGT(2)
!         ENDDO       
!         ACOOR_IN(2,1)=XG(1,2)
!         ACOOR_IN(2,2)=XG(2,2)
!         ACOOR_IN(2,3)=XG(3,2)
!         XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
!         XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!!
!         ACOOR_IN(3,1)=XG(2,2)
!         ACOOR_IN(3,2)=XG(1,2)
!         ACOOR_IN(3,3)=XG(3,2)
!         XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
!         XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!!
!         ACOOR_IN(4,1)=XG(2,2)
!         ACOOR_IN(4,2)=XG(3,2)
!         ACOOR_IN(4,3)=XG(1,2)
!         XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
!         XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        WGT(1)  = 0.166666666666667D0
        XG(1,1) = 0.659027622374092D0
        XG(2,1) = 0.231933368553031D0
        XG(3,1) = 0.109039009072877D0
!
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(2,3)=WGT(1)
        XABSC_WEIG_IN(3,3)=WGT(1)
        XABSC_WEIG_IN(4,3)=WGT(1)
        XABSC_WEIG_IN(5,3)=WGT(1)
        XABSC_WEIG_IN(6,3)=WGT(1)
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(2,1)=XG(1,1)
        ACOOR_IN(2,2)=XG(3,1)
        ACOOR_IN(2,3)=XG(2,1)
        XABSC_WEIG_IN(2,1)=XG(1,1)*X1+XG(3,1)*X2+XG(2,1)*X3
        XABSC_WEIG_IN(2,2)=XG(1,1)*Y1+XG(3,1)*Y2+XG(2,1)*Y3
!
        ACOOR_IN(3,1)=XG(2,1)
        ACOOR_IN(3,2)=XG(1,1)
        ACOOR_IN(3,3)=XG(3,1)
        XABSC_WEIG_IN(3,1)=XG(2,1)*X1+XG(1,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(3,2)=XG(2,1)*Y1+XG(1,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(4,1)=XG(2,1)
        ACOOR_IN(4,2)=XG(3,1)
        ACOOR_IN(4,3)=XG(1,1)
        XABSC_WEIG_IN(4,1)=XG(2,1)*X1+XG(3,1)*X2+XG(1,1)*X3
        XABSC_WEIG_IN(4,2)=XG(2,1)*Y1+XG(3,1)*Y2+XG(1,1)*Y3
!
        ACOOR_IN(5,1)=XG(3,1)
        ACOOR_IN(5,2)=XG(2,1)
        ACOOR_IN(5,3)=XG(1,1)
        XABSC_WEIG_IN(5,1)=XG(3,1)*X1+XG(2,1)*X2+XG(1,1)*X3
        XABSC_WEIG_IN(5,2)=XG(3,1)*Y1+XG(2,1)*Y2+XG(1,1)*Y3
!
        ACOOR_IN(6,1)=XG(3,1)
        ACOOR_IN(6,2)=XG(1,1)
        ACOOR_IN(6,3)=XG(2,1)
        XABSC_WEIG_IN(6,1)=XG(3,1)*X1+XG(1,1)*X2+XG(2,1)*X3
        XABSC_WEIG_IN(6,2)=XG(3,1)*Y1+XG(1,1)*Y2+XG(2,1)*Y3
!}}}
   CASE(4)
!{{{
        WGT(1)  = 0.375000000000000D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
        WGT(2)  = 0.104166666666667D0
        XG(1,2) = 0.736712498968435D0
        XG(2,2) = 0.237932366472434D0
        XG(3,2) = 0.025355134559132D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        XABSC_WEIG_IN(2,3)=WGT(2)
        XABSC_WEIG_IN(3,3)=WGT(2)
        XABSC_WEIG_IN(4,3)=WGT(2)
        XABSC_WEIG_IN(5,3)=WGT(2)
        XABSC_WEIG_IN(6,3)=WGT(2)
        XABSC_WEIG_IN(7,3)=WGT(2)
!
        ACOOR_IN(2,1)=XG(1,2)
        ACOOR_IN(2,2)=XG(2,2)
        ACOOR_IN(2,3)=XG(3,2)
        XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(3,1)=XG(1,2)
        ACOOR_IN(3,2)=XG(3,2)
        ACOOR_IN(3,3)=XG(2,2)
        XABSC_WEIG_IN(3,1)=XG(1,2)*X1+XG(3,2)*X2+XG(2,2)*X3
        XABSC_WEIG_IN(3,2)=XG(1,2)*Y1+XG(3,2)*Y2+XG(2,2)*Y3
!
        ACOOR_IN(4,1)=XG(2,2)
        ACOOR_IN(4,2)=XG(3,2)
        ACOOR_IN(4,3)=XG(1,2)
        XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        ACOOR_IN(5,1)=XG(2,2)
        ACOOR_IN(5,2)=XG(1,2)
        ACOOR_IN(5,3)=XG(3,2)
        XABSC_WEIG_IN(5,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(5,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(6,1)=XG(3,2)
        ACOOR_IN(6,2)=XG(1,2)
        ACOOR_IN(6,3)=XG(2,2)
        XABSC_WEIG_IN(6,1)=XG(3,2)*X1+XG(1,2)*X2+XG(2,2)*X3
        XABSC_WEIG_IN(6,2)=XG(3,2)*Y1+XG(1,2)*Y2+XG(2,2)*Y3
!
        ACOOR_IN(7,1)=XG(3,2)
        ACOOR_IN(7,2)=XG(2,2)
        ACOOR_IN(7,3)=XG(1,2)
        XABSC_WEIG_IN(7,1)=XG(3,2)*X1+XG(2,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(7,2)=XG(3,2)*Y1+XG(2,2)*Y2+XG(1,2)*Y3
!}}}
   CASE(5)
!{{{
!         WGT(1)  = 0.225000000000000D0
!         XG(1,1) = 0.333333333333333D0
!         XG(2,1) = 0.333333333333333D0
!         XG(3,1) = 0.333333333333333D0
!
!         WGT(2)  = 0.125939180544827D0
!         XG(1,2) = 0.797426985353087D0
!         XG(2,2) = 0.101286507323456D0
!         XG(3,2) = 0.101286507323456D0
!
!         WGT(3)  = 0.132394152788506D0
!         XG(1,3) = 0.470142064105115D0
!         XG(2,3) = 0.470142064105115D0
!         XG(3,3) = 0.059715871789770D0
!!
!         XABSC_WEIG_IN(1,3)=WGT(1)
!         ACOOR_IN(1,1)=XG(1,1)
!         ACOOR_IN(1,2)=XG(2,1)
!         ACOOR_IN(1,3)=XG(3,1)
!         XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
!         XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!!
!         DO I=2,4
!            XABSC_WEIG_IN(I,3)=WGT(2)
!         ENDDO
!         ACOOR_IN(2,1)=XG(1,2)
!         ACOOR_IN(2,2)=XG(2,2)
!         ACOOR_IN(2,3)=XG(3,2)
!         XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
!         XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!!
!         ACOOR_IN(3,1)=XG(2,2)
!         ACOOR_IN(3,2)=XG(1,2)
!         ACOOR_IN(3,3)=XG(3,2)
!         XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
!         XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!!
!         ACOOR_IN(4,1)=XG(2,2)
!         ACOOR_IN(4,2)=XG(3,2)
!         ACOOR_IN(4,3)=XG(1,2)
!         XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
!         XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!!
!         DO I=5,7
!            XABSC_WEIG_IN(I,3)=WGT(3)
!         ENDDO
!         ACOOR_IN(5,1)=XG(1,3)
!         ACOOR_IN(5,2)=XG(2,3)
!         ACOOR_IN(5,3)=XG(3,3)
!         XABSC_WEIG_IN(5,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
!         XABSC_WEIG_IN(5,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!!
!         ACOOR_IN(6,1)=XG(2,3)
!         ACOOR_IN(6,2)=XG(1,3)
!         ACOOR_IN(6,3)=XG(3,3)
!         XABSC_WEIG_IN(6,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
!         XABSC_WEIG_IN(6,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!!
!         ACOOR_IN(7,1)=XG(2,3)
!         ACOOR_IN(7,2)=XG(3,3)
!         ACOOR_IN(7,3)=XG(1,3)
!         XABSC_WEIG_IN(7,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
!         XABSC_WEIG_IN(7,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3

!******************************************

        WGT(1)  = 0.205950504760887D0
        XG(1,1) = 0.124949503233232D0
        XG(2,1) = 0.437525248383384D0
        XG(3,1) = 0.437525248383384D0
        WGT(2)  = 0.063691414286223D0
        XG(1,2) = 0.797112651860071D0
        XG(2,2) = 0.165409927389841D0
        XG(3,2) = 0.037477420750088D0
!
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(2,3)=WGT(1)
        XABSC_WEIG_IN(3,3)=WGT(1)
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(2,1)=XG(2,1)
        ACOOR_IN(2,2)=XG(1,1)
        ACOOR_IN(2,3)=XG(3,1)
        XABSC_WEIG_IN(2,1)=XG(2,1)*X1+XG(1,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(2,2)=XG(2,1)*Y1+XG(1,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(3,1)=XG(2,1)
        ACOOR_IN(3,2)=XG(3,1)
        ACOOR_IN(3,3)=XG(1,1)
        XABSC_WEIG_IN(3,1)=XG(2,1)*X1+XG(3,1)*X2+XG(1,1)*X3
        XABSC_WEIG_IN(3,2)=XG(2,1)*Y1+XG(3,1)*Y2+XG(1,1)*Y3
!
        XABSC_WEIG_IN(4,3)=WGT(2)
        XABSC_WEIG_IN(5,3)=WGT(2)
        XABSC_WEIG_IN(6,3)=WGT(2)
        XABSC_WEIG_IN(7,3)=WGT(2)
        XABSC_WEIG_IN(8,3)=WGT(2)
        XABSC_WEIG_IN(9,3)=WGT(2)
!
        ACOOR_IN(4,1)=XG(1,2)
        ACOOR_IN(4,2)=XG(2,2)
        ACOOR_IN(4,3)=XG(3,2)
        XABSC_WEIG_IN(4,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(4,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(5,1)=XG(1,2)
        ACOOR_IN(5,2)=XG(3,2)
        ACOOR_IN(5,3)=XG(2,2)
        XABSC_WEIG_IN(5,1)=XG(1,2)*X1+XG(3,2)*X2+XG(2,2)*X3
        XABSC_WEIG_IN(5,2)=XG(1,2)*Y1+XG(3,2)*Y2+XG(2,2)*Y3
!
        ACOOR_IN(6,1)=XG(2,2)
        ACOOR_IN(6,2)=XG(3,2)
        ACOOR_IN(6,3)=XG(1,2)
        XABSC_WEIG_IN(6,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(6,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        ACOOR_IN(7,1)=XG(2,2)
        ACOOR_IN(7,2)=XG(1,2)
        ACOOR_IN(7,3)=XG(3,2)
        XABSC_WEIG_IN(7,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(7,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(8,1)=XG(3,2)
        ACOOR_IN(8,2)=XG(1,2)
        ACOOR_IN(8,3)=XG(2,2)
        XABSC_WEIG_IN(8,1)=XG(3,2)*X1+XG(1,2)*X2+XG(2,2)*X3
        XABSC_WEIG_IN(8,2)=XG(3,2)*Y1+XG(1,2)*Y2+XG(2,2)*Y3
!
        ACOOR_IN(9,1)=XG(3,2)
        ACOOR_IN(9,2)=XG(2,2)
        ACOOR_IN(9,3)=XG(1,2)
        XABSC_WEIG_IN(9,1)=XG(3,2)*X1+XG(2,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(9,2)=XG(3,2)*Y1+XG(2,2)*Y2+XG(1,2)*Y3
!}}}
   CASE(6)
!{{{
        WGT(1)  = 0.050844906370207D0
        XG(1,1) = 0.873821971016996D0
        XG(2,1) = 0.063089014491502D0
        XG(3,1) = 0.063089014491502D0
        WGT(2)  = 0.116786275726379D0
        XG(1,2) = 0.501426509658179D0
        XG(2,2) = 0.249286745170910D0
        XG(3,2) = 0.249286745170911D0
        WGT(3)  = 0.082851075618374D0
        XG(1,3) = 0.636502499121399D0
        XG(2,3) = 0.310352451033785D0
        XG(3,3) = 0.053145049844816D0
!
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(2,3)=WGT(1)
        XABSC_WEIG_IN(3,3)=WGT(1)
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(2,1)=XG(2,1)
        ACOOR_IN(2,2)=XG(1,1)
        ACOOR_IN(2,3)=XG(3,1)
        XABSC_WEIG_IN(2,1)=XG(2,1)*X1+XG(1,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(2,2)=XG(2,1)*Y1+XG(1,1)*Y2+XG(3,1)*Y3
!
        ACOOR_IN(3,1)=XG(2,1)
        ACOOR_IN(3,2)=XG(3,1)
        ACOOR_IN(3,3)=XG(1,1)
        XABSC_WEIG_IN(3,1)=XG(2,1)*X1+XG(3,1)*X2+XG(1,1)*X3
        XABSC_WEIG_IN(3,2)=XG(2,1)*Y1+XG(3,1)*Y2+XG(1,1)*Y3
!
        XABSC_WEIG_IN(4,3)=WGT(2)
        XABSC_WEIG_IN(5,3)=WGT(2)
        XABSC_WEIG_IN(6,3)=WGT(2)
!
        ACOOR_IN(4,1)=XG(1,2)
        ACOOR_IN(4,2)=XG(2,2)
        ACOOR_IN(4,3)=XG(3,2)
        XABSC_WEIG_IN(4,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(4,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(5,1)=XG(2,2)
        ACOOR_IN(5,2)=XG(1,2)
        ACOOR_IN(5,3)=XG(3,2)
        XABSC_WEIG_IN(5,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(5,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(6,1)=XG(2,2)
        ACOOR_IN(6,2)=XG(3,2)
        ACOOR_IN(6,3)=XG(1,2)
        XABSC_WEIG_IN(6,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(6,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        XABSC_WEIG_IN(7,3)=WGT(3)
        XABSC_WEIG_IN(8,3)=WGT(3)
        XABSC_WEIG_IN(9,3)=WGT(3)
        XABSC_WEIG_IN(10,3)=WGT(3)
        XABSC_WEIG_IN(11,3)=WGT(3)
        XABSC_WEIG_IN(12,3)=WGT(3)
!
        ACOOR_IN(7,1)=XG(1,3)
        ACOOR_IN(7,2)=XG(2,3)
        ACOOR_IN(7,3)=XG(3,3)
        XABSC_WEIG_IN(7,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(7,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(8,1)=XG(1,3)
        ACOOR_IN(8,2)=XG(3,3)
        ACOOR_IN(8,3)=XG(2,3)
        XABSC_WEIG_IN(8,1)=XG(1,3)*X1+XG(3,3)*X2+XG(2,3)*X3
        XABSC_WEIG_IN(8,2)=XG(1,3)*Y1+XG(3,3)*Y2+XG(2,3)*Y3
!
        ACOOR_IN(9,1)=XG(2,3)
        ACOOR_IN(9,2)=XG(3,3)
        ACOOR_IN(9,3)=XG(1,3)
        XABSC_WEIG_IN(9,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(9,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3
!
        ACOOR_IN(10,1)=XG(2,3)
        ACOOR_IN(10,2)=XG(1,3)
        ACOOR_IN(10,3)=XG(3,3)
        XABSC_WEIG_IN(10,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(10,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(11,1)=XG(3,3)
        ACOOR_IN(11,2)=XG(1,3)
        ACOOR_IN(11,3)=XG(2,3)
        XABSC_WEIG_IN(11,1)=XG(3,3)*X1+XG(1,3)*X2+XG(2,3)*X3
        XABSC_WEIG_IN(11,2)=XG(3,3)*Y1+XG(1,3)*Y2+XG(2,3)*Y3
!
        ACOOR_IN(12,1)=XG(3,3)
        ACOOR_IN(12,2)=XG(2,3)
        ACOOR_IN(12,3)=XG(1,3)
        XABSC_WEIG_IN(12,1)=XG(3,3)*X1+XG(2,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(12,2)=XG(3,3)*Y1+XG(2,3)*Y2+XG(1,3)*Y3
!}}}
   CASE(7)
!{{{
        WGT(1) = -0.149570044467670D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
        WGT(2)  = 0.175615257433204D0
        XG(1,2) = 0.479308067841923D0
        XG(2,2) = 0.260345966079038D0
        XG(3,2) = 0.260345966079038D0
        WGT(3)   = 0.053347235608839D0
        XG(1,3) = 0.869739794195568D0
        XG(2,3) = 0.065130102902216D0
        XG(3,3) = 0.065130102902216D0
        WGT(4)   = 0.077113760890257D0
        XG(1,4) = 0.638444188569809D0
        XG(2,4) = 0.312865496004875D0
        XG(3,4) = 0.048690315425316D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        XABSC_WEIG_IN(2,3)=WGT(2)
        XABSC_WEIG_IN(3,3)=WGT(2)
        XABSC_WEIG_IN(4,3)=WGT(2)
!
        ACOOR_IN(2,1)=XG(1,2)
        ACOOR_IN(2,2)=XG(2,2)
        ACOOR_IN(2,3)=XG(3,2)
        XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(3,1)=XG(2,2)
        ACOOR_IN(3,2)=XG(1,2)
        ACOOR_IN(3,3)=XG(3,2)
        XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(4,1)=XG(2,2)
        ACOOR_IN(4,2)=XG(3,2)
        ACOOR_IN(4,3)=XG(1,2)
        XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        XABSC_WEIG_IN(5,3)=WGT(3)
        XABSC_WEIG_IN(6,3)=WGT(3)
        XABSC_WEIG_IN(7,3)=WGT(3)
!
        ACOOR_IN(5,1)=XG(1,3)
        ACOOR_IN(5,2)=XG(2,3)
        ACOOR_IN(5,3)=XG(3,3)
        XABSC_WEIG_IN(5,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(5,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(6,1)=XG(2,3)
        ACOOR_IN(6,2)=XG(1,3)
        ACOOR_IN(6,3)=XG(3,3)
        XABSC_WEIG_IN(6,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(6,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(7,1)=XG(2,3)
        ACOOR_IN(7,2)=XG(3,3)
        ACOOR_IN(7,3)=XG(1,3)
        XABSC_WEIG_IN(7,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(7,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3
!
        XABSC_WEIG_IN(8,3)=WGT(4)
        XABSC_WEIG_IN(9,3)=WGT(4)
        XABSC_WEIG_IN(10,3)=WGT(4)
        XABSC_WEIG_IN(11,3)=WGT(4)
        XABSC_WEIG_IN(12,3)=WGT(4)
        XABSC_WEIG_IN(13,3)=WGT(4)
!
        ACOOR_IN(8,1)=XG(1,4)
        ACOOR_IN(8,2)=XG(2,4)
        ACOOR_IN(8,3)=XG(3,4)
        XABSC_WEIG_IN(8,1)=XG(1,4)*X1+XG(2,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(8,2)=XG(1,4)*Y1+XG(2,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(9,1)=XG(1,4)
        ACOOR_IN(9,2)=XG(3,4)
        ACOOR_IN(9,3)=XG(2,4)
        XABSC_WEIG_IN(9,1)=XG(1,4)*X1+XG(3,4)*X2+XG(2,4)*X3
        XABSC_WEIG_IN(9,2)=XG(1,4)*Y1+XG(3,4)*Y2+XG(2,4)*Y3
!
        ACOOR_IN(10,1)=XG(2,4)
        ACOOR_IN(10,2)=XG(3,4)
        ACOOR_IN(10,3)=XG(1,4)
        XABSC_WEIG_IN(10,1)=XG(2,4)*X1+XG(3,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(10,2)=XG(2,4)*Y1+XG(3,4)*Y2+XG(1,4)*Y3
!
        ACOOR_IN(11,1)=XG(2,4)
        ACOOR_IN(11,2)=XG(1,4)
        ACOOR_IN(11,3)=XG(3,4)
        XABSC_WEIG_IN(11,1)=XG(2,4)*X1+XG(1,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(11,2)=XG(2,4)*Y1+XG(1,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(12,1)=XG(3,4)
        ACOOR_IN(12,2)=XG(1,4)
        ACOOR_IN(12,3)=XG(2,4)
        XABSC_WEIG_IN(12,1)=XG(3,4)*X1+XG(1,4)*X2+XG(2,4)*X3
        XABSC_WEIG_IN(12,2)=XG(3,4)*Y1+XG(1,4)*Y2+XG(2,4)*Y3
!
        ACOOR_IN(13,1)=XG(3,4)
        ACOOR_IN(13,2)=XG(2,4)
        ACOOR_IN(13,3)=XG(1,4)
        XABSC_WEIG_IN(13,1)=XG(3,4)*X1+XG(2,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(13,2)=XG(3,4)*Y1+XG(2,4)*Y2+XG(1,4)*Y3
!}}}
   CASE(8)
!{{{
        WGT(1)  = 0.144315607677787D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.095091634267285D0
        XG(1,2) = 0.081414823414554D0
        XG(2,2) = 0.459292588292723D0
        XG(3,2) = 0.459292588292723D0
!
        WGT(3)  = 0.103217370534718D0
        XG(1,3) = 0.658861384496480D0
        XG(2,3) = 0.170569307751760D0
        XG(3,3) = 0.170569307751760D0
!
        WGT(4)  = 0.032458497623198D0
        XG(1,4) = 0.898905543365938D0
        XG(2,4) = 0.050547228317031D0
        XG(3,4) = 0.050547228317031D0
!
        WGT(5)  = 0.027230314174435D0
        XG(1,5) = 0.008394777409958D0
        XG(2,5) = 0.263112829634638D0
        XG(3,5) = 0.728492392955404D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        XABSC_WEIG_IN(2,3)=WGT(2)
        XABSC_WEIG_IN(3,3)=WGT(2)
        XABSC_WEIG_IN(4,3)=WGT(2)

        ACOOR_IN(2,1)=XG(1,2)
        ACOOR_IN(2,2)=XG(2,2)
        ACOOR_IN(2,3)=XG(3,2)
        XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(3,1)=XG(2,2)
        ACOOR_IN(3,2)=XG(1,2)
        ACOOR_IN(3,3)=XG(3,2)
        XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(4,1)=XG(2,2)
        ACOOR_IN(4,2)=XG(3,2)
        ACOOR_IN(4,3)=XG(1,2)
        XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        XABSC_WEIG_IN(5,3)=WGT(3)
        XABSC_WEIG_IN(6,3)=WGT(3)
        XABSC_WEIG_IN(7,3)=WGT(3)
!
        ACOOR_IN(5,1)=XG(1,3)
        ACOOR_IN(5,2)=XG(2,3)
        ACOOR_IN(5,3)=XG(3,3)
        XABSC_WEIG_IN(5,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(5,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(6,1)=XG(2,3)
        ACOOR_IN(6,2)=XG(1,3)
        ACOOR_IN(6,3)=XG(3,3)
        XABSC_WEIG_IN(6,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(6,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(7,1)=XG(2,3)
        ACOOR_IN(7,2)=XG(3,3)
        ACOOR_IN(7,3)=XG(1,3)
        XABSC_WEIG_IN(7,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(7,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3
!
        XABSC_WEIG_IN(8,3)=WGT(4)
        XABSC_WEIG_IN(9,3)=WGT(4)
        XABSC_WEIG_IN(10,3)=WGT(4)
!
        ACOOR_IN(8,1)=XG(1,4)
        ACOOR_IN(8,2)=XG(2,4)
        ACOOR_IN(8,3)=XG(3,4)
        XABSC_WEIG_IN(8,1)=XG(1,4)*X1+XG(2,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(8,2)=XG(1,4)*Y1+XG(2,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(9,1)=XG(2,4)
        ACOOR_IN(9,2)=XG(1,4)
        ACOOR_IN(9,3)=XG(3,4)
        XABSC_WEIG_IN(9,1)=XG(2,4)*X1+XG(1,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(9,2)=XG(2,4)*Y1+XG(1,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(10,1)=XG(2,4)
        ACOOR_IN(10,2)=XG(3,4)
        ACOOR_IN(10,3)=XG(1,4)
        XABSC_WEIG_IN(10,1)=XG(2,4)*X1+XG(3,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(10,2)=XG(2,4)*Y1+XG(3,4)*Y2+XG(1,4)*Y3
!
        XABSC_WEIG_IN(11,3)=WGT(5)
        XABSC_WEIG_IN(12,3)=WGT(5)
        XABSC_WEIG_IN(13,3)=WGT(5)
        XABSC_WEIG_IN(14,3)=WGT(5)
        XABSC_WEIG_IN(15,3)=WGT(5)
        XABSC_WEIG_IN(16,3)=WGT(5)
!
        ACOOR_IN(11,1)=XG(1,5)
        ACOOR_IN(11,2)=XG(2,5)
        ACOOR_IN(11,3)=XG(3,5)
        XABSC_WEIG_IN(11,1)=XG(1,5)*X1+XG(2,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(11,2)=XG(1,5)*Y1+XG(2,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(12,1)=XG(1,5)
        ACOOR_IN(12,2)=XG(3,5)
        ACOOR_IN(12,3)=XG(2,5)
        XABSC_WEIG_IN(12,1)=XG(1,5)*X1+XG(3,5)*X2+XG(2,5)*X3
        XABSC_WEIG_IN(12,2)=XG(1,5)*Y1+XG(3,5)*Y2+XG(2,5)*Y3
!
        ACOOR_IN(13,1)=XG(2,5)
        ACOOR_IN(13,2)=XG(3,5)
        ACOOR_IN(13,3)=XG(1,5)
        XABSC_WEIG_IN(13,1)=XG(2,5)*X1+XG(3,5)*X2+XG(1,5)*X3
        XABSC_WEIG_IN(13,2)=XG(2,5)*Y1+XG(3,5)*Y2+XG(1,5)*Y3
!
        ACOOR_IN(14,1)=XG(2,5)
        ACOOR_IN(14,2)=XG(1,5)
        ACOOR_IN(14,3)=XG(3,5)
        XABSC_WEIG_IN(14,1)=XG(2,5)*X1+XG(1,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(14,2)=XG(2,5)*Y1+XG(1,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(15,1)=XG(3,5)
        ACOOR_IN(15,2)=XG(1,5)
        ACOOR_IN(15,3)=XG(2,5)
        XABSC_WEIG_IN(15,1)=XG(3,5)*X1+XG(1,5)*X2+XG(2,5)*X3
        XABSC_WEIG_IN(15,2)=XG(3,5)*Y1+XG(1,5)*Y2+XG(2,5)*Y3
!
        ACOOR_IN(16,1)=XG(3,5)
        ACOOR_IN(16,2)=XG(2,5)
        ACOOR_IN(16,3)=XG(1,5)
        XABSC_WEIG_IN(16,1)=XG(3,5)*X1+XG(2,5)*X2+XG(1,5)*X3
        XABSC_WEIG_IN(16,2)=XG(3,5)*Y1+XG(2,5)*Y2+XG(1,5)*Y3
!}}}
   CASE(9)
!{{{
        WGT(1)  = 0.097135796282799D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.031334700227139D0
        XG(1,2) = 0.020634961602525D0
        XG(2,2) = 0.489682519198738D0
        XG(3,2) = 0.489682519198738D0
!
        WGT(3)  = 0.077827541004774D0
        XG(1,3) = 0.125820817014127D0
        XG(2,3) = 0.437089591492937D0
        XG(3,3) = 0.437089591492937D0
!
        WGT(4)  = 0.079647738927210D0
        XG(1,4) = 0.623592928761935D0
        XG(2,4) = 0.188203535619033D0
        XG(3,4) = 0.188203535619033D0
!
        WGT(5)  = 0.025577675658698D0
        XG(1,5) = 0.910540973211095D0
        XG(2,5) = 0.044729513394453D0
        XG(3,5) = 0.044729513394453D0
!
        WGT(6)  = 0.043283539377289D0
        XG(1,6) = 0.036838412054736D0
        XG(2,6) = 0.221962989160766D0
        XG(3,6) = 0.741198598784498D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        XABSC_WEIG_IN(2,3)=WGT(2)
        XABSC_WEIG_IN(3,3)=WGT(2)
        XABSC_WEIG_IN(4,3)=WGT(2)
!
        ACOOR_IN(2,1)=XG(1,2)
        ACOOR_IN(2,2)=XG(2,2)
        ACOOR_IN(2,3)=XG(3,2)
        XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(3,1)=XG(2,2)
        ACOOR_IN(3,2)=XG(1,2)
        ACOOR_IN(3,3)=XG(3,2)
        XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(4,1)=XG(2,2)
        ACOOR_IN(4,2)=XG(3,2)
        ACOOR_IN(4,3)=XG(1,2)
        XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        XABSC_WEIG_IN(5,3)=WGT(3)
        XABSC_WEIG_IN(6,3)=WGT(3)
        XABSC_WEIG_IN(7,3)=WGT(3)
!
        ACOOR_IN(5,1)=XG(1,3)
        ACOOR_IN(5,2)=XG(2,3)
        ACOOR_IN(5,3)=XG(3,3)
        XABSC_WEIG_IN(5,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(5,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(6,1)=XG(2,3)
        ACOOR_IN(6,2)=XG(1,3)
        ACOOR_IN(6,3)=XG(3,3)
        XABSC_WEIG_IN(6,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(6,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(7,1)=XG(2,3)
        ACOOR_IN(7,2)=XG(3,3)
        ACOOR_IN(7,3)=XG(1,3)
        XABSC_WEIG_IN(7,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(7,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3
!
        XABSC_WEIG_IN(8,3)=WGT(4)
        XABSC_WEIG_IN(9,3)=WGT(4)
        XABSC_WEIG_IN(10,3)=WGT(4)
!
        ACOOR_IN(8,1)=XG(1,4)
        ACOOR_IN(8,2)=XG(2,4)
        ACOOR_IN(8,3)=XG(3,4)
        XABSC_WEIG_IN(8,1)=XG(1,4)*X1+XG(2,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(8,2)=XG(1,4)*Y1+XG(2,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(9,1)=XG(2,4)
        ACOOR_IN(9,2)=XG(1,4)
        ACOOR_IN(9,3)=XG(3,4)
        XABSC_WEIG_IN(9,1)=XG(2,4)*X1+XG(1,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(9,2)=XG(2,4)*Y1+XG(1,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(10,1)=XG(2,4)
        ACOOR_IN(10,2)=XG(3,4)
        ACOOR_IN(10,3)=XG(1,4)
        XABSC_WEIG_IN(10,1)=XG(2,4)*X1+XG(3,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(10,2)=XG(2,4)*Y1+XG(3,4)*Y2+XG(1,4)*Y3
!
        XABSC_WEIG_IN(11,3)=WGT(5)
        XABSC_WEIG_IN(12,3)=WGT(5)
        XABSC_WEIG_IN(13,3)=WGT(5)
!
        ACOOR_IN(11,1)=XG(1,5)
        ACOOR_IN(11,2)=XG(2,5)
        ACOOR_IN(11,3)=XG(3,5)
        XABSC_WEIG_IN(11,1)=XG(1,5)*X1+XG(2,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(11,2)=XG(1,5)*Y1+XG(2,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(12,1)=XG(2,5)
        ACOOR_IN(12,2)=XG(1,5)
        ACOOR_IN(12,3)=XG(3,5)
        XABSC_WEIG_IN(12,1)=XG(2,5)*X1+XG(1,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(12,2)=XG(2,5)*Y1+XG(1,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(13,1)=XG(2,5)
        ACOOR_IN(13,2)=XG(3,5)
        ACOOR_IN(13,3)=XG(1,5)
        XABSC_WEIG_IN(13,1)=XG(2,5)*X1+XG(3,5)*X2+XG(1,5)*X3
        XABSC_WEIG_IN(13,2)=XG(2,5)*Y1+XG(3,5)*Y2+XG(1,5)*Y3
!
        XABSC_WEIG_IN(14,3)=WGT(6)
        XABSC_WEIG_IN(15,3)=WGT(6)
        XABSC_WEIG_IN(16,3)=WGT(6)
        XABSC_WEIG_IN(17,3)=WGT(6)
        XABSC_WEIG_IN(18,3)=WGT(6)
        XABSC_WEIG_IN(19,3)=WGT(6)
!
        ACOOR_IN(14,1)=XG(1,6)
        ACOOR_IN(14,2)=XG(2,6)
        ACOOR_IN(14,3)=XG(3,6)
        XABSC_WEIG_IN(14,1)=XG(1,6)*X1+XG(2,6)*X2+XG(3,6)*X3
        XABSC_WEIG_IN(14,2)=XG(1,6)*Y1+XG(2,6)*Y2+XG(3,6)*Y3
!
        ACOOR_IN(15,1)=XG(1,6)
        ACOOR_IN(15,2)=XG(3,6)
        ACOOR_IN(15,3)=XG(2,6)
        XABSC_WEIG_IN(15,1)=XG(1,6)*X1+XG(3,6)*X2+XG(2,6)*X3
        XABSC_WEIG_IN(15,2)=XG(1,6)*Y1+XG(3,6)*Y2+XG(2,6)*Y3
!
        ACOOR_IN(16,1)=XG(2,6)
        ACOOR_IN(16,2)=XG(3,6)
        ACOOR_IN(16,3)=XG(1,6)
        XABSC_WEIG_IN(16,1)=XG(2,6)*X1+XG(3,6)*X2+XG(1,6)*X3
        XABSC_WEIG_IN(16,2)=XG(2,6)*Y1+XG(3,6)*Y2+XG(1,6)*Y3
!
        ACOOR_IN(17,1)=XG(2,6)
        ACOOR_IN(17,2)=XG(1,6)
        ACOOR_IN(17,3)=XG(3,6)
        XABSC_WEIG_IN(17,1)=XG(2,6)*X1+XG(1,6)*X2+XG(3,6)*X3
        XABSC_WEIG_IN(17,2)=XG(2,6)*Y1+XG(1,6)*Y2+XG(3,6)*Y3
!
        ACOOR_IN(18,1)=XG(3,6)
        ACOOR_IN(18,2)=XG(1,6)
        ACOOR_IN(18,3)=XG(2,6)
        XABSC_WEIG_IN(18,1)=XG(3,6)*X1+XG(1,6)*X2+XG(2,6)*X3
        XABSC_WEIG_IN(18,2)=XG(3,6)*Y1+XG(1,6)*Y2+XG(2,6)*Y3
!
        ACOOR_IN(19,1)=XG(3,6)
        ACOOR_IN(19,2)=XG(2,6)
        ACOOR_IN(19,3)=XG(1,6)
        XABSC_WEIG_IN(19,1)=XG(3,6)*X1+XG(2,6)*X2+XG(1,6)*X3
        XABSC_WEIG_IN(19,2)=XG(3,6)*Y1+XG(2,6)*Y2+XG(1,6)*Y3
!}}}
   CASE(10)
!{{{
        WGT(1)  = 0.090817990382754D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.036725957756467D0
        XG(1,2) = 0.028844733232685D0
        XG(2,2) = 0.485577633383657D0
        XG(3,2) = 0.485577633383657D0
!
        WGT(3)  = 0.045321059435528D0
        XG(1,3) = 0.781036849029926D0
        XG(2,3) = 0.109481575485037D0
        XG(3,3) = 0.109481575485037D0
!
        WGT(4)  = 0.072757916845420D0
        XG(1,4) = 0.141707219414880D0
        XG(2,4) = 0.307939838764121D0
        XG(3,4) = 0.550352941820999D0
!
        WGT(5)  = 0.028327242531057D0
        XG(1,5) = 0.025003534762686D0
        XG(2,5) = 0.246672560639903D0
        XG(3,5) = 0.728323904597411D0
!
        WGT(6)  = 0.009421666963733D0
        XG(1,6) = 0.009540815400299D0
        XG(2,6) = 0.066803251012200D0
        XG(3,6) = 0.923655933587500D0
!
        ACOOR_IN(1,1)=XG(1,1)
        ACOOR_IN(1,2)=XG(2,1)
        ACOOR_IN(1,3)=XG(3,1)
        XABSC_WEIG_IN(1,3)=WGT(1)
        XABSC_WEIG_IN(1,1)=XG(1,1)*X1+XG(2,1)*X2+XG(3,1)*X3
        XABSC_WEIG_IN(1,2)=XG(1,1)*Y1+XG(2,1)*Y2+XG(3,1)*Y3
!
        XABSC_WEIG_IN(2,3)=WGT(2)
        XABSC_WEIG_IN(3,3)=WGT(2)
        XABSC_WEIG_IN(4,3)=WGT(2)
!
        ACOOR_IN(2,1)=XG(1,2)
        ACOOR_IN(2,2)=XG(2,2)
        ACOOR_IN(2,3)=XG(3,2)
        XABSC_WEIG_IN(2,1)=XG(1,2)*X1+XG(2,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(2,2)=XG(1,2)*Y1+XG(2,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(3,1)=XG(2,2)
        ACOOR_IN(3,2)=XG(1,2)
        ACOOR_IN(3,3)=XG(3,2)
        XABSC_WEIG_IN(3,1)=XG(2,2)*X1+XG(1,2)*X2+XG(3,2)*X3
        XABSC_WEIG_IN(3,2)=XG(2,2)*Y1+XG(1,2)*Y2+XG(3,2)*Y3
!
        ACOOR_IN(4,1)=XG(2,2)
        ACOOR_IN(4,2)=XG(3,2)
        ACOOR_IN(4,3)=XG(1,2)
        XABSC_WEIG_IN(4,1)=XG(2,2)*X1+XG(3,2)*X2+XG(1,2)*X3
        XABSC_WEIG_IN(4,2)=XG(2,2)*Y1+XG(3,2)*Y2+XG(1,2)*Y3
!
        XABSC_WEIG_IN(5,3)=WGT(3)
        XABSC_WEIG_IN(6,3)=WGT(3)
        XABSC_WEIG_IN(7,3)=WGT(3)
!
        ACOOR_IN(5,1)=XG(1,3)
        ACOOR_IN(5,2)=XG(2,3)
        ACOOR_IN(5,3)=XG(3,3)
        XABSC_WEIG_IN(5,1)=XG(1,3)*X1+XG(2,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(5,2)=XG(1,3)*Y1+XG(2,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(6,1)=XG(2,3)
        ACOOR_IN(6,2)=XG(1,3)
        ACOOR_IN(6,3)=XG(3,3)
        XABSC_WEIG_IN(6,1)=XG(2,3)*X1+XG(1,3)*X2+XG(3,3)*X3
        XABSC_WEIG_IN(6,2)=XG(2,3)*Y1+XG(1,3)*Y2+XG(3,3)*Y3
!
        ACOOR_IN(7,1)=XG(2,3)
        ACOOR_IN(7,2)=XG(3,3)
        ACOOR_IN(7,3)=XG(1,3)
        XABSC_WEIG_IN(7,1)=XG(2,3)*X1+XG(3,3)*X2+XG(1,3)*X3
        XABSC_WEIG_IN(7,2)=XG(2,3)*Y1+XG(3,3)*Y2+XG(1,3)*Y3
!
        XABSC_WEIG_IN(8,3)=WGT(4)
        XABSC_WEIG_IN(9,3)=WGT(4)
        XABSC_WEIG_IN(10,3)=WGT(4)
        XABSC_WEIG_IN(11,3)=WGT(4)
        XABSC_WEIG_IN(12,3)=WGT(4)
        XABSC_WEIG_IN(13,3)=WGT(4)
!
        ACOOR_IN(8,1)=XG(1,4)
        ACOOR_IN(8,2)=XG(2,4)
        ACOOR_IN(8,3)=XG(3,4)
        XABSC_WEIG_IN(8,1)=XG(1,4)*X1+XG(2,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(8,2)=XG(1,4)*Y1+XG(2,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(9,1)=XG(1,4)
        ACOOR_IN(9,2)=XG(3,4)
        ACOOR_IN(9,3)=XG(2,4)
        XABSC_WEIG_IN(9,1)=XG(1,4)*X1+XG(3,4)*X2+XG(2,4)*X3
        XABSC_WEIG_IN(9,2)=XG(1,4)*Y1+XG(3,4)*Y2+XG(2,4)*Y3
!
        ACOOR_IN(10,1)=XG(2,4)
        ACOOR_IN(10,2)=XG(3,4)
        ACOOR_IN(10,3)=XG(1,4)
        XABSC_WEIG_IN(10,1)=XG(2,4)*X1+XG(3,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(10,2)=XG(2,4)*Y1+XG(3,4)*Y2+XG(1,4)*Y3
!
        ACOOR_IN(11,1)=XG(2,4)
        ACOOR_IN(11,2)=XG(1,4)
        ACOOR_IN(11,3)=XG(3,4)
        XABSC_WEIG_IN(11,1)=XG(2,4)*X1+XG(1,4)*X2+XG(3,4)*X3
        XABSC_WEIG_IN(11,2)=XG(2,4)*Y1+XG(1,4)*Y2+XG(3,4)*Y3
!
        ACOOR_IN(12,1)=XG(3,4)
        ACOOR_IN(12,2)=XG(1,4)
        ACOOR_IN(12,3)=XG(2,4)
        XABSC_WEIG_IN(12,1)=XG(3,4)*X1+XG(1,4)*X2+XG(2,4)*X3
        XABSC_WEIG_IN(12,2)=XG(3,4)*Y1+XG(1,4)*Y2+XG(2,4)*Y3
!
        ACOOR_IN(13,1)=XG(3,4)
        ACOOR_IN(13,2)=XG(2,4)
        ACOOR_IN(13,3)=XG(1,4)
        XABSC_WEIG_IN(13,1)=XG(3,4)*X1+XG(2,4)*X2+XG(1,4)*X3
        XABSC_WEIG_IN(13,2)=XG(3,4)*Y1+XG(2,4)*Y2+XG(1,4)*Y3
!
        XABSC_WEIG_IN(14,3)=WGT(5)
        XABSC_WEIG_IN(15,3)=WGT(5)
        XABSC_WEIG_IN(16,3)=WGT(5)
        XABSC_WEIG_IN(17,3)=WGT(5)
        XABSC_WEIG_IN(18,3)=WGT(5)
        XABSC_WEIG_IN(19,3)=WGT(5)
!
        ACOOR_IN(14,1)=XG(1,5)
        ACOOR_IN(14,2)=XG(2,5)
        ACOOR_IN(14,3)=XG(3,5)
        XABSC_WEIG_IN(14,1)=XG(1,5)*X1+XG(2,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(14,2)=XG(1,5)*Y1+XG(2,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(15,1)=XG(1,5)
        ACOOR_IN(15,2)=XG(3,5)
        ACOOR_IN(15,3)=XG(2,5)
        XABSC_WEIG_IN(15,1)=XG(1,5)*X1+XG(3,5)*X2+XG(2,5)*X3
        XABSC_WEIG_IN(15,2)=XG(1,5)*Y1+XG(3,5)*Y2+XG(2,5)*Y3
!
        ACOOR_IN(16,1)=XG(2,5)
        ACOOR_IN(16,2)=XG(3,5)
        ACOOR_IN(16,3)=XG(1,5)
        XABSC_WEIG_IN(16,1)=XG(2,5)*X1+XG(3,5)*X2+XG(1,5)*X3
        XABSC_WEIG_IN(16,2)=XG(2,5)*Y1+XG(3,5)*Y2+XG(1,5)*Y3
!
        ACOOR_IN(17,1)=XG(2,5)
        ACOOR_IN(17,2)=XG(1,5)
        ACOOR_IN(17,3)=XG(3,5)
        XABSC_WEIG_IN(17,1)=XG(2,5)*X1+XG(1,5)*X2+XG(3,5)*X3
        XABSC_WEIG_IN(17,2)=XG(2,5)*Y1+XG(1,5)*Y2+XG(3,5)*Y3
!
        ACOOR_IN(18,1)=XG(3,5)
        ACOOR_IN(18,2)=XG(1,5)
        ACOOR_IN(18,3)=XG(2,5)
        XABSC_WEIG_IN(18,1)=XG(3,5)*X1+XG(1,5)*X2+XG(2,5)*X3
        XABSC_WEIG_IN(18,2)=XG(3,5)*Y1+XG(1,5)*Y2+XG(2,5)*Y3
!
        ACOOR_IN(19,1)=XG(3,5)
        ACOOR_IN(19,2)=XG(2,5)
        ACOOR_IN(19,3)=XG(1,5)
        XABSC_WEIG_IN(19,1)=XG(3,5)*X1+XG(2,5)*X2+XG(1,5)*X3
        XABSC_WEIG_IN(19,2)=XG(3,5)*Y1+XG(2,5)*Y2+XG(1,5)*Y3
!
        XABSC_WEIG_IN(20,3)=WGT(6)
        XABSC_WEIG_IN(21,3)=WGT(6)
        XABSC_WEIG_IN(22,3)=WGT(6)
        XABSC_WEIG_IN(23,3)=WGT(6)
        XABSC_WEIG_IN(24,3)=WGT(6)
        XABSC_WEIG_IN(25,3)=WGT(6)
!
        ACOOR_IN(20,1)=XG(1,6)
        ACOOR_IN(20,2)=XG(2,6)
        ACOOR_IN(20,3)=XG(3,6)
        XABSC_WEIG_IN(20,1)=XG(1,6)*X1+XG(2,6)*X2+XG(3,6)*X3
        XABSC_WEIG_IN(20,2)=XG(1,6)*Y1+XG(2,6)*Y2+XG(3,6)*Y3
!
        ACOOR_IN(21,1)=XG(1,6)
        ACOOR_IN(21,2)=XG(3,6)
        ACOOR_IN(21,3)=XG(2,6)
        XABSC_WEIG_IN(21,1)=XG(1,6)*X1+XG(3,6)*X2+XG(2,6)*X3
        XABSC_WEIG_IN(21,2)=XG(1,6)*Y1+XG(3,6)*Y2+XG(2,6)*Y3
!
        ACOOR_IN(22,1)=XG(2,6)
        ACOOR_IN(22,2)=XG(3,6)
        ACOOR_IN(22,3)=XG(1,6)
        XABSC_WEIG_IN(22,1)=XG(2,6)*X1+XG(3,6)*X2+XG(1,6)*X3
        XABSC_WEIG_IN(22,2)=XG(2,6)*Y1+XG(3,6)*Y2+XG(1,6)*Y3
!
        ACOOR_IN(23,1)=XG(2,6)
        ACOOR_IN(23,2)=XG(1,6)
        ACOOR_IN(23,3)=XG(3,6)
        XABSC_WEIG_IN(23,1)=XG(2,6)*X1+XG(1,6)*X2+XG(3,6)*X3
        XABSC_WEIG_IN(23,2)=XG(2,6)*Y1+XG(1,6)*Y2+XG(3,6)*Y3
!
        ACOOR_IN(24,1)=XG(3,6)
        ACOOR_IN(24,2)=XG(1,6)
        ACOOR_IN(24,3)=XG(2,6)
        XABSC_WEIG_IN(24,1)=XG(3,6)*X1+XG(1,6)*X2+XG(2,6)*X3
        XABSC_WEIG_IN(24,2)=XG(3,6)*Y1+XG(1,6)*Y2+XG(2,6)*Y3
!
        ACOOR_IN(25,1)=XG(3,6)
        ACOOR_IN(25,2)=XG(2,6)
        ACOOR_IN(25,3)=XG(1,6)
        XABSC_WEIG_IN(25,1)=XG(3,6)*X1+XG(2,6)*X2+XG(1,6)*X3
        XABSC_WEIG_IN(25,2)=XG(3,6)*Y1+XG(2,6)*Y2+XG(1,6)*Y3
!}}}

   CASE DEFAULT
        WRITE(*,*) "ERROR, /GAUSSM/GAUTRIANGLE!"
        STOP
   END SELECT
!
   RETURN
!}}}
   END SUBROUTINE GAUTRIANGLE
!
   SUBROUTINE GAUTRIANGLE_OFE1(NGDOF,NN,AREA_COOR,GAUSS_WEI)
!{{{
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: NGDOF,NN
   REAL(8),DIMENSION(NN,3) :: AREA_COOR
   REAL(8),DIMENSION(NN)   :: GAUSS_WEI
!  ASSISTANT VARIABLES
   INTEGER :: I
   REAL(8) :: XG(3,79),WGT(79)
!
   SELECT CASE(NGDOF)
   CASE(1)
        WGT(1)  = 1.0D0
        XG(1,1) = 0.166666666666667D0
        XG(2,1) = 0.666666666666667D0
        XG(3,1) = 0.166666666666667D0
!
        GAUSS_WEI(:)=WGT(1)
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
   CASE(2)
!{{{
        WGT(1)  = 0.333333333333333D0
        XG(1,1) = 0.666666666666667D0
        XG(2,1) = 0.166666666666667D0
        XG(3,1) = 0.166666666666667D0
!
        GAUSS_WEI(:)=WGT(1)
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
!
        AREA_COOR(2,1)=XG(2,1)
        AREA_COOR(2,2)=XG(1,1)
        AREA_COOR(2,3)=XG(3,1)
!
        AREA_COOR(3,1)=XG(2,1)
        AREA_COOR(3,2)=XG(3,1)
        AREA_COOR(3,3)=XG(1,1)
!}}}
   CASE(3)
!{{{
!
        WGT(1)  = 0.166666666666667D0
        XG(1,1) = 0.659027622374092D0
        XG(2,1) = 0.231933368553031D0
        XG(3,1) = 0.109039009072877D0
!
        GAUSS_WEI(:)=WGT(1)
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
!
        AREA_COOR(2,1)=XG(1,1)
        AREA_COOR(2,2)=XG(3,1)
        AREA_COOR(2,3)=XG(2,1)
!
        AREA_COOR(3,1)=XG(2,1)
        AREA_COOR(3,2)=XG(1,1)
        AREA_COOR(3,3)=XG(3,1)
!
        AREA_COOR(4,1)=XG(2,1)
        AREA_COOR(4,2)=XG(3,1)
        AREA_COOR(4,3)=XG(1,1)
!
        AREA_COOR(5,1)=XG(3,1)
        AREA_COOR(5,2)=XG(2,1)
        AREA_COOR(5,3)=XG(1,1)
!
        AREA_COOR(6,1)=XG(3,1)
        AREA_COOR(6,2)=XG(1,1)
        AREA_COOR(6,3)=XG(2,1)
!}}}
   CASE(4)
!{{{
        WGT(1)  = 0.375000000000000D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
        WGT(2)  = 0.104166666666667D0
        XG(1,2) = 0.736712498968435D0
        XG(2,2) = 0.237932366472434D0
        XG(3,2) = 0.025355134559132D0
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
        GAUSS_WEI(1)  =WGT(1)
!
        GAUSS_WEI(2:7)=WGT(2)
!
        AREA_COOR(2,1)=XG(1,2)
        AREA_COOR(2,2)=XG(2,2)
        AREA_COOR(2,3)=XG(3,2)
!
        AREA_COOR(3,1)=XG(1,2)
        AREA_COOR(3,2)=XG(3,2)
        AREA_COOR(3,3)=XG(2,2)
!
        AREA_COOR(4,1)=XG(2,2)
        AREA_COOR(4,2)=XG(3,2)
        AREA_COOR(4,3)=XG(1,2)
!
        AREA_COOR(5,1)=XG(2,2)
        AREA_COOR(5,2)=XG(1,2)
        AREA_COOR(5,3)=XG(3,2)
!
        AREA_COOR(6,1)=XG(3,2)
        AREA_COOR(6,2)=XG(1,2)
        AREA_COOR(6,3)=XG(2,2)
!
        AREA_COOR(7,1)=XG(3,2)
        AREA_COOR(7,2)=XG(2,2)
        AREA_COOR(7,3)=XG(1,2)
!}}}
   CASE(5)
!{{{
        WGT(1)  = 0.205950504760887D0
        XG(1,1) = 0.124949503233232D0
        XG(2,1) = 0.437525248383384D0
        XG(3,1) = 0.437525248383384D0
        WGT(2)  = 0.063691414286223D0
        XG(1,2) = 0.797112651860071D0
        XG(2,2) = 0.165409927389841D0
        XG(3,2) = 0.037477420750088D0
!
        GAUSS_WEI(1:3)=WGT(1)
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
!
        AREA_COOR(2,1)=XG(2,1)
        AREA_COOR(2,2)=XG(1,1)
        AREA_COOR(2,3)=XG(3,1)
!
        AREA_COOR(3,1)=XG(2,1)
        AREA_COOR(3,2)=XG(3,1)
        AREA_COOR(3,3)=XG(1,1)
!
        GAUSS_WEI(4:9)=WGT(2)
!
        AREA_COOR(4,1)=XG(1,2)
        AREA_COOR(4,2)=XG(2,2)
        AREA_COOR(4,3)=XG(3,2)
!
        AREA_COOR(5,1)=XG(1,2)
        AREA_COOR(5,2)=XG(3,2)
        AREA_COOR(5,3)=XG(2,2)
 !
        AREA_COOR(6,1)=XG(2,2)
        AREA_COOR(6,2)=XG(3,2)
        AREA_COOR(6,3)=XG(1,2)
 !
        AREA_COOR(7,1)=XG(2,2)
        AREA_COOR(7,2)=XG(1,2)
        AREA_COOR(7,3)=XG(3,2)
!
        AREA_COOR(8,1)=XG(3,2)
        AREA_COOR(8,2)=XG(1,2)
        AREA_COOR(8,3)=XG(2,2)
!
        AREA_COOR(9,1)=XG(3,2)
        AREA_COOR(9,2)=XG(2,2)
        AREA_COOR(9,3)=XG(1,2)
!}}}
   CASE(6)
!{{{
        WGT(1)  = 0.050844906370207D0
        XG(1,1) = 0.873821971016996D0
        XG(2,1) = 0.063089014491502D0
        XG(3,1) = 0.063089014491502D0
        WGT(2)  = 0.116786275726379D0
        XG(1,2) = 0.501426509658179D0
        XG(2,2) = 0.249286745170910D0
        XG(3,2) = 0.249286745170911D0
        WGT(3)  = 0.082851075618374D0
        XG(1,3) = 0.636502499121399D0
        XG(2,3) = 0.310352451033785D0
        XG(3,3) = 0.053145049844816D0
!
        GAUSS_WEI(1:3)=WGT(1)
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
!
        AREA_COOR(2,1)=XG(2,1)
        AREA_COOR(2,2)=XG(1,1)
        AREA_COOR(2,3)=XG(3,1)
!
        AREA_COOR(3,1)=XG(2,1)
        AREA_COOR(3,2)=XG(3,1)
        AREA_COOR(3,3)=XG(1,1)
!
        GAUSS_WEI(4:6)=WGT(2)
!
        AREA_COOR(4,1)=XG(1,2)
        AREA_COOR(4,2)=XG(2,2)
        AREA_COOR(4,3)=XG(3,2)
!
        AREA_COOR(5,1)=XG(2,2)
        AREA_COOR(5,2)=XG(1,2)
        AREA_COOR(5,3)=XG(3,2)
!
        AREA_COOR(6,1)=XG(2,2)
        AREA_COOR(6,2)=XG(3,2)
        AREA_COOR(6,3)=XG(1,2)
!
        GAUSS_WEI(7:12)=WGT(3)
!
        AREA_COOR(7,1)=XG(1,3)
        AREA_COOR(7,2)=XG(2,3)
        AREA_COOR(7,3)=XG(3,3)
!
        AREA_COOR(8,1)=XG(1,3)
        AREA_COOR(8,2)=XG(3,3)
        AREA_COOR(8,3)=XG(2,3)
!
        AREA_COOR(9,1)=XG(2,3)
        AREA_COOR(9,2)=XG(3,3)
        AREA_COOR(9,3)=XG(1,3)
!
        AREA_COOR(10,1)=XG(2,3)
        AREA_COOR(10,2)=XG(1,3)
        AREA_COOR(10,3)=XG(3,3)
!
        AREA_COOR(11,1)=XG(3,3)
        AREA_COOR(11,2)=XG(1,3)
        AREA_COOR(11,3)=XG(2,3)
!
        AREA_COOR(12,1)=XG(3,3)
        AREA_COOR(12,2)=XG(2,3)
        AREA_COOR(12,3)=XG(1,3)
!}}}
   CASE(7)
!{{{
        WGT(1) = -0.149570044467670D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
        WGT(2)  = 0.175615257433204D0
        XG(1,2) = 0.479308067841923D0
        XG(2,2) = 0.260345966079038D0
        XG(3,2) = 0.260345966079038D0
        WGT(3)   = 0.053347235608839D0
        XG(1,3) = 0.869739794195568D0
        XG(2,3) = 0.065130102902216D0
        XG(3,3) = 0.065130102902216D0
        WGT(4)   = 0.077113760890257D0
        XG(1,4) = 0.638444188569809D0
        XG(2,4) = 0.312865496004875D0
        XG(3,4) = 0.048690315425316D0
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
        GAUSS_WEI(1)=WGT(1)
!
        GAUSS_WEI(2:4)=WGT(2)
!
        AREA_COOR(2,1)=XG(1,2)
        AREA_COOR(2,2)=XG(2,2)
        AREA_COOR(2,3)=XG(3,2)
!
        AREA_COOR(3,1)=XG(2,2)
        AREA_COOR(3,2)=XG(1,2)
        AREA_COOR(3,3)=XG(3,2)
 !
        AREA_COOR(4,1)=XG(2,2)
        AREA_COOR(4,2)=XG(3,2)
        AREA_COOR(4,3)=XG(1,2)
!
        GAUSS_WEI(5:7)=WGT(3)
!
        AREA_COOR(5,1)=XG(1,3)
        AREA_COOR(5,2)=XG(2,3)
        AREA_COOR(5,3)=XG(3,3)
 !
        AREA_COOR(6,1)=XG(2,3)
        AREA_COOR(6,2)=XG(1,3)
        AREA_COOR(6,3)=XG(3,3)
!
        AREA_COOR(7,1)=XG(2,3)
        AREA_COOR(7,2)=XG(3,3)
        AREA_COOR(7,3)=XG(1,3)
!
        GAUSS_WEI(8:13)=WGT(4)
!
        AREA_COOR(8,1)=XG(1,4)
        AREA_COOR(8,2)=XG(2,4)
        AREA_COOR(8,3)=XG(3,4)
!
        AREA_COOR(9,1)=XG(1,4)
        AREA_COOR(9,2)=XG(3,4)
        AREA_COOR(9,3)=XG(2,4)
!
        AREA_COOR(10,1)=XG(2,4)
        AREA_COOR(10,2)=XG(3,4)
        AREA_COOR(10,3)=XG(1,4)
!
        AREA_COOR(11,1)=XG(2,4)
        AREA_COOR(11,2)=XG(1,4)
        AREA_COOR(11,3)=XG(3,4)
!
        AREA_COOR(12,1)=XG(3,4)
        AREA_COOR(12,2)=XG(1,4)
        AREA_COOR(12,3)=XG(2,4)
!
        AREA_COOR(13,1)=XG(3,4)
        AREA_COOR(13,2)=XG(2,4)
        AREA_COOR(13,3)=XG(1,4)
!}}}
   CASE(8)
!{{{
        WGT(1)  = 0.144315607677787D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.095091634267285D0
        XG(1,2) = 0.081414823414554D0
        XG(2,2) = 0.459292588292723D0
        XG(3,2) = 0.459292588292723D0
!
        WGT(3)  = 0.103217370534718D0
        XG(1,3) = 0.658861384496480D0
        XG(2,3) = 0.170569307751760D0
        XG(3,3) = 0.170569307751760D0
!
        WGT(4)  = 0.032458497623198D0
        XG(1,4) = 0.898905543365938D0
        XG(2,4) = 0.050547228317031D0
        XG(3,4) = 0.050547228317031D0
!
        WGT(5)  = 0.027230314174435D0
        XG(1,5) = 0.008394777409958D0
        XG(2,5) = 0.263112829634638D0
        XG(3,5) = 0.728492392955404D0
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
        GAUSS_WEI(1)=WGT(1)
!
        GAUSS_WEI(2:4)=WGT(2)
!
        AREA_COOR(2,1)=XG(1,2)
        AREA_COOR(2,2)=XG(2,2)
        AREA_COOR(2,3)=XG(3,2)
!
        AREA_COOR(3,1)=XG(2,2)
        AREA_COOR(3,2)=XG(1,2)
        AREA_COOR(3,3)=XG(3,2)
!
        AREA_COOR(4,1)=XG(2,2)
        AREA_COOR(4,2)=XG(3,2)
        AREA_COOR(4,3)=XG(1,2)
!
        GAUSS_WEI(5:7)=WGT(3)
!
        AREA_COOR(5,1)=XG(1,3)
        AREA_COOR(5,2)=XG(2,3)
        AREA_COOR(5,3)=XG(3,3)
!
        AREA_COOR(6,1)=XG(2,3)
        AREA_COOR(6,2)=XG(1,3)
        AREA_COOR(6,3)=XG(3,3)
!
        AREA_COOR(7,1)=XG(2,3)
        AREA_COOR(7,2)=XG(3,3)
        AREA_COOR(7,3)=XG(1,3)
!
        GAUSS_WEI(8:10)=WGT(4)
!
        AREA_COOR(8,1)=XG(1,4)
        AREA_COOR(8,2)=XG(2,4)
        AREA_COOR(8,3)=XG(3,4)
 !
        AREA_COOR(9,1)=XG(2,4)
        AREA_COOR(9,2)=XG(1,4)
        AREA_COOR(9,3)=XG(3,4)
!
        AREA_COOR(10,1)=XG(2,4)
        AREA_COOR(10,2)=XG(3,4)
        AREA_COOR(10,3)=XG(1,4)
!
        GAUSS_WEI(11:16)=WGT(5)
!
        AREA_COOR(11,1)=XG(1,5)
        AREA_COOR(11,2)=XG(2,5)
        AREA_COOR(11,3)=XG(3,5)
!
        AREA_COOR(12,1)=XG(1,5)
        AREA_COOR(12,2)=XG(3,5)
        AREA_COOR(12,3)=XG(2,5)
!
        AREA_COOR(13,1)=XG(2,5)
        AREA_COOR(13,2)=XG(3,5)
        AREA_COOR(13,3)=XG(1,5)
!
        AREA_COOR(14,1)=XG(2,5)
        AREA_COOR(14,2)=XG(1,5)
        AREA_COOR(14,3)=XG(3,5)
!
        AREA_COOR(15,1)=XG(3,5)
        AREA_COOR(15,2)=XG(1,5)
        AREA_COOR(15,3)=XG(2,5)
!
        AREA_COOR(16,1)=XG(3,5)
        AREA_COOR(16,2)=XG(2,5)
        AREA_COOR(16,3)=XG(1,5)
!}}}
   CASE(9)
!{{{
        WGT(1)  = 0.097135796282799D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.031334700227139D0
        XG(1,2) = 0.020634961602525D0
        XG(2,2) = 0.489682519198738D0
        XG(3,2) = 0.489682519198738D0
!
        WGT(3)  = 0.077827541004774D0
        XG(1,3) = 0.125820817014127D0
        XG(2,3) = 0.437089591492937D0
        XG(3,3) = 0.437089591492937D0
!
        WGT(4)  = 0.079647738927210D0
        XG(1,4) = 0.623592928761935D0
        XG(2,4) = 0.188203535619033D0
        XG(3,4) = 0.188203535619033D0
!
        WGT(5)  = 0.025577675658698D0
        XG(1,5) = 0.910540973211095D0
        XG(2,5) = 0.044729513394453D0
        XG(3,5) = 0.044729513394453D0
!
        WGT(6)  = 0.043283539377289D0
        XG(1,6) = 0.036838412054736D0
        XG(2,6) = 0.221962989160766D0
        XG(3,6) = 0.741198598784498D0
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
        GAUSS_WEI(1)=WGT(1)
!
        GAUSS_WEI(2:4)=WGT(2)
!
        AREA_COOR(2,1)=XG(1,2)
        AREA_COOR(2,2)=XG(2,2)
        AREA_COOR(2,3)=XG(3,2)
!
        AREA_COOR(3,1)=XG(2,2)
        AREA_COOR(3,2)=XG(1,2)
        AREA_COOR(3,3)=XG(3,2)
!
        AREA_COOR(4,1)=XG(2,2)
        AREA_COOR(4,2)=XG(3,2)
        AREA_COOR(4,3)=XG(1,2)
!
        GAUSS_WEI(5:7)=WGT(3)
!
        AREA_COOR(5,1)=XG(1,3)
        AREA_COOR(5,2)=XG(2,3)
        AREA_COOR(5,3)=XG(3,3)
!
        AREA_COOR(6,1)=XG(2,3)
        AREA_COOR(6,2)=XG(1,3)
        AREA_COOR(6,3)=XG(3,3)
!
        AREA_COOR(7,1)=XG(2,3)
        AREA_COOR(7,2)=XG(3,3)
        AREA_COOR(7,3)=XG(1,3)
!
        GAUSS_WEI(8:10)=WGT(4)
!
        AREA_COOR(8,1)=XG(1,4)
        AREA_COOR(8,2)=XG(2,4)
        AREA_COOR(8,3)=XG(3,4)
!
        AREA_COOR(9,1)=XG(2,4)
        AREA_COOR(9,2)=XG(1,4)
        AREA_COOR(9,3)=XG(3,4)
!
        AREA_COOR(10,1)=XG(2,4)
        AREA_COOR(10,2)=XG(3,4)
        AREA_COOR(10,3)=XG(1,4)
!
        GAUSS_WEI(11:13)=WGT(5)
!
        AREA_COOR(11,1)=XG(1,5)
        AREA_COOR(11,2)=XG(2,5)
        AREA_COOR(11,3)=XG(3,5)
!
        AREA_COOR(12,1)=XG(2,5)
        AREA_COOR(12,2)=XG(1,5)
        AREA_COOR(12,3)=XG(3,5)
!
        AREA_COOR(13,1)=XG(2,5)
        AREA_COOR(13,2)=XG(3,5)
        AREA_COOR(13,3)=XG(1,5)
!
        GAUSS_WEI(14:19)=WGT(6)
!
        AREA_COOR(14,1)=XG(1,6)
        AREA_COOR(14,2)=XG(2,6)
        AREA_COOR(14,3)=XG(3,6)
!
        AREA_COOR(15,1)=XG(1,6)
        AREA_COOR(15,2)=XG(3,6)
        AREA_COOR(15,3)=XG(2,6)
!
        AREA_COOR(16,1)=XG(2,6)
        AREA_COOR(16,2)=XG(3,6)
        AREA_COOR(16,3)=XG(1,6)
!
        AREA_COOR(17,1)=XG(2,6)
        AREA_COOR(17,2)=XG(1,6)
        AREA_COOR(17,3)=XG(3,6)
!
        AREA_COOR(18,1)=XG(3,6)
        AREA_COOR(18,2)=XG(1,6)
        AREA_COOR(18,3)=XG(2,6)
!
        AREA_COOR(19,1)=XG(3,6)
        AREA_COOR(19,2)=XG(2,6)
        AREA_COOR(19,3)=XG(1,6)
!}}}
   CASE(10)
!{{{
        WGT(1)  = 0.090817990382754D0
        XG(1,1) = 0.333333333333333D0
        XG(2,1) = 0.333333333333333D0
        XG(3,1) = 0.333333333333333D0
!
        WGT(2)  = 0.036725957756467D0
        XG(1,2) = 0.028844733232685D0
        XG(2,2) = 0.485577633383657D0
        XG(3,2) = 0.485577633383657D0
!
        WGT(3)  = 0.045321059435528D0
        XG(1,3) = 0.781036849029926D0
        XG(2,3) = 0.109481575485037D0
        XG(3,3) = 0.109481575485037D0
!
        WGT(4)  = 0.072757916845420D0
        XG(1,4) = 0.141707219414880D0
        XG(2,4) = 0.307939838764121D0
        XG(3,4) = 0.550352941820999D0
!
        WGT(5)  = 0.028327242531057D0
        XG(1,5) = 0.025003534762686D0
        XG(2,5) = 0.246672560639903D0
        XG(3,5) = 0.728323904597411D0
!
        WGT(6)  = 0.009421666963733D0
        XG(1,6) = 0.009540815400299D0
        XG(2,6) = 0.066803251012200D0
        XG(3,6) = 0.923655933587500D0
!
        AREA_COOR(1,1)=XG(1,1)
        AREA_COOR(1,2)=XG(2,1)
        AREA_COOR(1,3)=XG(3,1)
        GAUSS_WEI(1)=WGT(1)
!
        GAUSS_WEI(2:4)=WGT(2)
!
        AREA_COOR(2,1)=XG(1,2)
        AREA_COOR(2,2)=XG(2,2)
        AREA_COOR(2,3)=XG(3,2)
!
        AREA_COOR(3,1)=XG(2,2)
        AREA_COOR(3,2)=XG(1,2)
        AREA_COOR(3,3)=XG(3,2)
!
        AREA_COOR(4,1)=XG(2,2)
        AREA_COOR(4,2)=XG(3,2)
        AREA_COOR(4,3)=XG(1,2)
!
        GAUSS_WEI(5:7)=WGT(3)
!
        AREA_COOR(5,1)=XG(1,3)
        AREA_COOR(5,2)=XG(2,3)
        AREA_COOR(5,3)=XG(3,3)
!
        AREA_COOR(6,1)=XG(2,3)
        AREA_COOR(6,2)=XG(1,3)
        AREA_COOR(6,3)=XG(3,3)
!
        AREA_COOR(7,1)=XG(2,3)
        AREA_COOR(7,2)=XG(3,3)
        AREA_COOR(7,3)=XG(1,3)
!
        GAUSS_WEI(8:13)=WGT(4)
!
        AREA_COOR(8,1)=XG(1,4)
        AREA_COOR(8,2)=XG(2,4)
        AREA_COOR(8,3)=XG(3,4)
!
        AREA_COOR(9,1)=XG(1,4)
        AREA_COOR(9,2)=XG(3,4)
        AREA_COOR(9,3)=XG(2,4)
!
        AREA_COOR(10,1)=XG(2,4)
        AREA_COOR(10,2)=XG(3,4)
        AREA_COOR(10,3)=XG(1,4)
!
        AREA_COOR(11,1)=XG(2,4)
        AREA_COOR(11,2)=XG(1,4)
        AREA_COOR(11,3)=XG(3,4)
!
        AREA_COOR(12,1)=XG(3,4)
        AREA_COOR(12,2)=XG(1,4)
        AREA_COOR(12,3)=XG(2,4)
!
        AREA_COOR(13,1)=XG(3,4)
        AREA_COOR(13,2)=XG(2,4)
        AREA_COOR(13,3)=XG(1,4)
!
        GAUSS_WEI(14:19)=WGT(5)
!
        AREA_COOR(14,1)=XG(1,5)
        AREA_COOR(14,2)=XG(2,5)
        AREA_COOR(14,3)=XG(3,5)
!
        AREA_COOR(15,1)=XG(1,5)
        AREA_COOR(15,2)=XG(3,5)
        AREA_COOR(15,3)=XG(2,5)
!
        AREA_COOR(16,1)=XG(2,5)
        AREA_COOR(16,2)=XG(3,5)
        AREA_COOR(16,3)=XG(1,5)
!
        AREA_COOR(17,1)=XG(2,5)
        AREA_COOR(17,2)=XG(1,5)
        AREA_COOR(17,3)=XG(3,5)
!
        AREA_COOR(18,1)=XG(3,5)
        AREA_COOR(18,2)=XG(1,5)
        AREA_COOR(18,3)=XG(2,5)
!
        AREA_COOR(19,1)=XG(3,5)
        AREA_COOR(19,2)=XG(2,5)
        AREA_COOR(19,3)=XG(1,5)
!
        GAUSS_WEI(20:25)=WGT(6)
!
        AREA_COOR(20,1)=XG(1,6)
        AREA_COOR(20,2)=XG(2,6)
        AREA_COOR(20,3)=XG(3,6)
!
        AREA_COOR(21,1)=XG(1,6)
        AREA_COOR(21,2)=XG(3,6)
        AREA_COOR(21,3)=XG(2,6)
!
        AREA_COOR(22,1)=XG(2,6)
        AREA_COOR(22,2)=XG(3,6)
        AREA_COOR(22,3)=XG(1,6)
!
        AREA_COOR(23,1)=XG(2,6)
        AREA_COOR(23,2)=XG(1,6)
        AREA_COOR(23,3)=XG(3,6)
!
        AREA_COOR(24,1)=XG(3,6)
        AREA_COOR(24,2)=XG(1,6)
        AREA_COOR(24,3)=XG(2,6)
!
        AREA_COOR(25,1)=XG(3,6)
        AREA_COOR(25,2)=XG(2,6)
        AREA_COOR(25,3)=XG(1,6)
!}}}
   CASE DEFAULT
        WRITE(*,*) "ERROR, /GAUSSM/GAUTRIANGLE!"
        STOP
   END SELECT
!
   RETURN
!}}}
   END SUBROUTINE GAUTRIANGLE_OFE1
!
!   SUBROUTINE NINTG_TETRA(NGDOF,AREA_COOR,GAUSS_WEI)
!!{{{
!   IMPLICIT NONE
!   INTEGER,INTENT(IN) :: NGDOF
!   REAL(8),DIMENSION(NITP_IN,4) :: AREA_COOR
!   REAL(8),DIMENSION(NITP_IN)   :: GAUSS_WEI
!!  ASSISTANT VARIABLES
!   REAL(8) :: XG(4,79),WGT(79)
!!
!   SELECT CASE(NGDOF)
!   CASE(4)
!!{{{
!!       6
!        WGT(1)  = 0.317460317460317450D-2
!        XG(1,1) = 0.5D0
!        XG(2,1) = 0.5D0
!        XG(3,1) = 0.0D0
!        XG(4,1) = 0.0D0
!!       4  
!        WGT(2)  = 0.147649707904967828D-1
!        XG(1,2) = 0.100526765225204467D0
!        XG(2,2) = 0.100526765225204467D0
!        XG(3,2) = 0.100526765225204467D0
!        XG(4,2) = 0.698419704324386603D0
!!       4
!        WGT(3)  = 0.221397911142651221D-1
!        XG(1,3) = 0.314372873493192195D0
!        XG(2,3) = 0.314372873493192195D0
!        XG(3,3) = 0.314372873493192195D0
!        XG(4,3) = 0.568813795204234229D-1
!!
!        CALL SET_VALUE(WGT(1),XG(:,1),AREA_COOR,GAUSS_WEI,1,6)
!        CALL SET_VALUE(WGT(2),XG(:,2),AREA_COOR,GAUSS_WEI,7,10)
!        CALL SET_VALUE(WGT(3),XG(:,3),AREA_COOR,GAUSS_WEI,11,14)
!!}}}
!   CASE(5)
!!{{{
!!       4
!        WGT(1)  = 0.602678571428571597D-2
!        XG(1,1) = 0.333333333333333333D0
!        XG(2,1) = 0.333333333333333333D0
!        XG(3,1) = 0.333333333333333333D0
!        XG(4,1) = 0.0D0
!!       1  
!        WGT(2)  = 0.302836780970891856D-1
!        XG(1,2) = 0.25D0
!        XG(2,2) = 0.25D0
!        XG(3,2) = 0.25D0
!        XG(4,2) = 0.25D0
!!       4
!        WGT(3)  = 0.116452490860289742D-1
!        XG(1,3) = 0.909090909090909091D-1
!        XG(2,3) = 0.909090909090909091D-1
!        XG(3,3) = 0.909090909090909091D-1
!        XG(4,3) = 0.727272727272727273D0
!!       6
!        WGT(4)  = 0.109491415613864534D-1
!        XG(1,4) = 0.665501535736642813D-1
!        XG(2,4) = 0.665501535736642813D-1
!        XG(3,4) = 0.433449846426335728D0
!        XG(4,4) = 0.433449846426335728D0
!!
!        CALL SET_VALUE(WGT(1),XG(:,1),AREA_COOR,GAUSS_WEI,1,4)
!        CALL SET_VALUE(WGT(2),XG(:,2),AREA_COOR,GAUSS_WEI,5,5)
!        CALL SET_VALUE(WGT(3),XG(:,3),AREA_COOR,GAUSS_WEI,6,9)
!        CALL SET_VALUE(WGT(4),XG(:,4),AREA_COOR,GAUSS_WEI,10,15)
!!}}}
!   CASE(6)
!!{{{
!!       4
!        WGT(1)  = 0.665379170969464506D-2
!        XG(1,1) = 0.214602871259151684D0
!        XG(2,1) = 0.214602871259151684D0
!        XG(3,1) = 0.214602871259151684D0
!        XG(4,1) = 0.356191386222544953D0
!!       4  
!        WGT(2)  = 0.167953517588677620D-2
!        XG(1,2) = 0.406739585346113397D-1
!        XG(2,2) = 0.406739585346113397D-1
!        XG(3,2) = 0.406739585346113397D-1
!        XG(4,2) = 0.877978124396165982D0
!!       4
!        WGT(3)  = 0.922619692394239843D-2
!        XG(1,3) = 0.322337890142275646D0
!        XG(2,3) = 0.322337890142275646D0
!        XG(3,3) = 0.322337890142275646D0
!        XG(4,3) = 0.329863295731730594D-1
!!       12
!        WGT(4)  = 0.803571428571428248D-2
!        XG(1,4) = 0.636610018750175299D-1
!        XG(2,4) = 0.636610018750175299D-1
!        XG(3,4) = 0.269672331458315867D0
!        XG(4,4) = 0.603005664791649076D0
!!
!        CALL SET_VALUE(WGT(1),XG(:,1),AREA_COOR,GAUSS_WEI,1,4)
!        CALL SET_VALUE(WGT(2),XG(:,2),AREA_COOR,GAUSS_WEI,5,8)
!        CALL SET_VALUE(WGT(3),XG(:,3),AREA_COOR,GAUSS_WEI,9,12)
!        CALL SET_VALUE(WGT(4),XG(:,4),AREA_COOR,GAUSS_WEI,13,24)
!!}}})
!   CASE DEFAULT
!        WRITE(*,*) "ERROR, /GAUSSM/GAUTRIANGLE!"
!        STOP
!   END SELECT
!!
!   RETURN
!!}}}
!   END SUBROUTINE NINTG_TETRA
!!
!   SUBROUTINE SET_VALUE(WEI,XG,AREA_COOR,GAUSS_WEI,S,E)
!!{{{
!   IMPLICIT NONE
!!  INPUTS
!   REAL(8) :: WEI
!   REAL(8),DIMENSION(4) :: XG
!   REAL(8),DIMENSION(NITP_IN,4) :: AREA_COOR
!   REAL(8),DIMENSION(NITP_IN)   :: GAUSS_WEI
!   INTEGER :: S,E
!!
!   SELECT CASE(E-S+1)
!   CASE(1)
!       GAUSS_WEI(S)=WEI
!!
!       AREA_COOR(S,1)=XG(1)
!       AREA_COOR(S,2)=XG(2)
!       AREA_COOR(S,3)=XG(3)
!       AREA_COOR(S,4)=XG(4)
!   CASE(4)
!       GAUSS_WEI(S:E)=WEI
!!
!       AREA_COOR(S,1)=XG(1)
!       AREA_COOR(S,2)=XG(2)
!       AREA_COOR(S,3)=XG(3)
!       AREA_COOR(S,4)=XG(4)
!!
!       AREA_COOR(S+1,1)=XG(1)
!       AREA_COOR(S+1,2)=XG(2)
!       AREA_COOR(S+1,3)=XG(4)
!       AREA_COOR(S+1,4)=XG(3)
!!
!       AREA_COOR(S+2,1)=XG(1)
!       AREA_COOR(S+2,2)=XG(4)
!       AREA_COOR(S+2,3)=XG(2)
!       AREA_COOR(S+2,4)=XG(3)
!!
!       AREA_COOR(S+3,1)=XG(4)
!       AREA_COOR(S+3,2)=XG(1)
!       AREA_COOR(S+3,3)=XG(2)
!       AREA_COOR(S+3,4)=XG(3)
!   CASE(6)
!       GAUSS_WEI(S:E)=WEI
!!
!       AREA_COOR(S,1)=XG(1)
!       AREA_COOR(S,2)=XG(2)
!       AREA_COOR(S,3)=XG(3)
!       AREA_COOR(S,4)=XG(4)
!!
!       AREA_COOR(S+1,1)=XG(1)
!       AREA_COOR(S+1,2)=XG(3)
!       AREA_COOR(S+1,3)=XG(2)
!       AREA_COOR(S+1,4)=XG(4)
!!
!       AREA_COOR(S+2,1)=XG(1)
!       AREA_COOR(S+2,2)=XG(3)
!       AREA_COOR(S+2,3)=XG(4)
!       AREA_COOR(S+2,4)=XG(2)
!!
!       AREA_COOR(S+3,1)=XG(3)
!       AREA_COOR(S+3,2)=XG(4)
!       AREA_COOR(S+3,3)=XG(1)
!       AREA_COOR(S+3,4)=XG(2)
!!
!       AREA_COOR(S+4,1)=XG(3)
!       AREA_COOR(S+4,2)=XG(1)
!       AREA_COOR(S+4,3)=XG(4)
!       AREA_COOR(S+4,4)=XG(2)
!!
!       AREA_COOR(S+5,1)=XG(3)
!       AREA_COOR(S+5,2)=XG(1)
!       AREA_COOR(S+5,3)=XG(2)
!       AREA_COOR(S+5,4)=XG(4)
!   CASE(12)
!       GAUSS_WEI(S:E)=WEI
!!
!       AREA_COOR(S,1)=XG(1)
!       AREA_COOR(S,2)=XG(2)
!       AREA_COOR(S,3)=XG(3)
!       AREA_COOR(S,4)=XG(4)
!!
!       AREA_COOR(S+1,1)=XG(1)
!       AREA_COOR(S+1,2)=XG(3)
!       AREA_COOR(S+1,3)=XG(2)
!       AREA_COOR(S+1,4)=XG(4)
!!
!       AREA_COOR(S+2,1)=XG(1)
!       AREA_COOR(S+2,2)=XG(3)
!       AREA_COOR(S+2,3)=XG(4)
!       AREA_COOR(S+2,4)=XG(2)
!!
!       AREA_COOR(S+3,1)=XG(1)
!       AREA_COOR(S+3,2)=XG(4)
!       AREA_COOR(S+3,3)=XG(2)
!       AREA_COOR(S+3,4)=XG(3)
!!
!       AREA_COOR(S+4,1)=XG(1)
!       AREA_COOR(S+4,2)=XG(4)
!       AREA_COOR(S+4,3)=XG(3)
!       AREA_COOR(S+4,4)=XG(2)
!!
!       AREA_COOR(S+5,1)=XG(1)
!       AREA_COOR(S+5,2)=XG(2)
!       AREA_COOR(S+5,3)=XG(4)
!       AREA_COOR(S+5,4)=XG(3)
!!
!       AREA_COOR(S+6,1)=XG(3)
!       AREA_COOR(S+6,2)=XG(1)
!       AREA_COOR(S+6,3)=XG(2)
!       AREA_COOR(S+6,4)=XG(4)
!!
!       AREA_COOR(S+7,1)=XG(3)
!       AREA_COOR(S+7,2)=XG(4)
!       AREA_COOR(S+7,3)=XG(1)
!       AREA_COOR(S+7,4)=XG(2)
!!
!       AREA_COOR(S+8,1)=XG(3)
!       AREA_COOR(S+8,2)=XG(1)
!       AREA_COOR(S+8,3)=XG(4)
!       AREA_COOR(S+8,4)=XG(2)
!!
!       AREA_COOR(S+9,1)=XG(4)
!       AREA_COOR(S+9,2)=XG(1)
!       AREA_COOR(S+9,3)=XG(2)
!       AREA_COOR(S+9,4)=XG(3)
!!
!       AREA_COOR(S+10,1)=XG(4)
!       AREA_COOR(S+10,2)=XG(1)
!       AREA_COOR(S+10,3)=XG(3)
!       AREA_COOR(S+10,4)=XG(2)
!!
!       AREA_COOR(S+11,1)=XG(4)
!       AREA_COOR(S+11,2)=XG(3)
!       AREA_COOR(S+11,3)=XG(2)
!       AREA_COOR(S+11,4)=XG(1)
!   CASE DEFAULT
!        WRITE(*,*) "ERROR, /GAUSSM/SET_VALUE!"
!        STOP
!   END SELECT
!!}}}
!   END SUBROUTINE SET_VALUE
!}}}
end module GAUSSM
