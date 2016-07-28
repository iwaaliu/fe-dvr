      subroutine lagrange_dvr_derivative(nfun,xbounds,pt,wt,df,ddf,ke_mat)
        use accuracy_real
        implicit NONE


     integer*4, INTENT(IN)   :: nfun
     real(KIND=idp), dimension(1:2), INTENT(IN)  :: xbounds
     real(KIND=idp), dimension(1:nfun), INTENT(OUT) :: pt 
     real(KIND=idp), dimension(1:nfun), INTENT(OUT) :: wt    
     real(KIND=idp), dimension(1:nfun,1:nfun), INTENT(OUT) :: ke_mat,df,ddf

     integer :: i,j,k
     real(KIND=idp), dimension(1:2) :: endpoints
     real(KIND=idp) :: xmin, xmax
     real(KIND=idp), dimension(1:nfun)           :: scratch  
     real(KIND=idp), dimension(1:nfun)           :: gauss_weights
     real(KIND=idp), dimension(1:nfun)           :: roots
     real(KIND=idp), dimension(1:nfun,1:nfun)    :: f
     real(KIND=idp) :: A,B,tsum

!-----------------------------------------!
! Endpoint of the interval for the region !
!-----------------------------------------!
      xmin = xbounds(1)
      xmax = xbounds(2)


!-----------------------------------!
! Interval for Legendre polynomial  !
!-----------------------------------!
      endpoints(1) = -1.0_idp
      endpoints(2) =  1.0_idp


!-------------------------------------------------!
! Find zeros of nfun'th-order Legendre quadrature !
! Gauss-Lobatto  quadrature is used here.         !
!-------------------------------------------------!
     call gaussq('legendre',nfun,0.d0,0.d0,2,endpoints,scratch,roots,gauss_weights)
   

! Set up rescaled position grids
      A = 0.5_idp * ABS(xmax - xmin)
      B = 0.5_idp * (xmax + xmin)
      do k = 1,nfun
         pt(k) = A * roots(k) + B
         wt(k) = A * gauss_weights(k) 
      end do

! Generate Lagrange interpolation polynomials and their first-
! and second-order derivatives.
     call lgngr(f,df,ddf,pt,pt,nfun,nfun,'all')


! Set up kinetic energy matrix (not normalized)
! 'j' is function index; 'i' is point index

      ke_mat = 0.0_idp

      do j = 1,nfun
      do i = 1,nfun
!! Version (1)
!!@@@ this is the version of using the 2nd order derivative and Bloch operators.
!!@@@         ke_mat(i,j)=ddf(i,j)*wt(i)+f(1,i)*df(1,j)-f(nfun,i)*df(nfun,j)
!!@@@ end of the version (1).

!! Version (2)
!!### this is the version of using the 1st order derivative and no Block operators.
           tsum = 0.0_idp
        do k=1,nfun
           tsum = tsum + wt(k)*df(k,i)*df(k,j) 
        end do
           ke_mat(i,j) = -tsum
!!### end of the version (2).

      end do
      end do
 

! Multiplying the factor of -1/2.
      ke_mat = -0.5_idp * ke_mat

      return
      end subroutine lagrange_dvr_derivative
