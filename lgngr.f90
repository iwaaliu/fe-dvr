!================================================!
! Finds Lagrange interpolation polynomials       !
! of function of x and their first and second    !
! derivatives on an arbitrary grid y.            !
!================================================!
!deck lgngr

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-06  Time: 11:11:54

SUBROUTINE lgngr(p,dp,ddp,x,y,nx,ny,drctv)

!***begin prologue     lgrply
!***date written       940504   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            lagrange polynomials at arbitrary points.
!***description
!***


!***references
!***routines called
!***end prologue       lgngr


IMPLICIT INTEGER (a-z)
REAL*8, INTENT(OUT)                      :: p(ny,nx)
REAL*8, INTENT(OUT)                      :: dp(ny,nx)
REAL*8, INTENT(OUT)                      :: ddp(ny,nx)
REAL*8, INTENT(IN)                       :: x(nx)
REAL*8, INTENT(IN)                       :: y(ny)
INTEGER, INTENT(IN)                      :: nx
INTEGER, INTENT(IN)                      :: ny
CHARACTER (LEN=*), INTENT(IN)            :: drctv


REAL*8 sn, ssn, fac

CHARACTER(LEN=80) :: title



COMMON /io/ inp, iout

!     generate polynomials and derivatives with respect to x

DO  i=1,ny
  zerfac = 0
  DO  j=1,nx
    fac =  y(i) - x(j)
    IF(ABS(fac) <= 1.d-10) THEN
      zerfac = j
    END IF
  END DO
  DO  j=1,nx
    p(i,j) = 1.d0
    DO  k=1,j-1
      p(i,j) = p(i,j)*( y(i) - x(k) ) / ( x(j) - x(k) )
    END DO
    DO  k=j+1,nx
      p(i,j) = p(i,j)*( y(i) - x(k) ) / ( x(j) - x(k) )
    END DO
    IF(drctv /= 'functions-only') THEN
      IF(ABS(p(i,j)) > 1.d-10) THEN
        sn = 0.d0
        ssn = 0.d0
        DO  k=1,j-1
          fac = 1.d0/( y(i) - x(k) )
          sn = sn + fac
          ssn = ssn + fac*fac
        END DO
        DO  k=j+1,nx
          fac = 1.d0/( y(i) - x(k) )
          sn = sn + fac
          ssn = ssn + fac*fac
        END DO
        dp(i,j) = sn*p(i,j)
        ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
      ELSE
        first=j
        second=zerfac
        IF(first > second) THEN
          first=zerfac
          second=j
        END IF
        sn = 1.d0
        ssn = 0.d0
        DO  k=1,first-1
          fac = 1.d0/( x(j) - x(k) )
          sn = sn*fac*( y(i) - x(k) )
          ssn = ssn + 1.d0/(y(i) - x(k))
        END DO
        DO  k=first+1,second-1
          fac = 1.d0/( x(j) - x(k) )
          sn = sn*fac*( y(i) - x(k) )
          ssn = ssn + 1.d0/( y(i) - x(k) )
        END DO
        DO  k=second+1,nx
          fac = 1.d0/( x(j) - x(k) )
          sn = sn*fac*( y(i) - x(k) )
          ssn = ssn + 1.d0/( y(i) - x(k) )
        END DO
        dp(i,j) = sn/( x(j) - x(zerfac) )
        ddp(i,j) = 2.d0*ssn*dp(i,j)
      END IF
    END IF
  END DO
  
END DO

return
END SUBROUTINE lgngr
