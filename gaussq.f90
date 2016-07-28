!========================================================!
! This subroutine generates the Gaussian quadratures.    !
!========================================================!
!deck @(#)gaussq.f 1.1 9/9/91

! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-02  Time: 11:45:58

 SUBROUTINE gaussq(ckind, n, alpha, beta, kpts, endpts, b, t, w)

!           this set of routines computes the nodes x(i) and weights
!        c(i) for gaussian-type quadrature rules with pre-assigned
!        nodes.  these are used when one wishes to approximate

!                 integral (from a to b)  f(x) w(x) dx

!                              n
!        by                   sum c  f(x )
!                             i=1  i    i

!        here w(x) is one of six possible non-negative weight
!        functions (listed below), and f(x) is the
!        function to be integrated.  gaussian quadrature is particularly
!        useful on infinite intervals (with appropriate weight
!        functions), since then other techniques often fail.

!           associated with each weight function w(x) is a set of
!        orthogonal polynomials.  the nodes x(i) are just the zeroes
!        of the proper n-th degree polynomial.

!     input parameters

!        kind     an integer between 0 and 6 giving the type of
!                 quadrature rule

!        kind = 0=  simpson's rule w(x) = 1 on (-1, 1) n must be odd.
!        kind = 1=  legendre quadrature, w(x) = 1 on (-1, 1)
!        kind = 2=  chebyshev quadrature of the first kind
!                   w(x) = 1/dsqrt(1 - x*x) on (-1, +1)
!        kind = 3=  chebyshev quadrature of the second kind
!                   w(x) = dsqrt(1 - x*x) on (-1, 1)
!        kind = 4=  hermite quadrature, w(x) = exp(-x*x) on
!                   (-infinity, +infinity)
!        kind = 5=  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
!                   beta on (-1, 1), alpha, beta .gt. -1.
!                   note= kind=2 and 3 are a special case of this.
!        kind = 6=  generalized laguerre quadrature, w(x) = exp(-x)*
!                   x**alpha on (0, +infinity), alpha .gt. -1

!        n        the number of points used for the quadrature rule
!        alpha    real*8 parameter used only for gauss-jacobi and gauss-
!                 laguerre quadrature (otherwise use 0.).
!        beta     real*8 parameter used only for gauss-jacobi quadrature--
!                 (otherwise use 0.).
!        kpts     (integer) normally 0, unless the left or right end-
!                 point (or both) of the interval is required to be a
!                 node (this is called gauss-radau or gauss-lobatto
!                 quadrature).  then kpts is the number of fixed
!                 endpoints (1 or 2).
!        endpts   real*8 array of length 2.  contains the values of
!                 any fixed endpoints, if kpts = 1 or 2.
!        b        real*8 scratch array of length n

!     output parameters (both arrays of length n)

!        t        will contain the desired nodes x(1),,,x(n)
!        w        will contain the desired weights c(1),,,c(n)

!     subroutines required

!        gbslve, class, and gbtql2 are provided. underflow may sometimes
!        occur, but it is harmless if the underflow interrupts are
!        turned off as they are on this machine.

!     accuracy

!        the routine was tested up to n = 512 for legendre quadrature,
!        up to n = 136 for hermite, up to n = 68 for laguerre, and up
!        to n = 10 or 20 in other cases.  in all but two instances,
!        comparison with tables in ref. 3 showed 12 or more significant
!        digits of accuracy.  the two exceptions were the weights for
!        hermite and laguerre quadrature, where underflow caused some
!        very small weights to be set to zero.  this is, of course,
!        completely harmless.

!     method

!           the coefficients of the three-term recurrence relation
!        for the corresponding set of orthogonal polynomials are
!        used to form a symmetric tridiagonal matrix, whose
!        eigenvalues (determined by the implicit ql-method with
!        shifts) are just the desired nodes.  the first components of
!        the orthonormalized eigenvectors, when properly scaled,
!        yield the weights.  this technique is much faster than using a
!        root-finder to locate the zeroes of the orthogonal polynomial.
!        for further details, see ref. 1.  ref. 2 contains details of
!        gauss-radau and gauss-lobatto quadrature only.

!     references

!        1.  golub, g. h., and welsch, j. h.,  calculation of gaussian
!            quadrature rules,  mathematics of computation 23 (april,
!            1969), pp. 221-230.
!        2.  golub, g. h.,  some modified matrix eigenvalue problems,
!            siam review 15 (april, 1973), pp. 318-334 (section 7).
!        3.  stroud and secrest, gaussian quadrature formulas, prentice-
!            hall, englewood cliffs, n.j., 1966.

!     ..................................................................


IMPLICIT REAL*8 (a-h,o-z)
CHARACTER (LEN=*), INTENT(IN)            :: ckind
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN OUT)                     :: alpha
REAL*8, INTENT(IN OUT)                     :: beta
INTEGER, INTENT(IN)                      :: kpts
REAL*8, INTENT(IN)                         :: endpts(2)
REAL*8, INTENT(IN OUT)                     :: b(n)
REAL*8, INTENT(OUT)                        :: t(n)
REAL*8, INTENT(OUT)                        :: w(n)
REAL*8  muzero


COMMON/io/inp,iout

IF (ckind == 'simpson') THEN
  kind=0
ELSE IF(ckind == 'legendre') THEN
  kind=1
ELSE IF(ckind == 'chebyshev-1') THEN
  kind=2
ELSE IF(ckind == 'chebyshev-2') THEN
  kind=3
ELSE IF(ckind == 'hermite') THEN
  kind=4
ELSE IF(ckind == 'jacobi') THEN
  kind=5
ELSE IF(ckind == 'laguerre') THEN
  kind=6
ELSE
  stop 
!  CALL lnkerr('error in quadrature type')
END IF
IF(kind == 0) THEN
  IF(2*(n/2) == n) THEN
   stop
!    CALL lnkerr('n must be odd for simpson rule')
  END IF
  IF(n <= 1) THEN
    t(1) = 0.d+00
    w(1) = 2.d+00
    RETURN
  END IF
  h = 2.d+00/(n-1)
  t(1) = -1.d+00
  t(n) = 1.d+00
  w(1) = h/3.d+00
  w(n) = h/3.d+00
  nm1 = n-1
  DO  i=2,nm1
    t(i) = t(i-1) + h
    w(i) = 4.d+00 - 2.d+00*(i-2*(i/2))
    w(i) = w(i)*h/3.d+00
  END DO
  RETURN
END IF

CALL class (kind, n, alpha, beta, b, t, muzero)

!           the matrix of coefficients is assumed to be symmetric.
!           the array t contains the diagonal elements, the array
!           b the off-diagonal elements.
!           make appropriate changes in the lower right 2 by 2
!           submatrix.

IF (kpts == 0)  GO TO 100
IF (kpts == 2)  GO TO  50

!           if kpts=1, only t(n) must be changed

t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
GO TO 100

!           if kpts=2, t(n) and b(n-1) must be recomputed

50 gam =gbslve(endpts(1), n, t, b)
t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
b(n-1) =  SQRT(t1)
t(n) = endpts(1) + gam*t1

!           note that the indices of the elements of b run from 1 to n-1
!           and thus the value of b(n) is arbitrary.
!           now compute the eigenvalues of the symmetric tridiagonal
!           matrix, which has been modified as necessary.
!           the method used is a ql-type method with origin shifting

100 w(1) = 1.0D0
DO  i = 2, n
  w(i) = 0.0D0
END DO

CALL gbtql2 (n, t, b, w, ierr)
DO  i = 1, n
  w(i) = muzero * w(i) * w(i)
END DO

RETURN
END SUBROUTINE gaussq


!deck @(#)class.f 1.1 9/9/91

! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-02  Time: 11:46:10

SUBROUTINE class(kind, n, alpha, beta, b, a, muzero)

!           this procedure supplies the coefficients a(j), b(j) of the
!        recurrence relation

!             b p (x) = (x - a ) p   (x) - b   p   (x)
!              j j            j   j-1       j-1 j-2

!        for the various classical (normalized) orthogonal polynomials,
!        and the zero-th moment

!             muzero = integral w(x) dx

!        of the given polynomial   weight function w(x).  since the
!        polynomials are orthonormalized, the tridiagonal matrix is
!        guaranteed to be symmetric.

!           the input parameter alpha is used only for laguerre and
!        jacobi polynomials, and the parameter beta is used only for
!        jacobi polynomials.  the laguerre and jacobi polynomials
!        require the gamma function.

!     ..................................................................


IMPLICIT REAL*8 (a-h,o-z)
INTEGER, INTENT(IN OUT)                  :: kind
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                         :: alpha
REAL*8, INTENT(IN)                         :: beta
REAL*8, INTENT(OUT)                        :: b(n)
REAL*8, INTENT(OUT)                        :: a(n)
REAL*8, INTENT(OUT)                      :: muzero


COMMON/io/inp,iout

DATA pi / 3.141592653589793D0  /

nm1 = n - 1
SELECT CASE ( kind )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO  20
  CASE (    3)
    GO TO  30
  CASE (    4)
    GO TO  40
  CASE (    5)
    GO TO  50
  CASE (    6)
    GO TO  60
END SELECT

!              kind = 1=  legendre polynomials p(x)
!              on (-1, +1), w(x) = 1.

10 muzero = 2.0D0
DO  i = 1, nm1
  a(i) = 0.0D0
  abi = i
  b(i) = abi/ SQRT(4.d0*abi*abi - 1.0D0  )
END DO
a(n) = 0.0D0
RETURN

!              kind = 2=  chebyshev polynomials of the first kind t(x)
!              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)

20 muzero = pi
DO  i = 1, nm1
  a(i) = 0.0D0
  b(i) = 0.5D0
END DO
b(1) =  SQRT(0.5D0  )
a(n) = 0.0D0
RETURN

!              kind = 3=  chebyshev polynomials of the second kind u(x)
!              on (-1, +1), w(x) = sqrt(1 - x*x)

30 muzero = pi/2.0D0
DO  i = 1, nm1
  a(i) = 0.0D0
  b(i) = 0.5D0
END DO
a(n) = 0.0D0
RETURN

!              kind = 4=  hermite polynomials h(x) on (-infinity,
!              +infinity), w(x) = exp(-x**2)

40 muzero =  SQRT(pi)
DO  i = 1, nm1
  a(i) = 0.0D0
  b(i) =  SQRT(i/2.0D0  )
END DO
a(n) = 0.0D0
RETURN

!              kind = 5=  jacobi polynomials p(alpha, beta)(x) on
!              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
!              beta greater than -1

50 ab = alpha + beta
abi = 2.0D0   + ab
muzero = 2.0D0   ** (ab + 1.0D0  ) * gamma(alpha + 1.0D0  ) * gamma(  &
    beta + 1.0D0  ) / gamma(abi)
a(1) = (beta - alpha)/abi
b(1) =  SQRT(4.0D0  *(1.0D0  + alpha)*(1.0D0   + beta)/((abi + 1.0D0  )*  &
    abi*abi))
a2b2 = beta*beta - alpha*alpha
DO  i = 2, nm1
  abi = 2.0D0  *i + ab
  a(i) = a2b2/((abi - 2.0D0  )*abi)
  b(i) =  SQRT (4.0D0  *i*(i + alpha)*(i + beta)*(i + ab)/  &
      ((abi*abi - 1)*abi*abi))
END DO
abi = 2.0D0  *n + ab
a(n) = a2b2/((abi - 2.0D0  )*abi)
RETURN

!              kind = 6=  laguerre polynomials l(alpha)(x) on
!              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
!              than -1.

60 muzero = gamma(alpha + 1.0D0)
DO  i = 1, nm1
  a(i) = 2.0D0  *i - 1.0D0   + alpha
  b(i) =  SQRT(i*(i + alpha))
END DO
a(n) = 2.0D0  *n - 1 + alpha
RETURN
END SUBROUTINE class



FUNCTION csevl(x,cs,n)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:15:24

!***begin prologue  csevl
!***date written   770401   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  c3a2
!***keywords  chebyshev,fnlib,special function
!***author  fullerton, w., (lanl)
!***purpose  evaluate the n-term chebyshev series cs at x.
!***description

! evaluate the n-term chebyshev series cs at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox
! and parker, chebyshev polynomials in numerical analysis, oxford press,
! page 56.

!       input arguments --
! x    value at which the series is to be evaluated.
! cs   array of n terms of a chebyshev series.  in eval-
!      uating cs, only half the first coefficient is summed.
! n    number of terms in array cs.
!***references  (none)
!***routines called  lnkerr
!***end prologue  csevl


IMPLICIT REAL*8 (a-h,o-z)
REAL*8, INTENT(IN)                         :: x
REAL*8, INTENT(IN)                         :: cs(*)
INTEGER, INTENT(IN)                      :: n



!***first executable statement  csevl
IF(n < 1) THEN
  CALL lnkerr( 'csevl   number of terms le 0')
END IF
IF(n > 1000) THEN
  CALL lnkerr ( 'csevl number of terms gt 1000')
END IF
IF (x < -1.0D0 .OR. x > 1.0D0) THEN
  CALL lnkerr( 'csevl x outside (-1,1,+1)')
END IF

b1=0.d0
b0=0.d0
twox=2.d0*x
DO  i=1,n
  b2=b1
  b1=b0
  ni=n+1-i
  b0=twox*b1-b2+cs(ni)
END DO

csevl = 0.5D0 * (b0-b2)

RETURN
END FUNCTION csevl




!deck @(#)gbslve.f 1.1 9/9/91

! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-02  Time: 11:46:27

FUNCTION gbslve(shift, n, a, b)

!       this procedure performs elimination to solve for the
!       n-th component of the solution delta to the equation

!             (jn - shift*identity) * delta  = en,

!       where en is the vector of all zeroes except for 1 in
!       the n-th position.

!       the matrix jn is symmetric tridiagonal, with diagonal
!       elements a(i), off-diagonal elements b(i).  this equation
!       must be solved to obtain the appropriate changes in the lower
!       2 by 2 submatrix of coefficients for orthogonal polynomials.



IMPLICIT REAL*8 (a-h,o-z)
REAL*8, INTENT(IN)                         :: shift
INTEGER, INTENT(IN)                        :: n
REAL*8, INTENT(IN)                         :: a(n)
REAL*8, INTENT(IN)                         :: b(n)


alpha = a(1) - shift
nm1 = n - 1
DO  i = 2, nm1
  alpha = a(i) - shift - b(i-1)**2/alpha
END DO
gbslve = 1.0D0  /alpha
RETURN
END FUNCTION gbslve




!deck @(#)gbtql2.f 1.1 9/9/91

! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-02  Time: 11:46:57

SUBROUTINE gbtql2(n, d, e, z, ierr)

!     this subroutine is a translation of the algol procedure imtql2,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).

!     this subroutine finds the eigenvalues and first components of the
!     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
!     method, and is adapted from the eispak routine imtql2

!     on input=

!        n is the order of the matrix;

!        d contains the diagonal elements of the input matrix;

!        e contains the subdiagonal elements of the input matrix
!          in its first n-1 positions.  e(n) is arbitrary;

!        z contains the first row of the identity matrix.

!      on output=

!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., ierr-1;

!        e has been destroyed;

!        z contains the first components of the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.  if an error exit is
!          made, z contains the eigenvectors associated with the stored
!          eigenvalues;

!        ierr is set to

!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.

!     ------------------------------------------------------------------


IMPLICIT REAL*8 (a-h,o-z)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN OUT)                     :: d(n)
REAL*8, INTENT(OUT)                        :: e(n)
REAL*8, INTENT(IN OUT)                     :: z(n)
INTEGER, INTENT(OUT)                     :: ierr
INTEGER :: i, j, k, l, m, ii, mml
REAL*8  machep
COMMON/io/inp,iout


!     ========== machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!                machep = 16.0d0**(-13) for long form arithmetic
!                on s360 ==========
machep=1.0E-14

ierr = 0
IF (n == 1) GO TO 1001

e(n) = 0.0D0
DO  l = 1, n
  j = 0
!     ========== look for small sub-diagonal element ==========
  105    DO  m = l, n
    IF (m == n) GO TO 120
    IF ( ABS(e(m)) <= machep * ( ABS(d(m)) +  ABS(d(m+1)))) GO TO 120
  END DO
  
  120    p = d(l)
  IF (m == l) CYCLE
  IF (j == 30) GO TO 1000
  j = j + 1
!     ========== form shift ==========
  g = (d(l+1) - p) / (2.0D0   * e(l))
  r =  SQRT(g*g+1.0D0  )
  g = d(m) - p + e(l) / (g +  SIGN(r, g))
  s = 1.0D0
  c = 1.0D0
  p = 0.0D0
  mml = m - l
!     ========== for i=m-1 step -1 until l do -- ==========
  DO  ii = 1, mml
    i = m - ii
    f = s * e(i)
    b = c * e(i)
    IF ( ABS(f) <  ABS(g)) GO TO 150
    c = g / f
    r =  SQRT(c*c+1.0D0  )
    e(i+1) = f * r
    s = 1.0D0   / r
    c = c * s
    GO TO 160
    150       s = f / g
    r =  SQRT(s*s+1.0D0  )
    e(i+1) = g * r
    c = 1.0D0   / r
    s = s * c
    160       g = d(i+1) - p
    r = (d(i) - g) * s + 2.0D0   * c * b
    p = s * r
    d(i+1) = g + p
    g = c * r - b
!     ========== form first component of vector ==========
    f = z(i+1)
    z(i+1) = s * z(i) + c * f
    z(i) = c * z(i) - s * f
    
  END DO
  
  d(l) = d(l) - p
  e(l) = g
  e(m) = 0.0D0
  GO TO 105
END DO
!     ========== order eigenvalues and eigenvectors ==========
DO  ii = 2, n
  i = ii - 1
  k = i
  p = d(i)
  
  DO  j = ii, n
    IF (d(j) >= p) CYCLE
    k = j
    p = d(j)
  END DO
  
  IF (k == i) CYCLE
  d(k) = d(i)
  d(i) = p
  
  p = z(i)
  z(i) = z(k)
  z(k) = p
  
END DO

GO TO 1001
!     ========== set error -- no convergence to an
!                eigenvalue after 30 iterations ==========
1000 ierr = l
1001 RETURN
!     ========== last card of gbtql2 ==========
END SUBROUTINE gbtql2



subroutine lnkerr(message)

CHARACTER (LEN=*), INTENT(IN)            :: message

write(*,*) message
stop

end subroutine lnkerr


FUNCTION gamma(x)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:18:39

!***begin prologue  gamma
!***date written   770601   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  c7a
!***keywords  gamma function,special function
!***author  fullerton, w., (lanl)
!***purpose  computes the gamma function.
!***description

! gamma computes the gamma function at x, where x is not 0, -1, -2, ....
! gamma and x are single precision.
!***references  (none)
!***routines called  csevl,gamlim,inits,r1mach,r9lgmc,lnkerr
!***end prologue  gamma

IMPLICIT REAL*8(a-h,o-z)
REAL*8, INTENT(IN)                         :: x
DIMENSION gcs(23)
DATA gcs   ( 1) / .008571195590989331D0/
DATA gcs   ( 2) / .004415381324841007D0/
DATA gcs   ( 3) / .05685043681599363D0/
DATA gcs   ( 4) /-.004219835396418561D0/
DATA gcs   ( 5) / .001326808181212460D0/
DATA gcs   ( 6) /-.0001893024529798880D0/
DATA gcs   ( 7) / .0000360692532744124D0/
DATA gcs   ( 8) /-.0000060567619044608D0/
DATA gcs   ( 9) / .0000010558295463022D0/
DATA gcs   (10) /-.0000001811967365542D0/
DATA gcs   (11) / .0000000311772496471D0/
DATA gcs   (12) /-.0000000053542196390D0/
DATA gcs   (13) / .0000000009193275519D0/
DATA gcs   (14) /-.0000000001577941280D0/
DATA gcs   (15) / .0000000000270798062D0/
DATA gcs   (16) /-.0000000000046468186D0/
DATA gcs   (17) / .0000000000007973350D0/
DATA gcs   (18) /-.0000000000001368078D0/
DATA gcs   (19) / .0000000000000234731D0/
DATA gcs   (20) /-.0000000000000040274D0/
DATA gcs   (21) / .0000000000000006910D0/
DATA gcs   (22) /-.0000000000000001185D0/
DATA gcs   (23) / .0000000000000000203D0/
DATA pi /3.14159265358979324D0/
! sq2pil is alog (sqrt (2.*pi) )
DATA sq2pil /0.91893853320467274D0/
DATA ngcs, xmin, xmax, dxrel /0, 3*0.0D0 /

! lanl dependent code removed 81.02.04

!***first executable statement  gamma
IF (ngcs /= 0) GO TO 10

! ---------------------------------------------------------------------
! initialize.  find legal bounds for x, and determine the number of
! terms in the series required to attain an accuracy ten times better
! than machine precision.

ngcs = inits (gcs, 23, 0.1D0*r1mach(3))

CALL gamlim (xmin, xmax)
dxrel = SQRT (r1mach(4))

! ---------------------------------------------------------------------
! finish initialization.  start evaluating gamma(x).

10   y = ABS(x)
IF (y > 10.0D0) GO TO 50

! compute gamma(x) for abs(x) .le. 10.0.  reduce interval and
! find gamma(1+y) for 0. .le. y .lt. 1. first of all.

n = x
IF (x < 0.d0) n = n - 1
y = x - FLOAT(n)
n = n - 1
gamma = 0.9375D0 + csevl(2.d0*y-1.d0, gcs, ngcs)
IF (n == 0) RETURN

IF (n > 0) GO TO 30

! compute gamma(x) for x .lt. 1.

n = -n
IF (x == 0.d0) CALL lnkerr ( 'gamma   x is 0')
IF (x < 0.d0 .AND. x+FLOAT(n-2) == 0.d0)  &
    CALL lnkerr (  'gamma   x is a negative integer') 
IF (x < (-0.5D0) .AND. ABS((x-AINT(x-0.5D0))/x) < dxrel) CALL  &
    lnkerr ( 'gamma   answer lt half precision because x too near negative integer')

DO  i=1,n
  gamma = gamma / (x+FLOAT(i-1))
END DO
RETURN

! gamma(x) for x .ge. 2.

30   DO  i=1,n
  gamma = (y+FLOAT(i))*gamma
END DO
RETURN

! compute gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).

50   IF (x > xmax) stop
!CALL lnkerr ( 'gamma   x so big gamma overflows')

gamma = 0.d0
IF (x < xmin) CALL lnkerr ( 'gamma   x so small gamma'// '  underflows')
IF (x < xmin) RETURN

gamma = EXP((y-0.5D0)*LOG(y) - y + sq2pil + r9lgmc(y) )
IF (x > 0.d0) RETURN

IF (ABS((x-AINT(x-0.5D0))/x) < dxrel) CALL lnkerr ( 'gamma '//  &
    'answer lt half precision, x too near negative integer')

sinpiy = SIN (pi*y)
IF (sinpiy == 0.d0) stop
!CALL lnkerr ( 'gamma   x is a negative '// 'integer')

gamma = -pi / (y*sinpiy*gamma)

RETURN
END FUNCTION gamma



FUNCTION r1mach(i)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:15:43

!***begin prologue  r1mach
!***date written   790101   (yymmdd)
!***revision date  860324   (yymmdd)
!***category no.  r1
!***keywords  machine constants
!***author  fox, p. a., (bell labs)
!           hall, a. d., (bell labs)
!           schryer, n. l., (bell labs)
!***purpose  returns single precision machine dependent constants
!***description

!     r1mach can be used to obtain machine-dependent parameters
!     for the local machine environment.  it is a function
!     subroutine with one (input) argument, and can be called
!     as follows, for example

!          a = r1mach(i)

!     where i=1,...,5.  the (output) value of a above is
!     determined by the (input) value of i.  the results for
!     various values of i are discussed below.

!  single-precision machine constants
!  r1mach(1) = b**(emin-1), the smallest positive magnitude.
!  r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
!  r1mach(3) = b**(-t), the smallest relative spacing.
!  r1mach(4) = b**(1-t), the largest relative spacing.
!  r1mach(5) = log10(b)
!***references  fox, p.a., hall, a.d., schryer, n.l, *framework for
!                 a portable library*, acm transactions on mathe-
!                 matical software, vol. 4, no. 2, june 1978,
!                 pp. 177-188.
!***routines called  lnkerr
!***end prologue  r1mach


INTEGER, INTENT(IN)                      :: i
INTEGER :: small(2)
INTEGER :: large(2)
INTEGER :: right(2)
INTEGER :: diver(2)
INTEGER :: LOG10(2)

REAL*8 rmach(5),r1mach

EQUIVALENCE (rmach(1),small(1))
EQUIVALENCE (rmach(2),large(1))
EQUIVALENCE (rmach(3),right(1))
EQUIVALENCE (rmach(4),diver(1))
EQUIVALENCE (rmach(5),LOG10(1))

!     machine constants for the burroughs 1700 system.

!     data rmach(1) / z400800000 /
!     data rmach(2) / z5ffffffff /
!     data rmach(3) / z4e9800000 /
!     data rmach(4) / z4ea800000 /
!     data rmach(5) / z500e730e8 /

!     machine constants for the burroughs 5700/6700/7700 systems.

!     data rmach(1) / o1771000000000000 /
!     data rmach(2) / o0777777777777777 /
!     data rmach(3) / o1311000000000000 /
!     data rmach(4) / o1301000000000000 /
!     data rmach(5) / o1157163034761675 /

!     machine constants for the cdc 6000/7000 series.

!     data rmach(1) / o"00564000000000000000" /
!     data rmach(2) / o"37767777777777777776" /
!     data rmach(3) / o"16414000000000000000" /
!     data rmach(4) / o"16424000000000000000" /
!     data rmach(5) / o"17164642023241175720" /


!     machine constants for sun

DATA rmach(1) /2.3d-308/
DATA rmach(2) /1.79d+308/
DATA rmach(3) /1.0d-14/
DATA rmach(4) /1.0d-14/

!     machine constants for the cray 1

!     data rmach(1) / 200034000000000000000b /
!     data rmach(2) / 577767777777777777776b /
!     data rmach(3) / 377224000000000000000b /
!     data rmach(4) / 377234000000000000000b /
!     data rmach(5) / 377774642023241175720b /

!     machine constants for the data general eclipse s/200

!     note - it may be appropriate to include the following card -
!     static rmach(5)

!     data small/20k,0/,large/77777k,177777k/
!     data right/35420k,0/,diver/36020k,0/
!     data log10/40423k,42023k/

!     machine constants for the harris 220

!     data small(1),small(2) / '20000000, '00000201 /
!     data large(1),large(2) / '37777777, '00000177 /
!     data right(1),right(2) / '20000000, '00000352 /
!     data diver(1),diver(2) / '20000000, '00000353 /
!     data log10(1),log10(2) / '23210115, '00000377 /

!     machine constants for the honeywell 600/6000 series.

!     data rmach(1) / o402400000000 /
!     data rmach(2) / o376777777777 /
!     data rmach(3) / o714400000000 /
!     data rmach(4) / o716400000000 /
!     data rmach(5) / o776464202324 /

!     machine constants for the hp 2100

!     3 word double precision with ftn4

!     data small(1), small(2) / 40000b,       1 /
!     data large(1), large(2) / 77777b, 177776b /
!     data right(1), right(2) / 40000b,    325b /
!     data diver(1), diver(2) / 40000b,    327b /
!     data log10(1), log10(2) / 46420b,  46777b /

!     machine constants for the hp 2100
!     4 word double precision with ftn4

!     data small(1), small(2) / 40000b,       1 /
!     data large91), large(2) / 77777b, 177776b /
!     data right(1), right(2) / 40000b,    325b /
!     data diver(1), diver(2) / 40000b,    327b /
!     data log10(1), log10(2) / 46420b,  46777b /

!     machine constants for the ibm 360/370 series,
!     the xerox sigma 5/7/9, the sel systems 85/86  and
!     the perkin elmer (interdata) 7/32.

!     data rmach(1) / z00100000 /
!     data rmach(2) / z7fffffff /
!     data rmach(3) / z3b100000 /
!     data rmach(4) / z3c100000 /
!     data rmach(5) / z41134413 /

!     machine constants for the pdp-10 (ka or ki processor).

!     data rmach(1) / "000400000000 /
!     data rmach(2) / "377777777777 /
!     data rmach(3) / "146400000000 /
!     data rmach(4) / "147400000000 /
!     data rmach(5) / "177464202324 /

!     machine constants for pdp-11 fortran supporting
!     32-bit integers (expressed in integer and octal).

!     data small(1) /    8388608 /
!     data large(1) / 2147483647 /
!     data right(1) /  880803840 /
!     data diver(1) /  889192448 /
!     data log10(1) / 1067065499 /

!     data rmach(1) / o00040000000 /
!     data rmach(2) / o17777777777 /
!     data rmach(3) / o06440000000 /
!     data rmach(4) / o06500000000 /
!     data rmach(5) / o07746420233 /

!     machine constants for pdp-11 fortran supporting
!     16-bit integers  (expressed in integer and octal).

!     data small(1),small(2) /   128,     0 /
!     data large(1),large(2) / 32767,    -1 /
!     data right(1),right(2) / 13440,     0 /
!     data diver(1),diver(2) / 13568,     0 /
!     data log10(1),log10(2) / 16282,  8347 /

!     data small(1),small(2) / o000200, o000000 /
!     data large(1),large(2) / o077777, o177777 /
!     data right(1),right(2) / o032200, o000000 /
!     data diver(1),diver(2) / o032400, o000000 /
!     data log10(1),log10(2) / o037632, o020233 /

!     machine constants for the univac 1100 series.

!     data rmach(1) / o000400000000 /
!     data rmach(2) / o377777777777 /
!     data rmach(3) / o146400000000 /
!     data rmach(4) / o147400000000 /
!     data rmach(5) / o177464202324 /

!     machine constants for the vax 11/780
!    (expressed in integer and hexadecimal)
!  ***the hex format below may not be suitable for unix systems***
!  *** the integer format should be ok for unix systems***

!     data small(1) /       128 /
!     data large(1) /    -32769 /
!     data right(1) /     13440 /
!     data diver(1) /     13568 /
!     data log10(1) / 547045274 /

!     data small(1) / z00000080 /
!     data large(1) / zffff7fff /
!     data right(1) / z00003480 /
!     data diver(1) / z00003500 /
!     data log10(1) / z209b3f9a /

!     machine constants for the z80 microprocessor

!     data small(1),small(2) /     0,    256/
!     data large(1),large(2) /    -1,   -129/
!     data right(1),right(2) /     0,  26880/
!     data diver(1),diver(2) /     0,  27136/
!     data log10(1),log10(2) /  8347,  32538/


!***first executable statement  r1mach
IF (i < 1  .OR.  i > 5) CALL lnkerr ( 'r1mach -- i out of bounds')

r1mach = rmach(i)
RETURN

END FUNCTION r1mach






FUNCTION r9lgmc(x)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:15:53

!***begin prologue  r9lgmc
!***date written   770801   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  c7e
!***keywords  correction factor,log gamma,special function
!***author  fullerton, w., (lanl)
!***purpose  computes the log gamma correction factor so that
!            log(gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x
!            + r9lgmc(x)
!***description

! compute the log gamma correction factor for x .ge. 10.0 so that
!  log (gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + r9lgmc(x)

! series for algm       on the interval  0.          to  1.00000d-02
!                                        with weighted error   3.40d-16
!                                         log weighted error  15.47
!                               significant figures required  14.39
!                                    decimal places required  15.86
!***references  (none)
!***routines called  csevl,inits,r1mach,lnkerr
!***end prologue  r9lgmc


IMPLICIT REAL*8 (a-h,o-z)
REAL*8, INTENT(IN)                         :: x

DIMENSION algmcs(6)
DATA algmcs( 1) /    .166638948045186D0 /
DATA algmcs( 2) /   -.0000138494817606D0 /
DATA algmcs( 3) /    .0000000098108256D0 /
DATA algmcs( 4) /   -.0000000000180912D0 /
DATA algmcs( 5) /    .0000000000000622D0 /
DATA algmcs( 6) /   -.0000000000000003D0 /
DATA nalgm, xbig, xmax / 0, 2*0.0D0 /
!***first executable statement  r9lgmc
IF (nalgm /= 0) GO TO 10
nalgm = inits (algmcs, 6, r1mach(3))
xbig = 1.0D0/SQRT(r1mach(3))
xmax = EXP (MIN(LOG(r1mach(2)/12.0D0), -LOG(12.0D0*r1mach(1))) )

10   IF (x < 10.0D0) CALL lnkerr ( 'r9lgmc  x must be ge 10')
IF (x >= xmax) GO TO 20

r9lgmc = 1.0D0/(12.0D0*x)
IF (x < xbig) r9lgmc = csevl (2.0D0*(10.d0/x)**2-1.d0, algmcs, nalgm)/x
RETURN

20   r9lgmc = 0.0D0
CALL lnkerr ( 'r9lgmc  x so big r9lgmc underflows')
RETURN

END FUNCTION r9lgmc



SUBROUTINE gamlim(xmin,xmax)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:15:30

!***begin prologue  gamlim
!***date written   770401   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  c7a,r2
!***keywords  gamma function,limits,special function
!***author  fullerton, w., (lanl)
!***purpose  computes the minimum and maximum bounds for x in gamma(x).
!***description

! calculate the minimum and maximum legal bounds for x in gamma(x).
! xmin and xmax are not the only bounds, but they are the only non-
! trivial ones to calculate.

!             output arguments --
! xmin   minimum legal value of x in gamma(x).  any smaller value of
!        x might result in underflow.
! xmax   maximum legal value of x in gamma(x).  any larger value will
!        cause overflow.
!***references  (none)
!***routines called  r1mach,lnkerr
!***end prologue  gamlim
!***first executable statement  gamlim


IMPLICIT REAL*8(a-h,o-z)
REAL*8, INTENT(OUT)                        :: xmin
REAL*8, INTENT(OUT)                        :: xmax

alnsml = LOG(r1mach(1))
xmin = -alnsml
DO  i=1,10
  xold = xmin
  xln = LOG(xmin)
  xmin = xmin - xmin*((xmin+0.5D0)*xln - xmin - 0.2258D0 + alnsml)  &
      / (xmin*xln + 0.5D0)
  IF (ABS(xmin-xold) < 0.005D0) GO TO 20
END DO
CALL lnkerr ( 'gamlim  unable to find xmin')

20   xmin = -xmin + 0.01D0

alnbig = LOG(r1mach(2))
xmax = alnbig
DO  i=1,10
  xold = xmax
  xln = LOG(xmax)
  xmax = xmax - xmax*((xmax-0.5D0)*xln - xmax + 0.9189D0 - alnbig)  &
      / (xmax*xln - 0.5D0)
  IF (ABS(xmax-xold) < 0.005D0) GO TO 40
END DO
CALL lnkerr ( 'gamlim  unable to find xmax')

40   xmax = xmax - 0.01D0
xmin = MAX (xmin, -xmax+1.d0)

RETURN
END SUBROUTINE gamlim




FUNCTION inits(os,nos,eta)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-12-18  Time: 10:15:35

!***begin prologue  inits
!***date written   770401   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  c3a2
!***keywords  initialize,orthogonal series,special function
!***author  fullerton, w., (lanl)
!***purpose  initializes an orthogonal series so that it defines the
!            number of terms to carry in the series to meet a specified
!            error.
!***description

! initialize the orthogonal series so that inits is the number of terms
! needed to insure the error is no larger than eta.  ordinarily, eta
! will be chosen to be one-tenth machine precision.

!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.
!***references  (none)
!***routines called  lnkerr
!***end prologue  inits


IMPLICIT REAL*8 (a-h,o-z)
REAL*8, INTENT(IN OUT)                     :: os(nos)
INTEGER, INTENT(IN)                      :: nos
REAL*8, INTENT(IN)                         :: eta


!***first executable statement  inits
IF (nos < 1) CALL lnkerr ( 'inits   number of coefficients lt 1')

ERR = 0.d0
DO  ii=1,nos
  i = nos + 1 - ii
  ERR = ERR + ABS(os(i))
  IF (ERR > eta) GO TO 20
END DO

20   IF (i == nos) CALL lnkerr ( 'inits   eta may be too small')
inits = i

RETURN
END FUNCTION inits

