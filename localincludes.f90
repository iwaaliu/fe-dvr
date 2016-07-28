module proppara
  use accuracy_real
      integer,parameter::in=20,m=200,imax=101
      integer,parameter::nm=in!-1
      real*8 xGaussPoint(in),wGaussPoint(in)
!      real*8 T(nm,nm),H(imax,nm,nm)
      real*8 dx,yy,x0
      integer ix,iy,i,j,k,kmax,nGaussPoint
      integer ij,ij2,i2,j2,ik,l,im,ijk
      character*256 LobattoFile,LegendreFile

      intrinsic CMPLX, CONJG, MIN, DBLE, IMAG
      external zgcoov, zgcrsv, zgccsv!,myzgcoov
      double precision tic, tac, clock
!*---  matrix data ...
!*---  BEWARE: these values must match those in zgmatv.f
      integer  nz, nmax, nzmax
      parameter( nmax = 9600, nzmax = 10000000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

!*---  arguments variables ...
      integer n, mmax, liwsp,lwork,iv,it
      INTEGER LWSP
!      parameter( iv=nm*imax-imax+1)
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = mmax+2 )
      integer iwsp(liwsp)
      double precision  tol, anorm, s1, s2,vnorm
      double precision dt,et,ur,ui,xr,E0,Eir,Exuv,times,tinit,tfinal
      real(KIND=idp) gr_rate,rr

      complex*16,allocatable:: vec(:),vec0(:), w(:),u(:,:),psi(:,:,:),psi0(:,:,:),coef0(:),vecl(:),a0(:)
      complex*16,allocatable::ut(:,:)
      real(KIND=idp),allocatable::v(:,:),v1(:,:),v2(:,:),pot(:),si(:),v12(:,:,:)
      real(KIND=idp),allocatable::F(:,:),f0(:,:),Eig(:),Ei(:),work(:),hc(:,:)
      real(KIND=idp),allocatable::f1(:,:),f2(:,:),vf12(:,:),f11(:,:),f22(:,:),df1(:),df2(:)

      integer,allocatable:: ipiv(:),kii(:)

      integer l1,l2,lC,l1p,l2p,lcp,index
      integer,allocatable,dimension(:,:)::llist
      real*8 ef,rima,rima1,rima2
      real*8 coef,cof3j,cof9j
      real*8,allocatable,dimension(:)::r1,r2
      real*8,allocatable,dimension(:,:,:)::f0l,DiagC,DiagF
      real*8,allocatable,dimension(:,:,:,:)::V0,VF
      complex*16,allocatable,dimension(:,:,:,:)::tranC,tranF

      complex*16 spec


      complex*16 wsp(lwsp)!,WSP0(LWSP)
      integer mprint,ntsteps
      integer:: nnz, itrace,iflag,info,nt=1000
      integer:: nl
      complex*16 ZERO, ONE, IU
      parameter( ZERO=(1.d-80,0.0d0), ONE=(1.0d0,0.0d0) )
!      parameter( IU=(1.d0,0.d0))!*(.0d0,1.0d0))
      
!      real(KIND=idp):: pi=3.141592632d0
      real*8 en,en0,en1
      intrinsic ABS

end module proppara


module mymataid

  contains
     subroutine myzgcoov ( x, y )
      implicit none
      complex*16 x(*), y(*)
!*
!*---  Computes y = A*x. A is passed via a fortran `common statement'.
!*---  A is assumed here to be under the COOrdinates storage format.
!*
      integer n, nz, nzmax
      parameter( nzmax = 10000000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )
      
      do j = 1,n
         y(j) = ZERO
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
    END subroutine myzgcoov
    
    
!*----------------------------------------------------------------------|
    subroutine myZGEXPV(  m, t, v, w, tol, anorm,&
         wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )
      !         call   myZGEXPV(im, dt,v,w, tol, 1*anorm,&
      !                    1.*wsp,1*lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      implicit none
      integer          n, m,  liwsp, itrace, iflag, iwsp(liwsp)
      INTEGER LWSP
      double precision t, tol, anorm
      complex*16,intent(in),dimension(:):: v
      complex*16,dimension(:)::wsp
      complex*16, intent(out),dimension(:):: w
      external         matvec

!*-----Purpose----------------------------------------------------------|
!*
!*---  ZGEXPV computes w = exp(t*A)*v
!*     for a Zomplex (i.e., complex double precision) matrix A 
!*
!*     It does not compute the matrix exponential in isolation but
!*     instead, it computes directly the action of the exponential
!*     operator on the operand vector. This way of doing so allows 
!*     for addressing large sparse problems. 
!*
!*     The method used is based on Krylov subspace projection
!*     techniques and the matrix under consideration interacts only
!*     via the external routine `matvec' performing the matrix-vector 
!*     product (matrix-free method).

      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,&
           mxreject = 0,&
           ideg     = 6,&
           delta    = 1.2d0,&
           gamma    = 0.9d0 )
      
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,&
           ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,&
           nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,&
           s_error, x_error, t_now, t_new, t_step, t_old,&
           xm, beta, break_tol, p1, p2, p3, eps, rndoff,&
           vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij
      
      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
    !  !$omp parallel default(private) 
      n=ubound(v,1)
      if(lwsp.ne.ubound(wsp,1)) print*, 'lwsp wrong'
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) then
         print*,'iflag,m,n(if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1'
         print*,'iflag,m,n(if ( liwsp.lt.m+2 ) iflag = -2'
         print*,'iflag,m,n( if ( m.ge.n .or. m.le.0 ) iflag = -3)'
         print*,iflag,m,n
         stop 'bad sizes (in input of myZGEXPV)'
      endif

      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
      break_tol = 1.0d-10


      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta 

      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

 100  if ( t_now.ge.t_out ) goto 500
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
!*---  Arnoldi loop ...
      j1v = iv + n
      do 200 j = 1,m

         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DZNRM2( n, wsp(j1v),1 )
!*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
!            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DZNRM2( n, wsp(j1v),1 )
!*
!*---  set 1 for the 2-corrected scheme ...
!*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
      ireject = 0
 401  continue
!*---  compute w = beta*V*exp(t_step*H)*e1 ...
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,&
                     wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else

         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue

      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif

      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.&
          (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
           ! return
         endif
         goto 401
      endif

      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )

      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )

      t_now = t_now + t_step

      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
      !!$OMP END PARALLEL         

    END subroutine myZGEXPV
!*----------------------------------------------------------------------|

end module mymataid
