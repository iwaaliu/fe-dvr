!-------------------------------------------------------------------------------!
! This code demonstrates the implementation of FE-DVR for one-electron atoms.   !
!-------------------------------------------------------------------------------!
! Language : Fortran 90                                                         !
! Parallized : NO                                                               !
! Initial version : November  27, 2006                                          !
! Last updated    : September 28, 2011                                          !
!                                                                               !
! Updated: includes expokit;1 Ch imag. TDSE tested.  Dec. 19, 2011 by iwaa      !
! Updated: includes anglib module,C-G,3j,6j,9j coef. Dec. 21, 2011 by iwaa      !
! Updated: includes laser field                      Dec. 22, 2011 by iwaa      !
! Updated: 2d without correlation(1/|r12|)           Dec. 25, 2011 by iwaa      !
! Updated: change to 2d code                         Dec. 26, 2011 by iwaa      !
! Updated: included laser field                      Jan. 10, 2012 by iwaa      !
! Last updated    :                                  Jan. 10, 2012 by iwaa      !
!-------------------------------------------------------------------------------!
      include "../../include/anglib.f90"


      program hydrogen_fedvr 
        use omp_lib
        use accuracy_real
        use globalmod_constants
        use globalmod_discrete
        use mymataid
        use proppara
        use anglib
        use laser
      implicit none
  
      integer                                   :: ierror
      real(KIND=idp) :: centrf,mycleb,my3j,my6j,my9j,rl

      write(idwrite,'(1x,a)') '===================================================='
      write(idwrite,'(1x,a)') '  Implimentation of FE-DVR in Hydrogen-like Atoms   '
      write(idwrite,'(1x,a)') '===================================================='
      write(idwrite,*)
      write(idwrite,*)
      print*, '# Define ZERO as (1.d-80,0.d00)!'

! Input parameters and then set up mesh points:
      call readin(ierror)
      if(ierror/=0) then
         write(idwrite,'(1x,a)') '# OOOoooops, something wrong !!!'
         write(idwrite,'(1x,a,i4)') '# ierror = ',ierror
      else
         write(idwrite,'(1x,a)') '# calling readin is down: :)'
      end if
      
      ! Diagonalize the Hamiltonian matrix:
      call hamiltonian_matrix
      print*,'prop...'
      call prop

!------------------------------------------------------------------------------!
      iu=(1.d0,0.d0)!*(0.d0,1.d0)  !sssssiiiiiiiiiiiiiiii 
!       lmax=20
      nl=lmax+1
      dlmax=dble(lmax)
      iv=(no_gp-1) * no_fe-1
     !  nsize=iv
      allocate(vec(iv),vec0(iv),w(iv),si(iv),vecl(0:lmax))
      allocate(V(0:lmax,iv),u(iv,iv),ut(iv,iv),pot(iv))
      allocate(V1(0:lmax,iv),V2(0:lmax,iv))
      allocate(F(0:lmax,0:lmax),f0(0:lmax,0:lmax),Eig(nL),Ei(0:lmax),hc(nl,nl))
      allocate(psi(0:lmax,iv,iv),psi0(0:lmax,iv,iv))

!      deallocate(Eig,Hc,ipiv)
!      allocate(r1(n),r2(n))
      allocate(coef0(0:2*lmax+3))
      allocate(f0l(0:2*lmax+3,nl,nl),DiagC(0:lmax,iv,iv))
      allocate(f1(nl,nl),f2(nl,nl),vf12(0:lmax,0:lmax))
      allocate(f11(nl,nl),f22(nl,nl),df1(0:lmax),df2(0:lmax))

      allocate(V0(0:lmax,0:lmax,iv,iv),tranC(0:lmax,0:lmax,iv,iv))
      allocate(VF(0:lmax,0:lmax,iv,iv),tranF(0:lmax,0:lmax,iv,iv))
      allocate(DiagF(0:lmax,iv,iv))

      allocate(llist(17,3))
      open(998,file='~/lib/l1l2Llist0-5.txt')
      read(998,*)
      do j=1,17
         read(998,*) index,llist(index,1:3)
         print*,index,llist(index,1:3)
      enddo
      print*, 'config reading done!'
      print*,2*lmax+3,llist(nl,1),llist(nl,2),min(2*lmax+3,llist(nl,1)+llist(nl,2))
      
      !$omp parallel default(none)  &
      !$OMP shared(llist,nl,lmax,iv,V0,f1,f2,f0l,r,coef0) private(i,j,L,ix,iy,L1,L2,Lc,L1p,L2p,Lcp)
      !$OMP do      
      do i=1,nl
         do j=1,nl
            L1 =llist(i,1)
            L2 =llist(i,2)
            LC =llist(i,3)
            l1p=llist(j,1)
            l2p=llist(j,2)
            LCp=llist(j,3)
            
            f1(i,j)=dsqrt(9.d0*(2*L1p+1)*(2*L2p+1)*(2*LCp+1))/pi*.25d0&
                 *mycleb(1,0,L1p,0,L1,0)*mycleb(0,0,L2p,0,L2,0)&
                 *mycleb(1,0,LCp,0,LC,0)*my9j(1,L1p,L1,0,L2p,L2,1,LCp,LC)

            f2(i,j)=dsqrt(9.d0*(2*L1p+1)*(2*L2p+1)*(2*LCp+1))/pi*.25d0&
                 *mycleb(0,0,L1p,0,L1,0)*mycleb(1,0,L2p,0,L2,0)&
                 *mycleb(1,0,LCp,0,LC,0)*my9j(0,L1p,L1,1,L2p,L2,1,LCp,LC)

            do L=0,min(2*lmax+3,llist(nl,1)+llist(nl,2)+3)  
               coef0(l)=dsqrt(((2*L+1.d0)**2)*(2*L1p+1.d0)*(2*L2p+1.d0)*(2*LCp+1.d0))
              ! print*, '30'

               f0l(L,i,j)=coef0(L)*mycleb(l,0,l1p,0,l1,0)&
                    *mycleb(L,0,L2p,0,L2,0)*mycleb(0,0,LCp,0,LC,0)&
                    *my9j(L,L1p,L1,L,L2p,L2,0,Lcp,Lc)             
               do iy=1,iv
                  do ix=1,iv
                     v0(i-1,j-1,ix,iy)=v0(i-1,j-1,ix,iy)&
                          +(-1)**L/dsqrt(2.d0*L+1.d0)*rl(r(ix+1),r(iy+1),l)*f0l(l,i,j)
                  end do
               end do
            enddo
         enddo
      enddo
      !$OMP end do
      !$OMP END PARALLEL
!      stop
      DO I=0,LMAX
         PRINT 2004,V0(0:LMAX,I,12,2)
      ENDDO
      print*,'i,j,l',i,j,l
   !$omp parallel default(none)  &
   !$OMP shared(nl,lmax,iv,V0) private(i)
   !$OMP do  
      DO I=0,LMAX
         PRINT 2004,i*1.,V0(0:LMAX-1,I,12,2)
      ENDDO
   !$OMP end do
   !$OMP END PARALLEL
      v0=0.d0*v0
! stop
      do i=1,nl
         do j=1,nl
            l1=llist(i,1)
            l2=llist(i,2)
            LC=llist(i,3)
            l1p=llist(j,1)
            l2p=llist(j,2)
            LCp=llist(j,3)
            
            do L=0,max(2*lmax+3,llist(nl,1)+llist(nl,2)+3)
               coef0(l)=(-1)**(LC+L2+L2p)*dsqrt((2*L1+1.d0)*(2*L1p+1.d0)*(2*L2+1.d0)*(2*L2p+1.d0))
               f0l(l,i,j)=coef0(l)*my3j(L1,L,L1P,0,0,0)&
                    *my3J(L2,L,L2P,0,0,0)*my6J(LC,L2P,L1P,L,L1,L2)

               do iy=1,iv
                  do ix=1,iv
                     v0(i-1,j-1,ix,iy)=v0(i-1,j-1,ix,iy)&
                          +rl(r(ix+1),r(iy+1),l)*f0l(l,i,j)
                  end do
               end do
            enddo
         enddo
      enddo
!~~~~~~~~~~~~~~~~~~~~~~~~
   !$omp parallel default(none)  &
   !$OMP shared(llist,nl,lmax,iv,V0,f1,f2,f0l,r) private(i,j,L,ix,iy,l1,l2,lc,L1p,L2p,lcp,coef0)
   !$OMP do  
      DO I=0,LMAX
         PRINT 2004,i*1.,V0(0:LMAX,I,12,2)
      ENDDO
   !$OMP end do
   !$OMP END PARALLEL

            print*, 'F1'
            print 2010,llist(1:nl,1)
            print 2010,llist(1:nl,2)
            print 2010,llist(1:nl,3)

            DO I=1,LMAX+1
               PRINT 2004,F1(1:nL,I)
            ENDDO

           print*, 'F2'
            print 2010,llist(1:nl,1)
            print 2010,llist(1:nl,2)
            print 2010,llist(1:nl,3)

            DO I=1,LMAX+1
               PRINT 2004,F2(1:nL,I)
            ENDDO
!STOP
      
      lwork = 3*nl-1
      ALLOCATE( work(lwork) )
      do ix=1,iv
         do iy=1,iv
            hc(1:nl,1:nl)=V0(0:lmax,0:lmax,ix,iy)
            LWORK = 3*nl-1
            call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)
            LWORK = MIN( 3*nl-1, INT( WORK( 1 ) ) )
           ! hc(1:nl,1:nl)=V0(0:lmax,0:lmax,ix,iy)
            !CALL DSYEV( 'Vectors', 'Upper', N, A, N, WW, WORK, LWORK, INFO )
           ! call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)
            diagC(0:lmax,ix,iy)=Eig(1:nl)
            TranC(0:lmax,0:lmax,ix,iy)=hc(1:nl,1:nl)
         enddo
      enddo


!      lwork = 3*nl-1
!      ALLOCATE( work(lwork) )

            hc=F1
            LWORK = 3*nl-1
            work=0*work
            call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)
!            LWORK = MIN( 3*nl-1, INT( WORK( 1 ) ) )
            df1(0:lmax)=4.*Pi*Eig(1:nl)/dsqrt(3.d0)
            F11(0:lmax,0:lmax)=hc(1:nl,1:nl)
            print*, 'tran F1'
            print 2010,llist(1:nl,1)
            print 2010,llist(1:nl,2)
            print 2010,llist(1:nl,3)

            print 2004,df1(0:lmax)
            DO I=0,LMAX
               PRINT 2004,F11(0:LMAX,I)
            ENDDO

            hc=F2
            LWORK = 3*nl-1
            work=0*work
            call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)
!            LWORK = MIN( 3*nl-1, INT( WORK( 1 ) ) )
            df2(0:lmax)=4.*Pi*Eig(1:nl)/dsqrt(3.d0)
            df2(0:lmax)=4.*Pi*Eig(1:nl)/dsqrt(3.d0)

            F22(0:lmax,0:lmax)=hc(1:nl,1:nl)

      print*, 'tran F2'
            print 2004, df2(0:lmax)
      DO I=0,LMAX
         PRINT 2004,F22(0:LMAX,I)
      ENDDO


      do iy=1,iv
         do ix=1,iv
            VF(0:lmax,0:lmax,ix,iy)=(r(ix+1)*F1(1:nl,1:nl)+r(iy+1)*F2(1:nl,1:nl))
         enddo
      enddo
      lwork = 3*nl-1
      do ix=1,iv
         do iy=1,iv
            hc(1:nl,1:nl)=VF(0:lmax,0:lmax,ix,iy)
            LWORK = 3*nl-1
            call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)
            LWORK = MIN( 3*nl-1, INT( WORK( 1 ) ) )
            diagF(0:lmax,ix,iy)=4.*pi/sqrt(3.d0)*Eig(1:nl)
            TranF(0:lmax,0:lmax,ix,iy)=hc(1:nl,1:nl)
         enddo
      enddo



!      deallocate(F1,F2)
!      allocate(f1(0:lmax,0:lmax),f2(0:lmax,0:lmax))
      
!      stop

            !CALL DSYEV( 'Vectors', 'Upper', N, A, N, WW, WORK, LWORK, INFO )
          !  call  DSYEV('V','U',nl,hc,nl,Eig,work,lwork,info)

          !  print*,eig
          !  stop
       centrf = 0.5_idp * dlang * (dlang + 1.0_idp)
!       cnuclear=1.d0
       print*,'size(r),iv',size(r),iv,dlang,centrf
       
       do i=1,iv
          rr=r(i+1)*r(i+1)
          do l=0,lmax
             index=L+1
             l1=llist(index,1)
             l2=llist(index,2)
             v1(L,i)= .5_idp*L1*(L1+1.0_idp) /rr
             v2(L,i)= .5_idp*L2*(L2+1.0_idp) /rr
          enddo
       enddo

       ij=0
       si=0.d0
       do i=1,no_fe
          do j=1,no_gp
             if(i.gt.1.and.j==1)ij=ij-1
             if(i==1.and.j==1) goto 2090
             if(i==no_fe.and.j==no_gp) goto 2091
             ij=ij+1
             si(ij)=si(ij)+( mat_reg(i)%wt(j) )
!             write(*,'(3i3,100f10.5)'),i,j,ij,si(ij),mat_reg(i)%wt(j)
2090         ijk=0  
          end do
       enddo
2091   si=sqrt(si)
       print*,ij,iv
       print 2001,si
  !     stop
      print*,'nz(or max(ijk))=',ijk
      allocate(kii(iv))
      ij=0
      do i=1, (no_gp-1) * no_fe-1!nsize
         do j=1, (no_gp-1) * no_fe-1!nsize
            if(abs(tot_ke_mat(i,j)).gt.1.e-8) then
               ij=ij+1
               a(ij)=one*tot_ke_mat(i,j)
               ia(ij)=i
               ja(ij)=j
               if(i==j) kii(i)=ij                          
            endif
         enddo
      enddo
      print*,'ij=',ij,(no_gp)**2*no_fe-(no_fe-1)-4*no_gp+2

      nz=ij
      allocate(a0(nz))
      print*,'nz(or max(ijk))=',nz
      a=a*(-iu)
      a0(1:nz)=a(1:nz)
!*---  compute the infinite norm of A ...
      n=no_fe*(no_gp-1)-1 !=nsize=iv

      do i = 1,ia(nz)
         wsp(i) = ZERO
      enddo

      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo

      anorm = wsp(1)
      do i = 2,ia(nz)
         if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
      enddo

!*---  set other input arguments ...
      tol = 1.0d-50
      im = min(n-1,50)

      itrace = 0
      u=0.d0
      do iy=1,iv
         do ix=1,iv
            u(ix,iy)=dexp(-(r(ix+1)*r(ix+1)+r(iy)*r(iy+1)+4.d0))!*(si(ix)*si(iy))
         enddo
      enddo
      psi(0,:,:)=u(:,:)
      goto 3012
      print*, 'Reading init data file'
      open(301,file='init')
      do ix=1,min(iv,89)
         do iy=1,min(iv,89)
            read(301,*) ur,ur,ur,ui
            psi(0,ix,iy)=ur*si(ix)*si(iy)+Zero
            psi(1,ix,iy)=ui*si(ix)*si(iy)+Zero
         enddo
         read(301,*) 
      enddo
      close(301)
      !u=psi(0,:,:)
3012  do iy=1,iv
         do ix=1,iv
            write(1000,*) r(ix+1),r(iy+1),abs(u(ix,iy))/(si(ix)*si(iy))
         enddo
         write(1000,*)
      enddo
      
      
3013  print*, 'reading done!'

!*---  compute w = exp(t*A)v with ZGEXPV ...
      vnorm=0.d0
      do l=0,1
         u=psi(l,:,:)
      do ix=1,iv
         vnorm=vnorm+sum(abs(conjg(u(ix,:))*u(ix,:)))
      enddo
      enddo
      psi=psi/sqrt(vnorm)
      psi0=psi
     ! psi=zero
     ! psi(0,:,:)=u(:,:)+ZERO
   !   psi(1,:,:)=u(:,:)+ZERO     
      print*,vnorm
      print*, 'initialized done!'


      vnorm=0.d0
      do l=0,lmax
         u=psi(l,:,:)+ZERO        
         do ix=1,iv
            vnorm=vnorm+sum(abs(conjg(u(ix,:))*u(ix,:)))
         enddo
         print*,vnorm
      enddo

      print*,'setup initial/final time:'
      dt=.1d-0
      tinit=0.d0
      tfinal=400.d0
      print*,'propagating...'  
      ntsteps=int((tfinal-tinit)/dt)+1
      print*,'dt,ntsteps:',dt,ntsteps
! propagation
      Eir=0.1d0
      Exuv=0.1d0
      do it=1,ntsteps
         print*,'time=',real(it*dt)
         times=it*dt
         if(times.gt.221) Eir=0.0d0
         if(abs(times-1*10.437).lt.(.5*20.87442)) then
            Exuv=0.10d0         
         else
            Exuv=0.d0
         endif
         et=Eir*dsin(0.057d0*times)*sin(pi*times/221)**2&
              +Exuv*dsin(0.903d0*(times))*sin(pi*(times-2*10.437)/20.87442)**2
         et=0.d0
         write(509,*) times,et
      !   et=*0.01d0

         do L=0,lmax
!Hx for exp(-iHx(t)*dt)
            a(kii(:))=a0(kii(:))+V1(L,:)*(-iu)
            wsp = wsp*0.d0*ZERO   
 !$omp parallel default(none)  &
 !$OMP shared(wsp,nz,a,ia) private(i)
 !$OMP do       
            do i = 1,nz
               wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
            enddo
 !$OMP end do
 !$OMP END PARALLEL 
            anorm = wsp(1)
            do i = 2,ia(nz)
               if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
            enddo

!   exp(-iHx(t)*dt/2)\psi
            u(:,:)=psi(L,:,:)!+ZERO

 !$omp parallel default(none)  &
 !$OMP shared(u,ut,iv,im,nz,a,ia,ja,dt,tol,anorm) private(wsp,iwsp,ix,itrace,iflag)
 !$OMP do
            do ix=1,iv
               itrace=0
               iflag=0
!               print*,'ix=',ix
               if(abs(sum(conjg(u(ix,:))*u(ix,:))).le.abs(zero)) then
                  ut(ix,:)=u(ix,:)
               else
                  call myZGEXPV(im,.5d0*dt,u(ix,:),ut(ix,:), tol,anorm,&
                       wsp,lwsp, iwsp,liwsp, myzgcoov, itrace, iflag )
               endif
              ! ut(ix,iv)=0.d0*zero !set boundary as 0
               u(ix,:)=ut(ix,:)
            enddo
 !$OMP end do
 !$OMP END PARALLEL         

!Hy for exp(-iHy(t)*dt)
            a(kii(:))=a0(kii(:))+V2(L,:)*(-iu)
            wsp = wsp*0.d0*ZERO   
 !$omp parallel default(none)  &
 !$OMP shared(wsp,nz,a,ia) private(i)
 !$OMP do 
            do i = 1,nz
               wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
            enddo
 !$OMP end do
 !$OMP END PARALLEL
            anorm = wsp(1)
            do i = 2,ia(nz)
               if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
            enddo
!   exp(-iHy(t)*dt/2)\psi
 !$omp parallel default(none)  &
 !$OMP shared(u,ut,iv,im,nz,a,ia,ja,dt,tol,anorm) private(wsp,iwsp,iy,itrace,iflag)
 !$OMP do
            do iy=1,iv
              ! vec=u(:,iy)!+ZERO
               itrace=0
               iflag=0
               if(abs(sum(conjg(u(:,iy))*u(:,iy))).le.abs(zero)) then
                  ut(:,iy)=u(:,iy)
               else
                  call myZGEXPV(im,.5d0*dt,u(:,iy),ut(:,iy), tol,anorm,&
                       wsp,lwsp, iwsp,liwsp, myzgcoov, itrace, iflag )
               endif
               ut(iv,iy)=0.d0*zero !set boundary as 0
               u(:,iy)=ut(:,iy)
            enddo
 !$OMP end do
 !$OMP END PARALLEL            
            psi(L,:,:)=u(:,:)
         enddo


!  1/|r1-r2|
        ! goto 3001 !skip  e-e interaction 1/|r1-r2|

         do ix=1,iv
            do iy=1,iv
               f=tranC(:,:,ix,iy)
               f0= transpose(f)
               vecl=psi(:,ix,iy)
               psi(:,ix,iy)=matmul(f0,vecl)
            enddo
         enddo
 !$omp parallel default(none)  &
 !$OMP shared(iv,psi,iu,diagC,dt) private(ix)
 !$OMP do
         do ix=1,iv
            psi(:,ix,:)=psi(:,ix,:)*cdexp(-iu*diagC(:,ix,:)*dt)
         enddo
 !$OMP end do
 !$OMP END PARALLEL        
         do iy=1,iv
            do ix=1,iv
               f=tranC(:,:,ix,iy)
               vecl=psi(:,ix,iy)
               psi(:,ix,iy)=matmul(f,vecl)
            enddo
         enddo

         
! e-field 
         if(abs(et).lt..1e-8) goto 3001 !ignore very weak field
         print*, 'field working'
         do ix=1,iv
            do iy=1,iv
               f=tranF(:,:,ix,iy)
               f0= transpose(f) 
               vecl=psi(:,ix,iy)
               psi(:,ix,iy)=matmul(f0,vecl)
            enddo
         enddo
 !$omp parallel default(none)  &
 !$OMP shared(iv,psi,iu,diagF,dt,et) private(ix)
 !$OMP do
         do ix=1,iv
            psi(:,ix,:)=psi(:,ix,:)*cdexp(-iu*et*diagF(:,ix,:)*dt)
 !$OMP end do
 !$OMP END PARALLEL          
         do iy=1,iv
            do ix=1,iv
               f=tranF(:,:,ix,iy)
               vecl=psi(:,ix,iy)
               psi(:,ix,iy)=matmul(f,vecl)
            enddo
         enddo
         

!  beginning of exp(-iH(t)*dt/2)\psi
3001     do L=0,lmax
!Hx for exp(-iHx(t)*dt)
            a(kii(:))=a0(kii(:))+V1(L,:)*(-iu)
            wsp = wsp*0.d0*ZERO   
 !$omp parallel default(none)  &
 !$OMP shared(wsp,nz,a,ia) private(i)
 !$OMP do       
            do i = 1,nz
               wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
            enddo
 !$OMP end do
 !$OMP END PARALLEL 
            anorm = wsp(1)
            do i = 2,ia(nz)
               if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
            enddo

!   exp(-iHx(t)*dt/2)\psi
            u(:,:)=psi(L,:,:)!+ZERO

 !$omp parallel default(none)  &
 !$OMP shared(u,ut,iv,im,nz,a,ia,ja,dt,tol,anorm) private(wsp,iwsp,ix,itrace,iflag)
 !$OMP do
            do ix=1,iv
               itrace=0
               iflag=0
!               print*,'ix=',ix
               if(abs(sum(conjg(u(ix,:))*u(ix,:))).le.abs(zero)) then
                  ut(ix,:)=u(ix,:)
               else
                  call myZGEXPV(im,.5d0*dt,u(ix,:),ut(ix,:), tol,anorm,&
                       wsp,lwsp, iwsp,liwsp, myzgcoov, itrace, iflag )
               endif
              ! ut(ix,iv)=0.d0*zero !set boundary as 0
               u(ix,:)=ut(ix,:)
            enddo
 !$OMP end do
 !$OMP END PARALLEL         

!Hy for exp(-iHy(t)*dt)
            a(kii(:))=a0(kii(:))+V2(L,:)*(-iu)
            wsp = wsp*0.d0*ZERO   
 !$omp parallel default(none)  &
 !$OMP shared(wsp,nz,a,ia) private(i)
 !$OMP do 
            do i = 1,nz
               wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
            enddo
 !$OMP end do
 !$OMP END PARALLEL
            anorm = wsp(1)
            do i = 2,ia(nz)
               if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
            enddo
!   exp(-iHy(t)*dt/2)\psi
 !$omp parallel default(none)  &
 !$OMP shared(u,ut,iv,im,nz,a,ia,ja,dt,tol,anorm) private(wsp,iwsp,iy,itrace,iflag)
 !$OMP do
            do iy=1,iv
              ! vec=u(:,iy)!+ZERO
               itrace=0
               iflag=0
               if(abs(sum(conjg(u(:,iy))*u(:,iy))).le.abs(zero)) then
                  ut(:,iy)=u(:,iy)
               else
                  call myZGEXPV(im,.5d0*dt,u(:,iy),ut(:,iy), tol,anorm,&
                       wsp,lwsp, iwsp,liwsp, myzgcoov, itrace, iflag )
               endif
               ut(iv,iy)=0.d0*zero !set boundary as 0
               u(:,iy)=ut(:,iy)
            enddo
 !$OMP end do
 !$OMP END PARALLEL            
            psi(L,:,:)=u(:,:)
         enddo
         
         
         if(abs(iu-1.d0).gt.1.d-5) then ! iu=(0.d0,1.d0) Imaginary time
            psi=(psi)
         else ! iu=(1.d0,0.d0) real time TDSE
! normalized
            vnorm=0.d0
            do l=0,lmax
               u=psi(l,:,:)
               do ix=1,iv
                  vnorm=vnorm+sum(abs(conjg(u(ix,:))*u(ix,:)))
               enddo
            enddo
            psi=(psi)/sqrt(vnorm)
         endif

!Energy
         en=0.d0
         en0=0.d0
         do L=0,lmax
            u=psi(l,:,:)  
!Hx for exp(-iHx(t)*dt)
            a(kii(:))=a0(kii(:))+V1(L,:)*(-iu)
            do ix=1,iv
               vec=u(ix,:)
               do i=1,nz
                  en=en+dble(conjg(vec(ia(i)))*vec(ja(i))*(iu*a(i)))*(-iu*iu)
               enddo
            enddo
!Hy for exp(-iHy(t)*dt)
            a(kii(:))=a0(kii(:))+V2(L,:)*(-iu)
            do i = 2,ia(nz)
               if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
            enddo
            do iy=1,iv
               vec=u(:,iy)
               do i=1,nz
                  en0=en0+dble(conjg(vec(ia(i)))*vec(ja(i))*(iu*a(i)))*(-iu*iu)
               enddo
            enddo
         enddo
! E(1/r12)
         en1=0.d0
         do ix=1,iv
            do iy=1,iv
               do l=0,lmax
                  do j=0,lmax
                     en1=en1+conjg(psi(j,ix,iy))*psi(l,ix,iy)*v0(j,l,ix,iy)
                  enddo
               enddo
            enddo
         enddo

         if(mod(it,1)==0) print 2004,times,en,en0,en1,en+en0+en1
         
!wavefunction snapshot
         if(mod(it,10)==1) then
            do ix=1,iv
               do iy=1,iv
                  write(1000+it/10,'(100E15.6)') r(ix+1),r(iy+1),abs(psi(0:lmax,ix,iy))/(si(ix)*si(iy))
               enddo
               write(1000+it/10,*)
            enddo
         endif

! density projection (Spectrum)
         spec=(0.d0,0.d0)
         do L=0,1
            do ix=1,iv
               spec=spec+sum(conjg(psi0(L,ix,:)*psi(L,ix,:)))
            enddo
         enddo
         write(516,'(3f15.6)') times,spec
         
      enddo
      
      
      write(*,*), 'End of Program! Done!',iu
      write(*,'(a6,f15.5)') 'rmax=',r(ubound(r,1))

200  Format(1x,100f3.1)
2000 Format(1x,100f6.3)
2001 Format(1x,100f10.3)
2002 Format(1x,100E15.6)
2003 Format(1x,i4,100f12.5)
2004 Format(1x,100f11.5)
2010 Format(1x,100i12)
    End program hydrogen_fedvr



      
      function mycleb(j1,m1,j2,m2,j,m)
        use anglib       
        integer j1,m1,j2,m2,j,m,l
        real*8 mycleb
        l=mod(m+j1-j2,2)
        !mycleb=(-1)**L*sqrt(2.d0*j+1)*cof3j( dble(j1),dble(j2), dble(j),dble(m1),dble(m2),dble(m) )
        mycleb=cleb(2*j1,2*m1,2*j2,2*m2,2*j,2*m)
        return
      end function mycleb

      double precision function my3j(j1,j2,j,m1,m2,m)
        use anglib      
        integer j1,j2,j,m1,m2,m
!        real*8 cof3j
        my3j=(-1)**(j2-j1-m)*cleb(2*j1,2*m1,2*j2,2*m2,2*j,2*m)/dsqrt(2.d0*j+1)
        return
      end function my3j

      double precision function my6j(j1,j2,j3,j4,j5,j6)
        use anglib      
        integer j1,j2,j3,j4,j5,j6
!        real*8 cof6j
        my6j=sixj( 2*j1, 2*j2, 2*j3, 2*j4, 2*j5, 2*j6)
        return
      end function my6j
      
      double precision function my9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
        use anglib      
        integer j1,j2,j3,j4,j5,j6,j7,j8,j9
!        real*8 cof9j
        my9j=ninej( 2*j1, 2*j2, 2*j3, 2*j4, 2*j5, 2*j6, 2*j7, 2*j8, 2*j9)
        return
      end function my9j

      double precision function rl(x,y,l) !r_<^l/r_>^{l+1}
        real*8 x,y,rmin,rmax!,rl0
        integer:: l
       !rl=0.d0
        rmin=min(x,y)
        rmax=1.0d0/max(x,y)
        rl=rmin*rmax
        rl=rl**L!(l*1.d0)
        rl=rl*rmax
        return
      end function rl
