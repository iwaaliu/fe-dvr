! ------------------------------------------------------------------!
! This subroutine generates the Hamiltonian matrix in FE-DVR bases. !
! ------------------------------------------------------------------!
      subroutine hamiltonian_matrix
      use accuracy_real
      use globalmod_constants
      use globalmod_discrete
   
      implicit NONE
   
      integer*4                                   :: i,j,m,k,ic,ip,ir,mp, nsize,ierr,lwork,jstop,ijk,ij,ik,jk
      real(KIND=idp)                              :: wtemp,wsum,wv,fact,centrf
      real(KIND=idp), dimension(:,:), allocatable :: temp
      real(KIND=idp), dimension(:), allocatable   :: energy,work,a
!      integer*4,dimension(:),allocatable          ::ia,ja
      character(len=1)                            :: jobz,uplo
   
      INTERFACE
      subroutine lagrange_dvr_derivative(nfun,xbounds,pt,wt,df,ddf,ke_mat)
          use accuracy_real
          implicit none
          integer, INTENT(IN)                                   :: nfun
          real(KIND=idp), dimension(1:2), INTENT(IN)            :: xbounds
          real(KIND=idp), dimension(1:nfun), INTENT(OUT)        :: pt
          real(KIND=idp), dimension(1:nfun), INTENT(OUT)        :: wt   
          real(KIND=idp), dimension(1:nfun,1:nfun), INTENT(OUT) :: ke_mat,df,ddf
      end subroutine lagrange_dvr_derivative
      END INTERFACE


!------------------------------------------------------------------------------!
! Set up DVR weights, derivatives, and kinetic energy matrix for each element. !
! Note that in this step, non-normalized bases are introduced yet.             !
! NO bridge functions are introduced in this step. Each element is             !
! independent from others.                                                     !
!------------------------------------------------------------------------------!
      do i = 1, no_fe
         call lagrange_dvr_derivative( no_gp, bounds(i,:), mat_reg(i)%pt, mat_reg(i)%wt,  & 
                                       mat_reg(i)%df, mat_reg(i)%ddf, mat_reg(i)%ke_mat )
      end do   



!-----------------------------------------------------------------------!
! Getting matrix elements in a given element. DVR bases are normalized, !
! but bridge functions are not full taken into account yet.             !
!-----------------------------------------------------------------------!
      DO i = 1, no_fe
        do k = 1,no_gp
          
          if(mat_reg(i)%wt(k) <= 0.0_idp)then
             write(idwrite,*) 'in subroutine : hamiltonian_matrix.f90'
             write(idwrite,*) 'weight factors <= 0 !!!'
             write(idwrite,*) 'i=',i,'  k=',k,' mat_reg(j,i)%wt(k)=',mat_reg(i)%wt(k)
             write(idwrite,*) 'I stopped !!!'
             return
          endif
          
             if((i>1).and.(k==1)) then
                mat_reg(i)%fac1(k) = 1.0_idp/SQRT( mat_reg(i)%wt(1)     & 
                                   + mat_reg(i-1)%wt(no_gp) )
             elseif((i<no_fe).and.(k==no_gp)) then
                mat_reg(i)%fac1(k) = 1.0_idp/SQRT( mat_reg(i)%wt(no_gp) &
                                   + mat_reg(i+1)%wt(1) )
             else
                mat_reg(i)%fac1(k) = 1.0_idp/SQRT( mat_reg(i)%wt(k) )
             end if  

             do m=1,no_gp
                  if((i>1).and.(m==1)) then
                     mat_reg(i)%fac2(m) = 1.0_idp/SQRT( mat_reg(i)%wt(1)     & 
                                        + mat_reg(i-1)%wt(no_gp) )  
                  elseif((i<no_fe).and.(m==no_gp )) then
                     mat_reg(i)%fac2(m) = 1.0_idp/SQRT(mat_reg(i)%wt(no_gp)  & 
                                        + mat_reg(i+1)%wt(1) )
                  else
                     mat_reg(i)%fac2(m) = 1.0_idp/SQRT( mat_reg(i)%wt(m) )
                  end if  
                  
      mat_reg(i)%ke_mat(k,m) = mat_reg(i)%ke_mat(k,m)              & 
                             * mat_reg(i)%fac1(k)*mat_reg(i)%fac2(m) 

             end do         
        end do
      END DO


!---------------------------------------------------------!
! The last step to get matrix elements of K.E. by taking  !
! the Bridge functions fully into account.                !
!---------------------------------------------------------!
      nsize = (no_gp-1) * no_fe-1

      ALLOCATE (tot_ke_mat(1:nsize,1:nsize), temp(1:nsize+1,1:nsize+1))
!      allocate(ia(nsize),ja(nsize))
!      ALLOCATE(A(10000000))
      tot_ke_mat = 0.0_idp
      temp = 0.0_idp

      ic = 0
      column: do ip = 1, no_fe
         do mp = 2, no_gp
            ic = ic + 1
            ir = 0
            row: do i = 1, no_fe
               do m = 2, no_gp
                  ir = ir + 1
                  if (i == ip-1 .and. m == no_gp) then
                     temp(ir,ic) = mat_reg(i+1)%ke_mat(1,mp)
                  else if (i == ip .and. (m < no_gp .or. mp < no_gp .or. i == no_fe)) then
                     temp(ir,ic) = mat_reg(i)%ke_mat(m,mp)
                  else if (i == ip .and. m == no_gp .and. mp == no_gp .and. i /= no_fe) then
                     temp(ir,ic) = mat_reg(i)%ke_mat(no_gp,mp) + mat_reg(i+1)%ke_mat(1,1)
                  else if (i == ip+1 .and. mp == no_gp) then
                     temp(ir,ic) = mat_reg(i)%ke_mat(m,1)
                  end if
               end do
            end do row
         end do
      end do column
      

      do ic = 1, nsize
         do ir = 1, nsize
            tot_ke_mat(ir,ic) = temp(ir,ic)
         end do
      end do




!         centrf = 0.5_idp * dlang * (dlang + 1.0_idp)
! only update the diagonal elements:
      do i = 1, nsize
         tot_ke_mat(i,i) = tot_ke_mat(i,i) - cnuclear / r(i+1)
      end do


!      do ic=1,20
!      write(*,'(100f10.3)'),tot_ke_mat(1:16,ic)
!      enddo


      return


      ALLOCATE( energy(nsize) )

!
! set jobz='N' for only eigen-values;
! set jobz='V' for both eigen-value and eigen-vectors.
!
      jobz = 'V'
      uplo = 'U'
      lwork = 3*nsize-1

      ALLOCATE( work(lwork) )

      call  DSYEV(jobz,uplo,nsize,tot_ke_mat,nsize,energy,work,lwork,ierr)

      if(ierr == 0) then
         write(idwrite,*)
         write(idwrite,*)
         write(idwrite,'(1x,a)') '# Eigenvalues of 1st min(40,nsize),analytical,deff:'
!         write(idwrite,'(1x,es22.15)') (energy(i),i=1,nsize)
         do i=1,min(40,nsize)
            write(idwrite,'(i3,10f20.10)')i, energy(i),- cnuclear**2*0.5/i/i,energy(i)/(-cnuclear**2*0.5/i/i)-1
         enddo
         write(idwrite,*)
         write(idwrite,*)
      else
         write(idwrite,*) '# Something wrong with the diagonalization of matrix !!!'
         write(idwrite,*) '# ierr = ',ierr
      end if
! stop

! Are you normalized to unit?
       wsum = DOT_PRODUCT( tot_ke_mat(:,1),tot_ke_mat(:,1) )
       write(idwrite,'(1x,a,f13.10)') '# Normalization factor (the lowest state) = ',wsum



! Getting the real wave functions, I need to make a transformation : 
      write(idwrite,'(1x,a)') '# wave function of the lowest state :'
      fact = 1.0_idp
      if( tot_ke_mat(1,1) < 0.0_idp ) fact = -1.0_idp
      
      k = 0
      do i = 1, no_fe
         jstop = no_gp
         if(i == no_fe) jstop = no_gp-1
         do j = 2, jstop
            k = k+1
            if(i /= no_fe.and.j == jstop) then
               wtemp = 1.0_idp/SQRT( mat_reg(i)%wt(j)      &
                    + mat_reg(i+1)%wt(1))  
            else
               wtemp = 1.0_idp/SQRT( mat_reg(i)%wt(j) )
            end if
            tot_ke_mat(k,1) = fact * wtemp * tot_ke_mat(k,1)
            
! this is only for the lowest s state: 
            wv = 2.0_idp * r(k+1) * EXP( -r(k+1) )
!      write(idwrite,'(1x,f10.6,4x,es22.15,4x,es22.15)') r(k+1), tot_ke_mat(k,1), wv
            
         end do
      end do
      
            print*, tot_ke_mat(1:5,1)

      DEALLOCATE( energy,work )
      
      return
    end subroutine hamiltonian_matrix
    
