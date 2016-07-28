!==================================================!
! This subroutine reads all the input parameters.  !
! Then r mesh points are set up.                   !
!==================================================! 
      subroutine readin(ierror)
        use accuracy_real
        use globalmod_discrete
        use globalmod_constants

        implicit NONE
        integer,                      intent(out) :: ierror
        integer                                   :: i,j,k
        real(KIND=idp)                            :: temp,rmax,rmin
        real(KIND=idp), dimension(1:2)            :: endpoints
        real(KIND=idp), dimension(:), allocatable :: scratch

        namelist/HYDROGEN/cnuclear,lang,lmax
        namelist/FEDVRDIS/no_gp,no_fe,bound_start,bound_end,bslope



      ierror = 0
! Getting 'no_gp' and 'no_fe'.
      read(idread,HYDROGEN)
         if(lang < 0) then
            write(idwrite,*) 'In subroutine : readin'
            write(idwrite,*) 'lang < 0 !!!'
            write(idwrite,*) 'lang = ',lang
            write(idwrite,*) 'I stopped !!!'
            ierror = 100
            return
         else
            dlang = dble(lang)
            write(idwrite,'(1x,a,f4.2)') '# Nuclear charge   = ',cnuclear
            write(idwrite,'(1x,a,i3)')   '# Angular momentum = ',lang
         end if


! Getting 'bound_start', 'bound_end', and 'bslope'.
      read(idread,FEDVRDIS)
         if( no_gp <= 0 .or. no_fe <= 0 ) then
             write(idwrite,*) 'In subroutine : readin'
             write(idwrite,*) 'no_gp <= 0 OR no_fe <=0 !!!'
             write(idwrite,*) 'no_gp = ',no_gp
             write(idwrite,*) 'no_fe = ',no_fe
             write(idwrite,*) 'I stopped !!!'
             ierror = 101
             return
         end if

         if(bound_start < 0.0_idp) then
             write(idwrite,*) 'In subroutine : readin'
             write(idwrite,*) 'bound_start < 0 !!!'
             write(idwrite,*) 'bound_start = ',bound_start
             write(idwrite,*) 'I stopped !!!'
             ierror = 102
             return
         end if

         if(bound_start >= bound_end) then
             write(idwrite,*) 'In subroutine : readin'
             write(idwrite,*) 'bound_start >= bound_end !!!'
             write(idwrite,*) 'bound_start = ',bound_start
             write(idwrite,*) 'bound_end   = ',bound_end
             write(idwrite,*) 'I stopped !!!'
             ierror = 103
             return
         end if


! Allocating arraies for each element.
! 'no_global' is the total No. of r mesh points
! (before throwing r_0 and r_max away).
      no_global = (no_gp-1) * no_fe + 1

      ALLOCATE( bounds(1:no_fe,1:2) )
      ALLOCATE( scratch(1:no_gp),r(no_global) )
      ALLOCATE( grids(1:no_gp),weights(1:no_gp) )
      ALLOCATE( mat_reg(1:no_fe) )
         do i = 1, no_fe
             ALLOCATE(  mat_reg(i)%ke_mat(no_gp,no_gp),      &
                        mat_reg(i)%eigvec_mat(no_gp,no_gp),  &
                        mat_reg(i)%pt(no_gp),                & 
                        mat_reg(i)%wt(no_gp),                &
                        mat_reg(i)%fac1(no_gp),              &
                        mat_reg(i)%fac2(no_gp),              &
                        mat_reg(i)%df(no_gp,no_gp),          &
                        mat_reg(i)%ddf(no_gp,no_gp),         &
                        mat_reg(i)%eigval_mat(no_gp))  
         end do
        

      endpoints(1) = -1.d0
      endpoints(2) = +1.d0

      call gaussq('legendre',no_gp,0.0_idp,0.0_idp,2,endpoints,scratch,grids,weights)

      write(idwrite,'(1x,a,i4)') '# No. of grid pints in each elements = ',no_gp
      write(idwrite,'(1x,a,i4)') '# No. of finite-elements             = ',no_fe
      write(idwrite,'(1x,a,i7)') '# No. of r mesh points (before throwing r_0 and r_max away) = ',no_global
      write(idwrite,*)
      write(idwrite,'(1x,a)') '# Grid poins     Weight factors in [-1,+1] :'
      write(idwrite,'(2x,f10.6,5x,f10.6)') ( grids(i),weights(i),i=1,no_gp )
      write(idwrite,*)


         temp = bound_end - bound_start
        
         i = 1
         bounds(i,1) = bound_start
         bounds(i,2) = bound_end
         do i = 2, no_fe
            bounds(i,1) = bounds(i-1,2)
!            bounds(i,2)=bounds(i,1)+temp
            bounds(i,2) = bounds(i,1) + temp * EXP( bslope*dble(i-1)/dble(no_fe) )
         end do

      
!--------------------------------
! Generating the r mesh points  !
!--------------------------------
         k=0
      do i = 1, no_fe
         rmin = bounds(i,1)
         rmax = bounds(i,2)
      do j = 1,no_gp
         k = k + 1
         r(k) = 0.5_idp * ( grids(j)*(rmax-rmin)+(rmax+rmin) )
             if(r(k) < 1.d-14) r(k) = 0.0_idp
      end do
         k = k - 1
      end do


      write(idwrite,'(1x,a)') '# Element                 r mesh points'
      write(idwrite,'(1x,a)') '# -------   -----------------------------------------------------------'

          j = 1
       do k = 1, no_fe
          write(idwrite,'(5x,i3,3x,5(f10.6,3x),5(/,11x,5(f10.6,3x)) )') &
     &          k,(r(i), i = j, j+no_gp-1)
          j = j + no_gp - 1
       end do

       write(idwrite,*)


      return
      end subroutine readin
