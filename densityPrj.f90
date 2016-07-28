real*8,allocatable,dimension(:)::x,xp
complex*16,allocatable,dimension(:)::p,d

real*8 dx,dp,emax,emin
real*8:: Pi=3.141592632d0,InvPi,dil
integer n,ix,ne,nl,j
complex*16:: i=(0.d0,1.d0)
InvPi=1./Pi
n=1400
print*,'input: n,ne,emax,emin,nl'
read(*,*)      n,ne,emax,emin,nl
print*,        n,ne,emax,emin,nl
allocate(x(n),d(n),xp(ne),p(ne))
!emax=1.d0
open(201,file='fort.516')
do j=0,nl,2
do ix=1,n
   read(201,'(3f15.6)')  x(ix),d(ix)
  ! read(201,*) x(ix),dx,d(ix)
   if(mod(ix,100)==1) print*,dil,ix,d(ix)
enddo
read(201,*)
dx=x(2)-x(1)

de=(emax-emin)/ne
if(abs(ne*de).lt.1.D-8) stop 'Error::Emax-Emin=0'
do ie=1,ne
   xp(ie)=emin+ie*de
   p(ie)=0.d0
   do ix=1,n
      p(ie)=p(ie)+cdexp(i*xp(ie)*x(ix))*(d(ix))*dx*InvPi
   enddo
enddo

do ie=1,ne
   write(209,'(100f10.5)') dil, xp(ie),p(ie)
enddo
   write(209,*)

enddo

end program
