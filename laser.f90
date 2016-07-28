MODULE  laser
! Laser field calculation
  use accuracy_real, only: idp
  use globalmod_constants
  implicit none
  integer, save               :: nostep, nthread
  integer, save               :: pulse_type
  integer, save               :: nion_block
  integer,parameter           :: nophoton_max = 5
  character(LEN=2), save      :: eunits
  real(KIND=idp), save        :: peakfield, phoene, timeduration, ramptime
  real(KIND=idp), save        :: opticycl_au, timestep_au
  integer, ALLOCATABLE        :: ion_block(:)

  private
  public ele_laser_field
  public nostep, nthread, pulse_type, nion_block
  public nophoton_max, eunits, peakfield, phoene, timeduration
  public ramptime, opticycl_au, timestep_au, ion_block

contains
  subroutine ele_laser_field (what_time, efield)
!==============================================================!
! Calculate the time-dependent electric field at a given time. ! 
! what_time and timeduration are given in units of             !
! laser-field periods.                                         !  
!==============================================================!
    use accuracy_real, only: idp
    use globalmod_constants, only: pi, twopi
    real(idp), intent(in)   :: what_time     ! time
    real(idp), intent(out)  :: efield        ! laser electric field
    real(idp)               :: temp
    real(idp)               :: eXUV,eIR
    real(idp)               :: tauXUV,delay,phased,FFXUV,OmegaXUV,&
                               irtimeduration
    irtimeduration=timeduration
    OmegaXUV=24.6/27.2
    delay=5.d0
    phased=0.d0
    tauxuv=10.d0
    

    select case (pulse_type)
    case (1)                            ! === for sin^2 profile === !
       temp = SIN(what_time * pi / timeduration)
       efield = peakfield * temp * temp * COS(twopi * what_time)
    case (2)        ! === Linear ramp-on and ramp-off profile === !
       if (what_time <= ramptime) then
          temp = what_time / ramptime
       else if (ramptime < what_time .AND. &
            what_time < timeduration - ramptime) then
          temp = 1.0_idp
       else if (what_time >= timeduration - ramptime) then
          temp = (timeduration - what_time) / ramptime
       end if
       efield = peakfield * temp * SIN(twopi * what_time)
    case (3)       ! === for cos^2 profile === !  

       temp = COS(what_time * pi / irtimeduration)
       eIR = peakfield * temp * temp * COS(twopi*phoene/27.2 * what_time)
       FFXUV=peakfield*0.1
       eXUV=FFXUV*Dcos(OmegaXUV*(what_time-delay)+PHASED)&
            *Dsin(pi*((what_time-delay)+TauXUV*0.5d0)/TauXUV)**2
       
       if((what_time-delay).lt.-TauXUV/2. .or. (what_time-delay).gt.TauXUV/2.) then
	  eXUV=0.
       endif
       efield=eIR+eXUV
       
    end select
  end subroutine ele_laser_field
end MODULE laser
