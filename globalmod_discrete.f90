      MODULE globalmod_discrete
          use accuracy_real
          implicit none
          integer*4                                      :: no_gp,no_fe,no_global
          real(KIND=idp)                                 :: bound_start,bound_end,bslope
          real(KIND=idp), dimension(:),   allocatable    :: grids,weights,r
          real(KIND=idp), dimension(:,:), allocatable    :: tot_ke_mat,bounds

          TYPE dvr_mat
               real(KIND=idp), dimension(:,:), pointer   :: ke_mat,eigvec_mat,df,ddf
               real(KIND=idp), dimension(:),   pointer   :: pt,wt,fac1,fac2,eigval_mat
          END TYPE dvr_mat

          TYPE (dvr_mat), dimension(:), allocatable      :: mat_reg 
      END MODULE globalmod_discrete
