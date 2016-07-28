      MODULE globalmod_constants
          use accuracy_real
          implicit none
          integer*4, parameter               ::  idread = 5, idwrite = 6
          real(KIND=idp)                     ::  epsilon_mach
          real(KIND=idp), parameter          ::  pi     = 3.141592653589793238_idp,   &
                                                 halfpi = 0.5_idp *pi             ,   &
                                                 quadpi = 0.25_idp*pi             ,   &
                                                 twopi  = 2.0_idp *pi 
          integer*4                          ::  lang,lmax
          real(KIND=idp)                     ::  cnuclear,dlang,dlmax
      END MODULE globalmod_constants
