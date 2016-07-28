!---  Utility timer routine (seconds).                                  
!---  uncomment the appropriate one and comment the others              
                                                                        
!---  SUN & SGI --                                                      
      double precision function clock() 
      real*4 etime, tm(2) 
      clock = etime( tm ) 
      END                                           
