      program testDO
              real :: a 
              integer :: i

              do a = 0.,5.,0.5  ! real
                  print*,'a =',a
              enddo
              do i = 0,5  ! integer
                  print*,i
              enddo
       end program

              
