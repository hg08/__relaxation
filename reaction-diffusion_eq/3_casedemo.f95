      program casedemo
              implicit none
              integer :: i
              integer, parameter :: low=3, high=5

              do i =1, 10

                  select case (i)
                  case (high+1:)
                      print*,i,"is greater than ",high
                  case (:low)
                      print*,i,"is less than ",low
                  case default
                      print*,i,"is nothing special."
                  end select
               end do
        end program casedemo
