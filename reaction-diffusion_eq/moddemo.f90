      program moddemo
              use fact
        implicit none
        integer :: n=0

        do while (n<1)
          print*,'Input a positive integer:'
          read*,n
        end do
        print*,n,'!=',factorial(n)
      end program moddemo

      !----------
      module fact
              !no varialbes in this module
      contains
        integer function factorial(n)
              implicit none
              integer,intent(in) :: n
              integer :: i,a
              a = 1
              do  i=1,n
                  a=a*i
              enddo 
              factorial = a
        end function factorial

        end module fact
