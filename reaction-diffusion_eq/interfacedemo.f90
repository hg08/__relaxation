      program interfacedemo
              implicit none

              interface !interface block
                      integer function factorial(n)
                              implicit none
                              integer,intent(in) :: n
                      end function factorial
              end interface 

        integer :: n=0

        do while (n<1)
          print*,'Input a positive integer:'
          read*,n
        end do
        print*,n,'!=',factorial(n)
      end program interfacedemo

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
