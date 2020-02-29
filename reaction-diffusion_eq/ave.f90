       program ave
          use ari
            implicit none
            real(kind=4) :: a,b,c
            print*,'(a,$)',"Input value of a,b,c :"
            read*, a, b, c

            print*,arithmetic_ave(a,b,c) 
       end program ave
       !--------------

       module ari
                !no variables
       contains
           real function arithmetic_ave(a,b,c) 
              implicit none
              real(kind=4) :: a,b,c,m
                   m=(a+b+c)/3.
                   arithmetic_ave=m 
           end function arithmetic_ave

        end module ari
