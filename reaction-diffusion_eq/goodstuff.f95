      module goodstuff
        implicit none

        interface apbxc 
          module procedure rapbxc,iapbxc
        end interface

      contains

        real function rapbxc(a,b,c) !actual function. one can define function in a module.
          real, intent(in):: a,b,c
          rapbxc = a+b*c
        end function rapbxc

        integer function iapbxc(a,b,c)
          integer,intent(in):: a,b,c
          iapbxc = a+b*c
        end function iapbxc

        end module goodstuff

        !-------------------
        program generic
          use goodstuff
          implicit none
          real :: a= 1.2,b=3.4,c=5.6
          integer :: i=2,j=4,k=6

          print*,apbxc(a,b,c) ! real arguments
          print*,apbxc(i,j,k) ! integer arguments

        end program generic
