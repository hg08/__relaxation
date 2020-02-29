      module coords
        implicit none

        type point 
          real :: x,y,z
        end type point

        interface operator (.distance.) !Defining new operator .distance.
          module procedure separation
        end interface 

      contains
        real function separation(a,b)
        !DO NOT NEED: real :: dist    ! the output,i.e., the result
         type (point),intent(in) :: a,b ! inten(in) is important here
        !WRONG: dist = sqrt((a%x-b%x)*(a%x-b%x)+ &
        !WRONG:   (a%y-b%y)*(a%y-b%y)+(a%z-b%z)*(a%z-b%z))
         separation = sqrt((a%x-b%x)*(a%x-b%x)+ &
           (a%y-b%y)*(a%y-b%y)+(a%z-b%z)*(a%z-b%z))
         end function separation
      end module coords  

      !----------------
      program distance 
        use coords
        real :: d
        type (point) :: p1,p2
        p1%x =3.2; p1%y=4.9; p1%z=5.0
        p2%x =5.3; p2%y=7.0; p2%z=4.9

        d = p1 .distance. p2
        print*,d

      end program distance
