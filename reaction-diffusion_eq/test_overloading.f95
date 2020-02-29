       module coordstuff
         implicit none

         type point
           real :: x,y,z
         end type point

         interface operator (+)
           module procedure pointplus
         end interface
         interface operator (-)
           module procedure pointminus
         end interface

       contains
         function pointplus(a,b)
           type (point) :: pointplus
           type (point),intent(in) :: a,b
           pointplus%x = a%x + b%x
           pointplus%y = a%y + b%y
           pointplus%z = a%z + b%z
         end function pointplus

         function pointminus(a,b)
           type (point) :: pointminus
           type (point),intent(in) :: a,b
           pointminus%x = a%x - b%x
           pointminus%y = a%y - b%y
           pointminus%z = a%z - b%z
         end function pointminus
       end module coordstuff
       
       program test_overloading
         use coordstuff
         type(point) :: p1,p2,p3
         p1%x= 1.; p1%y=0.5; p1%z = 0.7 
         p2%x= 0.; p2%y=1.0; p2%z =0.7
         
         p3= p2 - p1
         print*,p3
         p3= p2 + p1
         print*,p3

       end program test_overloading

