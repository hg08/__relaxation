       program aveInterfacedemo
            implicit none

            interface
                    real function ariAve(a,b,c)
                            implicit none
                            real,intent(in):: a,b,c
                    end function ariAve
                    real function geoAve(a,b,c)
                            implicit none
                            real,intent(in) :: a,b,c
                    end function geoAve
                    real function harAve(a,b,c)
                            implicit none
                            real,intent(in) :: a,b,c
                    end function harAve
            end interface

            real(kind=4) :: a,b,c
            print*,"Input the values of a,b,c :"
            read*, a, b, c

            print*,'The arithmetic average:', ariAve(a,b,c)
            print*,'The geometric average:', geoAve(a,b,c)
            print*,'The harmonic average:', harAve(a,b,c)

       end program aveInterfacedemo

       !-------------------
       ! Arithmetic average
       !-------------------
           real function ariAve(a,b,c) 
              implicit none
              integer, parameter :: n=3
              real(kind=4),intent(in) :: a,b,c
              real :: m
                   m=(a+b+c)/n
                   ariAve=m
           end function ariAve

        !------------------
        ! Geometric average
        !------------------
           real function geoAve(a,b,c) 
              implicit none
              integer, parameter :: n=3
              real(kind=4),intent(in) :: a,b,c
              real :: m
                   m=(a*b*c)**(1./n)
                   geoAve=m
           end function geoAve

        !------------------
        ! Harmonic average
        !------------------
           real function harAve(a,b,c) 
              implicit none
              integer, parameter :: n=3
              real(kind=4),intent(in) :: a,b,c
              real :: m
                   m=n/(1/a +1/b +1/c)
                   harAve=m
           end function harAve
 
