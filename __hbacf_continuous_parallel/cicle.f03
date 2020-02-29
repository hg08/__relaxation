program cicle
        use math_module
        implicit none
        integer :: i
        real,dimension(10) :: r
        real ::  a
        do i=1,10
           r(i)=i
           a=pi*r(i)**2
           print *, a
        enddo 
end program cicle
