      program typedemo1
       
        implicit none

        type person
          character(len=20) :: name
          integer :: birthyear
          integer,allocatable :: childyear(:)
        end type person  

        type(person) :: beatle(4)

        beatle(1)%name    = "Gang"
        beatle(1)%birthyear  = 1986
        allocate(beatle(1)%childyear(2))
        beatle(1)%childyear(1)=2017
        beatle(1)%childyear(2)=2020

        beatle(2)%name    = "Fang"
        beatle(2)%birthyear  = 1990
        allocate(beatle(2)%childyear(2))
        beatle(2)%childyear(1)=2017
        beatle(2)%childyear(2)=2020
         
        print*,beatle(1)%name, &
               beatle(1)%birthyear,&
               beatle(1)%childyear
        print*,beatle(2)%name, &
               beatle(2)%birthyear,&
               beatle(2)%childyear
       end program typedemo1
