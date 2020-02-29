      program pointertype
              implicit none

              type array2D
                real, pointer :: a(:,:) ! 2D array
                real,pointer :: b(:) ! 1D array
              end type
              
              type(array2D),allocatable :: T(:) ! 1D array of 2D-arrays-plus-1D-array. T(i) has two components: a(n,n) and b(n), where n is a function of i. 
              integer i,n,ng

              print '(a,$)', "No. of multigrid levels:"
              read*,ng
              allocate (T(ng))  ! allocate no. of grids
              do i = 1,ng
                 n = 2**(i+1)
                 print*,'Allocating grid',n,n
                 print*,'Allocating grid',n
                 allocate (T(i)%a(n,n))
                 allocate (T(i)%b(n))
                 call random_number(T(i)%a)
                 call random_number(T(i)%b)
                 print*,'The ',i,'th element of T%a:'
                 print*,T(i)%a
                 print*,'The ',i,'th element of T%b:'
                 print*,T(i)%b
              enddo 
              
       end program pointertype
