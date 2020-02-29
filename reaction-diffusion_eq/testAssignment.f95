      module coords3
        implicit none

        type point 
          real :: x,y,z
        end type point

        interface assignment (=)
          module procedure absvec
        end interface

      contains
        
        subroutine absvec(a,b)  !No output, it is called a subroutine!!
          real,intent(out) :: a
          type(point), intent(in) :: b
          a = sqrt(b%x**2 + b%y**2 + b%z**2)
        end subroutine absvec

      end module coords3

      !---------------------
      program testAssignment
        use coords3
        real :: d
        type(point) :: p
        
        print*, 'Input the three compoents of p bitte:'
        read*, p%x,p%y,p%z
        a = p  ! Assign p to a real number a
        print*, 'absolute value of p is:', a
      end program testAssignment
