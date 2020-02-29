      program secDeriv
        implicit none

        integer :: n,i
        real :: h
        real,allocatable :: y(:)

        print '(a,$)','Input grid spacing :'; read*,h
        print '(a,$)','Input number of points :'; read*,n

        allocate(y(n))
        !forall(i=1:n) y(i)=((i-1)*h)**4 ! y=x^4
        forall(i=1:n) y(i)=((i-1)*h)**2  ! (i-1)*h play the role of 'x'

        print*,'Second derivative =',d2(y,h)
        deallocate(y)

      contains
        function d2(f,h)
          real,intent(in) :: f(:),h
          real :: d2(size(f))
          integer :: i

          forall (i=2:size(f)-1) d2(i)=(f(i+1)+f(i-1)-2*f(i))/h**2
          d2(1)=0.
          d2(size(f))=0.
        end function d2

      end program SecDeriv
