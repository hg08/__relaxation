       program deriv1
               implicit none
               integer :: n,i
               real,allocatable :: y(:),dydx(:)
               real :: x,dx

               write(*,'(a,$)') 'Input no. of grid points:'; read*,n
               allocate (y(n),dydx(n)) ! allocate grid arrays
               dx = 10.0/(n-1)
               do i =1,n
                 x = (i-1)*dx
                 y(i)=cos(x)
               end do

               call derivative (y,n,dx,dydx)
               
               open(10,file='deriv1_out.dat')
               do i = 1,n
                 x = (x-1)*dx
                 write(10,*) i*1., dydx(i), -sin(x), -sin(x)-dydx(i)
               end do
               close(10)

              deallocate(y,dydx)


        contains
               subroutine derivative (a,np,h,aprime)
                    integer, intent(in) :: np
                    real, intent(in) :: a(np),h
                    real, intent(out) :: aprime(np)
                    integer :: i

                    do i =1, np-1
                      aprime(i) = (a(i+1)-a(i))/h
                    end do
                    aprime(np) =0.
                end subroutine derivative
        end program deriv1
