       program array_declarations
               implicit none

               real, dimension(5,5) :: a,b
               real :: c(3,5,7),d(-5:5), e(0:1) 
               real, allocatable:: f(:),g(:,:,:)
               integer n(3),i
               real :: sum1Darray

               write(*,'(a,$)') 'Input 3 array dimensions:' ! a : character
               read*,(n(i),i=1,3)
               allocate(f(n(1)),g(n(1),n(2),n(3))) ! f is a n-vecter
               print*,'The size of first array is:',size(f)
               print*,'The size of the second array is:',size(g)
               
               do i=1,n(1)
                 write(*,*)'The',i,'th elem. of arra for 1st Dim:'!a:character

                 read*,f(i)
               end do 

               print*,'The sum is',sum1Darray(f,n(1))        
               !main body of program

               deallocate (f,g)
        end program array_declarations

        !========
        real function sum1Darray (a,n)
                implicit none
                integer,intent(in):: n
                real,intent(in) :: a(n)
                integer i
                real :: sum=0

                do i=1,n
                  sum=sum +a(i)
                end do
                sum1Darray = sum
        end function sum1Darray
