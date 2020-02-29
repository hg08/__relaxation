      program declare
            implicit none

              character(len=15) :: names
              integer,parameter :: n=5
              real :: c= 5.4,ii
              integer :: i,m=39,sum_i
              real,dimension(1) :: a(1:10) ! a is 1D array
              real, allocatable :: b(:,:,:,:)
              do i =1,10
                a(i)=rand()-0.5
                print*,a(i)
              end do

              print*,'The nearest int. of',c,'is:',Nint(c)
              print*,'mod(',m,',',n,') is:', mod(m,n)
              do i = 1,10
                if (a(i) < 0)then
                  print*,'The first negative element is:',a(i)
                  exit
                end if
              end do
              sum_i=0
              do i = 12,124,2
                sum_i=sum_i+i
              end do
              print*,'The sum of even no. in [12, 124] is',sum_i
      end program declare
