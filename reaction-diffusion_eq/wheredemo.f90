      program wheredemo
              implicit none

              real a(5,5)
              call random_number(a)
              where (a >0.5) a=1. ! single-line where
              print*, a
              print*

              call random_number(a); a=a-0.5

              where (a<0.0) ! where block
                      a=0.
              elsewhere
                      a=sqrt(a)
              end where

              print*,a
      end program wheredemo
