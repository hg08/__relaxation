      !A type to contain several different variables at each grid point 
      program typedemo2
              implicit none

              ! slow-running code
              type gridpoint
                      real :: temperature, composition, vel(3)
              end type gridpoint

              type(gridpoint),allocatable:: grid(:,:)
              integer nx,ny

              read*,nx,ny
              allocate(grid(nx,ny))
              call random_number(grid%temperature)
              ! etc

              print*, grid%temperature
      end program typedemo2

