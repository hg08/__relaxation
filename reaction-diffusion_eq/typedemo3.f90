      module gridDefinition
              implicit none
              type grid
                integer nx,ny
                real dx,dy  !spacing
                real,allocatable,dimension(:,:) :: T
                real,allocatable,dimension(:,:,:) :: V
              end type grid
      end module gridDefinition

      program typedemo3
              use gridDefinition
              implicit none
              type(grid) :: stuff ! declaration of variable

              print '(a,$)',"Input nx  and ny :"
              read*,stuff%nx,stuff%ny
              call initialise_grid (stuff)
              print*, stuff%T
              print*, stuff%V

      contains

              subroutine initialise_grid(a)
              implicit none
              type(grid),intent(inout) :: a

              a%dx=1./a%dx; a%dy=a%dx
              allocate(a%T(a%nx,a%ny),a%V(2,a%nx,a%ny))
              call random_number (a%T); a%V =0.
              end subroutine initialise_grid

      end program typedemo3
