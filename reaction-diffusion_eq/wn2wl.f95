      !===================================================================
      ! this function convert wave number (in cm-1) into wavelength (in nm)
      !===================================================================

      module wn2wl 
        implicit none

        interface n2l 
          module procedure rn2l
        end interface

       contains

        real function rn2l(a) !actual function. one can define function in a module.
          real, intent(in):: a
          rn2l = 10000000.0/a
        end function rn2l

       end module wn2wl

      !=============
      !Main function
      !=============
      program generic  ! this function convert wave number (in cm-1) into wavelength (in nm)
         use wn2wl 
         implicit none
         real(kind=4) :: a 

         write(6,*)'What is the wavenumber of the radiation (in cm-1):'
         read(5,*) a

         print*, 'Wavelength:',rn2l(a),'nm' ! real arguments

      end program generic
