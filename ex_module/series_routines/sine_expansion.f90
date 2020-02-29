      PROGRAM sine_expansion
      USE series_routines, ONLY : expand_sine
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.14159265359
      REAL, DIMENSION(100) :: x, series
      INTEGER :: i
      DO i=1,100 !*** Define the range of x values to use
          x(i)=2.0*pi*(i-1)/ 99
      END DO
      series=expand_sine(x) !**** sin(x) up to five terms
      !**** Prints the results to the screen
      open(20,file="sine.dat")
      write(20,'(2e12.4)') (x(i), series(i), i = 1, 100)
      close(20)
      END PROGRAM sine_expansion
