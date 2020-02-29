      MODULE series_routines !***** example module *****
      IMPLICIT NONE

      CONTAINS
      FUNCTION expand_sine(x) !*** Finds series expansion
      REAL, DIMENSION(:) :: x
      REAL, DIMENSION(SIZE(x)) :: expand_sine
      expand_sine=x-x**3/factorial(3)+x**5/factorial(5)-&
      x**7/factorial(7)+x**9/factorial(9)
      END FUNCTION expand_sine
      ! *******************************************************
      FUNCTION factorial(n)
      !*** calculates factorials
      INTEGER, INTENT(IN) :: n
      REAL :: factorial, a
      INTEGER :: i
      a=1.0
      DO i = 1,n
          a = a*i
      END DO
      factorial=a
      END FUNCTION factorial

      END MODULE series_routines
