      MODULE mathsparse
      USE mathtypes
      TYPE(sparse),PRIVATE :: hidden_sparse
      CONTAINS
      SUBROUTINE prod(x, z)
      DOUBLE PRECISION, POINTER :: x(:), z(:)
      END SUBROUTINE prod
      END MODULE mathsparse
