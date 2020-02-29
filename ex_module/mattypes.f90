      MODULE mattypes
      TYPE sparse
      DOUBLE PRECISION, POINTER :: A(:)
      INTEGER, POINTER :: irow(:)    !indexing array
      INTEGER, POINTER :: jcol(:)    !indexing array
      INTEGER :: m,n
      INTEGER :: nnz
      !A: long vector with nonzero matrix entries
      !A is logically m times n
      !nnz: number of nonzeroes
      END TYPE sparse
      END MODULE mattypes
