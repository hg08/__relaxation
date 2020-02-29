      MODULE tools_0
      IMPLICIT NONE

      CONTAINS
        ! subroutine for reading the index file
        SUBROUTINE read_index_file(index_filename,np,ndx1,ndx2,ndx3)
            IMPLICIT NONE
            
            ! from upper function
            character(LEN=200),INTENT(INOUT) :: index_filename
            INTEGER, INTENT(IN) :: np
            INTEGER,allocatable,dimension(:),INTENT(INOUT) :: ndx1,ndx2,ndx3
            
            ! local
            INTEGER :: k 

            ! Initialization
            k=0

            index_filename=trim(index_filename)
            OPEN(10,file=index_filename)
            do k=1,np
                 read(10,*)ndx1(k),ndx2(k),ndx3(k)
            enddo
            CLOSE(10)
        END SUBROUTINE read_index_file
      END MODULE tools_0
