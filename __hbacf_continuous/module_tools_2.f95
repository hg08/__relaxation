      MODULE tools_2
      IMPLICIT NONE

      CONTAINS
        ! subroutine for smoothing data 
        SUBROUTINE moving_window_average(arr,n_neighbors) 
            IMPLICIT NONE
            ! from upper function
            real,allocatable,dimension (:),INTENT(INOUT) :: arr 
            INTEGER,INTENT(IN) :: n_neighbors

            ! local
            INTEGER :: i, j, new_size 
            real,allocatable,dimension (:) :: arr_extend 
            INTEGER :: width
            REAL :: mean

            ! Initialization
            i=0;j=0;new_size=0; width=0; mean=0.0

            width = n_neighbors*2 +1
            new_size = size(arr)+2*n_neighbors 
            ALLOCATE(arr_extend(new_size))
            do i= 1, n_neighbors
                arr_extend(i) = arr(1)
            end do
            do i = n_neighbors +1, n_neighbors+size(arr)
                arr_extend(i) = arr(i-n_neighbors)
            end do 
            do i = n_neighbors+size(arr)+1, new_size
                arr_extend(i) = arr(size(arr))
            end do 

            ! Calculate the averaged array (after smoothing)
            DO i=1,size(arr)
                mean = 0.0
                DO j= 1,width
                    mean = mean + arr_extend(i+j-1)
                ENDDO
                mean = mean/width
                arr(i) = mean
            ENDDO
            DEALLOCATE(arr_extend)
        END SUBROUTINE moving_window_average 
      END MODULE tools_2
