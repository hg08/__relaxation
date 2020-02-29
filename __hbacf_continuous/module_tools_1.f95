      MODULE tools_1
      IMPLICIT NONE

      CONTAINS
        ! subroutine for reading trajectory file
        SUBROUTINE read_trajectory(pos_filename,nat,nmo,nmo_start,ns,atom_type,x,y,z)
            IMPLICIT NONE
            ! from upper function
            character(LEN=200),INTENT(IN) :: pos_filename
            INTEGER, INTENT(IN) :: nat,nmo_start,nmo,ns
            real,allocatable,dimension (:,:),INTENT(INOUT) :: x,y,z
            character(LEN=3),INTENT(INOUT) :: atom_type  
            
            ! local
            INTEGER :: i, imovie, iatom
            i=0;imovie=0;iatom=0

            open(10,file=trim(pos_filename))     
            do i=1,(nat+2)*(nmo_start-1)
               read(10,*)!Neglect data of this line
            enddo
            do imovie=1,nmo
               read(10,*)!Neglect data of this line
               read(10,*)
               do iatom= 1,nat
                   read (10,*)atom_type,x(iatom,imovie),y(iatom,imovie),&
                              z(iatom,imovie)
               enddo

               do i=1,(nat+2)*(ns-1)
                   read(10,*) !Neglect (nat+2)*(ns-1) lines
               enddo

            enddo
            close(10)
            write(6,*) 'end of trajectory reading'
        END SUBROUTINE read_trajectory
      END MODULE tools_1
