SUBROUTINE sample(pos_filename,nmo_start,nmo_end,nat,ns,n_samples)
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) The subroutine sample.f95 reduce the size of the trajectory. 
  ! 0b) Date     Programmer     version
  !  20181225    Gang Huang     3c
  !
  ! Declaration
  USE parameter_shared
  USE atom_module
  USE wannier_center_module
  USE traj
  IMPLICIT NONE

!==========
!parameters
!==========
  !from main
  INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start  
  INTEGER,INTENT(IN) :: nmo_end  
  INTEGER,INTENT(IN) :: nat ! number of atoms
  character(LEN=*), INTENT(IN) :: pos_filename

  integer,parameter :: rank=8              
  INTEGER :: iatom, imovie
  INTEGER, INTENT(OUT) :: n_samples  !n_samples = INT(nmo/ns)
  
  ! Initialization
  iatom = 0
  imovie =0
  i =0
  n_samples = 0

!=======================
!read in trajectory file 
!=======================
  open(10,file=trim(pos_filename))
  ! READ the first two steps to get information of the system
  ! Now starting read data
  n_samples = sampling_number(nmo_start,nmo_end,ns)
  allocate(wannier_center_info(nat,n_samples))
  allocate(sampled_movie(n_samples))
  allocate(sampled_time(n_samples))
  allocate(sampled_energy(n_samples))
  to_skip: IF (nmo_start == 0) THEN
    CALL read_traj_wannier_center(10,nmo_start,nmo_end,ns,nat,n_samples)
  ELSE IF (nmo_start > 0 .AND. nmo_start < nmo_end) THEN
    CALL skip_lines(10, (nat+2)*(nmo_start-1)) ! use CALL to rerun a SUBROUTINE
    CALL read_traj_wannier_center(10,nmo_start,nmo_end,ns,nat,n_samples)
  ELSE
    WRITE(*,*) 'Please enter a right value of starting step'
  END IF to_skip

  close(10)
  write(6,*) 'End of trajectory reading.'
!=============
!write in file
!=============
  sampled_pos_filename = 'traj_pos_sampled.dat'
  open(10,file=sampled_pos_filename)
  do i =1,n_samples
    write (10,'(I8)') nat
    WRITE(10,100) 'i =',i-1,', time =',sampled_time(i),', E =',sampled_energy(i)
    100 FORMAT (1X,A3,I10,A8,F10.3,A5,F20.10)
    DO iatom = 1, nat
       WRITE(10,*) TRIM(wannier_center_info(iatom, i)%wannier_center_name), &
         wannier_center_info(iatom,i)%coord(1), &
         wannier_center_info(iatom,i)%coord(2), &
         wannier_center_info(iatom,i)%coord(3)
    ENDDO
  enddo
  write(6,*)'Sampled trajectory is written: ',sampled_pos_filename
  close(10)

  deallocate(sampled_movie, sampled_time,sampled_energy)

END SUBROUTINE sample
