!=============================================================================
! Modules: 
! All modules can be written in one fortran file. The name of the fortran file 
! (my_module.f95) and the names of the module (eg, atom_module) is different. 
! They must be distinguished!
!=============================================================================
MODULE atom_module
!
! Purpose:
!   To define the derived data type for atom
IMPLICIT NONE
TYPE :: atom
  CHARACTER(LEN=2) :: atom_name
  REAL :: mass
  REAL, DIMENSION(3) :: coord 
END TYPE atom

!Here we use Remi's info_atoms type
!Definition of a structure to store the informations about common elements (already used atoms) and new elements (to define when the program runs)
type :: info_atoms
   character(len=3) :: name!Symbol
   double precision :: mw!Weight
   integer :: amount!Number of atoms "name"
end type info_atoms

!Information about a specific bond
!Used to store efficiently the main information about an atom which makes a bond
type :: info_bond
   integer :: index !Index of the bonded atom. 
   integer :: third !Index of the third atom which will alow to determine the molecule framework
   double precision, dimension(3) :: cos !Cosine of the angle between the bond and the X,Y,Z axis (in the lab frame). It is used to
   !determine the enlongated length $\Delta r$
   double precision :: dist !Distance square between the studied atom and the bonded atom (Q: bonded atom is the H? A: Yes) 
end type info_bond  

! The array atom_info can be shared by subroutines  
TYPE(atom), ALLOCATABLE, DIMENSION(:,:) :: atom_info
END MODULE atom_module

MODULE types
! Derived data type to stroe array
TYPE :: real_array
  CHARACTER(2) :: atom_name
  INTEGER :: time_step
  REAL(8), DIMENSION(3) :: coord_array
  TYPE (real_array), POINTER :: p
END TYPE real_array
! Derived data type to store water molecules around a center atom
TYPE :: water_array
  INTEGER :: time_step
  CHARACTER(2) :: center_atom_name
  REAL(8), DIMENSION(9) :: coord_water_array
  REAL(8) :: d_RO
  REAL(8), DIMENSION(2) :: theta 
  TYPE (water_array), POINTER :: p
END TYPE water_array
END MODULE types

MODULE parameter_shared
! Purpose:
!   To declare data to share between routines.
IMPLICIT NONE
SAVE 
!character(LEN=50) :: filename, pos_filename
character(LEN=50) :: sampled_pos_filename
!INTEGER :: nat ! number of atoms
INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
REAL, ALLOCATABLE, DIMENSION(:) :: sampled_time, sampled_energy
!CHARACTER(LEN=2) :: str_of_center_atoms
!INTEGER :: num_of_kind_center_atoms 
END MODULE parameter_shared

!============
! Interfaces 
!============
MODULE interface_definitions
INTERFACE
  SUBROUTINE norm(v1,v2,length)
    !Purpose:
    ! to calculate the norm of the vector x2-x1, where x1 and x2 are arrays
   REAL(kind=8), dimension(3), intent(in)::v1,v2
   REAL(kind=8),intent(inout) :: length 
  END SUBROUTINE norm
END INTERFACE

INTERFACE
  SUBROUTINE orientation_cos(x1,y1,z1,x2,y2,z2, cos_1,cos_2,cos_3)
    !Purpose:
    ! to calculate the orientation cos of the vector x2-x1, where x1 and x2 are arrays
    IMPLICIT NONE
    !data dict
    REAL(kind=8), intent(in)::x1,x2,y1,y2,z1,z2
    REAL(kind=8),intent(inout) :: cos_1, cos_2, cos_3
    REAL(kind=8) :: length
  END SUBROUTINE orientation_cos
END INTERFACE

INTERFACE
  !=============================================================
  !Function which checks if an atom is inside the desired volume
  !=============================================================
  logical function f_in_vol(a,b,c,x,y,z)
    REAL(kind=8), intent(in) :: a,b,c,x,y,z
  end function f_in_vol
END INTERFACE


!========================================
! a function to convert integer to string
!========================================
interface
  function str(k)
    !   "Convert an integer to string."
    integer, intent(in) :: k
    character(len=20) :: str
  end function str
end interface
END MODULE interface_definitions

MODULE wannier_center_module
!
! Purpose:
!   To define the derived data type for a wannier center
IMPLICIT NONE
TYPE :: wannier_center
  CHARACTER(LEN=2) :: wannier_center_name
  INTEGER :: molecular_id
  REAL(kind=8) :: charge
  REAL(kind=8), DIMENSION(3) :: coord 
  REAL(kind=8), DIMENSION(3) :: image_coord 
END TYPE wannier_center

! The array atom_info can be shared by subroutines  
TYPE(wannier_center), ALLOCATABLE, DIMENSION(:,:) :: wannier_center_info
END MODULE wannier_center_module


MODULE traj
!
! Purpose: 
! To declear data related to the simulation and traj.
IMPLICIT NONE
!INTEGER :: n_samples  !n_samples = INT(nmo/ns)
!INTEGER :: nmo_start, nmo_end  ! To get the total number of moves
INTEGER :: i_sample, i_input,i ! dummy index
CONTAINS
  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
    !
    ! Purpose:
    !  To calculate the total numbers of samples one want to include in their analysis.
    ! Data dictionary
    INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves

    write(*,*) 'In Function sampling_number, mo_end:', nmo_end
    ! no. of samples = INT({no. of moves}/ns)
    positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
      write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
    ELSE IF (nmo_end < nmo_start) THEN
      write(*,*) 'Please note that starting step shoud not larger than  ending step.'
    ELSE IF (ns ==0) THEN
      sampling_number = nmo_end-(nmo_start-1)
    ELSE IF (nmo_end-(nmo_start-1) <= ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns + 1)
    ELSE IF (nmo_end-(nmo_start-1) > ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns)
    END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples)
    !
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    USE atom_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    !TYPE(atom), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: atom_info
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,120) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
      120 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
      write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
          atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
        !WRITE (*,*) & 
        !atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
        !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj
  !similar to read_traj, but not extracting sampled_movie(i_sample), sampled_time(i_sample) or sampled_energy(i_sample)
  SUBROUTINE read_traj_wannier_center(indx,nmo_start,nmo_end,ns,nat,n_samples)
    !
    ! Purpose:
    ! To read info from the ion+center file (format: ion+center.xyz)
    USE wannier_center_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,*)!Neglect data of this line 
      write(*,*) 'the step:', imovie
      ! Loop for ALL atoms(or wannier centers)
      inner: do iatom= 1,nat
        read (indx,*) wannier_center_info(iatom, i_sample)%wannier_center_name, wannier_center_info(iatom,i_sample)%coord(1), & 
          wannier_center_info(iatom,i_sample)%coord(2), wannier_center_info(iatom,i_sample)%coord(3)
        ! Initialize the image_coord attribute: ALL the image coordinates are set to be its own coordinates
        wannier_center_info(iatom,i_sample)%image_coord(1) = wannier_center_info(iatom,i_sample)%coord(1) 
        wannier_center_info(iatom,i_sample)%image_coord(2) = wannier_center_info(iatom,i_sample)%coord(2)
        wannier_center_info(iatom,i_sample)%image_coord(3) = wannier_center_info(iatom,i_sample)%coord(3)
        ! Initialize the charge
        charge: if(TRIM(wannier_center_info(iatom, i_sample)%wannier_center_name) == "O") then 
          wannier_center_info(iatom,i_sample)%charge = 6 
        else if(TRIM(wannier_center_info(iatom, i_sample)%wannier_center_name) == "H") then 
          wannier_center_info(iatom,i_sample)%charge = +1 
        else
          wannier_center_info(iatom,i_sample)%charge = -2
        end if charge
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj_wannier_center

  SUBROUTINE skip_lines(indx, i_input)
    !
    ! Purpose: 
    ! To skip lines when read data from the input
    IMPLICIT NONE
    INTEGER :: i
    INTEGER,INTENT(IN) :: i_input,indx
    do i=1,i_input
       read(indx,*) !Neglect (nat+2)*(ns-1) lines
    enddo    
  END SUBROUTINE skip_lines
END MODULE traj
