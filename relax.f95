!================================================================================================
!  20150213     
!  Autocorrelation of relative positon vector:
!  Auto correlation fuction of auto correlation fucntion of O-H vector, O-O vector, etc.
!  In this function, we calculate the relaxiation of O-H bonds.
!  We find that the auto correlation function is alway larger a positive value in a relatively 
!  long time period, idicating that the O-H bond has never breaking in this such a short time period.  
!================================================================================================

program  relaxiation 
  implicit none
  
  !Declaration
  integer, parameter :: nmovie_max=300000 ! max number of movie steps

  integer :: i,iatom,imovie,maxt,mt,natoms,nhmaxt,nmovie,nts
  integer :: nbonds
  real :: hmaxt
  integer, dimension(:,:), allocatable :: idx_bonds
  double precision,dimension(:,:), allocatable :: x,y,z,rx,ry,rz
  double precision, dimension(nmovie_max) :: corr
  character(len=3), dimension(:), allocatable :: atom_type

  !Initialization
  i=0 ; iatom=0 ; imovie=0 ; maxt=0 ; mt=0 ; natoms=0
  nhmaxt=0 ; nmovie=0 ; nts=0
  corr(:)=0.d0

  !data in BOX_VEL
  open(10,file='BOX_VEL_relax')
  read(10,*)nmovie          ! number of movie steps
  read(10,*)natoms          ! number of atoms per  molecule
  read(10,*)nbonds          !number of bonds to study

  if(nmovie.ge.nmovie_max)then
     write(6,*)'!!! nmovie > nmovie_max !!!'
     write(6,*)'!!! stop !!!'
     stop
  endif

  !Allocation
  allocate (x(natoms,nmovie),y(natoms,nmovie),z(natoms,nmovie),&
       atom_type(natoms))
  allocate (idx_bonds(2,nbonds),&
       rx(nbonds,nmovie),ry(nbonds,nmovie),rz(nbonds,nmovie),&
       )

  !Filling the idx_bonds array with the idx of the atoms
  !making the bonds
  do i=1,nbonds
     read(10,*)idx_bonds(1,i),idx_bonds(2,i)
  end do

  write(6,*)'number of movie steps:',nmovie
  write(6,*)'number of atoms in molecule:',natoms
  write(6,*)
  close(10)

  ! read in TRAJECTORY/VELOCITY file from CP2K 
  open(10,file='traj_pos.xyz')

  do imovie = 1,nmovie
     ! write(6,*)'imovie',imovie
     read(10,*)
     read(10,*)

     do iatom = 1,natoms
        read(10,*)atom_type(iatom),&
             x(iatom,imovie),y(iatom,imovie),z(iatom,imovie)
     enddo

  enddo
  close(10)
  write(6,*)'end of TRAJECTORY reading'

  ! _______________________________
  ! Order of atoms :  N(1) O(2) O(3) O(4) (NO3-), O(5) H(6) H(7) (Wat1), O(8) H(9) H(10) (Wat2), O(11) H(12) H(13) (Wat3), K(14)

  ! maxt = nmovie
  hmaxt = nmovie/2        ! Case one 
  nhmaxt = int(hmaxt)     ! We only consider thr first half of the correlation function. This is kind of FILTER!
  ! nhmaxt = maxt            ! Case two  (NO filter!!)


  !First: Calculation for ALL the times of ALL the studied bonds
  do i = 1,nbonds 
     rx(i,:) = x(idx_bonds(2,i),:)-x(idx_bonds(1,i),:)
     ry(i,:) = y(idx_bonds(2,i),:)-y(idx_bonds(1,i),:)
     rz(i,:) = z(idx_bonds(2,i),:)-z(idx_bonds(1,i),:)
  enddo

  
  do mt=0,maxt-1!Loop of Dt (=t-t0)
     do nts = 1,maxt-mt!Loop of t0
        do i = 1,nbonds ! O(5) in wat1
           corr(mt+1) = corr(mt+1) + (&
                rx(i,nts)*rx(i,nts+mt) +&
                ry(i,nts)*ry(i,nts+mt) +&
                rz(i,nts)*rz(i,nts+mt)&
                )   !Scalar product of vectors
        end do
     end do
     corr(mt+1) = corr(mt+1)/dble(maxt-mt)!Normalization by the number of time steps
  end do
  !Normalization by corr(1)=<r(0).r(0)>
  do mt=maxt,1,-1
     corr(mt)=corr(mt)/corr(1)
  enddo

  open(10,file='correl_O-H1.dat')
  do i = 1,nhmaxt
     write(10,*)i-1,corr(i)
  enddo
  write(6,*)'correlation results written in file correlation'


  !Close and deallocation
  close(10)
  deallocate (x,y,z,rx,ry,rz,atom_type,idx_bonds)

END program
