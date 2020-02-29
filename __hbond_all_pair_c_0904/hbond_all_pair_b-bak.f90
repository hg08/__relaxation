!============
! 2015/09/04
! huang gang
!============

!=========================================      
! If there is no Nitrogen atom in system!
!=========================================
      program hbond_all_pair
      implicit none

!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! For specific filename to analyzing data
      character(LEN=20) :: pos_filename        ! For specific trajectory file name to analyzing data
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775            ! max number of atoms
      integer,parameter :: nmovie_max=300000   ! max number of movie
      real,parameter :: rooc=29.16             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25             ! rOH (3.5**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real,parameter :: r_ohc_min=0.25         ! rOH (.5**2)
     ! real,parameter :: r_ONc=9.0             ! r_ON (3.0**2)
      real           :: r23
      !real,parameter :: cosphic=0.866          ! 1.732/2; phiC=pi/6.
      !real(kind=rk),parameter :: delta_t=0.005 ! fs
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW
      !real(kind=rk),allocatable,dimension (:,:,:,:) :: r23 !r12, r13, r23
      ! real(kind=rk),allocatable,dimension (:,:)     :: r_ON 
      ! real(kind=rk),allocatable,dimension (:,:,:,:)  :: cosphi, pm
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_OW,&
                                                    ndx_H
      call system_clock(begin_time,rat) 
!==================
!read data in input
!==================
      open(10,file='input_hbond_all_pair_b')
      read(10,*)filename
      read(10,*)pos_filename
      read(10,*)nmovie     !number of movie steps: 1
      read(10,*)natoms     !number of atoms per molecules
      if (nmovie >= nmovie_max)then
        write(6,*)'!! nmovie> nmovie_max !!'
        write(6,*)'!! stop !!'
        stop
      endif
      
      if (natoms >= nmax)then
        write(6,*)'!! natoms> nmax !!' 
        write(6,*)'!! stop !!'
        stop
      endif  
     ! write(6,*)'number of movie steps:',nmovie
     ! write(6,*)'number of atoms in molecule:',natoms
      close(10)

      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))

!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
         read(10,*)                  !Neglect data of this line
         read(10,*)                  !Neglect data of this line
         i=0
         ii=0
         do iatom= 1,natoms
            read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
            if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
            elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
            endif
         enddo
      enddo
      i_OW=i
      i_H=ii
      write(6,*)i_OW,i_H
      close(10)
!=======================
!read in trajectory file 
!=======================

!      open(10,file='traj_pos.xyz')     
!      do imovie=1,nmovie
!         read(10,*)                  !Neglect data of this line
!         read(10,*)                  !Neglect data of this line
!         do iatom= 1,natoms
!            read (10,*)atom_type(iatom),x(iatom,imovie),& 
!                       y(iatom,imovie),z(iatom,imovie)
!         enddo
!      enddo
!      close(10)
!      write(6,*)'end of trajectory reading'

!=======================================================
! calculate the number of O (H) atoms in water molecules
!=======================================================

      allocate(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      i=0
      ii=0
      jj=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_OW(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          else
          endif
      enddo 
     !TEST
     !write(6,*) jj
      
      deallocate(atom_type)

!===========================key part====================================
      !allocate(r23(1,i_OW,i_H))

      open(20,file=trim(filename)//'_list.dat')
          
          do ii=1,i_H
              m3=ndx_H(ii)
              do i1=1, i_OW-1      ! No O atom can not be bonded to itself 
                  m1=ndx_OW(i1)
                  do i2=i1+1, i_OW 
                      m2=ndx_OW(i2)
                      do j =1, 1
                           r23= (x(m2,j)-x(m3,j))**2+  &
                                (y(m2,j)-y(m3,j))**2+  &
                                (z(m2,j)-z(m3,j))**2
                           ! write(20,*) m2
                            if (r23<r_ohc .and.    &
                                r23>r_ohc_min ) then
                                write(20,*) m1,m2,m3
                            endif
                            if ((x(m1,j)-x(m3,j))**2+  &  ! Actually, we do not need to define the array r23
                               (y(m1,j)-y(m3,j))**2+  &
                               (z(m1,j)-z(m3,j))**2 <r_ohc &
                               .and.                       & 
                               (x(m1,j)-x(m3,j))**2+  &
                               (y(m1,j)-y(m3,j))**2+  &
                               (z(m1,j)-z(m3,j))**2 >r_ohc_min)then
                                write(20,*) m2,m1,m3
                            endif
                      enddo
                  enddo
              enddo
          enddo
      close(20)

      deallocate(ndx_OW,ndx_H,&
                x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      !STOP
      END
