!============
! 2018/01/18
! huang gang
!============

!===================================      
!We list the index for OH...I bonds.      
!===================================
      program hbond_all_pair
      implicit none
!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! For specific filename to analyzing data
      character(LEN=20) :: pos_filename        ! For specific trajectory file name to analyzing data
      integer,parameter :: rk=4              
      real,parameter :: r_oic=21.16             ! cutoff distance of rOI(4.6**2 ) 
      real,parameter :: r_oic_min=0.25          ! Minimum of rOI (0.5**2 )
      real,parameter :: r_hic=10.24            ! rHI (3.2**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real,parameter :: r_ohc_min=0.25         ! rOH (.5**2)
      real,parameter :: r_ONc=9.0             ! r_ON (3.0**2)
      real           :: r23,r13
      real,parameter :: cosphic=0.866          ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: delta_t=0.005 ! fs
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW,i_I
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_OW,&
                                                    ndx_H,ndx_I
      call system_clock(begin_time,rat) 
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajectory file:'
      read(5,*)pos_filename
      write(6,*)'What is the number of atoms in the system:'
      read(5,*)natoms     !number of atoms per molecules

      nmovie=1     !number of movie steps: 1
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
         jj=0
         do iatom= 1,natoms
            read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
            if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
            elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
            elseif(trim(atom_type(iatom)) .eq. 'I') then
                  jj=jj+1
            endif
         enddo
      enddo
      i_OW=i
      i_H=ii
      i_I=jj
      write(6,*)'Number of Ow','    Number of Hw', '   Number of I'
      write(6,*)i_OW,i_H,i_I
      close(10)
!=======================================================
! calculate the number of O (H) atoms in water molecules
!=======================================================

      allocate(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      allocate(ndx_I(i_I))    ! this should be put after i_I is defined
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
          elseif(trim(atom_type(iatom)) .eq. 'I') then
                 jj=jj+1
                 ndx_I(jj)=iatom
          else
          endif
      enddo 
     !TEST
     !write(6,*) jj
      
      deallocate(atom_type)

!========================
!Producing the O-O-H list
!========================      
      open(20,file=trim(filename)//'_OHI_list.dat')
          do i1=1, i_OW
              m1=ndx_OW(i1)
              do i2=1, i_H 
                  m2=ndx_H(i2)
                  do j =1, 1
                      do ii=1,i_I
                          m3=ndx_I(ii)
                          r13=(x(m1,j)-x(m3,j))**2+  &
                              (y(m1,j)-y(m3,j))**2+  &
                              (z(m1,j)-z(m3,j))**2 
                          r23=(x(m2,j)-x(m3,j))**2+  &
                              (y(m2,j)-y(m3,j))**2+  &
                              (z(m2,j)-z(m3,j))**2 
                          if (r13 <r_oic .and. r23<r_hic .and. & 
                              r13 >r_oic_min)then
                              write(20,*) m1,m2,m3
                          endif
                      enddo
                  enddo
              enddo
          enddo
      close(20)

      deallocate(ndx_OW,ndx_H,ndx_I,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
