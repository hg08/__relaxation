!============
! 2018/01/18
! huang gang
!============

!===============================================      
!Additioanlly,based on hbond_all_pair_b-version;
!input read from keyboard!
!===============================================
      program hbond_all_pair
      implicit none
!========================
!Parameters and variables
!========================
      character(LEN=20) :: filename            ! specific filename to analyz data
      character(LEN=20) :: pos_filename        ! specific trajectory filename to analyz data
      integer,parameter :: rk=4              
      !real,parameter :: rooc=29.16             ! cutoff distance of rOO (5.4**2 )
      !real,parameter :: rohc=12.25             ! rOH (3.5**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real           :: r23
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW,i_I
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_I,ndx_OW,ndx_H
      call system_clock(begin_time,rat) 

!==================
!Read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajectory file:'
      read(5,*)pos_filename
      write(6,*)'What is the number of atoms in the system:'
      read(5,*)natoms     !number of atoms per molecules

      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)!Neglect data of this line
          read(10,*)                  
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
      write(6,*)'Number of I','    Number of OW','    Number of HW'
      write(6,*)i_I,i_OW,i_H
      close(10)
!=======================================================
! Calculate the number of O (H) atoms in water molecules
!=======================================================
      allocate(ndx_I(i_I))      
      allocate(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      allocate(ndx_H(i_H))      ! this should be put after i_H is defined
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
      deallocate(atom_type)
!========================
!Producing the O-O-H list      
!========================
      open(20,file=trim(filename)//'_Chandra_list.dat')
          do i1=1, i_I
              m1=ndx_I(i1)
              do i2=i1+1, i_OW 
                  m2=ndx_OW(i2)
                  do j =1, 1
                      do ii=1,i_H
                          m3=ndx_H(ii)
                          r23= (x(m2,j)-x(m3,j))**2+  &
                               (y(m2,j)-y(m3,j))**2+  &
                               (z(m2,j)-z(m3,j))**2
                          ! write(20,*) m2 !test
                          if (r23<r_ohc) then
                              write(20,*) m1,m2,m3
                          endif
                      enddo
                  enddo
              enddo
          enddo
      close(20)
      deallocate(ndx_I,ndx_OW,ndx_H,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
