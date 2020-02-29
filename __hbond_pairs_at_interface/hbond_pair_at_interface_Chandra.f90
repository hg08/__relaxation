!============
! 2015/09/04
! huang gang
!============

!===============================================      
!Additioanlly,based on hbond_all_pair_b-version;
!input read from keyboard!
!===============================================
      program hbond_pair_at_interface
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
      real,parameter :: upperbound=24.78       ! x_maximun (Unit: Angstrom)
      real,parameter :: lowerbound=-0.93       ! x_minimun (Unit: Angstrom)
      real,parameter :: thickness=2            ! Thickness of interface (Unit: Angstrom)
      real           :: r23
      integer    :: begin_time,end_time,rat,j,nmovie,natoms,iatom,& 
          imovie,m1,m2,m3,i_H,i1,i2,ii,&
          iii1,i_OInterf1,i_OInterf2
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_H,&
          ndx_OInterf1,ndx_OInterf2
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
          ii=0
          iii1=0
          !iii2=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if(trim(atom_type(iatom)) .eq. 'O' .and.&
                    x(iatom,imovie) .le. upperbound .and.&
                    x(iatom,imovie) .lt. upperbound-thickness ) then
                  iii1=iii1+1
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
              endif
          enddo
      enddo
      i_H=ii
      i_OInterf1=iii1
      !i_OInterf2=iii2
      write(6,*)'Number of HW     ',&
                'Number of OIntrf1'
      write(6,*)i_H,i_OInterf1
      close(10)
!=======================================================
! Calculate the number of O (H) atoms in water molecules
!=======================================================
      allocate(ndx_H(i_H))      ! this should be put after i_H is defined
      if (i_OInterf1 .gt. 0) then
          allocate(ndx_OInterf1(i_OInterf1))! this should be put after i_OInterf1 is defined
      endif    
      do imovie=1,nmovie
      ii=0
      iii1=0
      do iatom=1,natoms
          if(trim(atom_type(iatom)) .eq. 'O' .and.&
              x(iatom,nmovie) .le. upperbound .and.&
              x(iatom,nmovie) .lt. upperbound-thickness ) then
              iii1=iii1+1
              ndx_OInterf1(iii1)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          else
          endif
      enddo 
      enddo
      deallocate(atom_type)
!====================================
!Producing the interfacial O-O-H list      
!====================================
      if (i_OInterf1 .gt.2) then
      open(20,file=trim(filename)//'_Interf1_Chandra_list.dat')
          do i1=1, i_OInterf1-1! No O atom can not be bonded to itself 
              m1=ndx_OInterf1(i1)
              do i2=i1+1, i_OInterf1 
                  m2=ndx_OInterf1(i2)
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
      endif
      deallocate(ndx_H,ndx_OInterf1,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
