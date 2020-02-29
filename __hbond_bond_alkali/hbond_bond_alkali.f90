!===========================
! 2015/10/27
! huang gang
!===========================
      program hbond_bond_alk
      implicit none
!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! specific filename to analyz data
      character(LEN=20) :: pos_filename        ! specific trajectory filename to analyz data
      integer,parameter :: rk=4              
      real,parameter :: rooc=12.25             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=6.0025            ! rOH (3.5**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real           :: r23,r13,r_shell
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_H,i_alk,i_N,i_I,&
                    i1,i2,ii,iii,i4,i5,i_OW,i_OW_shell_alk
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_OW,ndx_H,&
          ndx_OW_shell_alk,ndx_alk
      call system_clock(begin_time,rat) 
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajectory file:'
      read(5,*)pos_filename
      write(6,*)'What is the number of atoms in the system:'
      read(5,*)natoms!number of atoms per molecules
      write(6,*)'What is the distance r_LiO (r_shell):'
      read(5,*)r_shell

      nmovie=1!number of movie steps: 1
      r_shell =r_shell**2
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
          iii=0
          i4=0
          i5=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
              elseif (trim(atom_type(iatom)) .eq. 'Li' .or. &
                  trim(atom_type(iatom)) .eq. 'Na' .or.&
                  trim(atom_type(iatom)) .eq. 'K' ) then
                  iii=iii+1
              elseif(trim(atom_type(iatom)) .eq. 'I') then
                  i4=i4+1
              elseif(trim(atom_type(iatom)) .eq. 'N') then
                  i5=i5+1
              endif
          enddo
      enddo
      i_OW=i
      i_H=ii
      i_alk=iii
      i_I=i4
      i_N=i5
      write(6,*)'Number of Ow','    Number of alkali metal ions'
      write(6,*)i_OW,i_alk
      close(10)
!==================================
!Extracting ndx of alk cations
!and O (H) atoms in water molecules
!==================================
      allocate(ndx_OW(i_OW))! this should be put after i_OW is defined
      allocate(ndx_H(i_H))! this should be put after i_H is defined
      allocate(ndx_alk(i_alk))
      i=0
      ii=0
      iii=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_OW(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          elseif (trim(atom_type(iatom)) .eq. 'Li' .or.&
              trim(atom_type(iatom)) .eq. 'Na' .or.&
              trim(atom_type(iatom)) .eq. 'K' ) then
                 iii=iii+1
                 ndx_alk(iii)=iatom
          endif
      enddo 
     deallocate(atom_type)
!=========================================
!Print index of O in the shell of alk ions
!=========================================
      i=0
      do i1=1, i_OW-1! No O atom can not be bonded to itself 
          m1=ndx_OW(i1)
          do i2=1, i_alk 
              m2=ndx_alk(i2)
                  do j =1,1
                      r23= (x(m2,j)-x(m1,j))**2+  &
                           (y(m2,j)-y(m1,j))**2+  &
                           (z(m2,j)-z(m1,j))**2
                  if (r23<r_shell) then
                    i=i+1 
                  endif
                  enddo
          enddo
      enddo
      i_OW_shell_alk=i
      allocate(ndx_OW_shell_alk(i_OW_shell_alk)) 
      i=0
      do i1=1, i_OW-1
          m1=ndx_OW(i1)
          do i2=1, i_alk 
              m2=ndx_alk(i2)
              do j =1,1
                  r23= (x(m2,j)-x(m1,j))**2+  &
                       (y(m2,j)-y(m1,j))**2+  &
                       (z(m2,j)-z(m1,j))**2
                  if (r23<r_shell) then
                    i=i+1 
                    ndx_OW_shell_alk(i)=m1
                  endif
              enddo
          enddo
      enddo
!=====================================
!Print the O-O-H list in the alk shell     
!=====================================
      open(20,file=trim(filename)//'_alk_shell_list.dat')
      do i1=1, i_OW_shell_alk-1! No O atom can not be bonded to itself 
          m1=ndx_OW_shell_alk(i1)
          do i2=i1+1, i_OW_shell_alk 
              m2=ndx_OW_shell_alk(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      deallocate(ndx_OW,ndx_OW_shell_alk,ndx_alk,&
                 ndx_H,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
