!===========================
! 2015/10/27
! Author: Huang Gang
!===========================
!===========================================================================      
! The function 'hbond_bond_zno.f90' is different from 'hbond_shell2_hal.f90',
! because there are O atoms in zno, we can print '*_znoO-Ow-Hw_list
! .dat' and '*_OW-OW-HW_list.dat'.
!===========================================================================      
      program hbond_zno_wat
      implicit none
!=========================================      
! For specified k pairs of molecules，
! to calculate the average and the auto-
! correlation function of the H-Bond 
! population operator:
! <h>,<h_d>,<hh>, <hh_d>;      
! S_HB(t),S_d_HB(t),
! C_HB(t),C_d_HB(t).
! Here we list all the pairs of molecules 
! which may  form H bonds.      

! input filename 'input': 
     ! 1st line of input: sytem name
     ! 2nd line: posistion filename
     ! 3rd line: number of atoms
!========================
!parameters and variables
!========================
      character(LEN=200) :: filename             ! For specific filename to analyzing data
      character(LEN=200) :: pos_filename         ! For specific filename to analyzing data
      integer,parameter :: rk=4              
      real,parameter :: rooc=5**2                ! cutoff distance of rOO (Angstrom;to extract water molecules in NO3--shell )
      !real,parameter :: rohc=2.45**2            ! rOH (Angstrom)
      real(kind=rk) :: r_ohc                         ! rOH (1.2 A)
      real(kind=rk) :: r_ohc_min                     ! rOH lower bound(0.5 A)
      real(kind=rk) :: r_ON,r23
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_O,i_H,&
                    i1,i2,ii,i_OW,i_O_zno
      real(kind=rk),allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_OW,&
          ndx_O_zno, ndx_H
      call system_clock(begin_time,rat) 

!==================
!read data in input
!==================
      write(6,*)'What is the name of the syetem:' 
      read(5,*)filename   ! system name
      write(6,*)'What is the name of the trajectory file:' 
      read(5,*)pos_filename
      write(6,*)'What is the total number of atoms in the system:' 
      read(5,*)natoms     !number of atoms per molecules

      nmovie=1
      r_ohc=1.2
      r_ohc_min=0.5
      r_ohc=r_ohc**2
      r_ohc_min=r_ohc_min**2
      
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
      i_O=i
      i_H=ii

      close(10)
!============================

      allocate(ndx_O(i_O))    ! this should be put after i_O is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      i=0
      ii=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_O(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          else
          endif
      enddo 
      
      deallocate(atom_type)
!==================================================
! calculate the number of O atom in water molecules
!==================================================
      do j=1,1

      i_OW=0
      i_O_zno=0
      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          r_ON=z(m1,j) 
          if ((r_ON >9 .and. r_ON < 11) .or. &
             (r_ON >29 .and. r_ON < 32)) then
              i=i+1
          endif
          if (r_ON<29 .and. r_ON >11) then
              ii=ii+1    
              !write(6,*) i_O_zno
          endif
      enddo
      i_O_zno=i
      i_OW=ii
      write(6,*) "O in zno surface    :", i_O_zno
      write(6,*) "O in water molecules:",i_OW

      allocate(ndx_OW(i_OW))
      allocate(ndx_O_zno(i_O_zno))

      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          r_ON=z(m1,j)
          if ((r_ON >9 .and. r_ON < 11) .or. &
             (r_ON >29 .and. r_ON < 32)) then
              i=i+1
              ndx_O_zno(i)=m1
          endif
          if (r_ON<29 .and. r_ON >11) then
              ii=ii+1    
              ndx_OW(ii)=m1
              !write(6,*) i_O_zno
          endif
       enddo

      enddo !j-loop
!==============================
!allocate(r23(1,i_OW,i_OW,i_H))
!==============================
      do j =1,1
      !================================
      !print the list znoO-OH pairs
      !================================
      open(20,file=trim(filename)//'_znoO-OW-HW_list.dat')
          
          do i1=1, i_O_zno!i_O_zno: Nr. of O atom in ZnO surface 
              m1=ndx_O_zno(i1)
              do i2=1, i_O_zno 
                  m2=ndx_OW(i2)
                  do ii=1,i_H
                      m3=ndx_H(ii)
                          r23= (x(m2,j)-x(m3,j))**2+  &
                               (y(m2,j)-y(m3,j))**2+  &
                               (z(m2,j)-z(m3,j))**2
                          if (r23<r_ohc .and.    &
                              r23>r_ohc_min ) then
                              write(20,*) m1,m2,m3
                          endif
                  enddo
              enddo
          enddo
      close(20)
      
      enddo!j-loop

      deallocate(ndx_O,ndx_OW,ndx_O_zno,ndx_H,&
                x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
